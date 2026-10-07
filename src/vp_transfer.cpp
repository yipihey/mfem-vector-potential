// vp_transfer.cpp -- see vp_transfer.hpp.
#include "vp_transfer.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <sstream>

namespace vp
{

// ===========================================================================
// RemoteEvaluator
// ===========================================================================

RemoteEvaluator::RemoteEvaluator(ParMesh &M1, double bbox_rel_size_inc,
                                 double newt_tol)
   : mesh(&M1)
{
   MFEM_VERIFY(M1.GetNodes() != nullptr, "RemoteEvaluator: mesh has no nodes");
   finder.reset(new FindPointsGSLIB(M1.GetComm()));
   finder->Setup(M1, bbox_rel_size_inc, newt_tol);
   // Default boundary tolerance of MFEM is dist^2 <= 1e-8 (i.e. points up to
   // 1e-4 outside an element would be clamped to its border and accepted);
   // tighten so that such points go through the periodic-image retry.
   finder->SetDistanceToleranceForPointsFoundOnBoundary(1e-18);
   finder->SetDefaultInterpolationValue(0.0);
}

RemoteEvaluator::~RemoteEvaluator() { finder.reset(); }

namespace
{
// periodic image shifts, nearest (fewest non-zero components) first
std::vector<std::array<int, 3>> ImageShifts()
{
   std::vector<std::array<int, 3>> s;
   for (int nz = 1; nz <= 3; nz++)
   {
      for (int i = -1; i <= 1; i++)
         for (int j = -1; j <= 1; j++)
            for (int k = -1; k <= 1; k++)
            {
               const int c = (i != 0) + (j != 0) + (k != 0);
               if (c == nz) { s.push_back({i, j, k}); }
            }
   }
   return s;
}

long long GlobalSum(MPI_Comm comm, long long v)
{
   MPI_Allreduce(MPI_IN_PLACE, &v, 1, MPI_LONG_LONG, MPI_SUM, comm);
   return v;
}
}

long long RemoteEvaluator::Locate(const Vector &x)
{
   const double t0 = WallTime();
   MPI_Comm comm = mesh->GetComm();
   MFEM_VERIFY(x.Size() % 3 == 0, "Locate: need xyz triples");
   const int n = x.Size() / 3;
   Vector pts(3 * n);
   for (int i = 0; i < 3 * n; i++) { pts(i) = x(i) - std::floor(x(i)); }

   finder->FindPoints(pts, Ordering::byVDIM);
   std::vector<int> bad;
   {
      const Array<unsigned int> &code = finder->GetCode();
      for (int i = 0; i < n; i++) { if (code[i] == 2) { bad.push_back(i); } }
   }
   long long gbad = GlobalSum(comm, (long long)bad.size());
   long long nshift = 0;
   if (gbad > 0)
   {
      const auto shifts = ImageShifts();
      for (const auto &s : shifts)
      {
         Vector sub(3 * (int)bad.size());
         for (size_t j = 0; j < bad.size(); j++)
            for (int d = 0; d < 3; d++) { sub(3 * j + d) = pts(3 * bad[j] + d) + s[d]; }
         finder->FindPoints(sub, Ordering::byVDIM);
         const Array<unsigned int> &code = finder->GetCode();
         std::vector<int> still;
         for (size_t j = 0; j < bad.size(); j++)
         {
            if (code[j] == 2) { still.push_back(bad[j]); }
            else
            {
               for (int d = 0; d < 3; d++) { pts(3 * bad[j] + d) += s[d]; }
               nshift++;
            }
         }
         bad.swap(still);
         gbad = GlobalSum(comm, (long long)bad.size());
         if (gbad == 0) { break; }
      }
      // final search of the complete (shifted) point set: sets the internal
      // state used by Evaluate().
      finder->FindPoints(pts, Ordering::byVDIM);
   }
   nshifted = GlobalSum(comm, nshift);
   long long nb = 0, nbord = 0;
   {
      const Array<unsigned int> &code = finder->GetCode();
      for (int i = 0; i < n; i++)
      {
         if (code[i] == 2) { nb++; }
         else if (code[i] == 1) { nbord++; }
      }
   }
   nnotfound = GlobalSum(comm, nb);
   nborder = GlobalSum(comm, nbord);
   npts_local = n;
   set_id++;
   t_locate = WallTime() - t0;
   return set_id;
}

void RemoteEvaluator::Evaluate(const ParGridFunction &gf, Vector &vals)
{
   const double t0 = WallTime();
   const int vd = gf.VectorDim();
   finder->Interpolate(gf, vals, Ordering::byVDIM);
   MFEM_VERIFY(vals.Size() == npts_local * vd, "RemoteEvaluator: size mismatch");
   t_eval = WallTime() - t0;
}

// ===========================================================================
// Transfer
// ===========================================================================

const char *KindName(TransferKind k)
{
   switch (k)
   {
      case TransferKind::A_pt: return "A_pt";
      case TransferKind::A_int: return "A_int";
      case TransferKind::A_l2: return "A_l2";
      case TransferKind::B_pt: return "B_pt";
      case TransferKind::B_int: return "B_int";
      case TransferKind::B_l2: return "B_l2";
      case TransferKind::B_l2c: return "B_l2c";
   }
   return "?";
}

TransferKind ParseKind(const std::string &s)
{
   for (TransferKind k : {TransferKind::A_pt, TransferKind::A_int, TransferKind::A_l2,
                          TransferKind::B_pt, TransferKind::B_int, TransferKind::B_l2,
                          TransferKind::B_l2c})
   {
      if (s == KindName(k)) { return k; }
   }
   MFEM_ABORT("unknown transfer kind '" << s << "'");
   return TransferKind::A_pt;
}

bool IsARoute(TransferKind k)
{
   return k == TransferKind::A_pt || k == TransferKind::A_int || k == TransferKind::A_l2;
}

std::vector<TransferKind> ParseKindList(const std::string &csv)
{
   std::vector<TransferKind> v;
   if (csv == "all" || csv.empty())
   {
      return {TransferKind::A_pt, TransferKind::A_int, TransferKind::A_l2,
              TransferKind::B_pt, TransferKind::B_int, TransferKind::B_l2,
              TransferKind::B_l2c};
   }
   std::stringstream ss(csv);
   std::string tok;
   while (std::getline(ss, tok, ',')) { if (!tok.empty()) { v.push_back(ParseKind(tok)); } }
   return v;
}

int DefaultTransferNq(int p) { return p + 2; }

namespace
{
// Two-pass coefficient.  Collect pass: records the physical position of every
// point it is asked to evaluate (in call order) and returns 0.  Replay pass:
// returns the stored values in the same order.  The call order is that of
// vp::ProjectA/ProjectB (deterministic element/dof loops).
class PointCache : public VectorCoefficient
{
public:
   enum Mode { Collect, Replay } mode = Collect;
   std::vector<double> pts;      // byVDIM xyz
   const Vector *vals = nullptr;  // byVDIM
   size_t counter = 0;
   PointCache() : VectorCoefficient(3) {}
   using VectorCoefficient::Eval;
   void Eval(Vector &V, ElementTransformation &T, const IntegrationPoint &ip) override
   {
      V.SetSize(3);
      if (mode == Collect)
      {
         Vector x(3);
         T.Transform(ip, x);
         pts.push_back(x(0)); pts.push_back(x(1)); pts.push_back(x(2));
         V = 0.0;
      }
      else
      {
         MFEM_ASSERT(3 * (counter + 1) <= (size_t)vals->Size(), "PointCache: replay overrun");
         for (int c = 0; c < 3; c++) { V(c) = (*vals)(3 * counter + c); }
         counter++;
      }
   }
};

// Source values at the quadrature points of the target mesh.
class QuadCoef : public VectorCoefficient
{
public:
   const Vector *vals = nullptr;
   int np = 0;
   QuadCoef() : VectorCoefficient(3) {}
   using VectorCoefficient::Eval;
   void Eval(Vector &V, ElementTransformation &T, const IntegrationPoint &ip) override
   {
      V.SetSize(3);
      const long long idx = (long long)T.ElementNo * np + ip.index;
      for (int c = 0; c < 3; c++) { V(c) = (*vals)(3 * idx + c); }
   }
};

// Preconditioner z = P B P r with B the AMG V-cycle on Sa and P = I - w w^T/(w.w) the
// projection onto w^perp (w spans the kernel of Sa, S: D^T w = 0, i.e. the weights
// int phi_i of the L2 basis).  Without the projection the AMG coarse solve amplifies the
// roundoff component of the residual along the kernel.
class ProjectedPrec : public Solver
{
public:
   Solver &B;
   Vector w;
   real_t ww;
   MPI_Comm comm;
   mutable Vector rp;
   ProjectedPrec(Solver &B_, const Vector &w_, MPI_Comm c)
      : Solver(w_.Size()), B(B_), w(w_), comm(c), rp(w_.Size())
   { ww = InnerProduct(comm, w, w); }
   void SetOperator(const Operator &) override {}
   void Mult(const Vector &r, Vector &z) const override
   {
      rp = r;
      rp.Add(-InnerProduct(comm, r, w) / ww, w);
      B.Mult(rp, z);
      z.Add(-InnerProduct(comm, z, w) / ww, w);
   }
};

// Schur complement S = D M^{-1} D^T (matrix free; M^{-1} by inner CG).
class SchurOp : public Operator
{
public:
   const HypreParMatrix &D;
   CGSolver &Minv;
   mutable Vector w, z;
   mutable int inner_iters = 0;
   SchurOp(const HypreParMatrix &D_, CGSolver &Minv_)
      : Operator(D_.Height()), D(D_), Minv(Minv_), w(D_.Width()), z(D_.Width()) {}
   void Mult(const Vector &x, Vector &y) const override
   {
      D.MultTranspose(x, w);
      z = 0.0;
      Minv.Mult(w, z);
      inner_iters += Minv.GetNumIterations();
      D.Mult(z, y);
   }
};
}

struct Transfer::Impl
{
   TransferKind kind;
   RemoteEvaluator &ev;
   ParFiniteElementSpace &fes;
   TransferOptions opt;
   const Operators *ops;
   bool isA, isInt, isPt, isL2, isClean;
   int p, q, nq;
   long long set = -1;
   Vector ptsv;              // all evaluation points (byVDIM)
   Vector vals;

   // pt / int
   std::unique_ptr<ParGridFunction> tmp;
   PointCache pc;
   // L2
   QuadCoef qc;
   std::unique_ptr<ParLinearForm> lf;
   std::unique_ptr<ParBilinearForm> mass;
   OperatorHandle Mh;
   std::unique_ptr<OperatorJacobiSmoother> jac;
   std::unique_ptr<CGSolver> cg;
   Array<int> ess;
   const IntegrationRule *qir = nullptr;
   // clean
   std::unique_ptr<SchurOp> S;
   std::unique_ptr<CGSolver> cg_outer;
   std::unique_ptr<HypreParMatrix> Sa;     // D diag(M)^{-1} D^T (preconditioner matrix)
   std::unique_ptr<HypreBoomerAMG> amg;
   std::unique_ptr<ProjectedPrec> pprec;

   Impl(TransferKind k, RemoteEvaluator &e, const ParFiniteElementSpace &f,
        const Operators *o, const TransferOptions &op)
      : kind(k), ev(e), fes(const_cast<ParFiniteElementSpace &>(f)), opt(op), ops(o)
   {
      isA = IsARoute(k);
      isInt = (k == TransferKind::A_int || k == TransferKind::B_int);
      isPt = (k == TransferKind::A_pt || k == TransferKind::B_pt);
      isL2 = (k == TransferKind::A_l2 || k == TransferKind::B_l2 || k == TransferKind::B_l2c);
      isClean = (k == TransferKind::B_l2c);
      p = fes.GetMaxElementOrder();   // ND_p: p, RT_{p-1}: p
      const FiniteElementSpace *nfes = fes.GetParMesh()->GetNodalFESpace();
      q = nfes ? nfes->GetMaxElementOrder() : 1;
      nq = opt.nq > 0 ? opt.nq : DefaultTransferNq(p);
      if (isClean)
      {
         MFEM_VERIFY(ops != nullptr, "B_l2c needs the Operators (D_h) of the target mesh");
      }
   }
};

Transfer::Transfer(TransferKind k, RemoteEvaluator &src_eval,
                   const ParFiniteElementSpace &dst_space,
                   const Operators *dst_ops, const TransferOptions &opt)
   : kind(k), impl(new Impl(k, src_eval, dst_space, dst_ops, opt))
{
   Impl &I = *impl;
   MPI_Comm comm = I.fes.GetComm();
   const double t0 = WallTime();
   ParMesh &pm = *I.fes.GetParMesh();
   double tc = 0.0;

   if (I.isPt || I.isInt)
   {
      // collect pass
      const double c0 = WallTime();
      I.tmp.reset(new ParGridFunction(&I.fes));
      I.pc.mode = PointCache::Collect;
      const ProjMode mode = I.isInt ? ProjMode::Integrated : ProjMode::Pointwise;
      if (I.isA) { ProjectA(*I.tmp, I.pc, mode, I.nq); }
      else { ProjectB(*I.tmp, I.pc, mode, I.nq); }
      I.ptsv.SetSize((int)I.pc.pts.size());
      std::copy(I.pc.pts.begin(), I.pc.pts.end(), I.ptsv.GetData());
      std::vector<double>().swap(I.pc.pts);
      tc = WallTime() - c0;
   }
   else
   {
      // L2: points = quadrature points of the target mesh
      const double c0 = WallTime();
      const int order = opt.l2_order > 0 ? opt.l2_order : 2 * I.p + 2 * I.q + 2;
      I.qir = &IntRules.Get(Geometry::CUBE, order);
      const int np = I.qir->GetNPoints();
      I.qc.np = np;
      const int ne = pm.GetNE();
      I.ptsv.SetSize(3 * np * ne);
      Vector x(3);
      for (int e = 0; e < ne; e++)
      {
         ElementTransformation *T = pm.GetElementTransformation(e);
         for (int k = 0; k < np; k++)
         {
            T->Transform(I.qir->IntPoint(k), x);
            for (int d = 0; d < 3; d++) { I.ptsv(3 * ((long long)e * np + k) + d) = x(d); }
         }
      }
      tc = WallTime() - c0;

      // right hand side and (PA) mass operator
      I.lf.reset(new ParLinearForm(&I.fes));
      I.lf->AddDomainIntegrator(new VectorFEDomainLFIntegrator(I.qc, I.qir));
      I.mass.reset(new ParBilinearForm(&I.fes));
      auto *mi = new VectorFEMassIntegrator;
      mi->SetIntRule(I.qir);   // same rule as the rhs: exact identity on the same mesh
      I.mass->AddDomainIntegrator(mi);
      I.mass->SetAssemblyLevel(AssemblyLevel::PARTIAL);
      I.mass->Assemble();
      I.mass->FormSystemMatrix(I.ess, I.Mh);
      I.jac.reset(new OperatorJacobiSmoother(*I.mass, I.ess));
      I.cg.reset(new CGSolver(comm));
      I.cg->SetOperator(*I.Mh);
      I.cg->SetPreconditioner(*I.jac);
      I.cg->SetAbsTol(0.0);
      I.cg->SetMaxIter(opt.max_iter);
      I.cg->SetPrintLevel(-1);
      I.cg->SetRelTol(opt.l2_rtol);
   }

   if (I.isClean && I.opt.clean_amg)
   {
      // preconditioner of the outer CG: AMG on  Sa = D diag(M)^{-1} D^T
      Vector dg(I.fes.GetTrueVSize());
      I.mass->AssembleDiagonal(dg);
      std::unique_ptr<HypreParMatrix> Dt(I.ops->D->Transpose());
      Vector dinv(dg.Size());
      for (int i = 0; i < dg.Size(); i++) { dinv(i) = 1.0 / dg(i); }
      Dt->ScaleRows(dinv);
      I.Sa.reset(ParMult(I.ops->D.get(), Dt.get()));
      I.Sa->CopyRowStarts();
      I.Sa->CopyColStarts();
      I.amg.reset(new HypreBoomerAMG(*I.Sa));
      I.amg->SetPrintLevel(0);
      // kernel vector of S: w_i = int phi_i over the L2 basis
      ParLinearForm wl(I.ops->sp.L2.get());
      ConstantCoefficient one(1.0);
      wl.AddDomainIntegrator(new DomainLFIntegrator(one));
      wl.Assemble();
      Vector w(I.ops->D->Height());
      wl.ParallelAssemble(w);
      I.pprec.reset(new ProjectedPrec(*I.amg, w, comm));
   }

   stats.t_collect = tc;
   const double l0 = WallTime();
   I.set = I.ev.Locate(I.ptsv);
   stats.t_locate = WallTime() - l0;
   stats.npoints = GlobalSum(comm, (long long)I.ptsv.Size() / 3);
   stats.nnotfound = I.ev.NumNotFound();
   stats.nshifted = I.ev.NumShifted();
   stats.t_setup = WallTime() - t0;
}

Transfer::~Transfer() = default;

void Transfer::Apply(const ParGridFunction &src, ParGridFunction &dst)
{
   Impl &I = *impl;
   MPI_Comm comm = I.fes.GetComm();
   const double t0 = WallTime();
   stats.t_eval = stats.t_dofs = stats.t_solve = 0.0;
   stats.iters = stats.iters_l2 = stats.inner_iters = 0;
   if (I.ev.CurrentSet() != I.set)
   {
      I.set = I.ev.Locate(I.ptsv);   // the evaluator was used for another plan
   }
   double t = WallTime();
   I.ev.Evaluate(src, I.vals);
   stats.t_eval = WallTime() - t;

   if (I.isPt || I.isInt)
   {
      t = WallTime();
      I.pc.mode = PointCache::Replay;
      I.pc.vals = &I.vals;
      I.pc.counter = 0;
      const ProjMode mode = I.isInt ? ProjMode::Integrated : ProjMode::Pointwise;
      if (I.isA) { ProjectA(dst, I.pc, mode, I.nq); }
      else { ProjectB(dst, I.pc, mode, I.nq); }
      MFEM_VERIFY(I.pc.counter * 3 == (size_t)I.vals.Size(),
                  "Transfer: collect/replay point count mismatch");
      stats.t_dofs = WallTime() - t;
   }
   else
   {
      t = WallTime();
      I.qc.vals = &I.vals;
      I.lf->Assemble();
      Vector rhs(I.fes.GetTrueVSize());
      I.lf->ParallelAssemble(rhs);
      stats.t_dofs = WallTime() - t;
      t = WallTime();
      Vector x(I.fes.GetTrueVSize());
      x = 0.0;
      I.cg->SetRelTol(I.opt.l2_rtol);
      I.cg->Mult(rhs, x);
      stats.iters_l2 = I.cg->GetNumIterations();
      stats.iters = stats.iters_l2;
      MFEM_VERIFY(I.cg->GetConverged(), "Transfer: L2 mass CG did not converge");
      if (I.isClean)
      {
         MFEM_VERIFY(I.ops != nullptr && I.ops->D != nullptr, "missing D_h");
         const HypreParMatrix &D = *I.ops->D;
         // div before cleaning
         Vector r(D.Height());
         D.Mult(x, r);
         {
            real_t mx = r.Normlinf();
            MPI_Allreduce(MPI_IN_PLACE, &mx, 1, MPITypeMap<real_t>::mpi_type, MPI_MAX, comm);
            stats.div_before_clean = mx;
         }
         // inner solver on M_RT (tight), outer CG on S = D M^{-1} D^T
         CGSolver &inner = *I.cg;
         inner.SetRelTol(I.opt.clean_inner_rtol);
         if (!I.S)
         {
            I.S.reset(new SchurOp(D, inner));
            I.cg_outer.reset(new CGSolver(comm));
            I.cg_outer->SetOperator(*I.S);
            if (I.pprec) { I.cg_outer->SetPreconditioner(*I.pprec); }
            I.cg_outer->SetAbsTol(0.0);
            I.cg_outer->SetMaxIter(I.opt.max_iter);
            I.cg_outer->SetPrintLevel(-1);
         }
         I.cg_outer->SetRelTol(I.opt.clean_outer_rtol);
         // Defect correction around the (AMG-preconditioned) CG on S: the CG
         // criterion is the preconditioned residual norm, the Euclidean norm of
         // D x is only reduced to ~1e-9..1e-11 by one solve, so repeat with the
         // true residual until ||D x||_2 <= clean_outer_rtol ||D x_0||_2 (or an
         // absolute floor ~ roundoff of D x).
         const real_t dscale = MaxAbsEntry(D) * std::sqrt(InnerProduct(comm, x, x));
         const real_t r0n = std::sqrt(InnerProduct(comm, r, r));
         const real_t target = std::max(I.opt.clean_outer_rtol * r0n, 1e-13 * dscale);
         I.S->inner_iters = 0;
         Vector lam(D.Height()), w(D.Width()), z(D.Width());
         int rounds = 0;
         int outer_total = 0;
         for (;;)
         {
            D.Mult(x, r);
            const real_t rn = std::sqrt(InnerProduct(comm, r, r));
            if (rn <= target || rounds >= 6) { break; }
            I.cg_outer->SetAbsTol(0.0);
            lam = 0.0;
            I.cg_outer->Mult(r, lam);
            outer_total += I.cg_outer->GetNumIterations();
            D.MultTranspose(lam, w);
            z = 0.0;
            inner.Mult(w, z);
            I.S->inner_iters += inner.GetNumIterations();
            x -= z;
            rounds++;
         }
         stats.iters = outer_total;
         stats.clean_rounds = rounds;
         stats.inner_iters = I.S->inner_iters;
         inner.SetRelTol(I.opt.l2_rtol);
      }
      dst.SetFromTrueDofs(x);
      stats.t_solve = WallTime() - t;
   }
   stats.t_apply = WallTime() - t0;
}

// -- one-shot wrappers --------------------------------------------------------

namespace
{
TransferStats OneShot(TransferKind k, RemoteEvaluator &ev, const ParGridFunction &src,
                      ParGridFunction &dst, const Operators *ops, const TransferOptions &o)
{
   Transfer T(k, ev, *dst.ParFESpace(), ops, o);
   T.Apply(src, dst);
   return T.Stats();
}
}

TransferStats TransferA_Point(RemoteEvaluator &ev, const ParGridFunction &a1,
                              ParGridFunction &a2, const TransferOptions &o)
{ return OneShot(TransferKind::A_pt, ev, a1, a2, nullptr, o); }
TransferStats TransferA_Int(RemoteEvaluator &ev, const ParGridFunction &a1,
                            ParGridFunction &a2, const TransferOptions &o)
{ return OneShot(TransferKind::A_int, ev, a1, a2, nullptr, o); }
TransferStats TransferA_L2(RemoteEvaluator &ev, const ParGridFunction &a1,
                           ParGridFunction &a2, const TransferOptions &o)
{ return OneShot(TransferKind::A_l2, ev, a1, a2, nullptr, o); }
TransferStats TransferB_Point(RemoteEvaluator &ev, const ParGridFunction &b1,
                              ParGridFunction &b2, const TransferOptions &o)
{ return OneShot(TransferKind::B_pt, ev, b1, b2, nullptr, o); }
TransferStats TransferB_Int(RemoteEvaluator &ev, const ParGridFunction &b1,
                            ParGridFunction &b2, const TransferOptions &o)
{ return OneShot(TransferKind::B_int, ev, b1, b2, nullptr, o); }
TransferStats TransferB_L2(RemoteEvaluator &ev, const ParGridFunction &b1,
                           ParGridFunction &b2, const TransferOptions &o)
{ return OneShot(TransferKind::B_l2, ev, b1, b2, nullptr, o); }
TransferStats TransferB_L2Clean(RemoteEvaluator &ev, const ParGridFunction &b1,
                                ParGridFunction &b2, const Operators &ops2,
                                const TransferOptions &o)
{ return OneShot(TransferKind::B_l2c, ev, b1, b2, &ops2, o); }

void ComposeB(const Operators &ops, const Vector &b0_t, const Vector &a_t,
              Vector &b_t)
{
   b_t.SetSize(ops.C->Height());
   ops.C->Mult(a_t, b_t);
   b_t += b0_t;
}

void ComposeB(const Spaces &sp, const Operators &ops, VectorCoefficient &B0,
              const Vector &a_t, Vector &b_t)
{
   ParGridFunction b0(sp.RT.get());
   ProjectB(b0, B0, ProjMode::Integrated);
   Vector b0t;
   b0.GetTrueDofs(b0t);
   ComposeB(ops, b0t, a_t, b_t);
}

// ===========================================================================
// MeshDiag
// ===========================================================================

struct MeshDiag::Impl
{
   const Spaces &sp;
   const Operators &ops;
   const std::vector<SliceFamily> &slices;
   std::string field;
   FieldParams fp;
   Vector B0eff, b0_t;
   int ord = 0, ord_ex = 0;
   std::unique_ptr<VectorFunctionCoefficient> Bcoef;
   Field F;
   // gauge diagnostics (lazily built)
   bool gauge_built = false;
   std::unique_ptr<ParBilinearForm> Mnd, K;
   OperatorHandle Mh;
   std::unique_ptr<HypreParMatrix> Kmat;
   std::unique_ptr<HypreBoomerAMG> amg;
   std::unique_ptr<CGSolver> cgK;
   Array<int> ess;
   Impl(const Spaces &s, const Operators &o, const std::vector<SliceFamily> &sl)
      : sp(s), ops(o), slices(sl) {}

   void BuildGauge()
   {
      if (gauge_built) { return; }
      gauge_built = true;
      const FiniteElementSpace *nfes = sp.pmesh->GetNodalFESpace();
      const int q = nfes ? nfes->GetMaxElementOrder() : 1;
      const IntegrationRule &ir = IntRules.Get(Geometry::CUBE, 2 * sp.p + 2 * q + 2);
      Mnd.reset(new ParBilinearForm(sp.ND.get()));
      auto *vi = new VectorFEMassIntegrator;
      vi->SetIntRule(&ir);
      Mnd->AddDomainIntegrator(vi);
      Mnd->SetAssemblyLevel(AssemblyLevel::PARTIAL);
      Mnd->Assemble();
      Mnd->FormSystemMatrix(ess, Mh);
      K.reset(new ParBilinearForm(sp.H1.get()));
      auto *di = new DiffusionIntegrator;
      di->SetIntRule(&ir);
      K->AddDomainIntegrator(di);
      K->Assemble();
      K->Finalize();
      Kmat.reset(K->ParallelAssemble());
      amg.reset(new HypreBoomerAMG(*Kmat));
      amg->SetPrintLevel(0);
      cgK.reset(new CGSolver(sp.pmesh->GetComm()));
      cgK->SetOperator(*Kmat);
      cgK->SetPreconditioner(*amg);
      cgK->SetRelTol(1e-10);
      cgK->SetAbsTol(0.0);
      cgK->SetMaxIter(500);
      cgK->SetPrintLevel(-1);
   }
};

MeshDiag::MeshDiag(const Spaces &sp, const Operators &ops,
                   const std::vector<SliceFamily> &slices, const std::string &field,
                   const FieldParams &fp, const Vector &B0eff, const Vector &b0_t)
   : impl(new Impl(sp, ops, slices))
{
   Impl &I = *impl;
   I.field = field;
   I.fp = fp;
   I.B0eff = B0eff;
   I.b0_t = b0_t;
   I.F = MakeField(field, fp);
   const FiniteElementSpace *nfes = sp.pmesh->GetNodalFESpace();
   const int q = nfes ? nfes->GetMaxElementOrder() : 1;
   I.ord = 2 * sp.p + 2 * q + 2;
   I.ord_ex = I.ord;
   Vector B0 = B0eff;
   auto curl = I.F.curlA;
   I.Bcoef.reset(new VectorFunctionCoefficient(3, [curl, B0](const Vector &x, Vector &B)
   {
      curl(x, B);
      B += B0;
   }));
   normB_ex = L2NormVec(*sp.pmesh, *I.Bcoef, I.ord_ex);
   E_ex = vp::ExactEnergy(field, fp) + 0.5 * (B0eff * B0eff);
   hmean = MeshSize(*sp.pmesh).hmean;
}

MeshDiag::~MeshDiag() = default;

real_t MeshDiag::RTDistance(const Vector &b1_t, const Vector &b2_t)
{
   Vector d(b1_t);
   d -= b2_t;
   ParGridFunction g(impl->sp.RT.get());
   g.SetFromTrueDofs(d);
   return std::sqrt(InnerGF(g, g, impl->ord));
}

StateDiag MeshDiag::Evaluate(const Vector *a_t, const Vector &b_t, bool want_grad_frac)
{
   Impl &I = *impl;
   const Spaces &sp = I.sp;
   MPI_Comm comm = sp.pmesh->GetComm();
   StateDiag r;
   const bool prof = (std::getenv("VP_PROF") != nullptr) && Mpi::Root();
   double tp = WallTime();
   auto P = [&](const char *w)
   {
      if (prof) { const double t = WallTime(); std::cout << "    [diag] " << w << " " << t - tp << std::endl; tp = t; }
   };
   ParGridFunction bgf(sp.RT.get()), bcgf(sp.RT.get()), agf(sp.ND.get());
   bgf.SetFromTrueDofs(b_t);
   if (a_t)
   {
      agf.SetFromTrueDofs(*a_t);
      Vector bc(b_t);
      bc -= I.b0_t;
      bcgf.SetFromTrueDofs(bc);
   }
   // one pass over the quadrature points: ||b-B||^2, ||b||^2, (a,bc), ||a||^2, ||bc||^2
   {
      ParMesh &pm = *sp.pmesh;
      const IntegrationRule &ir = IntRules.Get(Geometry::CUBE, I.ord);
      DenseMatrix bv, bcv, av;
      Vector Bex(3);
      real_t s[5] = {0, 0, 0, 0, 0};   // err2, bb, a.bc, aa, bcbc
      for (int e = 0; e < pm.GetNE(); e++)
      {
         ElementTransformation *T = pm.GetElementTransformation(e);
         bgf.GetVectorValues(*T, ir, bv);
         if (a_t)
         {
            agf.GetVectorValues(*T, ir, av);
            bcgf.GetVectorValues(*T, ir, bcv);
         }
         for (int k = 0; k < ir.GetNPoints(); k++)
         {
            const IntegrationPoint &ip = ir.IntPoint(k);
            T->SetIntPoint(&ip);
            const real_t w = ip.weight * T->Weight();
            I.Bcoef->Eval(Bex, *T, ip);
            real_t e2 = 0, b2 = 0;
            for (int d = 0; d < 3; d++)
            {
               const real_t df = bv(d, k) - Bex(d);
               e2 += df * df;
               b2 += bv(d, k) * bv(d, k);
            }
            s[0] += w * e2;
            s[1] += w * b2;
            if (a_t)
            {
               real_t abc = 0, aa = 0, cc = 0;
               for (int d = 0; d < 3; d++)
               {
                  abc += av(d, k) * bcv(d, k);
                  aa += av(d, k) * av(d, k);
                  cc += bcv(d, k) * bcv(d, k);
               }
               s[2] += w * abc; s[3] += w * aa; s[4] += w * cc;
            }
         }
      }
      MPI_Allreduce(MPI_IN_PLACE, s, 5, MPITypeMap<real_t>::mpi_type, MPI_SUM, comm);
      r.errB_L2 = std::sqrt(s[0]);
      r.errB_rel = r.errB_L2 / normB_ex;
      r.energy = 0.5 * s[1];
      r.energy_rel_ex = (r.energy - E_ex) / E_ex;
      if (a_t)
      {
         r.helicity = s[2];
         r.a_L2 = std::sqrt(s[3]);
         r.hel_scale = std::sqrt(s[3]) * std::sqrt(s[4]);
      }
   }
   P("pass");
   const DivDiag dd = DivergenceDiag(I.ops, b_t);
   r.div_l2 = dd.l2;
   r.div_max = dd.maxabs;
   r.div_rel_l2 = dd.l2 / (normB_ex / hmean);
   r.div_rel_max = dd.maxabs / (normB_ex / hmean);
   P("div");
   real_t fmax = 0.0;
   for (const auto &fam : I.slices)
   {
      fmax = std::max(fmax, (real_t)std::fabs(SliceFlux(bgf, fam) - I.B0eff(fam.dir)));
   }
   r.flux_err = fmax;
   if (a_t)
   {
      I.BuildGauge();
      P("buildgauge");
      Vector Ma(a_t->Size()), g(sp.H1->GetTrueVSize());
      I.Mh->Mult(*a_t, Ma);
      I.ops.G->MultTranspose(Ma, g);
      r.coulomb_res = std::sqrt(InnerProduct(comm, g, g));
      if (want_grad_frac)
      {
         Vector phi(g.Size());
         phi = 0.0;
         I.cgK->Mult(g, phi);
         const real_t gp = InnerProduct(comm, phi, g);   // ||G phi||_M^2
         const real_t aa = InnerProduct(comm, *a_t, Ma); // ||a||_M^2
         r.grad_frac = std::sqrt(std::max(gp, (real_t)0.0) / aa);
         P("gradfrac");
      }
   }
   return r;
}

} // namespace vp
