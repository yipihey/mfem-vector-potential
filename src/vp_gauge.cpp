#include "vp_gauge.hpp"

#include <cmath>
#include <random>

namespace vp
{

/// B r = P0 AMG(P0 r), P0 = projection onto mean-free vectors: keeps the preconditioner symmetric on the
/// complement of the constant null space (the AMG coarse solve of a singular matrix can otherwise inject
/// arbitrary constants and stall CG; seen at p=1 where the whole hierarchy is tiny).
struct CoulombGauge::MeanFreePrec : public Solver
{
   HypreBoomerAMG &amg;
   MPI_Comm comm;
   long long nglob;
   mutable Vector tmp;
   MeanFreePrec(HypreBoomerAMG &a, MPI_Comm c, long long n) : Solver(a.Height()), amg(a), comm(c), nglob(n) {}
   void Strip(Vector &v) const
   {
      double s = 0.0;
      for (int i = 0; i < v.Size(); i++) { s += v(i); }
      MPI_Allreduce(MPI_IN_PLACE, &s, 1, MPI_DOUBLE, MPI_SUM, comm);
      s /= (double)nglob;
      for (int i = 0; i < v.Size(); i++) { v(i) -= s; }
   }
   void Mult(const Vector &x, Vector &y) const override
   {
      tmp = x;
      Strip(tmp);
      amg.Mult(tmp, y);
      Strip(y);
   }
   void SetOperator(const Operator &) override {}
};

CoulombGauge::CoulombGauge(const Spaces &s, const Operators &o) : sp(s), ops(o)
{
   const double t0 = WallTime();
   const FiniteElementSpace *nfes = sp.pmesh->GetNodalFESpace();
   const int q = nfes ? nfes->GetMaxElementOrder() : 1;
   ir_order = 2 * sp.p + 2 * q + 2;
   const IntegrationRule &ir = IntRules.Get(Geometry::CUBE, ir_order);
   // PA mass operator for the (cheap, low-memory) applications ...
   Mnd.reset(new ParBilinearForm(sp.ND.get()));
   auto *vi = new VectorFEMassIntegrator;
   vi->SetIntRule(&ir);
   Mnd->AddDomainIntegrator(vi);
   Mnd->SetAssemblyLevel(AssemblyLevel::PARTIAL);
   Mnd->Assemble();
   Array<int> ess;
   Mnd->FormSystemMatrix(ess, Mh);
   // ... and a temporary fully assembled copy for L = G^T M G
   {
      ParBilinearForm Mfull(sp.ND.get());
      auto *vi2 = new VectorFEMassIntegrator;
      vi2->SetIntRule(&ir);
      Mfull.AddDomainIntegrator(vi2);
      Mfull.Assemble();
      Mfull.Finalize();
      std::unique_ptr<HypreParMatrix> Mp(Mfull.ParallelAssemble());
      L.reset(RAP(Mp.get(), ops.G.get()));   // G^T M G
   }
   amg.reset(new HypreBoomerAMG(*L));
   amg->SetPrintLevel(0);
   cg.reset(new CGSolver(sp.pmesh->GetComm()));
   cg->SetOperator(*L);
   prec.reset(new MeanFreePrec(*amg, sp.pmesh->GetComm(), (long long)sp.H1->GlobalTrueVSize()));
   cg->SetPreconditioner(*prec);
   cg->SetAbsTol(0.0);
   cg->SetMaxIter(500);
   cg->SetPrintLevel(-1);
   t_setup = WallTime() - t0;
}

CoulombGauge::~CoulombGauge() = default;

void CoulombGauge::RemoveMean(Vector &v) const
{
   MPI_Comm comm = sp.pmesh->GetComm();
   double s = 0.0;
   for (int i = 0; i < v.Size(); i++) { s += v(i); }
   MPI_Allreduce(MPI_IN_PLACE, &s, 1, MPI_DOUBLE, MPI_SUM, comm);
   s /= (double)sp.H1->GlobalTrueVSize();
   for (int i = 0; i < v.Size(); i++) { v(i) -= s; }
}

void CoulombGauge::Rhs(const Vector &a_t, Vector &Ma, Vector &g) const
{
   Ma.SetSize(a_t.Size());
   Mh->Mult(a_t, Ma);
   g.SetSize(sp.H1->GetTrueVSize());
   ops.G->MultTranspose(Ma, g);
}

int CoulombGauge::SolveL(Vector &g, Vector &chi, real_t rtol) const
{
   RemoveMean(g);
   chi.SetSize(g.Size());
   chi = 0.0;
   cg->SetRelTol(rtol);
   cg->Mult(g, chi);
   if (cg->GetNumIterations() >= 500 && Mpi::Root())
   {
      std::cout << "CoulombGauge: CG did not converge (iters " << cg->GetNumIterations()
                << ", final rel norm " << cg->GetFinalRelNorm() << ")" << std::endl;
   }
   RemoveMean(chi);
   last_iters = cg->GetNumIterations();
   return last_iters;
}

real_t CoulombGauge::MassNormSq(const Vector &a_t) const
{
   Vector Ma(a_t.Size());
   Mh->Mult(a_t, Ma);
   return InnerProduct(sp.pmesh->GetComm(), a_t, Ma);
}

real_t CoulombGauge::CoulombResidual(const Vector &a_t) const
{
   Vector Ma, g;
   Rhs(a_t, Ma, g);
   return std::sqrt(InnerProduct(sp.pmesh->GetComm(), g, g));
}

real_t CoulombGauge::GradFraction(const Vector &a_t) const
{
   Vector Ma, g, chi;
   Rhs(a_t, Ma, g);
   const real_t aa = InnerProduct(sp.pmesh->GetComm(), a_t, Ma);
   Vector g0(g);
   SolveL(g, chi, 1e-10);
   const real_t gp = InnerProduct(sp.pmesh->GetComm(), chi, g0);   // ||G chi||_M^2
   return std::sqrt(std::max(gp, (real_t)0.0) / aa);
}

GaugeFixStats CoulombGauge::Apply(Vector &a_t, real_t rtol, bool diagnostics) const
{
   MPI_Comm comm = sp.pmesh->GetComm();
   GaugeFixStats st;
   Vector Ma, g, chi, Gchi(a_t.Size());
   const double t0 = WallTime();
   Rhs(a_t, Ma, g);
   const real_t aa = InnerProduct(comm, a_t, Ma);
   Vector g0(g);
   st.iters = SolveL(g, chi, rtol);
   ops.G->Mult(chi, Gchi);
   a_t -= Gchi;
   st.time = WallTime() - t0;
   st.resid_before = std::sqrt(InnerProduct(comm, g0, g0));
   if (diagnostics)
   {
      Vector Mg(a_t.Size());
      Mh->Mult(Gchi, Mg);
      st.grad_fraction_before = std::sqrt(std::max(InnerProduct(comm, Gchi, Mg), (real_t)0.0) / aa);
      st.resid_after = CoulombResidual(a_t);
      st.grad_fraction_after = GradFraction(a_t);
   }
   return st;
}

GaugeFixStats CoulombGauge::JacobiGaugeSmooth(Vector &a_t, int nsweeps, bool diagnostics) const
{
   MPI_Comm comm = sp.pmesh->GetComm();
   GaugeFixStats st;
   const int nh = sp.H1->GetTrueVSize();
   if (dinv.Size() != nh)
   {
      Vector d(nh);
      L->GetDiag(d);
      dinv.SetSize(nh);
      for (int i = 0; i < nh; i++) { dinv(i) = 1.0 / d(i); }
      // power iteration for lambda_max(D^{-1} L)
      Vector x(nh), y(nh);
      std::mt19937 gen(12345);
      std::uniform_real_distribution<double> U(-1.0, 1.0);
      for (int i = 0; i < nh; i++) { x(i) = U(gen); }
      double lam = 1.0;
      for (int k = 0; k < 30; k++)
      {
         const double nx = std::sqrt(InnerProduct(comm, x, x));
         x /= nx;
         L->Mult(x, y);
         for (int i = 0; i < nh; i++) { y(i) *= dinv(i); }
         lam = InnerProduct(comm, x, y);
         x = y;
      }
      lmax = 1.05 * lam;
   }
   if (diagnostics) { st.grad_fraction_before = GradFraction(a_t); }
   Vector Ma, g, chi(nh), Gchi(a_t.Size());
   const double t0 = WallTime();
   Rhs(a_t, Ma, g);
   RemoveMean(g);
   Vector g0(g);
   // Chebyshev iteration (Saad, Alg. 12.1) for L chi = g, preconditioner D^{-1}, x0 = 0
   {
      const double lmin = lmax / 30.0;
      const double theta = 0.5 * (lmax + lmin), delta = 0.5 * (lmax - lmin), sigma = theta / delta;
      double rho = 1.0 / sigma;
      Vector r(g), z(nh), d(nh), Ad(nh);
      for (int i = 0; i < nh; i++) { z(i) = dinv(i) * r(i); d(i) = z(i) / theta; }
      chi = d;
      for (int k = 1; k < nsweeps; k++)
      {
         L->Mult(d, Ad);
         r -= Ad;
         for (int i = 0; i < nh; i++) { z(i) = dinv(i) * r(i); }
         const double rho_new = 1.0 / (2.0 * sigma - rho);
         for (int i = 0; i < nh; i++) { d(i) = rho_new * rho * d(i) + 2.0 * rho_new / delta * z(i); }
         chi += d;
         rho = rho_new;
      }
   }
   ops.G->Mult(chi, Gchi);
   a_t -= Gchi;
   st.time = WallTime() - t0;
   st.iters = nsweeps;
   st.resid_before = std::sqrt(InnerProduct(comm, g0, g0));
   if (diagnostics)
   {
      st.resid_after = CoulombResidual(a_t);
      st.grad_fraction_after = GradFraction(a_t);
   }
   return st;
}

real_t CoulombGauge::StiffnessMismatch() const
{
   const IntegrationRule &ir = IntRules.Get(Geometry::CUBE, ir_order);
   ParBilinearForm K(sp.H1.get());
   auto *di = new DiffusionIntegrator;
   di->SetIntRule(&ir);
   K.AddDomainIntegrator(di);
   K.Assemble();
   K.Finalize();
   std::unique_ptr<HypreParMatrix> Km(K.ParallelAssemble());
   std::unique_ptr<HypreParMatrix> Diff(Add(1.0, *L, -1.0, *Km));
   return MaxAbsEntry(*Diff) / MaxAbsEntry(*Km);
}

} // namespace vp
