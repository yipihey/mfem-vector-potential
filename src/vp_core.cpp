// vp_core.cpp -- see vp_core.hpp.
#include "vp_core.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <random>
#include <sstream>

namespace vp
{

// ===========================================================================
// 1. Mesh
// ===========================================================================

namespace
{
// Undeformed node coordinates keyed by mesh address.  Filled by
// MakePeriodicBox (overwriting stale entries from a freed mesh that lived at
// the same address) or, as a fallback, at the first DeformMesh call.
std::map<const ParMesh *, Vector> &X0Registry()
{
   static std::map<const ParMesh *, Vector> reg;
   return reg;
}
}

ParMesh *MakePeriodicBox(MPI_Comm comm, int N, int geom_order)
{
   MFEM_VERIFY(N >= 2, "MakePeriodicBox: N must be >= 2");
   MFEM_VERIFY(geom_order >= 1, "MakePeriodicBox: geometric order >= 1");

   Mesh m0 = Mesh::MakeCartesian3D(N, N, N, Element::HEXAHEDRON, 1.0, 1.0, 1.0);
   std::vector<Vector> translations(3);
   for (int i = 0; i < 3; i++)
   {
      translations[i].SetSize(3);
      translations[i] = 0.0;
      translations[i](i) = 1.0;
   }
   Mesh m = Mesh::MakePeriodic(m0, m0.CreatePeriodicVertexMapping(translations));

   const long long N3 = (long long)N * N * N;
   MFEM_VERIFY(m.GetNE() == N3, "periodic box: wrong number of elements");
   MFEM_VERIFY(m.GetNumFaces() == 3 * N3,
               "periodic box: #faces = " << m.GetNumFaces() << " != 3 N^3 = "
               << 3 * N3 << " (mesh is not fully periodic)");

   // Periodic meshes need discontinuous (L2) nodes.
   m.SetCurvature(geom_order, /*discont=*/true, 3, Ordering::byNODES);

   ParMesh *pm = new ParMesh(comm, m);
   X0Registry()[pm] = *pm->GetNodes();
   return pm;
}

namespace
{
// Normal (not normalized) of face f at its centre, and the centre point.
void FaceCentreNormal(ParMesh &pm, int f, Vector &x, Vector &nor)
{
   ElementTransformation *T = pm.GetFaceTransformation(f);
   IntegrationPoint ip;
   ip.Set2(0.5, 0.5);
   T->SetIntPoint(&ip);
   x.SetSize(3);
   T->Transform(ip, x);
   nor.SetSize(3);
   CalcOrtho(T->Jacobian(), nor);
}
}

std::vector<SliceFamily> RecordSlices(ParMesh &pm, int N,
                                      const std::vector<int> &js)
{
   std::vector<SliceFamily> fams;
   for (int dir = 0; dir < 3; dir++)
   {
      for (int j : js) { fams.emplace_back(); fams.back().dir = dir; fams.back().j = j; }
   }
   Vector x, nor;
   for (int f = 0; f < pm.GetNumFaces(); f++)
   {
      FaceCentreNormal(pm, f, x, nor);
      const real_t nn = nor.Norml2();
      for (auto &fam : fams)
      {
         const int d = fam.dir;
         if (std::fabs(nor(d)) < 0.999 * nn) { continue; }
         real_t fr = x(d) * N - fam.j;            // distance to plane in cells
         fr -= N * std::floor(fr / N + 0.5);      // periodic wrap
         if (std::fabs(fr) < 1e-8)
         {
            fam.faces.Append(f);
            fam.sign.Append(nor(d) > 0 ? 1 : -1);
         }
      }
   }
   for (auto &fam : fams)
   {
      long long n = fam.faces.Size();
      // count each face once globally: only faces whose tdofs we own are
      // used by SliceFlux; for the diagnostic count, the global number of
      // (possibly duplicated shared) entries is an upper bound.
      MPI_Allreduce(MPI_IN_PLACE, &n, 1, MPI_LONG_LONG, MPI_SUM, pm.GetComm());
      fam.nfaces_global = n;
   }
   return fams;
}

real_t SliceFlux(const ParGridFunction &B, const SliceFamily &fam)
{
   const ParFiniteElementSpace *fes = B.ParFESpace();
   Array<int> dofs;
   real_t flux = 0.0;
   static std::map<int, Vector> wcache; // 1D GL weights by n
   for (int i = 0; i < fam.faces.Size(); i++)
   {
      fes->GetFaceDofs(fam.faces[i], dofs);
      const int nd = dofs.Size();
      const int n = (int)std::lround(std::sqrt((real_t)nd));
      MFEM_VERIFY(n * n == nd, "SliceFlux: face dof count is not a square");
      // ownership: count the face once globally
      if (fes->GetLocalTDofNumber(FiniteElementSpace::DecodeDof(dofs[0])) < 0)
      {
         continue;
      }
      auto it = wcache.find(n);
      if (it == wcache.end())
      {
         const IntegrationRule &ir = IntRules.Get(Geometry::SEGMENT, 2 * n - 1);
         MFEM_VERIFY(ir.GetNPoints() == n, "unexpected GL rule size");
         Vector w(n);
         for (int k = 0; k < n; k++) { w(k) = ir.IntPoint(k).weight; }
         it = wcache.emplace(n, w).first;
      }
      const Vector &w = it->second;
      real_t s = 0.0;
      for (int k = 0; k < nd; k++)
      {
         const int dd = dofs[k];
         const int ld = FiniteElementSpace::DecodeDof(dd);
         const real_t sg = (dd < 0) ? -1.0 : 1.0;
         s += w(k % n) * w(k / n) * sg * B(ld);
      }
      flux += fam.sign[i] * s;
   }
   MPI_Allreduce(MPI_IN_PLACE, &flux, 1, MPITypeMap<real_t>::mpi_type, MPI_SUM,
                 fes->GetComm());
   return flux;
}

real_t SliceFluxQuadrature(const ParGridFunction &B, const SliceFamily &fam,
                           int qorder)
{
   ParMesh *pm = B.ParFESpace()->GetParMesh();
   const IntegrationRule &ir = IntRules.Get(Geometry::SQUARE, qorder);
   Vector v(3), nor(3);
   real_t flux = 0.0;
   for (int i = 0; i < fam.faces.Size(); i++)
   {
      FaceElementTransformations *FT =
         pm->GetFaceElementTransformations(fam.faces[i]);
      real_t s = 0.0;
      for (int k = 0; k < ir.GetNPoints(); k++)
      {
         const IntegrationPoint &ip = ir.IntPoint(k);
         FT->SetAllIntPoints(&ip);
         B.GetVectorValue(FT->Elem1No, FT->GetElement1IntPoint(), v);
         CalcOrtho(FT->Face->Jacobian(), nor);
         s += ip.weight * (v * nor);
      }
      flux += fam.sign[i] * s;
   }
   MPI_Allreduce(MPI_IN_PLACE, &flux, 1, MPITypeMap<real_t>::mpi_type, MPI_SUM,
                 pm->GetComm());
   return flux;
}

// ===========================================================================
// 2. Deformation
// ===========================================================================

void DisplacementField(int variant, const real_t *x, real_t *f)
{
   auto f1 = [&](real_t k, real_t *o)
   {
      o[0] = std::sin(k * x[1]) * std::sin(k * x[2]) / k;
      o[1] = std::sin(k * x[2]) * std::sin(k * x[0]) / k;
      o[2] = std::sin(k * x[0]) * std::sin(k * x[1]) / k;
   };
   real_t a[3], b[3];
   switch (variant)
   {
      case 1: f1(kTwoPi, f); break;
      case 2: f1(2 * kTwoPi, f); break;
      case 3:
         f1(kTwoPi, a); f1(2 * kTwoPi, b);
         for (int i = 0; i < 3; i++) { f[i] = a[i] + 0.5 * b[i]; }
         break;
      default: MFEM_ABORT("unknown deformation variant " << variant);
   }
}

void DeformMesh(ParMesh &pm, real_t eps, int variant)
{
   GridFunction *nodes = pm.GetNodes();
   MFEM_VERIFY(nodes, "DeformMesh: mesh has no nodes (call SetCurvature)");
   auto &x0map = X0Registry();
   auto it = x0map.find(&pm);
   if (it == x0map.end()) { it = x0map.emplace(&pm, *nodes).first; }
   MFEM_VERIFY(it->second.Size() == nodes->Size(),
               "DeformMesh: cached undeformed nodes have the wrong size");
   const Vector &x0 = it->second;
   const FiniteElementSpace *nfes = nodes->FESpace();
   const int nd = nfes->GetNDofs();
   MFEM_VERIFY(nfes->GetVDim() == 3, "nodes must be 3D");
   const bool bynodes = (nfes->GetOrdering() == Ordering::byNODES);
   for (int i = 0; i < nd; i++)
   {
      real_t x[3], f[3];
      for (int c = 0; c < 3; c++) { x[c] = x0(bynodes ? c * nd + i : 3 * i + c); }
      DisplacementField(variant, x, f);
      for (int c = 0; c < 3; c++)
      {
         (*nodes)(bynodes ? c * nd + i : 3 * i + c) = x[c] + eps * f[c];
      }
   }
   pm.NodesUpdated();
}

JacobianStatsResult JacobianStats(ParMesh &pm)
{
   JacobianStatsResult r;
   const FiniteElementSpace *nfes = pm.GetNodalFESpace();
   const int q = nfes ? nfes->GetMaxElementOrder() : 1;
   r.quad_order = 2 * q + 2;
   const IntegrationRule &ir = IntRules.Get(Geometry::CUBE, r.quad_order);
   // element corners
   std::vector<IntegrationPoint> corners(8);
   const IntegrationRule *cv = Geometries.GetVertices(Geometry::CUBE);
   for (int c = 0; c < 8; c++) { corners[c] = cv->IntPoint(c); }

   // detJ is normalised by the reference-cell volume 1/NE_global (the box has
   // volume 1 and every deformation is volume-preserving as a map of the
   // torus), so the undeformed mesh has detJ = 1.
   long long ne_glob = pm.GetNE();
   MPI_Allreduce(MPI_IN_PLACE, &ne_glob, 1, MPI_LONG_LONG, MPI_SUM, pm.GetComm());
   const real_t nrm = (real_t)ne_glob;
   real_t mn = 1e300, mx = -1e300, rmin = 1e300, rmax = -1e300;
   long long neg = 0;
   for (int e = 0; e < pm.GetNE(); e++)
   {
      ElementTransformation *T = pm.GetElementTransformation(e);
      real_t emin = 1e300, emax = -1e300;
      auto visit = [&](const IntegrationPoint &ip)
      {
         T->SetIntPoint(&ip);
         const real_t w = nrm * T->Jacobian().Det();
         emin = std::min(emin, w);
         emax = std::max(emax, w);
         if (w <= 0.0) { neg++; }
      };
      for (int k = 0; k < ir.GetNPoints(); k++) { visit(ir.IntPoint(k)); }
      for (auto &c : corners) { visit(c); }
      mn = std::min(mn, emin);
      mx = std::max(mx, emax);
      const real_t ratio = (emax > 0) ? emin / emax : -1.0;
      rmin = std::min(rmin, ratio);
      rmax = std::max(rmax, ratio);
   }
   MPI_Comm comm = pm.GetComm();
   MPI_Allreduce(&mn, &r.min_det, 1, MPITypeMap<real_t>::mpi_type, MPI_MIN, comm);
   MPI_Allreduce(&mx, &r.max_det, 1, MPITypeMap<real_t>::mpi_type, MPI_MAX, comm);
   MPI_Allreduce(&rmin, &r.min_ratio, 1, MPITypeMap<real_t>::mpi_type, MPI_MIN, comm);
   MPI_Allreduce(&rmax, &r.max_ratio, 1, MPITypeMap<real_t>::mpi_type, MPI_MAX, comm);
   MPI_Allreduce(&neg, &r.n_neg, 1, MPI_LONG_LONG, MPI_SUM, comm);
   return r;
}

MeshSizeResult MeshSize(ParMesh &pm)
{
   real_t vol = 0.0, hmin = 1e300, hmax = -1e300;
   for (int e = 0; e < pm.GetNE(); e++)
   {
      vol += pm.GetElementVolume(e);
      const real_t h = pm.GetElementSize(e, 0);
      hmin = std::min(hmin, h);
      hmax = std::max(hmax, h);
   }
   long long ne = pm.GetNE();
   MPI_Comm comm = pm.GetComm();
   MPI_Allreduce(MPI_IN_PLACE, &vol, 1, MPITypeMap<real_t>::mpi_type, MPI_SUM, comm);
   MPI_Allreduce(MPI_IN_PLACE, &ne, 1, MPI_LONG_LONG, MPI_SUM, comm);
   MPI_Allreduce(MPI_IN_PLACE, &hmin, 1, MPITypeMap<real_t>::mpi_type, MPI_MIN, comm);
   MPI_Allreduce(MPI_IN_PLACE, &hmax, 1, MPITypeMap<real_t>::mpi_type, MPI_MAX, comm);
   MeshSizeResult r;
   r.hmean = std::cbrt(vol / ne);
   r.hmin = hmin;
   r.hmax = hmax;
   r.ne = ne;
   return r;
}

// ===========================================================================
// 3. Spaces and operators
// ===========================================================================

Spaces::Spaces(ParMesh &pm, int p_) : p(p_), pmesh(&pm)
{
   MFEM_VERIFY(p >= 1, "Spaces: p >= 1");
   h1c.reset(new H1_FECollection(p, 3));
   ndc.reset(new ND_FECollection(p, 3));
   rtc.reset(new RT_FECollection(p - 1, 3));
   l2c.reset(new L2_FECollection(p - 1, 3));
   H1.reset(new ParFiniteElementSpace(&pm, h1c.get()));
   ND.reset(new ParFiniteElementSpace(&pm, ndc.get()));
   RT.reset(new ParFiniteElementSpace(&pm, rtc.get()));
   L2.reset(new ParFiniteElementSpace(&pm, l2c.get()));
}

Operators::Operators(const Spaces &s) : sp(s)
{
   const double t0 = WallTime();
   {
      ParDiscreteLinearOperator op(sp.H1.get(), sp.ND.get());
      op.AddDomainInterpolator(new GradientInterpolator);
      op.Assemble();
      op.Finalize();
      G.reset(op.ParallelAssemble());
   }
   {
      ParDiscreteLinearOperator op(sp.ND.get(), sp.RT.get());
      op.AddDomainInterpolator(new CurlInterpolator);
      op.Assemble();
      op.Finalize();
      C.reset(op.ParallelAssemble());
   }
   {
      ParDiscreteLinearOperator op(sp.RT.get(), sp.L2.get());
      op.AddDomainInterpolator(new DivergenceInterpolator);
      op.Assemble();
      op.Finalize();
      D.reset(op.ParallelAssemble());
   }
   t_assemble = WallTime() - t0;
}

double Operators::ApplyCurlTimed(const Vector &a_t, Vector &b_t, int nrep_min,
                                 double tmin) const
{
   MPI_Comm comm = sp.pmesh->GetComm();
   b_t.SetSize(C->Height());
   C->Mult(a_t, b_t); // warm up
   int nrep = 0;
   MPI_Barrier(comm);
   const double t0 = WallTime();
   double el = 0.0;
   do
   {
      for (int i = 0; i < nrep_min; i++) { C->Mult(a_t, b_t); }
      nrep += nrep_min;
      el = WallTime() - t0;
      // make the decision collective
      MPI_Allreduce(MPI_IN_PLACE, &el, 1, MPI_DOUBLE, MPI_MAX, comm);
   }
   while (el < tmin);
   return el / nrep;
}

real_t MaxAbsEntry(const HypreParMatrix &A)
{
   SparseMatrix diag, offd;
   HYPRE_BigInt *cmap = nullptr;
   A.GetDiag(diag);
   A.GetOffd(offd, cmap);
   real_t m = std::max(diag.NumNonZeroElems() ? diag.MaxNorm() : 0.0,
                       offd.NumNonZeroElems() ? offd.MaxNorm() : 0.0);
   MPI_Allreduce(MPI_IN_PLACE, &m, 1, MPITypeMap<real_t>::mpi_type, MPI_MAX,
                 A.GetComm());
   return m;
}

real_t Operators::MaxAbsDC() const
{
   std::unique_ptr<HypreParMatrix> P(ParMult(D.get(), C.get()));
   return MaxAbsEntry(*P);
}

real_t Operators::MaxAbsCG() const
{
   std::unique_ptr<HypreParMatrix> P(ParMult(C.get(), G.get()));
   return MaxAbsEntry(*P);
}

// ===========================================================================
// 4. Fields
// ===========================================================================

namespace
{
// A = amp (sin ky cos kz, sin kz cos kx, sin kx cos ky), exact curl:
//   curl_x = -k (sin kx sin ky + cos kz cos kx)
//   curl_y = -k (sin ky sin kz + cos kx cos ky)
//   curl_z = -k (sin kz sin kx + cos ky cos kz)
void ModA(real_t k, real_t amp, const Vector &x, Vector &A, bool add)
{
   const real_t sx = std::sin(k * x(0)), cx = std::cos(k * x(0));
   const real_t sy = std::sin(k * x(1)), cy = std::cos(k * x(1));
   const real_t sz = std::sin(k * x(2)), cz = std::cos(k * x(2));
   const real_t v[3] = {sy * cz, sz * cx, sx * cy};
   for (int i = 0; i < 3; i++) { A(i) = (add ? A(i) : 0.0) + amp * v[i]; }
}
void ModCurl(real_t k, real_t amp, const Vector &x, Vector &B, bool add)
{
   const real_t sx = std::sin(k * x(0)), cx = std::cos(k * x(0));
   const real_t sy = std::sin(k * x(1)), cy = std::cos(k * x(1));
   const real_t sz = std::sin(k * x(2)), cz = std::cos(k * x(2));
   const real_t v[3] = {-k * (sx * sy + cz * cx), -k * (sy * sz + cx * cy),
                        -k * (sz * sx + cy * cz)
                       };
   for (int i = 0; i < 3; i++) { B(i) = (add ? B(i) : 0.0) + amp * v[i]; }
}
}

Field MakeField(const std::string &name, const FieldParams &fp)
{
   Field f;
   f.name = name;
   const real_t k = kTwoPi;
   if (name == "abc")
   {
      // A = (a sin kz + c cos ky, b sin kx + a cos kz, c sin ky + b cos kx).
      // curl A = k A exactly for ANY a,b,c (Beltrami / ABC flow).
      const real_t a = fp.a, b = fp.b, c = fp.c;
      f.A = [=](const Vector &x, Vector &A)
      {
         A.SetSize(3);
         A(0) = a * std::sin(k * x(2)) + c * std::cos(k * x(1));
         A(1) = b * std::sin(k * x(0)) + a * std::cos(k * x(2));
         A(2) = c * std::sin(k * x(1)) + b * std::cos(k * x(0));
      };
      f.curlA = [=](const Vector &x, Vector &B)
      {
         B.SetSize(3);
         B(0) = k * (a * std::sin(k * x(2)) + c * std::cos(k * x(1)));
         B(1) = k * (b * std::sin(k * x(0)) + a * std::cos(k * x(2)));
         B(2) = k * (c * std::sin(k * x(1)) + b * std::cos(k * x(0)));
      };
   }
   else if (name == "mod")
   {
      f.A = [=](const Vector &x, Vector &A) { A.SetSize(3); ModA(k, 1.0, x, A, false); };
      f.curlA = [=](const Vector &x, Vector &B) { B.SetSize(3); ModCurl(k, 1.0, x, B, false); };
   }
   else if (name == "mod2")
   {
      const real_t amp = fp.mod2_amp;
      f.A = [=](const Vector &x, Vector &A)
      {
         A.SetSize(3);
         ModA(k, 1.0, x, A, false);
         ModA(2 * k, amp, x, A, true);
      };
      f.curlA = [=](const Vector &x, Vector &B)
      {
         B.SetSize(3);
         ModCurl(k, 1.0, x, B, false);
         ModCurl(2 * k, amp, x, B, true);
      };
   }
   else
   {
      MFEM_ABORT("unknown field '" << name << "' (abc, mod, mod2)");
   }
   return f;
}

// Exact helicity (b0 = 0):
//  abc : A.curlA = k|A|^2, <|A|^2> = a^2+b^2+c^2 (cross terms average to 0)
//        => H = k (a^2+b^2+c^2).
//  mod : every term in A.curlA contains an odd trigonometric factor in some
//        coordinate (e.g. sin kx sin ky sin ky cos kz) so H = 0 exactly;
//  mod2: the cross terms between wavenumbers k and 2k are orthogonal => H = 0.
real_t ExactHelicity(const std::string &name, const FieldParams &fp)
{
   if (name == "abc") { return kTwoPi * (fp.a * fp.a + fp.b * fp.b + fp.c * fp.c); }
   return 0.0;
}

// Exact energy 1/2 int |curl A|^2:
//  abc : 1/2 k^2 (a^2+b^2+c^2)
//  mod : <|curl A|^2> = k^2 * 3/2 => E = 3 k^2/4
//  mod2: orthogonal modes; second mode has k -> 2k and amplitude amp:
//        E = 3 k^2/4 (1 + 4 amp^2).
real_t ExactEnergy(const std::string &name, const FieldParams &fp)
{
   const real_t k = kTwoPi;
   if (name == "abc") { return 0.5 * k * k * (fp.a * fp.a + fp.b * fp.b + fp.c * fp.c); }
   if (name == "mod") { return 0.75 * k * k; }
   if (name == "mod2") { return 0.75 * k * k * (1.0 + 4.0 * fp.mod2_amp * fp.mod2_amp); }
   MFEM_ABORT("unknown field");
   return 0.0;
}

real_t Gauge::chi(const Vector &x) const
{
   return g * std::sin(k * x(0)) * std::sin(k * x(1)) * std::sin(k * x(2));
}
void Gauge::grad(const Vector &x, Vector &v) const
{
   v.SetSize(3);
   const real_t sx = std::sin(k * x(0)), cx = std::cos(k * x(0));
   const real_t sy = std::sin(k * x(1)), cy = std::cos(k * x(1));
   const real_t sz = std::sin(k * x(2)), cz = std::cos(k * x(2));
   v(0) = g * k * cx * sy * sz;
   v(1) = g * k * sx * cy * sz;
   v(2) = g * k * sx * sy * cz;
}

real_t SelfTestCurl(const FieldParams &fp, real_t tol, bool *pass)
{
   std::mt19937 gen(12345);
   std::uniform_real_distribution<real_t> U(0.0, 1.0);
   const real_t h = 1e-5;
   real_t worst = 0.0;
   for (const char *nm : {"abc", "mod", "mod2"})
   {
      Field F = MakeField(nm, fp);
      for (int t = 0; t < 20; t++)
      {
         Vector x(3), xp(3), xm(3), Ap(3), Am(3), Bex(3);
         for (int i = 0; i < 3; i++) { x(i) = U(gen); }
         // d[i][j] = d A_i / d x_j
         real_t d[3][3];
         for (int j = 0; j < 3; j++)
         {
            xp = x; xm = x; xp(j) += h; xm(j) -= h;
            F.A(xp, Ap); F.A(xm, Am);
            for (int i = 0; i < 3; i++) { d[i][j] = (Ap(i) - Am(i)) / (2 * h); }
         }
         const real_t fd[3] = {d[2][1] - d[1][2], d[0][2] - d[2][0],
                               d[1][0] - d[0][1]
                              };
         F.curlA(x, Bex);
         for (int i = 0; i < 3; i++) { worst = std::max(worst, std::fabs(fd[i] - Bex(i))); }
      }
   }
   if (pass) { *pass = (worst < tol); }
   return worst;
}

// ---------------------------------------------------------------------------
// Projection
// ---------------------------------------------------------------------------

namespace
{

// Coefficient evaluating the i-th reference basis function of a vector FE.
class BasisCoef : public VectorCoefficient
{
   const FiniteElement &fe;
   int idx;
   mutable DenseMatrix sh;
public:
   BasisCoef(const FiniteElement &fe_, int i)
      : VectorCoefficient(3), fe(fe_), idx(i), sh(fe_.GetDof(), 3) {}
   void Eval(Vector &V, ElementTransformation &, const IntegrationPoint &ip) override
   {
      fe.CalcVShape(ip, sh);
      V.SetSize(3);
      for (int d = 0; d < 3; d++) { V(d) = sh(idx, d); }
   }
};

// Integrated dofs of a vector field on ONE element, computed with an n-point
// Gauss-Legendre rule per direction (nq points), i.e. (almost) exact
// integrals for smooth fields.  The dofs are indexed like the standard
// (Gauss-Legendre open basis) element, with the same orientation signs as
// MFEM's own integrated-dof implementation (ND_HexahedronElement::
// ProjectIntegrated / RT_HexahedronElement::ProjectIntegrated), which we do
// NOT call because (i) they are protected, (ii) they use a quadrature rule of
// order p only (not exact for non-polynomial fields), and (iii) the MFEM
// interpolators (curl, grad) assume the standard dof representation.
//   ND: d_idx = sign * int_{sub-edge} v.(J e_c) dxi   (sub-edge between
//       consecutive Gauss-Lobatto points; closed points in the other 2 dirs)
//   RT: d_idx = sign * int_{sub-face} (adj(J) v)_c dxi1 dxi2  (flux through
//       the sub-face between Gauss-Lobatto points at the closed position)
// Dof ordering inside the loops follows the MFEM dof_map construction.
void IntegratedDofs(bool nd, const FiniteElement &sfe, VectorCoefficient &vc,
                    ElementTransformation &T, int nq, Vector &d)
{
   const auto *te = dynamic_cast<const TensorBasisElement *>(&sfe);
   MFEM_VERIFY(te, "IntegratedDofs: not a tensor-product element");
   const Array<int> &dm = te->GetDofMap();
   const int P = sfe.GetOrder();
   const real_t *cp = poly1d.ClosedPoints(P, BasisType::GaussLobatto);
   const IntegrationRule &ir1 = IntRules.Get(Geometry::SEGMENT, 2 * nq - 1);
   const int n1 = ir1.GetNPoints();
   d.SetSize(sfe.GetDof());
   Vector v(3);
   IntegrationPoint ip;
   int o = 0;
   for (int c = 0; c < 3; c++)
   {
      int im, jm, km;
      if (nd)
      {
         im = (c == 0) ? P : P + 1;
         jm = (c == 1) ? P : P + 1;
         km = (c == 2) ? P : P + 1;
      }
      else
      {
         im = (c == 0) ? P + 1 : P;
         jm = (c == 1) ? P + 1 : P;
         km = (c == 2) ? P + 1 : P;
      }
      for (int k = 0; k < km; k++)
         for (int j = 0; j < jm; j++)
            for (int i = 0; i < im; i++)
            {
               int idx = dm[o++];
               real_t sg = 1.0;
               if (idx < 0) { idx = -1 - idx; sg = -1.0; }
               const int ijk[3] = {i, j, k};
               real_t val = 0.0, meas = 1.0;
               if (nd)
               {
                  const real_t h = cp[ijk[c] + 1] - cp[ijk[c]];
                  meas = h;
                  for (int q = 0; q < n1; q++)
                  {
                     real_t pt[3] = {cp[i], cp[j], cp[k]};
                     pt[c] += h * ir1.IntPoint(q).x;
                     ip.Set3(pt[0], pt[1], pt[2]);
                     T.SetIntPoint(&ip);
                     vc.Eval(v, T, ip);
                     const DenseMatrix &J = T.Jacobian();
                     real_t dot = 0.0;
                     for (int a = 0; a < 3; a++) { dot += v(a) * J(a, c); }
                     val += ir1.IntPoint(q).weight * dot;
                  }
               }
               else
               {
                  // the two in-plane (open, sub-interval) directions
                  const int a1 = (c + 1) % 3, a2 = (c + 2) % 3;
                  const real_t h1 = cp[ijk[a1] + 1] - cp[ijk[a1]];
                  const real_t h2 = cp[ijk[a2] + 1] - cp[ijk[a2]];
                  meas = h1 * h2;
                  for (int q1 = 0; q1 < n1; q1++)
                     for (int q2 = 0; q2 < n1; q2++)
                     {
                        real_t pt[3] = {cp[i], cp[j], cp[k]};
                        pt[c] = cp[ijk[c]];
                        pt[a1] += h1 * ir1.IntPoint(q1).x;
                        pt[a2] += h2 * ir1.IntPoint(q2).x;
                        ip.Set3(pt[0], pt[1], pt[2]);
                        T.SetIntPoint(&ip);
                        vc.Eval(v, T, ip);
                        const DenseMatrix &AJ = T.AdjugateJacobian();
                        real_t dot = 0.0;
                        for (int b = 0; b < 3; b++) { dot += AJ(c, b) * v(b); }
                        val += ir1.IntPoint(q1).weight * ir1.IntPoint(q2).weight * dot;
                     }
               }
               d(idx) = sg * val * meas;
            }
   }
}

struct LocalConv
{
   DenseMatrix Minv;                    // (std dofs) = Minv * (integrated dofs)
};

LocalConv &GetConv(bool nd, const FiniteElement &sfe)
{
   static std::map<std::pair<int, int>, std::unique_ptr<LocalConv>> cache;
   const auto key = std::make_pair(nd ? 0 : 1, sfe.GetOrder());
   auto it = cache.find(key);
   if (it != cache.end()) { return *it->second; }

   std::unique_ptr<LocalConv> lc(new LocalConv);
   const int n = sfe.GetDof();

   // reference (identity) transformation of the unit cube
   H1_FECollection h1(1, 3);
   IsoparametricTransformation T;
   T.SetFE(h1.FiniteElementForGeometry(Geometry::CUBE));
   T.GetPointMat().SetSize(3, 8);
   const IntegrationRule *cv = Geometries.GetVertices(Geometry::CUBE);
   for (int v = 0; v < 8; v++)
   {
      T.GetPointMat()(0, v) = cv->IntPoint(v).x;
      T.GetPointMat()(1, v) = cv->IntPoint(v).y;
      T.GetPointMat()(2, v) = cv->IntPoint(v).z;
   }

   DenseMatrix M(n, n);
   Vector d(n);
   for (int i = 0; i < n; i++)
   {
      BasisCoef bc(sfe, i);
      IntegratedDofs(nd, sfe, bc, T, sfe.GetOrder() + 2, d);   // exact: polynomial
      for (int j = 0; j < n; j++) { M(j, i) = d(j); }
   }
   DenseMatrixInverse inv(M);
   inv.GetInverseMatrix(lc->Minv);
   auto &ref = *lc;
   cache[key] = std::move(lc);
   return ref;
}

void ProjectIntegratedImpl(ParGridFunction &g, VectorCoefficient &vc, bool nd,
                           int nq)
{
   ParFiniteElementSpace *fes = g.ParFESpace();
   const FiniteElement &sfe =
      *fes->FEColl()->FiniteElementForGeometry(Geometry::CUBE);
   LocalConv &lc = GetConv(nd, sfe);
   if (nq <= 0)
   {
      const FiniteElementSpace *nfes = fes->GetParMesh()->GetNodalFESpace();
      const int q = nfes ? nfes->GetMaxElementOrder() : 1;
      nq = sfe.GetOrder() + q + 5;
   }
   const int n = sfe.GetDof();
   Vector d(n), s(n);
   Array<int> vdofs;
   for (int e = 0; e < fes->GetNE(); e++)
   {
      fes->GetElementVDofs(e, vdofs);
      ElementTransformation *T = fes->GetElementTransformation(e);
      IntegratedDofs(nd, sfe, vc, *T, nq, d);
      lc.Minv.Mult(d, s);
      g.SetSubVector(vdofs, s);
   }
   // make the shared dofs consistent
   Vector t;
   g.GetTrueDofs(t);
   g.SetFromTrueDofs(t);
}

} // namespace

void ProjectA(ParGridFunction &a, VectorCoefficient &vc, ProjMode mode, int nq)
{
   if (mode == ProjMode::Pointwise)
   {
      a.ProjectCoefficient(vc);
      Vector t; a.GetTrueDofs(t); a.SetFromTrueDofs(t);
   }
   else { ProjectIntegratedImpl(a, vc, true, nq); }
}

void ProjectB(ParGridFunction &b, VectorCoefficient &vc, ProjMode mode, int nq)
{
   if (mode == ProjMode::Pointwise)
   {
      b.ProjectCoefficient(vc);
      Vector t; b.GetTrueDofs(t); b.SetFromTrueDofs(t);
   }
   else { ProjectIntegratedImpl(b, vc, false, nq); }
}

// ===========================================================================
// 5. Diagnostics
// ===========================================================================

namespace
{
void FillRules(const IntegrationRule *irs[], int order)
{
   for (int i = 0; i < Geometry::NumGeom; i++) { irs[i] = &IntRules.Get(i, order); }
}
}

real_t L2ErrorVec(const ParGridFunction &u, VectorCoefficient &ex, int order)
{
   const IntegrationRule *irs[Geometry::NumGeom];
   FillRules(irs, order);
   return const_cast<ParGridFunction &>(u).ComputeL2Error(ex, irs);
}

real_t L2NormVec(ParMesh &pm, VectorCoefficient &vc, int order)
{
   return std::sqrt(InnerCoef(pm, vc, vc, order));
}

real_t InnerCoef(ParMesh &pm, VectorCoefficient &u, VectorCoefficient &v,
                 int order)
{
   const IntegrationRule &ir = IntRules.Get(Geometry::CUBE, order);
   Vector a(3), b(3);
   real_t s = 0.0;
   for (int e = 0; e < pm.GetNE(); e++)
   {
      ElementTransformation *T = pm.GetElementTransformation(e);
      for (int k = 0; k < ir.GetNPoints(); k++)
      {
         const IntegrationPoint &ip = ir.IntPoint(k);
         T->SetIntPoint(&ip);
         u.Eval(a, *T, ip);
         v.Eval(b, *T, ip);
         s += ip.weight * T->Weight() * (a * b);
      }
   }
   MPI_Allreduce(MPI_IN_PLACE, &s, 1, MPITypeMap<real_t>::mpi_type, MPI_SUM,
                 pm.GetComm());
   return s;
}

real_t InnerGF(const ParGridFunction &u, const ParGridFunction &v, int order)
{
   ParMesh *pm = u.ParFESpace()->GetParMesh();
   const IntegrationRule &ir = IntRules.Get(Geometry::CUBE, order);
   DenseMatrix uv, vv;
   real_t s = 0.0;
   for (int e = 0; e < pm->GetNE(); e++)
   {
      ElementTransformation *T = pm->GetElementTransformation(e);
      u.GetVectorValues(*T, ir, uv);
      v.GetVectorValues(*T, ir, vv);
      for (int k = 0; k < ir.GetNPoints(); k++)
      {
         const IntegrationPoint &ip = ir.IntPoint(k);
         T->SetIntPoint(&ip);
         real_t dot = 0.0;
         for (int d = 0; d < 3; d++) { dot += uv(d, k) * vv(d, k); }
         s += ip.weight * T->Weight() * dot;
      }
   }
   MPI_Allreduce(MPI_IN_PLACE, &s, 1, MPITypeMap<real_t>::mpi_type, MPI_SUM,
                 pm->GetComm());
   return s;
}

DivDiag DivergenceDiag(const Operators &ops, const Vector &b_t)
{
   const Spaces &sp = ops.sp;
   Vector d_t(ops.D->Height());
   ops.D->Mult(b_t, d_t);
   DivDiag r;
   real_t mx = d_t.Normlinf();
   MPI_Allreduce(&mx, &r.maxabs, 1, MPITypeMap<real_t>::mpi_type, MPI_MAX,
                 sp.pmesh->GetComm());
   ParGridFunction dg(sp.L2.get());
   dg.SetFromTrueDofs(d_t);
   ConstantCoefficient zero(0.0);
   const FiniteElementSpace *nfes = sp.pmesh->GetNodalFESpace();
   const int q = nfes ? nfes->GetMaxElementOrder() : 1;
   const IntegrationRule *irs[Geometry::NumGeom];
   FillRules(irs, 2 * sp.p + 2 * q + 2);
   r.l2 = dg.ComputeL2Error(zero, irs);
   return r;
}

namespace
{
int DefaultMassOrder(const Spaces &sp)
{
   const FiniteElementSpace *nfes = sp.pmesh->GetNodalFESpace();
   const int q = nfes ? nfes->GetMaxElementOrder() : 1;
   return 2 * sp.p + 2 * q;
}
}

real_t EnergyMass(const Spaces &sp, const Vector &b_t, int qorder)
{
   if (qorder < 0) { qorder = DefaultMassOrder(sp); }
   ParBilinearForm m(sp.RT.get());
   auto *vi = new VectorFEMassIntegrator;
   vi->SetIntRule(&IntRules.Get(Geometry::CUBE, qorder));
   m.AddDomainIntegrator(vi);
   m.Assemble();
   m.Finalize();
   std::unique_ptr<HypreParMatrix> M(m.ParallelAssemble());
   Vector Mb(b_t.Size());
   M->Mult(b_t, Mb);
   return 0.5 * InnerProduct(sp.pmesh->GetComm(), b_t, Mb);
}

real_t HelicityMass(const Spaces &sp, const Vector &a_t, const Vector &b_t,
                    int qorder)
{
   if (qorder < 0) { qorder = DefaultMassOrder(sp); }
   // trial = ND (a), test = RT (b): matrix rows = RT, cols = ND.
   ParMixedBilinearForm m(sp.ND.get(), sp.RT.get());
   auto *vi = new VectorFEMassIntegrator;
   vi->SetIntRule(&IntRules.Get(Geometry::CUBE, qorder));
   m.AddDomainIntegrator(vi);
   m.Assemble();
   m.Finalize();
   std::unique_ptr<HypreParMatrix> M(m.ParallelAssemble());
   Vector Ma(M->Height());
   M->Mult(a_t, Ma);
   return InnerProduct(sp.pmesh->GetComm(), b_t, Ma);
}

int MassCGIterations(const Spaces &sp, bool pa, real_t rtol)
{
   ParBilinearForm m(sp.ND.get());
   m.AddDomainIntegrator(new VectorFEMassIntegrator);
   if (pa) { m.SetAssemblyLevel(AssemblyLevel::PARTIAL); }
   m.Assemble();
   Array<int> ess;
   OperatorHandle Ah;
   if (pa)
   {
      m.FormSystemMatrix(ess, Ah);
   }
   else
   {
      m.Finalize();
      Ah.Reset(m.ParallelAssemble());
   }
   std::unique_ptr<Solver> prec;
   std::unique_ptr<HypreSmoother> hs;
   if (pa) { prec.reset(new OperatorJacobiSmoother(m, ess)); }
   else
   {
      hs.reset(new HypreSmoother(*Ah.As<HypreParMatrix>(), HypreSmoother::Jacobi));
   }
   const int n = sp.ND->GetTrueVSize();
   Vector rhs(n), x(n);
   std::mt19937 gen(777 + sp.pmesh->GetMyRank());
   std::uniform_real_distribution<real_t> U(-1.0, 1.0);
   for (int i = 0; i < n; i++) { rhs(i) = U(gen); }
   x = 0.0;
   CGSolver cg(sp.pmesh->GetComm());
   cg.SetOperator(*Ah);
   if (pa) { cg.SetPreconditioner(*prec); } else { cg.SetPreconditioner(*hs); }
   cg.SetRelTol(rtol);
   cg.SetAbsTol(0.0);
   cg.SetMaxIter(5000);
   cg.SetPrintLevel(-1);
   cg.Mult(rhs, x);
   return cg.GetNumIterations();
}

double WallTime() { return MPI_Wtime(); }

double RssMB(MPI_Comm comm, double *peak_mb)
{
   double rss = 0.0, hwm = 0.0;
   std::ifstream f("/proc/self/status");
   std::string line;
   while (std::getline(f, line))
   {
      long kb;
      if (line.compare(0, 6, "VmRSS:") == 0 && sscanf(line.c_str() + 6, "%ld", &kb) == 1)
      {
         rss = kb / 1024.0;
      }
      else if (line.compare(0, 6, "VmHWM:") == 0 && sscanf(line.c_str() + 6, "%ld", &kb) == 1)
      {
         hwm = kb / 1024.0;
      }
   }
   MPI_Allreduce(MPI_IN_PLACE, &rss, 1, MPI_DOUBLE, MPI_MAX, comm);
   MPI_Allreduce(MPI_IN_PLACE, &hwm, 1, MPI_DOUBLE, MPI_MAX, comm);
   if (peak_mb) { *peak_mb = hwm; }
   return rss;
}

// ===========================================================================
// 6. CSV
// ===========================================================================

void CsvRow::SetStr(const std::string &key, const std::string &s)
{
   for (auto &kv : items)
   {
      if (kv.first == key) { kv.second = s; return; }
   }
   items.emplace_back(key, s);
}
void CsvRow::Set(const std::string &key, double v)
{
   char buf[64];
   if (std::isnan(v)) { SetStr(key, "nan"); return; }
   snprintf(buf, sizeof(buf), "%.16e", v);
   SetStr(key, buf);
}
void CsvRow::Set(const std::string &key, long long v)
{
   SetStr(key, std::to_string(v));
}
void CsvRow::Set(const std::string &key, const std::string &s) { SetStr(key, s); }

void CsvAppend(MPI_Comm comm, const std::string &path, const CsvRow &row)
{
   int rank;
   MPI_Comm_rank(comm, &rank);
   if (rank != 0) { return; }
   std::string header, line;
   for (size_t i = 0; i < row.Items().size(); i++)
   {
      if (i) { header += ","; line += ","; }
      header += row.Items()[i].first;
      line += row.Items()[i].second;
   }
   std::string existing;
   {
      std::ifstream f(path);
      if (f) { std::getline(f, existing); }
   }
   if (!existing.empty() && existing != header)
   {
      MFEM_ABORT("CsvAppend: header mismatch in " << path << "\n  file: "
                 << existing << "\n  new : " << header);
   }
   std::ofstream out(path, std::ios::app);
   MFEM_VERIFY(out, "cannot open " << path);
   if (existing.empty()) { out << header << "\n"; }
   out << line << "\n";
}

} // namespace vp
