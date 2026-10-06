// vp_core.hpp -- core library for the vector-potential (A-based) magnetic
// representation study.  See PLAN.md and docs/design_note.md.
//
// Notation: p = polynomial order of ND (A); B lives in RT_{p-1}, div B in
// L2_{p-1}, gauge functions in H1_p.  q = geometric (mesh node) order.
//
// All operators act on TRUE dofs (HypreParMatrix from ParallelAssemble).
// Everything is MPI-parallel and works with a single rank.
#pragma once

#include "mfem.hpp"

#include <functional>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace vp
{

using namespace mfem;

constexpr real_t kTwoPi = 6.283185307179586476925286766559;

// ---------------------------------------------------------------------------
// 1. Mesh construction and topological slice families
// ---------------------------------------------------------------------------

/// Fully periodic unit box with N^3 hexahedra, geometric order geom_order,
/// nodes in a discontinuous (L2) space.  Aborts if the mesh is not fully
/// periodic (#faces != 3 N^3) or has the wrong number of elements.
ParMesh *MakePeriodicBox(MPI_Comm comm, int N, int geom_order);

/// A family of mesh faces forming the topological image of the plane
/// x_dir = j/N.  face indices are ParMesh local face numbers; sign[i] is the
/// sign of the i-th component (i = dir) of the face normal (as used by the RT
/// face dofs) on the UNDEFORMED mesh.  Topology never changes under
/// deformation so the lists remain valid.
struct SliceFamily
{
   int dir = 0;
   int j = 0;
   Array<int> faces;
   Array<int> sign;
   long long nfaces_global = 0;   ///< (diagnostic) should be N^2
};

/// Record the slice families {x_dir = j/N} for dir=0,1,2 and each j in js.
/// Must be called on the undeformed mesh.
std::vector<SliceFamily> RecordSlices(ParMesh &pmesh, int N,
                                      const std::vector<int> &js);

/// Net flux of the RT field B through the face family (global, MPI reduced).
/// Computed exactly from the RT face dofs (Gauss-Legendre weights): for
/// RT_{p-1} the p^2 face dofs are point values of the reference flux density
/// at a tensor Gauss-Legendre grid, so the sum with GL weights integrates
/// B.n exactly.
real_t SliceFlux(const ParGridFunction &B_rt, const SliceFamily &fam);

/// Serial/one-rank cross-check of SliceFlux using face quadrature of B.n via
/// GetFaceElementTransformations (only valid when there are no shared faces).
real_t SliceFluxQuadrature(const ParGridFunction &B_rt, const SliceFamily &fam,
                           int qorder);

// ---------------------------------------------------------------------------
// 2. Mesh deformation and quality
// ---------------------------------------------------------------------------

/// nodes = x0 + eps f_variant(x0).  The undeformed nodes are cached (keyed by
/// the mesh pointer) at the first call so repeated calls never accumulate.
/// variant: 1: f1 = (sin2piy sin2piz, sin2piz sin2pix, sin2pix sin2piy)/(2pi)
///          2: f2 = same with 4pi inside, divided by 4pi
///          3: f1 + 0.5 f2
void DeformMesh(ParMesh &pmesh, real_t eps, int variant);

/// Evaluate the displacement field f_variant at x (for tests/analysis).
void DisplacementField(int variant, const real_t *x, real_t *f);

struct JacobianStatsResult
{
   real_t min_det = 0, max_det = 0;       ///< global over quad points+corners
   long long n_neg = 0;                   ///< #(points with detJ <= 0)
   real_t min_ratio = 0, max_ratio = 0;   ///< over elements of min/max detJ
   int quad_order = 0;
};
/// Jacobian statistics with a quadrature rule of order 2 q + 2 (plus the 8
/// element corners), globally reduced.  detJ is normalised by the undeformed
/// cell volume 1/NE (so it equals 1 on the undeformed box).
JacobianStatsResult JacobianStats(ParMesh &pmesh);

struct MeshSizeResult { real_t hmean, hmin, hmax; long long ne; };
/// hmean = (Vol/NE)^(1/3); hmin/hmax over elements of GetElementSize(.,0)
/// (det(J)^(1/3) relative to the reference element, so hmin<=hmean<=hmax).
MeshSizeResult MeshSize(ParMesh &pmesh);

// ---------------------------------------------------------------------------
// 3. Spaces and discrete operators
// ---------------------------------------------------------------------------

struct Spaces
{
   int p;
   ParMesh *pmesh;
   std::unique_ptr<H1_FECollection> h1c;
   std::unique_ptr<ND_FECollection> ndc;
   std::unique_ptr<RT_FECollection> rtc;
   std::unique_ptr<L2_FECollection> l2c;
   std::unique_ptr<ParFiniteElementSpace> H1, ND, RT, L2;

   Spaces(ParMesh &pm, int p_);
   long long NdofsH1() const { return H1->GlobalTrueVSize(); }
   long long NdofsND() const { return ND->GlobalTrueVSize(); }
   long long NdofsRT() const { return RT->GlobalTrueVSize(); }
   long long NdofsL2() const { return L2->GlobalTrueVSize(); }
};

struct Operators
{
   const Spaces &sp;
   std::unique_ptr<HypreParMatrix> G, C, D;   // H1->ND, ND->RT, RT->L2
   double t_assemble = 0.0;                   // wall time to build G, C, D

   explicit Operators(const Spaces &s);

   /// b_t = C_h a_t (true dofs)
   void ApplyCurl(const Vector &a_t, Vector &b_t) const { C->Mult(a_t, b_t); }
   /// Applies C_h repeatedly (at least nrep_min times and >= tmin seconds),
   /// returns the mean wall time per application (max over ranks).
   double ApplyCurlTimed(const Vector &a_t, Vector &b_t, int nrep_min = 5,
                         double tmin = 0.2) const;

   /// max |entry| of D_h C_h (sparse product) and C_h G_h; global.
   real_t MaxAbsDC() const;
   real_t MaxAbsCG() const;
};

/// max |entry| of a HypreParMatrix (global over ranks).
real_t MaxAbsEntry(const HypreParMatrix &A);

// ---------------------------------------------------------------------------
// 4. Analytic fields (all 1-periodic, k = 2 pi)
// ---------------------------------------------------------------------------

struct FieldParams
{
   // ABC coefficients: any nonzero values give a Beltrami field curl A = k A.
   real_t a = 1.0, b = 0.8, c = 0.6;
   real_t mod2_amp = 0.3;     // amplitude of the 4pi part of "mod2"
};

/// A named analytic vector potential with its exact curl.
struct Field
{
   std::string name;
   std::function<void(const Vector &, Vector &)> A, curlA;
};

/// "abc", "mod", "mod2".  Aborts on unknown name.
Field MakeField(const std::string &name, const FieldParams &fp = FieldParams());

/// Exact helicity int A.B over the unit torus for the named field with
/// b0 = 0 (analytic closed form, see comments in vp_core.cpp).
real_t ExactHelicity(const std::string &name, const FieldParams &fp);
/// Exact magnetic energy 1/2 int |curl A|^2 (closed form).
real_t ExactEnergy(const std::string &name, const FieldParams &fp);

/// Gauge function chi = g sin(k x) sin(k y) sin(k z) with k = kfac*pi (kfac
/// = 2: smooth, 6: rough).  Provides chi and grad chi.
struct Gauge
{
   real_t g, k;
   Gauge(real_t g_, real_t kfac_pi) : g(g_), k(kfac_pi * M_PI) {}
   real_t chi(const Vector &x) const;
   void grad(const Vector &x, Vector &v) const;
};

/// Finite-difference self test of the analytic curls (20 random points,
/// central differences, tolerance tol).  Returns the max deviation found.
real_t SelfTestCurl(const FieldParams &fp, real_t tol, bool *pass);

enum class ProjMode { Pointwise, Integrated };

/// Project A (ND) / B (RT).
///  Pointwise: standard ParGridFunction::ProjectCoefficient (interpolatory
///    dofs: point values of the tangential/normal component).
///  Integrated: dofs are the line integrals of A.t over the sub-edges between
///    consecutive Gauss-Lobatto points (ND) / surface integrals of B.n over
///    the sub-faces (RT) -- the same functionals as MFEM's IntegratedGLL basis
///    (ND_HexahedronElement::ProjectIntegrated, RT_HexahedronElement::
///    ProjectIntegrated), but evaluated with an nq-point Gauss rule per
///    direction (nq<=0: p+q+5) instead of MFEM's order-p rule, so that they are
///    exact to roundoff for smooth fields.  MFEM's own routines are protected,
///    use a low-order rule, and its integrated bases are not supported by the
///    curl/gradient interpolators (they assume nodal dofs), so we compute
///    the integrated dofs d with our own routine and convert them to the
///    standard basis of the SAME function space with the local reference
///    matrix s = M^{-1} d, M_ji = L_j^{int}(phi_i^{std}).  The result commutes
///    with the discrete operators: C_h Pi_ND^int A = Pi_RT^int curl A and
///    G_h Pi_H1 chi = Pi_ND^int grad chi to roundoff, on any (curved) mesh.
void ProjectA(ParGridFunction &a, VectorCoefficient &vc, ProjMode mode,
              int nq = -1);
void ProjectB(ParGridFunction &b, VectorCoefficient &vc, ProjMode mode,
              int nq = -1);

// ---------------------------------------------------------------------------
// 5. Diagnostics (all MPI-reduced)
// ---------------------------------------------------------------------------

/// L2 error of a vector-valued grid function vs a coefficient (rule order).
real_t L2ErrorVec(const ParGridFunction &u, VectorCoefficient &ex, int order);
/// L2 norm of a vector coefficient over the mesh.
real_t L2NormVec(ParMesh &pm, VectorCoefficient &vc, int order);
/// int u.v over the mesh for two coefficients.
real_t InnerCoef(ParMesh &pm, VectorCoefficient &u, VectorCoefficient &v,
                 int order);
/// int u_h . v_h for two grid functions (any vector FE spaces).
real_t InnerGF(const ParGridFunction &u, const ParGridFunction &v, int order);

struct DivDiag
{
   real_t l2 = 0;       ///< ||D_h b||_{L2}
   real_t maxabs = 0;   ///< max |(D_h b)_i| over true dofs
};
DivDiag DivergenceDiag(const Operators &ops, const Vector &b_t);

/// Energy 1/2 b^T M_RT b with M_RT fully assembled (diagnostic; memory heavy
/// for large p).
real_t EnergyMass(const Spaces &sp, const Vector &b_t, int qorder = -1);
/// Helicity a^T M_mix b with M_mix: ND (trial) x RT (test) mixed mass matrix
/// (assembled), i.e. int A_h . B_h.
real_t HelicityMass(const Spaces &sp, const Vector &a_t, const Vector &b_t,
                    int qorder = -1);

/// Jacobi-preconditioned CG on the ND mass matrix (partial assembly if
/// pa==true) for a random rhs; returns iteration count (rel tol 1e-12).
int MassCGIterations(const Spaces &sp, bool pa = true, real_t rtol = 1e-12);

double WallTime();                 ///< MPI_Wtime
/// Resident set size (VmRSS) in MB from /proc/self/status, max over ranks;
/// peak (VmHWM) optionally returned.
double RssMB(MPI_Comm comm, double *peak_mb = nullptr);

// ---------------------------------------------------------------------------
// 6. CSV writer
// ---------------------------------------------------------------------------

class CsvRow
{
public:
   /// Add or replace a column (keeps first-insertion order).
   void Set(const std::string &key, double v);
   void Set(const std::string &key, long long v);
   void Set(const std::string &key, int v) { Set(key, (long long)v); }
   void Set(const std::string &key, const std::string &s);
   void Set(const std::string &key, const char *s) { Set(key, std::string(s)); }
   const std::vector<std::pair<std::string, std::string>> &Items() const
   { return items; }
private:
   std::vector<std::pair<std::string, std::string>> items;
   void SetStr(const std::string &key, const std::string &s);
};

/// Append the row to the CSV file (rank 0 of comm only); writes the header if
/// the file is new/empty; aborts if an existing header differs.
void CsvAppend(MPI_Comm comm, const std::string &path, const CsvRow &row);

} // namespace vp
