// vp_transfer.hpp -- mesh-to-mesh transfer operators (A-route: ND, B-route: RT)
// and diagnostics on the target mesh.  See docs/design_note.md section 5.
//
//   RemoteEvaluator   evaluates a ParGridFunction (ND/RT/H1) of a source mesh
//                     M1 at arbitrary physical points (any rank), via GSLIB
//                     FindPoints with periodic wrapping into the unit cell.
//   Transfer          one operator = one (kind, source mesh, target space)
//                     plan; Apply() can be called repeatedly (point location
//                     is done once).
//   TransferA_*/B_*   one-shot convenience wrappers.
//   MeshDiag          diagnostics of a state (a, b) on a mesh.
#pragma once

#include "vp_core.hpp"

#include <cmath>
#include <memory>
#include <string>
#include <vector>

namespace vp
{

// ---------------------------------------------------------------------------
// RemoteEvaluator
// ---------------------------------------------------------------------------

/// Evaluates grid functions of the source mesh at physical points.
///
/// Periodic handling: query points x are first wrapped into [0,1)^3
/// (x - floor(x)).  Because deformed periodic meshes are stored as a
/// deformed *copy of the unit cell* (the cells next to the periodic boundary
/// stick out of [0,1)^3 by up to eps*max|f|), a wrapped point can still be
/// outside every element of the unit-cell representation; such points are
/// searched again at the periodic images x + s, s in {-1,0,1}^3 (only for the
/// few points that were not found).  Points are only accepted if inside an
/// element (code 0) or on an element boundary within 1e-9 (bdr tolerance of
/// GSLIB is tightened from its default 1e-4 to this, otherwise points that
/// are slightly outside the unit cell image would be clamped to the border
/// element).
class RemoteEvaluator
{
public:
   explicit RemoteEvaluator(ParMesh &M1, double bbox_rel_size_inc = 0.1,
                            double newt_tol = 1e-12);
   ~RemoteEvaluator();
   RemoteEvaluator(const RemoteEvaluator &) = delete;
   RemoteEvaluator &operator=(const RemoteEvaluator &) = delete;

   /// Locate the points x (ordering byVDIM: x0 y0 z0 x1 y1 z1 ...; these are
   /// the points of THIS rank, any location).  Returns an id of the located
   /// set; the evaluator remembers only the last set.  Not-found points are
   /// counted (global) in NumNotFound(); they evaluate to 0.
   long long Locate(const Vector &x);
   long long CurrentSet() const { return set_id; }
   long long NumPoints() const { return npts_local; }
   long long NumNotFound() const { return nnotfound; }   ///< global, last Locate
   long long NumShifted() const { return nshifted; }     ///< global: needed an image shift
   long long NumOnBoundary() const { return nborder; }   ///< global: code 1
   double TimeLocate() const { return t_locate; }
   double TimeEval() const { return t_eval; }

   /// Values of gf at the located points, ordering byVDIM
   /// (npts x vdim, vdim = gf.VectorDim()).
   void Evaluate(const ParGridFunction &gf, Vector &vals);

   ParMesh &Mesh() { return *mesh; }

private:
   ParMesh *mesh;
   std::unique_ptr<FindPointsGSLIB> finder;
   long long set_id = 0, npts_local = 0, nnotfound = 0, nshifted = 0, nborder = 0;
   double t_locate = 0, t_eval = 0;
};

// ---------------------------------------------------------------------------
// Transfer operators
// ---------------------------------------------------------------------------

enum class TransferKind { A_pt, A_int, A_l2, B_pt, B_int, B_l2, B_l2c };

const char *KindName(TransferKind k);
TransferKind ParseKind(const std::string &s);          ///< "A_pt", ... aborts if unknown
bool IsARoute(TransferKind k);
std::vector<TransferKind> ParseKindList(const std::string &csv);  ///< comma list or "all"

struct TransferOptions
{
   /// Gauss points per direction for the *_int functionals.  <=0: default
   /// DefaultTransferNq(p) = p + 2.
   int nq = -1;
   /// Quadrature order of the rhs of the L2 projections (<=0: 2p+2q+2).
   int l2_order = -1;
   real_t l2_rtol = 1e-12;          ///< CG rel tol of the L2-projection mass solve
   real_t clean_outer_rtol = 1e-12; ///< outer CG on the Schur complement
   real_t clean_inner_rtol = 1e-14; ///< inner CG on M_RT
   bool clean_amg = true;           ///< precondition the outer CG with AMG on D diag(M)^{-1} D^T
   int max_iter = 5000;
};

int DefaultTransferNq(int p);

struct TransferStats
{
   double t_setup = 0;     ///< plan construction: point collection + location + mass setup
   double t_locate = 0;    ///<   of which FindPoints (+ image retries)
   double t_collect = 0;   ///<   of which collecting the physical points
   double t_eval = 0;      ///< last Apply: evaluation of the source (interpolation)
   double t_dofs = 0;      ///< last Apply: dof functional / rhs assembly
   double t_solve = 0;     ///< last Apply: CG solves (L2, clean)
   double t_apply = 0;     ///< last Apply total
   int iters = 0;          ///< CG iterations (L2: mass solve; clean: outer, after the L2 solve)
   int iters_l2 = 0;       ///< (clean) iterations of the L2 mass solve
   int inner_iters = 0;    ///< (clean) total inner M_RT CG iterations
   int clean_rounds = 0;   ///< (clean) defect-correction rounds
   long long npoints = 0;  ///< number of evaluation points (this rank summed globally)
   long long nnotfound = 0;
   long long nshifted = 0;
   double div_before_clean = std::nan("");   ///< (clean) max|D_h b| before cleaning
   double TotalOneShot() const { return t_setup + t_apply; }
};

/// A transfer plan from a source mesh (via its RemoteEvaluator) to the space
/// of `dst` (ND for the A ops, RT for the B ops).  The points at which the
/// source has to be evaluated are collected at construction (same loops as
/// vp::ProjectA/ProjectB), located once, and Apply() can be repeated for
/// different source functions.  `dst_ops` is only needed for B_l2c (D_h).
class Transfer
{
public:
   Transfer(TransferKind kind, RemoteEvaluator &src_eval,
            const ParFiniteElementSpace &dst_space,
            const Operators *dst_ops = nullptr,
            const TransferOptions &opt = TransferOptions());
   ~Transfer();
   /// dst <- T(src).  src lives on the source mesh (ND for A ops, RT for B ops).
   void Apply(const ParGridFunction &src, ParGridFunction &dst);
   const TransferStats &Stats() const { return stats; }
   TransferKind Kind() const { return kind; }
   struct Impl;
private:
   TransferKind kind;
   TransferStats stats;
   std::unique_ptr<Impl> impl;
};

/// One-shot wrappers (build a plan, apply once; the returned stats include
/// the setup time).  dst lives on M2 (ND resp. RT).
TransferStats TransferA_Point(RemoteEvaluator &ev, const ParGridFunction &a1,
                              ParGridFunction &a2, const TransferOptions &o = {});
TransferStats TransferA_Int(RemoteEvaluator &ev, const ParGridFunction &a1,
                            ParGridFunction &a2, const TransferOptions &o = {});
TransferStats TransferA_L2(RemoteEvaluator &ev, const ParGridFunction &a1,
                           ParGridFunction &a2, const TransferOptions &o = {});
TransferStats TransferB_Point(RemoteEvaluator &ev, const ParGridFunction &b1,
                              ParGridFunction &b2, const TransferOptions &o = {});
TransferStats TransferB_Int(RemoteEvaluator &ev, const ParGridFunction &b1,
                            ParGridFunction &b2, const TransferOptions &o = {});
TransferStats TransferB_L2(RemoteEvaluator &ev, const ParGridFunction &b1,
                           ParGridFunction &b2, const TransferOptions &o = {});
/// B_L2 followed by projection onto ker(D_h) (needs D_h of M2).
TransferStats TransferB_L2Clean(RemoteEvaluator &ev, const ParGridFunction &b1,
                                ParGridFunction &b2, const Operators &ops2,
                                const TransferOptions &o = {});

/// b = Pi_RT^int(B0) + C_h a   (true dofs; b0_t = projected mean field).
void ComposeB(const Operators &ops, const Vector &b0_t, const Vector &a_t,
              Vector &b_t);
/// Same, projecting the (constant) mean field B0 with integrated dofs.
void ComposeB(const Spaces &sp, const Operators &ops, VectorCoefficient &B0,
              const Vector &a_t, Vector &b_t);

// ---------------------------------------------------------------------------
// Diagnostics on a mesh
// ---------------------------------------------------------------------------

struct StateDiag
{
   real_t errB_L2 = NAN, errB_rel = NAN;      ///< vs the analytic field B0 + curl A
   real_t div_l2 = NAN, div_max = NAN;        ///< ||D_h b||, max|D_h b|  (raw)
   real_t div_rel_l2 = NAN, div_rel_max = NAN;///< scaled by ||B||/h
   real_t energy = NAN;                       ///< 1/2 int |B_h|^2
   real_t energy_rel_ex = NAN;                ///< (E_h - E_exact)/E_exact
   real_t helicity = NAN;                     ///< int A_h . (C_h a) (A-route state only)
   real_t hel_scale = NAN;                    ///< ||a||_L2 ||C_h a||_L2
   real_t flux_err = NAN;                     ///< max_family |flux(b) - B0_dir|
   real_t a_L2 = NAN;                         ///< ||a||_L2
   real_t coulomb_res = NAN;                  ///< ||G_h^T M_ND a||_2 (weak divergence of A)
   real_t grad_frac = NAN;                    ///< ||P_grad a||_L2 / ||a||_L2 (gradient part of A)
};

/// Diagnostics context of one mesh: caches the exact norms, the PA ND mass
/// operator, the H1 Laplacian (for the gradient fraction) ...
class MeshDiag
{
public:
   MeshDiag(const Spaces &sp, const Operators &ops,
            const std::vector<SliceFamily> &slices, const std::string &field,
            const FieldParams &fp, const Vector &B0eff, const Vector &b0_t);
   ~MeshDiag();
   /// a_t may be nullptr (B-route state).  Helicity / ||a|| / Coulomb /
   /// gradient fraction are then NaN.
   StateDiag Evaluate(const Vector *a_t, const Vector &b_t, bool want_grad_frac = true);
   /// ||u - v||_L2 for two RT true dof vectors.
   real_t RTDistance(const Vector &b1_t, const Vector &b2_t);
   real_t ExactNormB() const { return normB_ex; }
   real_t ExactEnergyTotal() const { return E_ex; }
   real_t Hmean() const { return hmean; }
private:
   struct Impl;
   std::unique_ptr<Impl> impl;
   real_t normB_ex = 0, E_ex = 0, hmean = 0;
};

} // namespace vp
