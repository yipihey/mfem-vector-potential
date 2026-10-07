// vp_gauge.hpp -- Coulomb (discrete minimum-norm) gauge fixing of the ND vector potential.
//
//   a  <-  a - G_h chi,   L chi = G_h^T M_ND a,   L = G_h^T M_ND G_h.
//
// L is the H1_p stiffness matrix (it equals the ParBilinearForm DiffusionIntegrator matrix
// to roundoff, see StiffnessMismatch()); it is singular (constants), the solve is done with
// AMG-preconditioned CG, the right-hand side and the solution are orthogonalised against the
// constant vector.  After Apply():  G_h^T M_ND a = 0 (to the solve tolerance), i.e. a is the
// M_ND-orthogonal-to-range(G_h) representative of its gauge class, and C_h a (hence B) is
// unchanged (C_h G_h = 0).
//
// All vectors are TRUE-dof vectors.  Everything works with any number of MPI ranks.
#pragma once

#include "vp_core.hpp"

#include <memory>

namespace vp
{

struct GaugeFixStats
{
   int iters = 0;                   ///< CG iterations of the gauge solve (jacobi variant: #sweeps)
   double time = 0.0;               ///< wall time of the gauge fix itself (rhs, solve, update); excl. diagnostics below
   real_t resid_before = 0, resid_after = 0;   ///< ||G_h^T M_ND a||_2 (Coulomb residual, same scale as StateDiag::coulomb_res)
   real_t grad_fraction_before = 0, grad_fraction_after = 0;   ///< ||G chi||_M/||a||_M, chi = L^{-1} G^T M a
};

class CoulombGauge
{
public:
   CoulombGauge(const Spaces &sp, const Operators &ops);
   ~CoulombGauge();

   /// Exact Coulomb projection of a_t (in place).  rtol: rel. tolerance of the (preconditioned) CG.
   /// If diagnostics==false the *_after entries and the grad_fraction_before are left at 0 (saves
   /// a second solve); resid_before is always set.
   GaugeFixStats Apply(Vector &a_t, real_t rtol = 1e-12, bool diagnostics = true) const;

   /// ||G chi||_M/||a||_M with chi = L^{-1} G^T M a (rel. tol 1e-10).
   real_t GradFraction(const Vector &a_t) const;

   /// EXPERIMENTAL local variant: a <- a - G chi~, chi~ = nsweeps steps of Chebyshev-accelerated Jacobi
   /// (target interval [lmax/30, 1.05 lmax] of D^{-1} L, zero initial guess) for L chi = G^T M a; only
   /// mat-vecs, no global solve.  It damps only the rough (high-frequency) part of the gradient content.
   /// Diagnostics (before/after fractions) are computed with exact solves, outside the timing.
   GaugeFixStats JacobiGaugeSmooth(Vector &a_t, int nsweeps, bool diagnostics = true) const;

   /// ||G^T M_ND a||_2.
   real_t CoulombResidual(const Vector &a_t) const;

   /// max|L - K| / max|K| with K the H1 DiffusionIntegrator stiffness matrix (same quadrature rule).
   /// Builds K on demand (diagnostic).  Should be ~1e-15.
   real_t StiffnessMismatch() const;

   /// ||M_ND a||-type helper: a^T M a.
   real_t MassNormSq(const Vector &a_t) const;

   const HypreParMatrix &Lmat() const { return *L; }
   int LastIterations() const { return last_iters; }
   double SetupTime() const { return t_setup; }

private:
   const Spaces &sp;
   const Operators &ops;
   int ir_order = 0;
   std::unique_ptr<ParBilinearForm> Mnd;
   OperatorHandle Mh;
   std::unique_ptr<HypreParMatrix> L;
   std::unique_ptr<HypreBoomerAMG> amg;
   struct MeanFreePrec;                        ///< AMG wrapped by projections onto mean-free vectors
   std::unique_ptr<Solver> prec;
   mutable std::unique_ptr<CGSolver> cg;
   mutable Vector dinv;           ///< 1/diag(L)
   mutable real_t lmax = 0;       ///< estimate of lambda_max(D^{-1} L)
   double t_setup = 0.0;
   mutable int last_iters = 0;

   void RemoveMean(Vector &v) const;
   /// chi = L^{-1} g (g is modified: mean removed).  Returns CG iterations.
   int SolveL(Vector &g, Vector &chi, real_t rtol) const;
   void Rhs(const Vector &a_t, Vector &Ma, Vector &g) const;
};

} // namespace vp
