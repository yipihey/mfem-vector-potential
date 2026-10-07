// exp5_common.hpp -- helpers of E5 (smooth ALE cycle with rezones):
//   * closed-form gradient of the displacement field f (variants 1,2,3) and the
//     material map Phi_c(X) = X + c f(X), c = eps s(t);
//   * Newton inverse of Phi_c (1e-14, analytic Jacobian, backtracking);
//   * the Cauchy push-forward of the initial field (B = F B0(X)/det F,
//     A = F^{-T} A0(X)) evaluated at PHYSICAL points x (the label X = Phi_c^{-1}(x)
//     is obtained by Newton, independent of how often the mesh was rezoned);
//   * moving a mesh from its material labels.
#pragma once
#include "exp_transfer_common.hpp"

namespace vp
{
namespace e5
{

/// d f_i / d x_j of DisplacementField(variant): G[3 i + j].
inline void GradF(int variant, const real_t *x, real_t *G)
{
   struct Mode { real_t k, w; };
   std::vector<Mode> modes;
   switch (variant)
   {
      case 1: modes = {{kTwoPi, 1.0}}; break;
      case 2: modes = {{2 * kTwoPi, 1.0}}; break;
      case 3: modes = {{kTwoPi, 1.0}, {2 * kTwoPi, 0.5}}; break;
      default: MFEM_ABORT("GradF: unknown variant " << variant);
   }
   for (int i = 0; i < 9; i++) { G[i] = 0.0; }
   for (const Mode &m : modes)
   {
      const real_t sx = std::sin(m.k * x[0]), sy = std::sin(m.k * x[1]), sz = std::sin(m.k * x[2]);
      const real_t cx = std::cos(m.k * x[0]), cy = std::cos(m.k * x[1]), cz = std::cos(m.k * x[2]);
      // f0 = sy sz / k, f1 = sz sx / k, f2 = sx sy / k
      G[1] += m.w * cy * sz;  G[2] += m.w * sy * cz;
      G[3] += m.w * sz * cx;  G[5] += m.w * cz * sx;
      G[6] += m.w * cx * sy;  G[7] += m.w * sx * cy;
   }
}

inline real_t SProfile(real_t t, real_t T) { const real_t u = std::sin(M_PI * t / T); return u * u; }

inline void Det3Inv(const real_t *F, real_t &det, real_t *Finv)
{
   det = F[0] * (F[4] * F[8] - F[5] * F[7]) - F[1] * (F[3] * F[8] - F[5] * F[6]) +
         F[2] * (F[3] * F[7] - F[4] * F[6]);
   const real_t id = 1.0 / det;
   Finv[0] = (F[4] * F[8] - F[5] * F[7]) * id;
   Finv[1] = (F[2] * F[7] - F[1] * F[8]) * id;
   Finv[2] = (F[1] * F[5] - F[2] * F[4]) * id;
   Finv[3] = (F[5] * F[6] - F[3] * F[8]) * id;
   Finv[4] = (F[0] * F[8] - F[2] * F[6]) * id;
   Finv[5] = (F[2] * F[3] - F[0] * F[5]) * id;
   Finv[6] = (F[3] * F[7] - F[4] * F[6]) * id;
   Finv[7] = (F[1] * F[6] - F[0] * F[7]) * id;
   Finv[8] = (F[0] * F[4] - F[1] * F[3]) * id;
}

/// Phi_c(X) - y  (residual of the inverse problem); returns the max norm.
inline real_t PhiResidual(int variant, real_t c, const real_t *X, const real_t *y, real_t *r)
{
   real_t f[3];
   DisplacementField(variant, X, f);
   real_t m = 0.0;
   for (int i = 0; i < 3; i++) { r[i] = X[i] + c * f[i] - y[i]; m = std::max(m, std::fabs(r[i])); }
   return m;
}

/// Solve X + c f(X) = y for X (Newton from X = y, analytic Jacobian, backtracking).
/// Returns the number of iterations; *resid_out = final max-norm residual.
inline int InverseMap(int variant, real_t c, const real_t *y, real_t *X,
                      real_t tol = 1e-14, real_t *resid_out = nullptr, int maxit = 80)
{
   for (int i = 0; i < 3; i++) { X[i] = y[i]; }
   real_t r[3];
   real_t nr = PhiResidual(variant, c, X, y, r);
   int it = 0;
   while (nr > tol)
   {
      MFEM_VERIFY(it < maxit, "InverseMap: Newton did not converge, residual " << nr << " c=" << c);
      real_t G[9], J[9], Ji[9], det;
      GradF(variant, X, G);
      for (int i = 0; i < 9; i++) { J[i] = c * G[i]; }
      J[0] += 1.0; J[4] += 1.0; J[8] += 1.0;
      Det3Inv(J, det, Ji);
      real_t d[3];
      for (int i = 0; i < 3; i++)
      {
         d[i] = -(Ji[3 * i] * r[0] + Ji[3 * i + 1] * r[1] + Ji[3 * i + 2] * r[2]);
      }
      real_t lam = 1.0, Xn[3], rn[3], nn = nr;
      for (int bt = 0; bt < 30; bt++)
      {
         for (int i = 0; i < 3; i++) { Xn[i] = X[i] + lam * d[i]; }
         nn = PhiResidual(variant, c, Xn, y, rn);
         if (nn < nr || nn <= tol) { break; }
         lam *= 0.5;
      }
      for (int i = 0; i < 3; i++) { X[i] = Xn[i]; r[i] = rn[i]; }
      nr = nn;
      it++;
   }
   if (resid_out) { *resid_out = nr; }
   return it;
}

/// Cauchy push-forward of the initial field at physical points.
///   B(x) = F B0(X)/det F,   A(x) = F^{-T} A0(X),   X = Phi_c^{-1}(x),  F = I + c grad f(X).
struct PushForward
{
   Field F;
   Vector B0;      ///< constant mean field added to curl A
   int variant = 1;
   real_t c = 0.0;
   mutable long long newton_iters = 0, npts = 0;

   void Eval(const real_t *x, real_t *B, real_t *A) const
   {
      real_t X[3], G[9], Fm[9], Fi[9], det;
      newton_iters += InverseMap(variant, c, x, X);
      npts++;
      GradF(variant, X, G);
      for (int i = 0; i < 9; i++) { Fm[i] = c * G[i]; }
      Fm[0] += 1.0; Fm[4] += 1.0; Fm[8] += 1.0;
      Det3Inv(Fm, det, Fi);
      Vector Xv(3), v(3);
      for (int i = 0; i < 3; i++) { Xv(i) = X[i]; }
      if (B)
      {
         F.curlA(Xv, v);
         for (int i = 0; i < 3; i++) { v(i) += B0(i); }
         for (int i = 0; i < 3; i++)
         {
            B[i] = (Fm[3 * i] * v(0) + Fm[3 * i + 1] * v(1) + Fm[3 * i + 2] * v(2)) / det;
         }
      }
      if (A)
      {
         F.A(Xv, v);
         for (int i = 0; i < 3; i++) { A[i] = Fi[i] * v(0) + Fi[3 + i] * v(1) + Fi[6 + i] * v(2); }
      }
   }
};

/// VectorCoefficient wrapper (B or A push-forward at the physical point of the integration point).
class PushForwardCoef : public VectorCoefficient
{
public:
   PushForwardCoef(const PushForward &pf_, bool wantA_) : VectorCoefficient(3), pf(pf_), wantA(wantA_) {}
   using VectorCoefficient::Eval;
   void Eval(Vector &V, ElementTransformation &T, const IntegrationPoint &ip) override
   {
      Vector x(3);
      T.Transform(ip, x);
      V.SetSize(3);
      if (wantA) { pf.Eval(x.GetData(), nullptr, V.GetData()); }
      else { pf.Eval(x.GetData(), V.GetData(), nullptr); }
   }
private:
   const PushForward &pf;
   bool wantA;
};

// ---- mesh helpers (node vector is byNODES or byVDIM, like vp::DeformMesh) --------------

/// nodes = labels + c f(labels) (per node dof), then NodesUpdated.
inline void MoveNodes(ParMesh &pm, const Vector &lab, real_t c, int variant)
{
   GridFunction *nodes = pm.GetNodes();
   const FiniteElementSpace *nfes = nodes->FESpace();
   const int nd = nfes->GetNDofs();
   MFEM_VERIFY(lab.Size() == nodes->Size(), "MoveNodes: label size mismatch");
   const bool bynodes = (nfes->GetOrdering() == Ordering::byNODES);
   for (int i = 0; i < nd; i++)
   {
      real_t x[3], f[3];
      for (int d = 0; d < 3; d++) { x[d] = lab(bynodes ? d * nd + i : 3 * i + d); }
      DisplacementField(variant, x, f);
      for (int d = 0; d < 3; d++) { (*nodes)(bynodes ? d * nd + i : 3 * i + d) = x[d] + c * f[d]; }
   }
   pm.NodesUpdated();
}

inline void SetNodes(ParMesh &pm, const Vector &X)
{
   GridFunction *nodes = pm.GetNodes();
   MFEM_VERIFY(X.Size() == nodes->Size(), "SetNodes: size mismatch");
   *nodes = X;
   pm.NodesUpdated();
}

/// Labels of the nodes y (node vector layout of `pm`): lab = Phi_c^{-1}(y), Newton 1e-14.
/// Returns the max residual |Phi_c(lab) - y|_inf (verified < tol_check) and the mean #iterations.
inline real_t ComputeLabels(const ParMesh &pm, int variant, real_t c, const Vector &y, Vector &lab,
                            double *mean_iters = nullptr)
{
   const FiniteElementSpace *nfes = pm.GetNodes()->FESpace();
   const int nd = nfes->GetNDofs();
   const bool bynodes = (nfes->GetOrdering() == Ordering::byNODES);
   lab.SetSize(y.Size());
   real_t rmax = 0.0;
   long long its = 0;
   for (int i = 0; i < nd; i++)
   {
      real_t yy[3], X[3], res;
      for (int d = 0; d < 3; d++) { yy[d] = y(bynodes ? d * nd + i : 3 * i + d); }
      its += InverseMap(variant, c, yy, X, 1e-14, &res);
      rmax = std::max(rmax, res);
      for (int d = 0; d < 3; d++) { lab(bynodes ? d * nd + i : 3 * i + d) = X[d]; }
   }
   if (mean_iters) { *mean_iters = (double)its / std::max(1, nd); }
   return rmax;
}

} // namespace e5
} // namespace vp
