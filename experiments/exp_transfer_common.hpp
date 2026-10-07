// exp_transfer_common.hpp -- helpers shared by exp3_remap and exp4_repeat.
#pragma once
#include "vp_transfer.hpp"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>

namespace vp
{

/// Gauss points per direction for the integrated dofs of the ANALYTIC fields (vs 8 in vp_core's default:
/// 6 gives < 1e-9 relative error for k <= 4 pi, N >= 8).
constexpr int kExactNq = 6;

/// Mesh description: N^3 periodic box, optionally deformed x' = x + eps f_fv(x).
struct MeshSpec
{
   int N = 8;
   real_t eps = 0.0;
   int fv = 1;
   std::string Str() const
   {
      std::ostringstream s;
      s << "N" << N << (eps != 0.0 ? "_def" : "_uni");
      if (eps != 0.0) { s << eps << "f" << fv; }
      return s.str();
   }
};

/// "uniform" | "deformed:eps:fv" | "fine:N2"  (N1 = default N)
inline MeshSpec ParseMeshSpec(const std::string &s, int N_default)
{
   MeshSpec m;
   m.N = N_default;
   std::vector<std::string> tok;
   std::stringstream ss(s);
   std::string t;
   while (std::getline(ss, t, ':')) { tok.push_back(t); }
   MFEM_VERIFY(!tok.empty(), "empty mesh spec");
   if (tok[0] == "uniform") { m.eps = 0.0; }
   else if (tok[0] == "deformed")
   {
      MFEM_VERIFY(tok.size() >= 2, "deformed:eps[:fv]");
      m.eps = std::stod(tok[1]);
      m.fv = tok.size() > 2 ? std::stoi(tok[2]) : 1;
   }
   else if (tok[0] == "fine")
   {
      MFEM_VERIFY(tok.size() >= 2, "fine:N2");
      m.N = std::stoi(tok[1]);
      m.eps = 0.0;
   }
   else { MFEM_ABORT("unknown mesh spec '" << s << "'"); }
   return m;
}

/// A mesh with its spaces, discrete operators and slice families.
struct MeshCase
{
   MeshSpec spec;
   std::unique_ptr<ParMesh> pm;
   std::vector<SliceFamily> slices;
   std::unique_ptr<Spaces> sp;
   std::unique_ptr<Operators> ops;
   JacobianStatsResult js;
   MeshSizeResult msz;
   MeshCase(MPI_Comm comm, const MeshSpec &s, int p, int q) : spec(s)
   {
      pm.reset(MakePeriodicBox(comm, s.N, q));
      slices = RecordSlices(*pm, s.N, {0, s.N / 2});
      if (s.eps != 0.0) { DeformMesh(*pm, s.eps, s.fv); }
      js = JacobianStats(*pm);
      msz = MeshSize(*pm);
      MFEM_VERIFY(js.n_neg == 0, "mesh " << s.Str() << " has non-positive Jacobians");
      sp.reset(new Spaces(*pm, p));
      ops.reset(new Operators(*sp));
   }
};

/// Exact data of the test problem.
struct Problem
{
   std::string field = "abc";
   FieldParams fp;
   bool b0on = false;
   Vector B0eff;
   Field F;
   real_t gauge_g = 0.0, gauge_kfac = 2.0;
   Problem(const std::string &f, bool b0, real_t g, real_t gk)
      : field(f), b0on(b0), gauge_g(g), gauge_kfac(gk)
   {
      F = MakeField(field, fp);
      B0eff.SetSize(3);
      B0eff = 0.0;
      if (b0on) { B0eff(0) = 0.3; B0eff(1) = -0.2; B0eff(2) = 0.5; }
   }
   /// A coefficient (optionally with gauge term grad chi)
   std::unique_ptr<VectorFunctionCoefficient> ACoef(bool with_gauge) const
   {
      auto Af = F.A;
      real_t g = (with_gauge ? gauge_g : 0.0), gk = gauge_kfac;
      return std::unique_ptr<VectorFunctionCoefficient>(new VectorFunctionCoefficient(
         3, [Af, g, gk](const Vector &x, Vector &A)
      {
         Af(x, A);
         if (g != 0.0) { Gauge gg(g, gk); Vector v; gg.grad(x, v); A += v; }
      }));
   }
   std::unique_ptr<VectorFunctionCoefficient> BCoef() const   // exact total B
   {
      auto cf = F.curlA;
      Vector b0 = B0eff;
      return std::unique_ptr<VectorFunctionCoefficient>(new VectorFunctionCoefficient(
         3, [cf, b0](const Vector &x, Vector &B) { cf(x, B); B += b0; }));
   }
   real_t ExactE() const { return vp::ExactEnergy(field, fp) + 0.5 * (B0eff * B0eff); }
   real_t ExactH() const { return vp::ExactHelicity(field, fp); }
};

/// Everything attached to one mesh for the state handling.
struct MeshState
{
   MeshCase &mc;
   Vector b0_t;                      // Pi_RT^int(B0)
   std::unique_ptr<MeshDiag> diag;
   MeshState(MeshCase &m, const Problem &pr) : mc(m)
   {
      ParGridFunction b0(mc.sp->RT.get());
      VectorConstantCoefficient B0c(pr.B0eff);
      // constant field: integrand adj(J) B0 is a polynomial of degree <= 2q per direction
      ProjectB(b0, B0c, ProjMode::Integrated, mc.sp->pmesh->GetNodalFESpace()->GetMaxElementOrder() + 2);
      b0.GetTrueDofs(b0_t);
      diag.reset(new MeshDiag(*mc.sp, *mc.ops, mc.slices, pr.field, pr.fp, pr.B0eff, b0_t));
   }
   /// a_t <- Pi_ND^int (exact A [+gauge])
   void ProjectExactA(const Problem &pr, bool with_gauge, Vector &a_t) const
   {
      auto Ac = pr.ACoef(with_gauge);
      ParGridFunction a(mc.sp->ND.get());
      ProjectA(a, *Ac, ProjMode::Integrated, kExactNq);
      a.GetTrueDofs(a_t);
   }
   /// b_t <- Pi_RT^int (exact B)
   void ProjectExactB(const Problem &pr, Vector &b_t) const
   {
      auto Bc = pr.BCoef();
      ParGridFunction b(mc.sp->RT.get());
      ProjectB(b, *Bc, ProjMode::Integrated, kExactNq);
      b.GetTrueDofs(b_t);
   }
   void Compose(const Vector &a_t, Vector &b_t) const { ComposeB(*mc.ops, b0_t, a_t, b_t); }
};

inline void SetDiag(CsvRow &r, const std::string &pre, const StateDiag &d)
{
   r.Set(pre + "errB_L2", (double)d.errB_L2);
   r.Set(pre + "errB_rel", (double)d.errB_rel);
   r.Set(pre + "div_L2", (double)d.div_l2);
   r.Set(pre + "div_max", (double)d.div_max);
   r.Set(pre + "div_rel_L2", (double)d.div_rel_l2);
   r.Set(pre + "div_rel_max", (double)d.div_rel_max);
   r.Set(pre + "energy", (double)d.energy);
   r.Set(pre + "energy_rel_ex", (double)d.energy_rel_ex);
   r.Set(pre + "helicity", (double)d.helicity);
   r.Set(pre + "hel_scale", (double)d.hel_scale);
   r.Set(pre + "flux_err", (double)d.flux_err);
   r.Set(pre + "a_L2", (double)d.a_L2);
   r.Set(pre + "coulomb", (double)d.coulomb_res);
   r.Set(pre + "grad_frac", (double)d.grad_frac);
}

} // namespace vp
