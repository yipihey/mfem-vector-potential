// exp_common.hpp -- shared driver of exp1_static (undeformed) and exp2_distort
// (deformed mesh).  One run = one CSV row.
#pragma once
#include "vp_core.hpp"
#include <cmath>
#include <iostream>
#include <iomanip>

namespace vp
{

inline int RunStatic(int argc, char *argv[], bool distort, const char *expname)
{
   Mpi::Init(argc, argv);
   Hypre::Init();
   MPI_Comm comm = MPI_COMM_WORLD;
   const int rank = Mpi::WorldRank();

   int N = 8, p = 2, q = -1, fvariant = 1;
   real_t eps = 0.0;
   std::string field = "abc", proj = "pt", csv;
   bool b0on = false, pa_cg = true;
   real_t gauge_g = 0.0, gauge_kfac = 2.0;
   OptionsParser args(argc, argv);
   args.AddOption(&N, "-N", "--N", "Cells per direction.");
   args.AddOption(&p, "-p", "--order", "ND order p (B in RT_{p-1}).");
   args.AddOption(&q, "-q", "--geom-order", "Geometric order (default: p).");
   args.AddOption(&field, "-field", "--field", "abc | mod | mod2.");
   args.AddOption(&b0on, "-b0", "--b0", "-no-b0", "--no-b0",
                  "Add the mean field B0=(0.3,-0.2,0.5).");
   args.AddOption(&proj, "-proj", "--proj", "pt (pointwise) | int (integrated).");
   args.AddOption(&csv, "-csv", "--csv", "Append one row to this CSV file.");
   args.AddOption(&gauge_g, "-gauge", "--gauge", "Gauge amplitude g (A += grad chi).");
   args.AddOption(&gauge_kfac, "-gk", "--gauge-k", "Gauge wavenumber in units of pi (2 smooth, 6 rough).");
   args.AddOption(&pa_cg, "-pacg", "--pa-cg", "-no-pacg", "--no-pa-cg", "Partial assembly for the ND-mass CG.");
   if (distort)
   {
      args.AddOption(&eps, "-eps", "--eps", "Deformation amplitude.");
      args.AddOption(&fvariant, "-fvariant", "--fvariant", "Displacement variant 1,2,3.");
   }
   args.Parse();
   if (!args.Good()) { if (rank == 0) { args.PrintUsage(std::cout); } return 1; }
   if (rank == 0) { args.PrintOptions(std::cout); }
   if (q < 0) { q = p; }
   const ProjMode mode = (proj == "int") ? ProjMode::Integrated : ProjMode::Pointwise;
   MFEM_VERIFY(proj == "int" || proj == "pt", "-proj must be pt or int");

   const double T0 = WallTime();
   FieldParams fp;
   const real_t nan = std::nan("");

   CsvRow row;
   // ---- declare all columns (NaN by default) so that invalid runs still
   // produce a row with the identical header ----
   row.Set("exp", expname);
   row.Set("N", N); row.Set("p", p); row.Set("q", q);
   row.Set("field", field); row.Set("b0", (int)b0on); row.Set("proj", proj);
   row.Set("gauge", gauge_g); row.Set("gauge_k", gauge_kfac);
   row.Set("eps", eps); row.Set("fvariant", distort ? fvariant : 0);
   row.Set("np", Mpi::WorldSize());
   row.Set("valid", 1);
   for (const char *k : {"h", "hmin", "hmax", "min_detJ", "max_detJ", "min_ratio", "max_ratio",
                         "n_neg_detJ", "ndofs_h1", "ndofs_nd", "ndofs_rt", "ndofs_l2",
                         "max_DC", "max_CG", "errA_L2", "errB_L2", "normB_L2", "errB_rel",
                         "errA_rel", "div_L2", "div_max", "div_rel_L2", "div_rel_max",
                         "b0proj_err_L2", "b0proj_div_max", "b0proj_div_L2",
                         "energy_h", "energy_exact", "energy_rel_err", "energy_mass_diff",
                         "helicity_h", "helicity_exact", "helicity_abs_err",
                         "flux_err_max", "flux_curlpart_max",
                         "flux_err_x0", "flux_err_xmid", "flux_err_y0", "flux_err_ymid",
                         "flux_err_z0", "flux_err_zmid",
                         "cg_iters", "t_mesh", "t_ops", "t_curl", "t_proj", "t_diag",
                         "t_total", "rss_mb", "peak_rss_mb"})
   {
      row.Set(k, nan);
   }

   // ---- mesh ----
   double t = WallTime();
   std::unique_ptr<ParMesh> pm(MakePeriodicBox(comm, N, q));
   const std::vector<int> js = {0, N / 2};
   std::vector<SliceFamily> slices = RecordSlices(*pm, N, js);
   if (distort && eps != 0.0) { DeformMesh(*pm, eps, fvariant); }
   const auto js_stats = JacobianStats(*pm);
   const auto msz = MeshSize(*pm);
   row.Set("t_mesh", WallTime() - t);
   row.Set("h", msz.hmean); row.Set("hmin", msz.hmin); row.Set("hmax", msz.hmax);
   row.Set("min_detJ", js_stats.min_det); row.Set("max_detJ", js_stats.max_det);
   row.Set("min_ratio", js_stats.min_ratio); row.Set("max_ratio", js_stats.max_ratio);
   row.Set("n_neg_detJ", (long long)js_stats.n_neg);
   if (rank == 0)
   {
      std::cout << "mesh: N=" << N << " q=" << q << " eps=" << eps << " fvariant=" << fvariant
                << "  detJ in [" << js_stats.min_det << ", " << js_stats.max_det
                << "], #nonpositive=" << js_stats.n_neg << "  min/max ratio in ["
                << js_stats.min_ratio << ", " << js_stats.max_ratio << "]" << std::endl;
   }
   if (js_stats.n_neg > 0)
   {
      row.Set("valid", 0);
      row.Set("t_total", WallTime() - T0);
      if (rank == 0) { std::cout << "INVALID mesh (non-positive Jacobian): no diagnostics." << std::endl; }
      if (!csv.empty()) { CsvAppend(comm, csv, row); }
      return 0;
   }

   // ---- spaces and operators ----
   Spaces sp(*pm, p);
   Operators ops(sp);
   row.Set("t_ops", ops.t_assemble);
   row.Set("ndofs_h1", (long long)sp.NdofsH1()); row.Set("ndofs_nd", (long long)sp.NdofsND());
   row.Set("ndofs_rt", (long long)sp.NdofsRT()); row.Set("ndofs_l2", (long long)sp.NdofsL2());
   const real_t dc = ops.MaxAbsDC(), cg = ops.MaxAbsCG();
   row.Set("max_DC", dc); row.Set("max_CG", cg);
   if (rank == 0)
   {
      std::cout << "dofs: H1=" << sp.NdofsH1() << " ND=" << sp.NdofsND() << " RT=" << sp.NdofsRT()
                << " L2=" << sp.NdofsL2() << "   max|D_h C_h|=" << dc << "  max|C_h G_h|=" << cg << std::endl;
   }

   // ---- coefficients ----
   Field F = MakeField(field, fp);
   Vector B0v(3);
   B0v(0) = 0.3; B0v(1) = -0.2; B0v(2) = 0.5;
   Vector B0eff(3);
   B0eff = 0.0;
   if (b0on) { B0eff = B0v; }
   Gauge gg(gauge_g, gauge_kfac);
   VectorFunctionCoefficient Acoef(3, [&](const Vector &x, Vector &A)
   {
      F.A(x, A);
      if (gauge_g != 0.0) { Vector g; gg.grad(x, g); A += g; }
   });
   VectorFunctionCoefficient Bcoef(3, [&](const Vector &x, Vector &B)
   {
      F.curlA(x, B);
      B += B0eff;
   });
   VectorFunctionCoefficient Bcurl(3, F.curlA);
   VectorConstantCoefficient B0c(B0v), B0effc(B0eff);
   const int ord = 2 * p + 2 * q + 2;       // error-norm quadrature order
   const int ord_ex = 2 * p + 2 * q + 6;    // exact-solution quadrature order

   // ---- projections ----
   t = WallTime();
   ParGridFunction a(sp.ND.get()), b0gf(sp.RT.get()), bgf(sp.RT.get()), bcgf(sp.RT.get());
   ProjectA(a, Acoef, mode);
   Vector at, b0t, bt(sp.RT->GetTrueVSize()), bct(sp.RT->GetTrueVSize());
   a.GetTrueDofs(at);
   // exactness check of the B0 projection (always with the nominal B0)
   {
      ParGridFunction chk(sp.RT.get());
      ProjectB(chk, B0c, mode);
      Vector chkt; chk.GetTrueDofs(chkt);
      const auto dd0 = DivergenceDiag(ops, chkt);
      row.Set("b0proj_err_L2", L2ErrorVec(chk, B0c, ord));
      row.Set("b0proj_div_max", dd0.maxabs);
      row.Set("b0proj_div_L2", dd0.l2);
   }
   ProjectB(b0gf, B0effc, mode);
   b0gf.GetTrueDofs(b0t);
   row.Set("t_proj", WallTime() - t);

   // ---- b = Pi_RT(B0) + C_h a ----
   ops.ApplyCurl(at, bct);
   bt = b0t;
   bt += bct;
   bgf.SetFromTrueDofs(bt);
   bcgf.SetFromTrueDofs(bct);
   Vector tmp;
   row.Set("t_curl", ops.ApplyCurlTimed(at, tmp));

   // ---- errors, divergence ----
   t = WallTime();
   const real_t errA = L2ErrorVec(a, Acoef, ord);
   const real_t errB = L2ErrorVec(bgf, Bcoef, ord);
   const real_t normA = L2NormVec(*pm, Acoef, ord_ex);
   const real_t normB = L2NormVec(*pm, Bcoef, ord_ex);
   row.Set("errA_L2", errA); row.Set("errB_L2", errB); row.Set("normB_L2", normB);
   row.Set("errA_rel", errA / normA); row.Set("errB_rel", errB / normB);
   const DivDiag dd = DivergenceDiag(ops, bt);
   row.Set("div_L2", dd.l2); row.Set("div_max", dd.maxabs);
   row.Set("div_rel_L2", dd.l2 / (normB / msz.hmean));
   row.Set("div_rel_max", dd.maxabs / (normB / msz.hmean));

   // ---- energy (matrix-free; matrix based cross-check for small sizes) ----
   const real_t Eh = 0.5 * InnerGF(bgf, bgf, ord);
   const real_t Eex = ExactEnergy(field, fp) + 0.5 * (B0eff * B0eff);
   row.Set("energy_h", Eh); row.Set("energy_exact", Eex);
   row.Set("energy_rel_err", (Eh - Eex) / Eex);
   if (sp.NdofsRT() <= 200000)
   {
      row.Set("energy_mass_diff", (EnergyMass(sp, bt, ord) - Eh) / Eex);
   }

   // ---- helicity (fluctuating part: A_h . (C_h a)) ----
   const real_t Hh = InnerGF(a, bcgf, ord);
   const real_t Hex = ExactHelicity(field, fp);   // gauge-invariant (periodic gauge)
   row.Set("helicity_h", Hh); row.Set("helicity_exact", Hex);
   row.Set("helicity_abs_err", std::fabs(Hh - Hex));

   // ---- slice fluxes ----
   real_t fmax = 0.0, fcmax = 0.0;
   for (auto &fam : slices)
   {
      const real_t fl = SliceFlux(bgf, fam);
      const real_t e = std::fabs(fl - B0eff(fam.dir));
      fmax = std::max(fmax, e);
      fcmax = std::max(fcmax, std::fabs(SliceFlux(bcgf, fam)));
      const char dn[3] = {'x', 'y', 'z'};
      row.Set(std::string("flux_err_") + dn[fam.dir] + (fam.j == 0 ? "0" : "mid"), e);
   }
   row.Set("flux_err_max", fmax); row.Set("flux_curlpart_max", fcmax);
   row.Set("t_diag", WallTime() - t);

   // ---- conditioning proxy ----
   row.Set("cg_iters", MassCGIterations(sp, pa_cg));

   double peak;
   const double rss = RssMB(comm, &peak);
   row.Set("rss_mb", rss); row.Set("peak_rss_mb", peak);
   row.Set("t_total", WallTime() - T0);

   if (rank == 0)
   {
      std::cout << std::scientific << std::setprecision(4)
                << "errA_L2=" << errA << " errB_L2=" << errB << " (rel " << errB / normB << ")"
                << "\ndiv: L2=" << dd.l2 << " max=" << dd.maxabs
                << " rel_L2=" << dd.l2 / (normB / msz.hmean)
                << "\nenergy rel err=" << (Eh - Eex) / Eex << "  helicity h=" << Hh << " exact=" << Hex
                << "\nflux err max=" << fmax << " (curl part max " << fcmax << ")"
                << std::defaultfloat << "\nt_assemble=" << ops.t_assemble
                << " s  rss=" << rss << " MB  t_total=" << WallTime() - T0 << " s" << std::endl;
   }
   if (!csv.empty()) { CsvAppend(comm, csv, row); }
   return 0;
}

} // namespace vp
