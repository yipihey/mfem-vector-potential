// exp3_remap.cpp -- E3: one remap M1 -> M2 with the A-route (ND transfer + B = B0 + C_h a)
// and the B-route (RT transfer).  One CSV row per transfer operator.
//
//   mpirun -np 1 bin/exp3_remap -N 8 -p 2 -field abc -b0 -eps1 0.3 -m2 uniform -ops all -csv out.csv
#include "exp_transfer_common.hpp"
#include <cstdlib>

using namespace vp;

int main(int argc, char *argv[])
{
   Mpi::Init(argc, argv);
   Hypre::Init();
   MPI_Comm comm = MPI_COMM_WORLD;
   const int rank = Mpi::WorldRank();

   int N = 8, p = 2, q = -1, fv1 = 1, nq = -1;
   real_t eps1 = 0.3, gauge_g = 0.0, gauge_k = 2.0;
   std::string field = "abc", m2s = "uniform", opss = "all", csv, tag, gfixstr = "none";
   bool b0on = false, gradfrac = true;
   OptionsParser args(argc, argv);
   args.AddOption(&N, "-N", "--N", "Cells per direction of M1.");
   args.AddOption(&p, "-p", "--order", "ND order p (B in RT_{p-1}).");
   args.AddOption(&q, "-q", "--geom-order", "Geometric order (default p).");
   args.AddOption(&field, "-field", "--field", "abc | mod | mod2.");
   args.AddOption(&b0on, "-b0", "--b0", "-no-b0", "--no-b0", "Add mean field B0=(0.3,-0.2,0.5).");
   args.AddOption(&eps1, "-eps1", "--eps1", "Deformation amplitude of M1.");
   args.AddOption(&fv1, "-fv1", "--fv1", "Displacement variant of M1.");
   args.AddOption(&m2s, "-m2", "--m2", "uniform | deformed:eps2:fv2 | fine:N2");
   args.AddOption(&gauge_g, "-gauge", "--gauge", "Gauge amplitude g: A1 += grad chi on M1 before the transfer.");
   args.AddOption(&gauge_k, "-gk", "--gauge-k", "Gauge wavenumber in units of pi (2 smooth, 6 rough).");
   args.AddOption(&gfixstr, "-gfix", "--gfix", "Gauge fix applied AFTER the A-route transfer: none | coulomb | jacobi:k "
                  "(named -gfix because -gauge is the gauge-perturbation amplitude here).");
   args.AddOption(&opss, "-ops", "--ops", "comma list of A_pt,A_int,A_l2,B_pt,B_int,B_l2,B_l2c or all");
   args.AddOption(&nq, "-nq", "--nq", "Gauss points per direction for the *_int transfers (default p+2).");
   args.AddOption(&gradfrac, "-gradfrac", "--gradfrac", "-no-gradfrac", "--no-gradfrac", "Compute the gradient fraction of A.");
   args.AddOption(&csv, "-csv", "--csv", "Append rows to this CSV file.");
   args.AddOption(&tag, "-tag", "--tag", "Free label stored in the CSV (series name).");
   args.Parse();
   if (!args.Good()) { if (rank == 0) { args.PrintUsage(std::cout); } return 1; }
   if (rank == 0) { args.PrintOptions(std::cout); }
   if (q < 0) { q = p; }
   const double T0 = WallTime();
   double Tlast = T0;
   auto Tick = [&](const char *what)
   {
      const double t = WallTime();
      if (rank == 0) { std::cout << "  [t] " << what << " " << t - Tlast << " s" << std::endl; }
      Tlast = t;
   };

   MeshSpec s1; s1.N = N; s1.eps = eps1; s1.fv = fv1;
   const MeshSpec s2 = ParseMeshSpec(m2s, N);
   Problem pr(field, b0on, gauge_g, gauge_k);
   const bool with_gauge = (gauge_g != 0.0);

   MeshCase M1(comm, s1, p, q), M2(comm, s2, p, q);
   Tick("meshes+ops");
   const GaugeSpec gspec = GaugeSpec::Parse(gfixstr);
   std::unique_ptr<CoulombGauge> cg2;
   if (gspec.On()) { cg2.reset(new CoulombGauge(*M2.sp, *M2.ops)); Tick("gauge setup"); }
   MeshState S1(M1, pr), S2(M2, pr);
   Tick("mesh state (b0, diag)");
   if (rank == 0)
   {
      std::cout << "M1 " << s1.Str() << " detJ [" << M1.js.min_det << "," << M1.js.max_det << "] ND="
                << M1.sp->NdofsND() << " RT=" << M1.sp->NdofsRT() << "\nM2 " << s2.Str() << " detJ ["
                << M2.js.min_det << "," << M2.js.max_det << "] ND=" << M2.sp->NdofsND()
                << " RT=" << M2.sp->NdofsRT() << std::endl;
   }

   // ---- states on M1 ----
   Vector a1, a1ng, b1, tmp;
   S1.ProjectExactA(pr, with_gauge, a1);
   S1.ProjectExactA(pr, false, a1ng);
   S1.Compose(a1, b1);
   Vector b1chk;
   S1.Compose(a1ng, b1chk);
   const real_t b1_gauge_dev = S1.diag->RTDistance(b1, b1chk);   // should be roundoff
   const StateDiag d1 = S1.diag->Evaluate(&a1, b1, gradfrac);
   Tick("M1 projection+diag");
   // ---- reference (floors) on M2: projection of the exact field ----
   Vector a2ref, b2refA;
   S2.ProjectExactA(pr, false, a2ref);
   S2.Compose(a2ref, b2refA);
   const StateDiag dref_A = S2.diag->Evaluate(&a2ref, b2refA, gradfrac);
   // reference of the B-route: Pi_RT^int(B_exact) = B0 + C_h Pi_ND^int A (commuting diagram), so dref_B == dref_A
   const StateDiag dref_B = dref_A;
   Tick("M2 reference projections+diag");
   // M1 -> exact quantities
   const real_t E1 = d1.energy, H1 = d1.helicity;

   TransferOptions topt;
   topt.nq = nq;
   RemoteEvaluator ev(*M1.pm);
   Tick("evaluator setup");

   for (TransferKind kind : ParseKindList(opss))
   {
      const double Tk = WallTime();
      const bool isA = IsARoute(kind);
      ParGridFunction a1gf(M1.sp->ND.get()), a1nggf(M1.sp->ND.get()), b1gf(M1.sp->RT.get());
      a1gf.SetFromTrueDofs(a1);
      a1nggf.SetFromTrueDofs(a1ng);
      b1gf.SetFromTrueDofs(b1);
      Transfer T(kind, ev, isA ? *M2.sp->ND : *M2.sp->RT, M2.ops.get(), topt);
      ParGridFunction dst(isA ? M2.sp->ND.get() : M2.sp->RT.get());
      T.Apply(isA ? a1gf : b1gf, dst);
      TransferStats st = T.Stats();   // copy (the gauge-free apply below overwrites)
      Vector x2, b2, a2;
      dst.GetTrueDofs(x2);
      StateDiag d2;
      real_t pollution = std::nan(""), pollution_a = std::nan("");
      if (isA)
      {
         a2 = x2;
         S2.Compose(a2, b2);
         d2 = S2.diag->Evaluate(&a2, b2, gradfrac);
         if (with_gauge)
         {
            ParGridFunction dst2(M2.sp->ND.get());
            T.Apply(a1nggf, dst2);
            Vector a2ng, b2ng;
            dst2.GetTrueDofs(a2ng);
            S2.Compose(a2ng, b2ng);
            Vector zero(b2.Size());
            zero = 0.0;
            const real_t nb = S2.diag->RTDistance(b2, zero);
            pollution = S2.diag->RTDistance(b2, b2ng) / nb;
            // a-difference measured in the ND L2 norm relative to ||a2ng||
            Vector da(a2);
            da -= a2ng;
            ParGridFunction g1(M2.sp->ND.get()), g2(M2.sp->ND.get());
            g1.SetFromTrueDofs(da);
            g2.SetFromTrueDofs(a2ng);
            pollution_a = std::sqrt(InnerGF(g1, g1, 2 * p + 2 * q + 2) /
                                    InnerGF(g2, g2, 2 * p + 2 * q + 2));
         }
      }
      else
      {
         b2 = x2;
         d2 = S2.diag->Evaluate(nullptr, b2, false);
      }
      // optional gauge fix after the transfer (A-route only): B must not change
      GaugeFixStats gst;
      StateDiag d2pre = d2;
      real_t g_dB = std::nan(""), g_dH = std::nan("");
      const bool gon = isA && gspec.On();
      if (gon)
      {
         const Vector b2pre = b2;
         gst = gspec.Apply(*cg2, a2, true);
         S2.Compose(a2, b2);
         d2 = S2.diag->Evaluate(&a2, b2, gradfrac);
         Vector zero(b2.Size());
         zero = 0.0;
         g_dB = S2.diag->RTDistance(b2pre, b2) / S2.diag->RTDistance(b2pre, zero);
         g_dH = (d2.helicity - d2pre.helicity) / d2pre.hel_scale;
      }
      const double t_total = WallTime() - Tk;
      Tick(KindName(kind));

      CsvRow r;
      r.Set("exp", "exp3_remap"); r.Set("tag", tag);
      r.Set("op", KindName(kind)); r.Set("route", isA ? "A" : "B");
      r.Set("p", p); r.Set("q", q); r.Set("field", field); r.Set("b0", (int)b0on);
      r.Set("N1", s1.N); r.Set("eps1", (double)s1.eps); r.Set("fv1", s1.fv);
      r.Set("N2", s2.N); r.Set("eps2", (double)s2.eps); r.Set("fv2", s2.fv);
      r.Set("m2", m2s);
      r.Set("gauge", (double)gauge_g); r.Set("gauge_k", (double)gauge_k);
      r.Set("nq", T.Kind() == TransferKind::A_int || T.Kind() == TransferKind::B_int
                  ? (topt.nq > 0 ? topt.nq : DefaultTransferNq(p)) : 0);
      r.Set("np", Mpi::WorldSize());
      r.Set("min_detJ1", (double)M1.js.min_det); r.Set("min_detJ2", (double)M2.js.min_det);
      r.Set("ndofs_nd1", (long long)M1.sp->NdofsND()); r.Set("ndofs_rt1", (long long)M1.sp->NdofsRT());
      r.Set("ndofs_nd2", (long long)M2.sp->NdofsND()); r.Set("ndofs_rt2", (long long)M2.sp->NdofsRT());
      r.Set("b1_gauge_dev", (double)b1_gauge_dev);
      SetDiag(r, "m1_", d1);
      SetDiag(r, "ref_A_", dref_A);
      SetDiag(r, "ref_B_", dref_B);
      SetDiag(r, "m2_", d2);
      r.Set("m2_dE_rel1", (double)((d2.energy - E1) / E1));
      r.Set("m2_dH_rel1", isA ? (double)((d2.helicity - H1) / std::fabs(H1)) : std::nan(""));
      r.Set("m2_dH_scaled", isA ? (double)((d2.helicity - H1) / d1.hel_scale) : std::nan(""));
      r.Set("m2_aratio", isA ? (double)(d2.a_L2 / d1.a_L2) : std::nan(""));
      r.Set("energy_exact", (double)pr.ExactE());
      r.Set("helicity_exact", (double)pr.ExactH());
      r.Set("pollution_B", (double)pollution);
      r.Set("pollution_A", (double)pollution_a);
      r.Set("gfix", isA ? gspec.str : std::string("none"));
      r.Set("gfix_iters", gon ? gst.iters : 0);
      r.Set("gfix_time", gon ? gst.time : 0.0);
      r.Set("gfix_setup", gon ? cg2->SetupTime() : 0.0);
      r.Set("gfix_resid_before", gon ? (double)gst.resid_before : std::nan(""));
      r.Set("gfix_resid_after", gon ? (double)gst.resid_after : std::nan(""));
      r.Set("gfix_gf_before", gon ? (double)gst.grad_fraction_before : std::nan(""));
      r.Set("gfix_gf_after", gon ? (double)gst.grad_fraction_after : std::nan(""));
      r.Set("gfix_dB_rel", (double)g_dB);
      r.Set("gfix_dH_scaled", (double)g_dH);
      r.Set("pre_grad_frac", isA ? (double)d2pre.grad_frac : std::nan(""));
      r.Set("pre_helicity", isA ? (double)d2pre.helicity : std::nan(""));
      r.Set("t_setup", st.t_setup); r.Set("t_locate", st.t_locate); r.Set("t_collect", st.t_collect);
      r.Set("t_eval", st.t_eval); r.Set("t_dofs", st.t_dofs); r.Set("t_solve", st.t_solve);
      r.Set("t_apply", st.t_apply); r.Set("t_oneshot", st.t_setup + st.t_apply);
      r.Set("t_total_op", t_total);
      r.Set("iters", st.iters); r.Set("iters_l2", st.iters_l2); r.Set("inner_iters", st.inner_iters);
      r.Set("div_before_clean", st.div_before_clean);
      r.Set("npoints", st.npoints); r.Set("nnotfound", st.nnotfound); r.Set("nshifted", st.nshifted);
      double peak;
      r.Set("rss_mb", RssMB(comm, &peak)); r.Set("peak_rss_mb", peak);
      if (!csv.empty()) { CsvAppend(comm, csv, r); }
      if (rank == 0)
      {
         std::cout << std::scientific << std::setprecision(3) << std::left << std::setw(6) << KindName(kind)
                   << " errB_rel=" << d2.errB_rel << " (ref " << (isA ? dref_A.errB_rel : dref_B.errB_rel)
                   << ") div_rel=" << d2.div_rel_l2 << " dE/E1=" << (d2.energy - E1) / E1
                   << " t_oneshot=" << st.t_setup + st.t_apply << " iters=" << st.iters
                   << " notfound=" << st.nnotfound;
         if (gon)
         {
            std::cout << " | gfix " << gspec.str << ": iters " << gst.iters << " time " << gst.time << " gradfrac "
                      << d2pre.grad_frac << "->" << d2.grad_frac << " dB " << g_dB << " dH " << g_dH;
         }
         if (with_gauge && isA) { std::cout << " pollutionB=" << pollution << " A=" << pollution_a; }
         std::cout << std::defaultfloat << std::endl;
      }
   }
   if (rank == 0) { std::cout << "total wall " << WallTime() - T0 << " s" << std::endl; }
   return 0;
}
