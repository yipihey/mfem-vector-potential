// exp7_perf.cpp -- E7: cost of ONE remap M1 (f1, eps 0.3) -> M2 (uniform), field abc, b0 on.
//
//   mpirun -np K bin/exp7_perf -N 16 -p 2 -op A_int -gauge coulomb -reps 2 -csv results/exp7_perf.csv
//   -op none : baseline, only mesh + spaces + operator construction (both meshes).
//
// All stages are timed with MPI_Wtime between barriers, max over ranks.  Memory is VmRSS /
// VmHWM (max over ranks, and sum over ranks).  Diagnostics (error vs the exact field) are
// computed after the last memory snapshot and are not part of any timing.
#include "exp_transfer_common.hpp"
#include <fstream>

using namespace vp;

namespace
{
MPI_Comm comm;

double MaxT(double t)
{
   MPI_Allreduce(MPI_IN_PLACE, &t, 1, MPI_DOUBLE, MPI_MAX, comm);
   return t;
}
double Tic() { MPI_Barrier(comm); return MPI_Wtime(); }
double Toc(double t0) { return MaxT(MPI_Wtime() - t0); }

void Mem(double &rss_max, double &rss_sum, double &hwm_max, double &hwm_sum)
{
   double rss = 0, hwm = 0;
   std::ifstream f("/proc/self/status");
   std::string line;
   long kb;
   while (std::getline(f, line))
   {
      if (line.compare(0, 6, "VmRSS:") == 0 && sscanf(line.c_str() + 6, "%ld", &kb) == 1) { rss = kb / 1024.0; }
      else if (line.compare(0, 6, "VmHWM:") == 0 && sscanf(line.c_str() + 6, "%ld", &kb) == 1) { hwm = kb / 1024.0; }
   }
   rss_max = rss; hwm_max = hwm; rss_sum = rss; hwm_sum = hwm;
   MPI_Allreduce(MPI_IN_PLACE, &rss_max, 1, MPI_DOUBLE, MPI_MAX, comm);
   MPI_Allreduce(MPI_IN_PLACE, &hwm_max, 1, MPI_DOUBLE, MPI_MAX, comm);
   MPI_Allreduce(MPI_IN_PLACE, &rss_sum, 1, MPI_DOUBLE, MPI_SUM, comm);
   MPI_Allreduce(MPI_IN_PLACE, &hwm_sum, 1, MPI_DOUBLE, MPI_SUM, comm);
}
long long GMax(long long v)
{
   MPI_Allreduce(MPI_IN_PLACE, &v, 1, MPI_LONG_LONG, MPI_MAX, comm);
   return v;
}
} // namespace

int main(int argc, char *argv[])
{
   Mpi::Init(argc, argv);
   Hypre::Init();
   comm = MPI_COMM_WORLD;
   const int rank = Mpi::WorldRank();

   int N = 8, p = 2, q = -1, reps = 2, nq = -1;
   real_t epsB = 0.3;
   std::string field = "abc", opstr = "A_int", csv, tag, gaugestr = "none";
   OptionsParser args(argc, argv);
   args.AddOption(&N, "-N", "--N", "Cells per direction.");
   args.AddOption(&p, "-p", "--order", "ND order p.");
   args.AddOption(&q, "-q", "--geom-order", "Geometric order (default p).");
   args.AddOption(&field, "-field", "--field", "abc | mod | mod2.");
   args.AddOption(&epsB, "-eps", "--eps", "Deformation of the source mesh M1 (f1).");
   args.AddOption(&opstr, "-op", "--op", "A_pt|A_int|A_l2|B_pt|B_int|B_l2|B_l2c|none.");
   args.AddOption(&nq, "-nq", "--nq", "Gauss points per direction for *_int.");
   args.AddOption(&gaugestr, "-gauge", "--gauge", "none | coulomb | jacobi:k (A-route only).");
   args.AddOption(&reps, "-reps", "--reps", "Timed repetitions of the apply/gauge/curl stages (min kept).");
   args.AddOption(&csv, "-csv", "--csv", "Append the row to this CSV file.");
   args.AddOption(&tag, "-tag", "--tag", "Series label.");
   args.Parse();
   if (!args.Good()) { if (rank == 0) { args.PrintUsage(std::cout); } return 1; }
   if (q < 0) { q = p; }
   const bool baseline = (opstr == "none");
   const TransferKind kind = baseline ? TransferKind::A_int : ParseKind(opstr);
   const bool isA = IsARoute(kind);
   const GaugeSpec gspec = GaugeSpec::Parse(isA ? gaugestr : "none");
   const int np = Mpi::WorldSize();

   CsvRow row;
   row.Set("exp", "exp7_perf"); row.Set("tag", tag);
   row.Set("op", baseline ? "none" : KindName(kind));
   row.Set("gauge", gspec.str);
   row.Set("opname", baseline ? std::string("none") : (std::string(KindName(kind)) + (gspec.On() ? "+" + gspec.str : "")));
   row.Set("N", N); row.Set("p", p); row.Set("q", q); row.Set("np", np);
   row.Set("field", field); row.Set("eps", (double)epsB); row.Set("reps", reps);

   double m_rss, m_rsum, m_hwm, m_hsum;
   Mem(m_rss, m_rsum, m_hwm, m_hsum);
   row.Set("rss_start_mb", m_rss);

   // ---------------- 1. mesh + spaces + operators (both meshes) ----------------
   double t0 = Tic();
   std::unique_ptr<ParMesh> pm1(MakePeriodicBox(comm, N, q));
   std::unique_ptr<ParMesh> pm2(MakePeriodicBox(comm, N, q));
   const double t_mesh_gen = Toc(t0);
   t0 = Tic();
   DeformMesh(*pm1, epsB, 1);
   const double t_deform = Toc(t0);
   t0 = Tic();
   std::unique_ptr<Spaces> sp1(new Spaces(*pm1, p));
   std::unique_ptr<Spaces> sp2(new Spaces(*pm2, p));
   const double t_spaces = Toc(t0);
   t0 = Tic();
   std::unique_ptr<Operators> ops1(new Operators(*sp1));
   std::unique_ptr<Operators> ops2(new Operators(*sp2));
   const double t_ops = Toc(t0);
   const double t_ops2_internal = MaxT(ops2->t_assemble);
   row.Set("t_mesh", t_mesh_gen + t_deform); row.Set("t_spaces", t_spaces); row.Set("t_ops", t_ops);
   row.Set("t_build", t_mesh_gen + t_deform + t_spaces + t_ops);
   row.Set("t_ops_M2_internal", t_ops2_internal);
   Mem(m_rss, m_rsum, m_hwm, m_hsum);
   row.Set("rss_build_mb", m_rss); row.Set("rss_build_sum_mb", m_rsum);
   row.Set("peak_build_mb", m_hwm); row.Set("peak_build_sum_mb", m_hsum);

   // DOF counts etc
   row.Set("ne", (long long)pm2->GetGlobalNE());
   row.Set("nd", sp2->NdofsND()); row.Set("rt", sp2->NdofsRT()); row.Set("l2", sp2->NdofsL2()); row.Set("h1", sp2->NdofsH1());
   row.Set("loc_nd_max", GMax(sp2->ND->GetTrueVSize())); row.Set("loc_rt_max", GMax(sp2->RT->GetTrueVSize()));
   row.Set("loc_ne_max", GMax(pm2->GetNE()));
   row.Set("h", 1.0 / N);
   row.Set("nnz_C", (long long)ops2->C->NNZ()); row.Set("nnz_D", (long long)ops2->D->NNZ());
   row.Set("nnz_G", (long long)ops2->G->NNZ());
   row.Set("nnz_C_per_nd", (double)ops2->C->NNZ() / sp2->NdofsND());
   row.Set("nnz_D_per_rt", (double)ops2->D->NNZ() / sp2->NdofsRT());

   if (baseline)
   {
      row.Set("mem_note", "baseline");
      if (!csv.empty()) { CsvAppend(comm, csv, row); }
      if (rank == 0)
      {
         std::cout << "baseline N=" << N << " p=" << p << " np=" << np << " build " << t_mesh_gen + t_deform + t_spaces + t_ops
                   << " s  (mesh " << t_mesh_gen + t_deform << " spaces " << t_spaces << " ops " << t_ops << ")  ND=" << sp2->NdofsND()
                   << " rss " << m_rss << " MB peak " << m_hwm << std::endl;
      }
      return 0;
   }

   // ---------------- source state on M1 (initial condition, not timed) ----------------
   Problem pr(field, true, 0.0, 2.0);
   Vector a1_t, b1_t, b0_1;
   {
      auto Ac = pr.ACoef(false);
      ParGridFunction a(sp1->ND.get());
      ProjectA(a, *Ac, ProjMode::Integrated, kExactNq);
      a.GetTrueDofs(a1_t);
      ParGridFunction b0(sp1->RT.get());
      VectorConstantCoefficient B0c(pr.B0eff);
      ProjectB(b0, B0c, ProjMode::Integrated, q + 2);
      b0.GetTrueDofs(b0_1);
      ComposeB(*ops1, b0_1, a1_t, b1_t);
   }
   ParGridFunction src(isA ? sp1->ND.get() : sp1->RT.get());
   src.SetFromTrueDofs(isA ? a1_t : b1_t);
   ParGridFunction dst(isA ? sp2->ND.get() : sp2->RT.get());

   // ---------------- 2. RemoteEvaluator (GSLIB setup on M1) ----------------
   t0 = Tic();
   std::unique_ptr<RemoteEvaluator> ev(new RemoteEvaluator(*pm1));
   const double t_evsetup = Toc(t0);
   row.Set("t_ev_setup", t_evsetup);

   // ---------------- 3. transfer plan (points collection, FindPoints, mass/AMG setup) ----------------
   TransferOptions topt;
   topt.nq = nq;
   t0 = Tic();
   std::unique_ptr<Transfer> plan(new Transfer(kind, *ev, isA ? *sp2->ND : *sp2->RT, ops2.get(), topt));
   const double t_plan = Toc(t0);
   const TransferStats &st0 = plan->Stats();
   const double t_collect = MaxT(st0.t_collect), t_locate = MaxT(st0.t_locate);
   row.Set("t_plan", t_plan);
   row.Set("t_collect", t_collect); row.Set("t_findpoints", t_locate);
   row.Set("t_plan_other", std::max(0.0, t_plan - t_collect - t_locate));   // mass operator, rhs, AMG (non-lazy part)
   row.Set("npoints", (long long)st0.npoints);
   row.Set("nnotfound", (long long)st0.nnotfound); row.Set("nshifted", (long long)st0.nshifted);
   Mem(m_rss, m_rsum, m_hwm, m_hsum);
   row.Set("rss_plan_mb", m_rss);

   // ---------------- 4. Apply: first call, then reps (min over reps by t_apply) ----------------
   double a_first = 0;
   {
      const double ta = Tic();
      plan->Apply(src, dst);
      a_first = Toc(ta);
   }
   row.Set("t_apply_first", a_first);
   double best = 1e300, b_eval = 0, b_dofs = 0, b_solve = 0;
   int its = 0, its_l2 = 0, inner = 0, rounds = 0;
   double div_before = std::nan("");
   for (int r = 0; r < reps; r++)
   {
      const double ta = Tic();
      plan->Apply(src, dst);
      const double tt = Toc(ta);
      const TransferStats &s = plan->Stats();
      if (tt < best)
      {
         best = tt; b_eval = MaxT(s.t_eval); b_dofs = MaxT(s.t_dofs); b_solve = MaxT(s.t_solve);
         its = s.iters; its_l2 = s.iters_l2; inner = s.inner_iters; rounds = s.clean_rounds;
         div_before = s.div_before_clean;
      }
      else { MaxT(0); MaxT(0); MaxT(0); }
   }
   row.Set("t_apply", best); row.Set("t_eval", b_eval); row.Set("t_dofs", b_dofs); row.Set("t_solve", b_solve);
   row.Set("iters", its); row.Set("iters_l2", its_l2); row.Set("inner_iters", inner); row.Set("clean_rounds", rounds);
   Vector x_t;
   dst.GetTrueDofs(x_t);   // result of the last apply (identical every time)

   // ---------------- 5. gauge fix (A-route) ----------------
   double t_gsetup = 0, t_gfirst = 0, t_gauge = 1e300;
   int g_iters = 0;
   if (isA && gspec.On())
   {
      t0 = Tic();
      CoulombGauge cg(*sp2, *ops2);
      t_gsetup = Toc(t0);
      Vector a_work;
      {
         a_work = x_t;
         const double tg = Tic();
         GaugeFixStats gs = gspec.Apply(cg, a_work, false);
         t_gfirst = Toc(tg);
         g_iters = gs.iters;
      }
      for (int r = 0; r < reps; r++)
      {
         a_work = x_t;
         const double tg = Tic();
         GaugeFixStats gs = gspec.Apply(cg, a_work, false);
         const double tt = Toc(tg);
         if (tt < t_gauge) { t_gauge = tt; g_iters = gs.iters; }
      }
      x_t = a_work;     // gauge-fixed a (used below)
      Mem(m_rss, m_rsum, m_hwm, m_hsum);
      row.Set("rss_gauge_mb", m_rss);
   }
   else { t_gauge = 0.0; }
   row.Set("t_gauge_setup", t_gsetup); row.Set("t_gauge_first", t_gfirst); row.Set("t_gauge", t_gauge);
   row.Set("gauge_iters", g_iters);

   // ---------------- 6. curl apply and composition ----------------
   double t_curl = 0, t_comp = 0, t_b0 = 0;
   Vector b_t, b0_2;
   if (isA)
   {
      Vector tmp;
      t_curl = 1e300;
      for (int r = 0; r < reps; r++) { t_curl = std::min(t_curl, ops2->ApplyCurlTimed(x_t, tmp, 10, 0.0)); }
      t0 = Tic();
      {
         ParGridFunction b0(sp2->RT.get());
         VectorConstantCoefficient B0c(pr.B0eff);
         ProjectB(b0, B0c, ProjMode::Integrated, q + 2);
         b0.GetTrueDofs(b0_2);
      }
      t_b0 = Toc(t0);
      t_comp = 1e300;
      for (int r = 0; r < reps; r++)
      {
         const double tc = Tic();
         for (int i = 0; i < 10; i++) { ComposeB(*ops2, b0_2, x_t, b_t); }
         t_comp = std::min(t_comp, Toc(tc) / 10.0);
      }
   }
   else { b_t = x_t; }
   row.Set("t_curl", t_curl); row.Set("t_b0proj", t_b0); row.Set("t_compose", t_comp);   // compose = b0 + C a (incl. curl)
   row.Set("t_curl_per_nd", isA ? t_curl / sp2->NdofsND() : 0.0);

   // ---------------- totals ----------------
   const double t_remap_steady = b_eval + b_dofs + b_solve + t_gauge + (isA ? t_comp : 0.0);
   const double t_remap_oneshot = t_evsetup + t_plan + t_remap_steady;
   row.Set("t_remap_steady", t_remap_steady);        // evaluate + dof functionals/solves + gauge + compose, plan reused
   row.Set("t_remap_oneshot", t_remap_oneshot);      // + RemoteEvaluator setup + FindPoints + plan
   row.Set("t_total_with_build", t_remap_oneshot + t_mesh_gen + t_deform + t_spaces + t_ops);
   Mem(m_rss, m_rsum, m_hwm, m_hsum);
   row.Set("rss_end_mb", m_rss); row.Set("rss_end_sum_mb", m_rsum);
   row.Set("peak_mb", m_hwm); row.Set("peak_sum_mb", m_hsum);
   row.Set("peak_mb_per_nd_kB", m_hsum * 1024.0 / sp2->NdofsND());

   double errB_rel = 0;
   // ---------------- diagnostics (outside the timings and the memory snapshots) ----------------
   {
      auto Bc = pr.BCoef();
      ParGridFunction bgf(sp2->RT.get());
      bgf.SetFromTrueDofs(b_t);
      const real_t err = L2ErrorVec(bgf, *Bc, q + 3);
      const real_t nrm = L2NormVec(*pm2, *Bc, q + 3);
      errB_rel = err / nrm;
      row.Set("errB_rel", errB_rel);
      Vector dv;
      const DivDiag dd = DivergenceDiag(*ops2, b_t);
      row.Set("div_max", (double)dd.maxabs);
      row.Set("div_before_clean", div_before);
   }
   if (!csv.empty()) { CsvAppend(comm, csv, row); }
   if (rank == 0)
   {
      std::cout << std::setprecision(4) << row.Items()[2].second << " N=" << N << " p=" << p << " np=" << np << " ND=" << sp2->NdofsND()
                << "  build " << t_mesh_gen + t_deform + t_spaces + t_ops << "  evsetup " << t_evsetup << "  plan " << t_plan
                << " (collect " << t_collect << " find " << t_locate << ")  eval " << b_eval << "  dofs " << b_dofs << "  solve " << b_solve
                << " (its " << its << "/" << inner << ")  gauge " << t_gauge << " (" << g_iters << ")  curl " << t_curl
                << "  steady " << t_remap_steady << "  oneshot " << t_remap_oneshot << "  peak " << m_hwm << " MB/rank, " << m_hsum
                << " MB total  errB_rel " << errB_rel << std::endl;
   }
   return 0;
}
