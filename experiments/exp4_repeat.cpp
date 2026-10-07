// exp4_repeat.cpp -- E4: repeated remaps, ping-pong M_A (uniform) <-> M_B (deformed).
// State carried as A (A-route ops) or as the total B (B-route ops).  One CSV row per
// (op, logged step).  Transfer plans (point location, mass operators) are built once.
//
//   mpirun -np 1 bin/exp4_repeat -N 8 -p 2 -field abc -epsB 0.3 -n 100 -every 5 -ops all -csv out.csv
#include "exp_transfer_common.hpp"

using namespace vp;

int main(int argc, char *argv[])
{
   Mpi::Init(argc, argv);
   Hypre::Init();
   MPI_Comm comm = MPI_COMM_WORLD;
   const int rank = Mpi::WorldRank();

   int N = 8, p = 2, q = -1, fvB = 1, fvA = 1, nsteps = 100, every = 5, nq = -1;
   real_t epsA = 0.0, epsB = 0.3;
   std::string field = "abc", opss = "all", csv, tag;
   bool b0on = false, gradfrac = true;
   OptionsParser args(argc, argv);
   args.AddOption(&N, "-N", "--N", "Cells per direction.");
   args.AddOption(&p, "-p", "--order", "ND order p.");
   args.AddOption(&q, "-q", "--geom-order", "Geometric order (default p).");
   args.AddOption(&field, "-field", "--field", "abc | mod | mod2.");
   args.AddOption(&b0on, "-b0", "--b0", "-no-b0", "--no-b0", "Mean field B0.");
   args.AddOption(&epsA, "-epsA", "--epsA", "Deformation of mesh A (0: uniform).");
   args.AddOption(&fvA, "-fvA", "--fvA", "Variant of mesh A.");
   args.AddOption(&epsB, "-epsB", "--epsB", "Deformation of mesh B.");
   args.AddOption(&fvB, "-fvB", "--fvB", "Variant of mesh B.");
   args.AddOption(&nsteps, "-n", "--n", "Number of transfers (A->B->A ... counts each).");
   args.AddOption(&every, "-every", "--every", "Diagnostics every k steps.");
   args.AddOption(&opss, "-ops", "--ops", "comma list or all.");
   args.AddOption(&nq, "-nq", "--nq", "Gauss points per direction for *_int.");
   args.AddOption(&gradfrac, "-gradfrac", "--gradfrac", "-no-gradfrac", "--no-gradfrac", "Gradient fraction of A.");
   args.AddOption(&csv, "-csv", "--csv", "Append rows to this CSV file.");
   args.AddOption(&tag, "-tag", "--tag", "Series label.");
   args.Parse();
   if (!args.Good()) { if (rank == 0) { args.PrintUsage(std::cout); } return 1; }
   if (rank == 0) { args.PrintOptions(std::cout); }
   if (q < 0) { q = p; }
   const double T0 = WallTime();

   MeshSpec sA; sA.N = N; sA.eps = epsA; sA.fv = fvA;
   MeshSpec sB; sB.N = N; sB.eps = epsB; sB.fv = fvB;
   Problem pr(field, b0on, 0.0, 2.0);
   MeshCase MA(comm, sA, p, q), MB(comm, sB, p, q);
   MeshState SA(MA, pr), SB(MB, pr);
   MeshState *S[2] = {&SA, &SB};
   MeshCase *M[2] = {&MA, &MB};

   // exact-projection references on both meshes (A-route state; b = B0 + C a)
   Vector aref[2], bref[2];
   StateDiag dref[2];
   for (int m = 0; m < 2; m++)
   {
      S[m]->ProjectExactA(pr, false, aref[m]);
      S[m]->Compose(aref[m], bref[m]);
      dref[m] = S[m]->diag->Evaluate(&aref[m], bref[m], gradfrac);
   }
   if (rank == 0)
   {
      std::cout << "N=" << N << " p=" << p << " q=" << q << "  A: " << sA.Str() << "  B: " << sB.Str()
                << "  ND=" << MA.sp->NdofsND() << " RT=" << MA.sp->NdofsRT() << "  n=" << nsteps << std::endl;
   }

   TransferOptions topt;
   topt.nq = nq;
   RemoteEvaluator evA(*MA.pm), evB(*MB.pm);   // evaluators of the SOURCE meshes
   RemoteEvaluator *ev[2] = {&evA, &evB};

   for (TransferKind kind : ParseKindList(opss))
   {
      const double Tk = WallTime();
      const bool isA = IsARoute(kind);
      // plan[m]: source mesh m -> target mesh 1-m
      std::unique_ptr<Transfer> plan[2];
      double t_setup_total = 0.0;
      for (int m = 0; m < 2; m++)
      {
         MeshCase &dst = *M[1 - m];
         plan[m].reset(new Transfer(kind, *ev[m], isA ? *dst.sp->ND : *dst.sp->RT, dst.ops.get(), topt));
         t_setup_total += plan[m]->Stats().t_setup;
      }
      // state on mesh A at step 0
      Vector a_t = aref[0], b_t = bref[0];
      if (!isA)
      {
         // B-route: initial B = Pi_RT^int(exact B) (identical to B0 + C a_ref to roundoff)
         b_t = bref[0];
      }
      ParGridFunction src_gf[2] = {ParGridFunction(isA ? MA.sp->ND.get() : MA.sp->RT.get()),
                                   ParGridFunction(isA ? MB.sp->ND.get() : MB.sp->RT.get())};
      ParGridFunction dst_gf[2] = {ParGridFunction(isA ? MA.sp->ND.get() : MA.sp->RT.get()),
                                   ParGridFunction(isA ? MB.sp->ND.get() : MB.sp->RT.get())};
      // first discrete state seen on each mesh (reference for ||B_n - B_first||)
      Vector bfirst[2];
      bfirst[0] = b_t;
      StateDiag d0;
      double t_cum = t_setup_total, t_step = 0.0, t_diag_cum = 0.0;
      int mesh = 0;           // 0: A, 1: B
      long long it_total = 0, inner_total = 0;
      int it_last = 0, inner_last = 0;
      long long nnf = plan[0]->Stats().nnotfound + plan[1]->Stats().nnotfound;
      real_t E0 = 0, H0 = 0, a0 = 0;

      auto log = [&](int step)
      {
         const double td = WallTime();
         const StateDiag d = S[mesh]->diag->Evaluate(isA ? &a_t : nullptr, b_t, gradfrac);
         if (step == 0) { d0 = d; E0 = d.energy; H0 = d.helicity; a0 = d.a_L2; }
         real_t dfirst = std::nan("");
         if (!bfirst[mesh].Size()) { bfirst[mesh] = b_t; }
         {
            Vector zero(b_t.Size());
            zero = 0.0;
            dfirst = S[mesh]->diag->RTDistance(b_t, bfirst[mesh]) / S[mesh]->diag->RTDistance(bfirst[mesh], zero);
         }
         const StateDiag &r = dref[mesh];
         CsvRow row;
         row.Set("exp", "exp4_repeat"); row.Set("tag", tag);
         row.Set("op", KindName(kind)); row.Set("route", isA ? "A" : "B");
         row.Set("p", p); row.Set("q", q); row.Set("N", N); row.Set("field", field); row.Set("b0", (int)b0on);
         row.Set("epsA", (double)epsA); row.Set("epsB", (double)epsB); row.Set("fvB", fvB);
         row.Set("nsteps", nsteps); row.Set("nq", (kind == TransferKind::A_int || kind == TransferKind::B_int)
                                                  ? (nq > 0 ? nq : DefaultTransferNq(p)) : 0);
         row.Set("np", Mpi::WorldSize());
         row.Set("step", step); row.Set("mesh", mesh == 0 ? "A" : "B");
         SetDiag(row, "", d);
         row.Set("E_over_E0", (double)(d.energy / E0));
         row.Set("E_over_Eref", (double)(d.energy / r.energy));
         row.Set("dB_first", (double)dfirst);
         row.Set("errB_over_ref", (double)(d.errB_rel / r.errB_rel));
         row.Set("H_over_H0", isA ? (double)(d.helicity / H0) : std::nan(""));
         row.Set("H_over_Href", isA ? (double)(d.helicity / r.helicity) : std::nan(""));
         row.Set("H_dev_scaled", isA ? (double)((d.helicity - r.helicity) / r.hel_scale) : std::nan(""));
         row.Set("a_over_aref", isA ? (double)(d.a_L2 / r.a_L2) : std::nan(""));
         row.Set("a_over_a0", isA ? (double)(d.a_L2 / a0) : std::nan(""));
         row.Set("t_setup", t_setup_total);
         row.Set("t_cum", t_cum - t_setup_total);
         row.Set("t_step_last", t_step);
         row.Set("iters_last", it_last); row.Set("inner_last", inner_last);
         row.Set("iters_total", it_total); row.Set("inner_total", inner_total);
         row.Set("nnotfound", nnf);
         if (!csv.empty()) { CsvAppend(comm, csv, row); }
         t_diag_cum += WallTime() - td;
         return d;
      };

      log(0);
      for (int step = 1; step <= nsteps; step++)
      {
         const double ts = WallTime();
         const int src = mesh, dstm = 1 - mesh;
         if (isA)
         {
            src_gf[src].SetFromTrueDofs(a_t);
            plan[src]->Apply(src_gf[src], dst_gf[dstm]);
            dst_gf[dstm].GetTrueDofs(a_t);
            S[dstm]->Compose(a_t, b_t);
         }
         else
         {
            src_gf[src].SetFromTrueDofs(b_t);
            plan[src]->Apply(src_gf[src], dst_gf[dstm]);
            dst_gf[dstm].GetTrueDofs(b_t);
         }
         mesh = dstm;
         t_step = WallTime() - ts;
         t_cum += t_step;
         const auto &st = plan[src]->Stats();
         it_last = st.iters; inner_last = st.inner_iters;
         it_total += st.iters; inner_total += st.inner_iters;
         if (!bfirst[mesh].Size()) { bfirst[mesh] = b_t; }
         if (step % every == 0 || step == nsteps || step <= 2)
         {
            const StateDiag d = log(step);
            if (rank == 0 && (step % (10 * every) == 0 || step == nsteps || step == 1))
            {
               std::cout << std::scientific << std::setprecision(3) << KindName(kind) << " step " << std::setw(4)
                         << step << " " << (mesh ? "B" : "A") << "  E/E0=" << d.energy / E0
                         << " errB=" << d.errB_rel << " div_rel=" << d.div_rel_l2
                         << " |a|=" << d.a_L2 << " H/H0=" << (isA ? d.helicity / H0 : 0.0)
                         << " t_cum=" << t_cum - t_setup_total << std::defaultfloat << std::endl;
            }
         }
      }
      if (rank == 0)
      {
         std::cout << KindName(kind) << " done: setup " << t_setup_total << " s, transfers " << t_cum - t_setup_total
                   << " s (" << (t_cum - t_setup_total) / nsteps << " s/remap), diagnostics " << t_diag_cum
                   << " s, wall " << WallTime() - Tk << " s, not found " << nnf << std::endl;
      }
   }
   if (rank == 0) { std::cout << "total wall " << WallTime() - T0 << " s" << std::endl; }
   return 0;
}
