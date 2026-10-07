// exp5_ale.cpp -- E5: smooth ALE cycle with rezones on the periodic torus.
//
// Material motion Phi_t(X) = X + eps s(t) f(X), s = sin^2(pi t/T), f = vp::DisplacementField(variant).
// The field is frozen in: the true-dof vectors a (A-route) / b (B-route) are CONSTANT while the nodes
// move (Lagrangian phase); at prescribed rezone times the state is transferred from the deformed
// mesh to a fresh uniform mesh with a vp transfer operator, and the nodes then continue to move from
// the material labels X_u = Phi_{t_r}^{-1}(y).  The exact field is the Cauchy push-forward of the
// initial field, evaluated at physical points by a Newton inverse of Phi_t.
//
//   mpirun -np 1 bin/exp5_ale -N 8 -p 2 -q 2 -eps 0.3 -K 100 -rezones 5 -op A_int -gauge coulomb -csv out.csv
//   mpirun -np 1 bin/exp5_ale -selftest
#include "exp5_common.hpp"
#include "vp_gauge.hpp"

#include <random>

using namespace vp;
using namespace vp::e5;

namespace
{

struct Cfg
{
   int N = 8, p = 2, q = -1, K = 100, nrez = 0, variant = 1, every = 1;
   real_t eps = 0.3, T = 1.0;
   std::string op = "B_l2c", gauge = "none", field = "abc", tag, csv;
   bool b0 = false;
   bool verbose = true;
};

struct StepRec
{
   int step = 0, rezone = 0;
   double t = 0, s = 0, minDetJ = 0, errB = 0, errA = NAN, E_B = 0, E_ex = 0, H = NAN, flux_err = 0,
          div_scaled = 0, div_max = 0, normB = 0;
};

struct Result
{
   std::vector<StepRec> rec;
   Vector a0, b0v, aT, bT;
   double final_dB = NAN, final_E_ratio = NAN, final_H_ratio = NAN, final_flux_err = NAN,
          final_div = NAN, final_errB = NAN, final_errA = NAN, t_transfer_total = 0, a_change_inf = NAN,
          b_change_inf = NAN, max_label_resid = 0, max_dc = 0;
   long long nnotfound = 0;
};

/// One mesh of the ping-pong pair with everything attached to it.
struct Lag
{
   std::unique_ptr<MeshCase> mc;
   Vector X0, lab;
   std::unique_ptr<CoulombGauge> cg;
   std::unique_ptr<ParGridFunction> a_gf, b_gf;
   Lag(MPI_Comm comm, int N, int p, int q)
   {
      MeshSpec ms;
      ms.N = N;
      ms.eps = 0.0;
      mc.reset(new MeshCase(comm, ms, p, q));
      X0 = *mc->pm->GetNodes();
      lab = X0;
      a_gf.reset(new ParGridFunction(mc->sp->ND.get()));
      b_gf.reset(new ParGridFunction(mc->sp->RT.get()));
   }
   ParMesh &pm() { return *mc->pm; }
   const CoulombGauge &Gauge()
   {
      if (!cg) { cg.reset(new CoulombGauge(*mc->sp, *mc->ops)); }
      return *cg;
   }
};

struct Meas
{
   double errB = 0, normB = 0, errA = NAN, normA = NAN, E_B = 0, H = NAN, aa = NAN, div2 = 0;
};

/// One quadrature sweep over the moving mesh: errors vs the push-forward, energy, helicity.
Meas Sweep(Lag &m, const PushForward &pf, const Vector *a_t, const Vector &b_t, int order)
{
   ParMesh &pm = m.pm();
   m.b_gf->SetFromTrueDofs(b_t);
   // physical div b_h at quadrature points from the L2 point values of the REFERENCE (uniform-mesh) D_h:
   // value = div_ref/detJ_0 with detJ_0 = N^-3, so div_phys = value * detJ_0/detJ_current.
   Vector dvec(m.mc->ops->D->Height());
   m.mc->ops->D->Mult(b_t, dvec);
   ParGridFunction dgf(m.mc->sp->L2.get());
   dgf.SetFromTrueDofs(dvec);
   const double detJ0 = 1.0 / ((double)pm.GetNE());
   const bool wantA = (a_t != nullptr);
   if (wantA) { m.a_gf->SetFromTrueDofs(*a_t); }
   const IntegrationRule &ir = IntRules.Get(Geometry::CUBE, order);
   DenseMatrix bv, av;
   Meas r;
   double s[8] = {0, 0, 0, 0, 0, 0, 0, 0};   // errB2, nB2, errA2, nA2, EB, H, aa
   real_t x[3], B[3], A[3];
   Vector xv(3);
   for (int e = 0; e < pm.GetNE(); e++)
   {
      ElementTransformation *T = pm.GetElementTransformation(e);
      m.b_gf->GetVectorValues(*T, ir, bv);
      if (wantA) { m.a_gf->GetVectorValues(*T, ir, av); }
      for (int k = 0; k < ir.GetNPoints(); k++)
      {
         const IntegrationPoint &ip = ir.IntPoint(k);
         T->SetIntPoint(&ip);
         T->Transform(ip, xv);
         for (int d = 0; d < 3; d++) { x[d] = xv(d); }
         pf.Eval(x, B, wantA ? A : nullptr);
         const double w = ip.weight * T->Weight();
         double e2 = 0, n2 = 0, bb = 0;
         for (int d = 0; d < 3; d++)
         {
            const double db = bv(d, k) - B[d];
            e2 += db * db; n2 += B[d] * B[d]; bb += bv(d, k) * bv(d, k);
         }
         s[0] += w * e2; s[1] += w * n2; s[4] += 0.5 * w * bb;
         const double dv = dgf.GetValue(e, ip) * detJ0 / T->Weight();
         s[7] += w * dv * dv;
         if (wantA)
         {
            double ea = 0, na = 0, ah = 0, aa = 0;
            for (int d = 0; d < 3; d++)
            {
               const double da = av(d, k) - A[d];
               ea += da * da; na += A[d] * A[d]; ah += av(d, k) * bv(d, k); aa += av(d, k) * av(d, k);
            }
            s[2] += w * ea; s[3] += w * na; s[5] += w * ah; s[6] += w * aa;
         }
      }
   }
   MPI_Allreduce(MPI_IN_PLACE, s, 8, MPI_DOUBLE, MPI_SUM, pm.GetComm());
   r.errB = std::sqrt(s[0]); r.normB = std::sqrt(s[1]); r.E_B = s[4]; r.div2 = s[7];
   if (wantA) { r.errA = std::sqrt(s[2]); r.normA = std::sqrt(s[3]); r.H = s[5]; r.aa = std::sqrt(s[6]); }
   return r;
}

double SliceFluxErr(Lag &m, const ParGridFunction &bgf, const Vector &B0, double flux[3])
{
   double fmax = 0.0;
   for (int d = 0; d < 3; d++) { flux[d] = NAN; }
   for (const auto &fam : m.mc->slices)
   {
      const double fl = SliceFlux(bgf, fam);
      fmax = std::max(fmax, std::fabs(fl - B0(fam.dir)));
      if (fam.j == 0) { flux[fam.dir] = fl; }
   }
   return fmax;
}

const char *kCols[] = {"kind", "step", "t", "s", "rezone"};

void Run(const Cfg &cfg0, Result &res, MPI_Comm comm)
{
   Cfg cfg = cfg0;
   const int rank = Mpi::WorldRank();
   if (cfg.q < 0) { cfg.q = cfg.p; }
   const int p = cfg.p, q = cfg.q, N = cfg.N, K = cfg.K;
   const TransferKind kind = ParseKind(cfg.op);
   const bool isA = IsARoute(kind);
   // coulomb: Coulomb projection of the transferred a on the new uniform mesh (as specified);
   // coulomb_pre (extra): Coulomb projection of a on the DEFORMED source mesh right before the transfer
   // (a CoulombGauge is built for the current geometry each time)
   const bool coulomb = (cfg.gauge == "coulomb");
   const bool coulomb_pre = (cfg.gauge == "coulomb_pre");
   MFEM_VERIFY(cfg.gauge == "none" || coulomb || coulomb_pre, "-gauge must be none|coulomb|coulomb_pre");
   MFEM_VERIFY(isA || cfg.gauge == "none", "gauge fixing only applies to A-route operators");
   MFEM_VERIFY(!(isA && cfg.b0 && cfg.nrez > 0),
               "A-route with a mean field B0 and rezones is not implemented: the deformed mean field has a "
               "fluctuating part whose potential would need a curl inversion (see results/exp5_summary.md)");
   MFEM_VERIFY(K > 0 && K % std::max(1, cfg.nrez) == 0, "K must be a positive multiple of #rezones");
   const double Tw0 = WallTime();
   const real_t dt = cfg.T / K;

   Problem pr(cfg.field, cfg.b0, 0.0, 2.0);
   PushForward pf;
   pf.F = pr.F;
   pf.B0 = pr.B0eff;
   pf.variant = cfg.variant;
   // quadrature order of all diagnostics (cheaper than 2p+2q; changes errB by < 1e-3 relative, see summary)
   const int qorder = 2 * p + q + 1;

   Lag M0(comm, N, p, q), M1(comm, N, p, q);
   Lag *M[2] = {&M0, &M1};
   int cur = 0;

   // initial state on the uniform mesh 0: a = Pi_ND^int A_exact, b = b0 + C a
   Vector b0_t;
   {
      ParGridFunction b0(M0.mc->sp->RT.get());
      VectorConstantCoefficient B0c(pr.B0eff);
      ProjectB(b0, B0c, ProjMode::Integrated, q + 2);
      b0.GetTrueDofs(b0_t);
   }
   Vector a_t, b_t;
   {
      auto Ac = pr.ACoef(false);
      ParGridFunction a(M0.mc->sp->ND.get());
      ProjectA(a, *Ac, ProjMode::Integrated, kExactNq);
      a.GetTrueDofs(a_t);
      ComposeB(*M0.mc->ops, b0_t, a_t, b_t);
   }
   res.a0 = a_t;
   res.b0v = b_t;
   const Vector a_init = a_t, b_init = b_t;

   TransferOptions topt;
   std::vector<int> rez_steps;
   for (int r = 1; r <= cfg.nrez; r++) { rez_steps.push_back((int)std::llround((double)r * K / cfg.nrez)); }
   auto is_rez = [&](int k) { return std::find(rez_steps.begin(), rez_steps.end(), k) != rez_steps.end(); };

   double E0 = NAN, H0 = NAN, t_cum_transfer = 0.0;
   long long nnf_total = 0;
   const double kN = N;   // 1/h

   auto make_row = [&](const std::string &kindname, int step, int rezone, const StepRec &r, double t_rez,
                       double t_setup, double t_apply, double t_gauge, int iters, int inner,
                       const GaugeFixStats *gs, const double *flux) -> CsvRow
   {
      CsvRow row;
      row.Set("exp", "exp5_ale"); row.Set("tag", cfg.tag);
      row.Set("p", p); row.Set("q", q); row.Set("N", N); row.Set("eps", (double)cfg.eps);
      row.Set("K", K); row.Set("nrez", cfg.nrez);
      row.Set("op", cfg.op); row.Set("route", isA ? "A" : "B"); row.Set("gauge", cfg.gauge);
      row.Set("b0", (int)cfg.b0); row.Set("field", cfg.field); row.Set("variant", cfg.variant);
      row.Set("kind", kindname); row.Set("step", step); row.Set("t", r.t); row.Set("s", r.s);
      row.Set("rezone", rezone);
      row.Set("minDetJ", r.minDetJ);
      row.Set("errB_rel", r.errB); row.Set("errA_rel", r.errA);
      row.Set("normB_ex", r.normB);
      row.Set("E_B", r.E_B); row.Set("E_ex", r.E_ex);
      row.Set("E_err_rel", (r.E_B - r.E_ex) / r.E_ex);
      row.Set("E_over_E0", r.E_B / E0);
      row.Set("E_over_Eex", r.E_B / r.E_ex);
      row.Set("H", r.H); row.Set("H_over_H0", r.H / H0);
      row.Set("flux_err", r.flux_err);
      row.Set("flux_x", flux ? flux[0] : (double)NAN);
      row.Set("flux_y", flux ? flux[1] : (double)NAN);
      row.Set("flux_z", flux ? flux[2] : (double)NAN);
      row.Set("div_scaled", r.div_scaled); row.Set("div_max", r.div_max);
      row.Set("t_rezone", t_rez); row.Set("t_setup", t_setup); row.Set("t_apply", t_apply);
      row.Set("t_gauge", t_gauge);
      row.Set("t_cum_transfer", t_cum_transfer); row.Set("t_cum_wall", WallTime() - Tw0);
      row.Set("iters", iters); row.Set("inner_iters", inner);
      row.Set("gauge_iters", gs ? gs->iters : -1);
      row.Set("gauge_resid_before", gs ? (double)gs->resid_before : (double)NAN);
      row.Set("gauge_resid_after", gs ? (double)gs->resid_after : (double)NAN);
      row.Set("gauge_grad_before", gs ? (double)gs->grad_fraction_before : (double)NAN);
      row.Set("gauge_grad_after", gs ? (double)gs->grad_fraction_after : (double)NAN);
      row.Set("nnotfound", nnf_total);
      for (const char *c : {"final_dB", "final_E_ratio", "final_H_ratio", "final_flux_err", "final_div",
                            "final_errB", "a_change_inf", "t_transfer_total"})
      {
         row.Set(c, (double)NAN);
      }
      return row;
   };

   // measure the state on the current mesh (time t, nodes already at Phi_t(labels))
   auto measure = [&](int step, int k_rez, double t, double s) -> std::pair<StepRec, double *>
   {
      static double flux[3];
      Lag &m = *M[cur];
      pf.c = cfg.eps * s;
      StepRec r;
      r.step = step; r.rezone = k_rez; r.t = t; r.s = s;
      r.minDetJ = JacobianStats(m.pm()).min_det;
      const Meas ms = Sweep(m, pf, isA ? &a_t : nullptr, b_t, qorder);
      r.errB = ms.errB / ms.normB;
      r.errA = isA ? ms.errA / ms.normA : (double)NAN;
      r.normB = ms.normB;
      r.E_B = ms.E_B;
      r.E_ex = 0.5 * ms.normB * ms.normB;
      r.H = (isA && !cfg.b0) ? ms.H : (double)NAN;
      r.flux_err = SliceFluxErr(m, *m.b_gf, pr.B0eff, flux);
      // max |D_h b| over dofs with D_h of the (uniform-mesh) reference geometry (= reference divergence x N^3)
      Vector dvec(m.mc->ops->D->Height());
      m.mc->ops->D->Mult(b_t, dvec);
      double dmax = dvec.Normlinf();
      MPI_Allreduce(MPI_IN_PLACE, &dmax, 1, MPI_DOUBLE, MPI_MAX, comm);
      r.div_scaled = std::sqrt(ms.div2) / (ms.normB * kN);   // pointwise div b_h on the CURRENT mesh
      r.div_max = dmax;
      return {r, flux};
   };

   auto record = [&](const char *kindname, const StepRec &r, const double *flux, int rezone, double t_rez,
                     double t_setup, double t_apply, double t_gauge, int iters, int inner,
                     const GaugeFixStats *gs)
   {
      res.rec.push_back(r);
      if (!cfg.csv.empty())
      {
         CsvAppend(comm, cfg.csv, make_row(kindname, r.step, rezone, r, t_rez, t_setup, t_apply, t_gauge,
                                           iters, inner, gs, flux));
      }
   };

   // t = 0
   {
      auto mr = measure(0, 0, 0.0, 0.0);
      E0 = mr.first.E_B;
      H0 = mr.first.H;
      record("step", mr.first, mr.second, 0, 0, 0, 0, 0, 0, 0, nullptr);
      if (rank == 0 && cfg.verbose)
      {
         std::cout << "t=0: errB=" << mr.first.errB << " errA=" << mr.first.errA << " E_B=" << E0 << " E_ex="
                   << mr.first.E_ex << " H=" << H0 << " minDetJ=" << mr.first.minDetJ << std::endl;
      }
   }

   for (int k = 1; k <= K; k++)
   {
      const double t = (k == K) ? cfg.T : k * dt;
      const double s = (k == K) ? 0.0 : SProfile(t, cfg.T);
      Lag &m = *M[cur];
      MoveNodes(m.pm(), m.lab, cfg.eps * s, cfg.variant);   // Lagrangian step: only the nodes move
      const bool rz = is_rez(k);
      if (k % cfg.every == 0 || rz || k == K)
      {
         auto mr = measure(k, 0, t, s);
         record("step", mr.first, mr.second, 0, 0, 0, 0, 0, 0, 0, nullptr);
         if (rank == 0 && cfg.verbose && (k % 10 == 0 || rz))
         {
            std::cout << "k=" << k << " t=" << t << " errB=" << mr.first.errB << " E/Eex=" << mr.first.E_B / mr.first.E_ex
                      << " H/H0=" << mr.first.H / H0 << " div=" << mr.first.div_scaled << " minDetJ="
                      << mr.first.minDetJ << std::endl;
         }
      }
      if (rz)
      {
         const double tr0 = WallTime();
         Lag &S = *M[cur];
         Lag &D = *M[1 - cur];
         // target: fresh uniform mesh with material labels X_u = Phi_t^{-1}(y)
         SetNodes(D.pm(), D.X0);
         double mi = 0;
         const real_t lres = ComputeLabels(D.pm(), cfg.variant, cfg.eps * s, D.X0, D.lab, &mi);
         res.max_label_resid = std::max(res.max_label_resid, (double)lres);
         MFEM_VERIFY(lres < 1e-13, "label Newton residual " << lres);
         GaugeFixStats gs;
         double t_gauge = 0.0;
         if (coulomb_pre)
         {
            const double tg0 = WallTime();
            CoulombGauge cgs(*S.mc->sp, *S.mc->ops);   // M_ND of the current (deformed) geometry
            gs = cgs.Apply(a_t, 1e-12, true);
            t_gauge = WallTime() - tg0 - 0.0;
         }
         RemoteEvaluator ev(S.pm());
         Transfer plan(kind, ev, isA ? *D.mc->sp->ND : *D.mc->sp->RT, D.mc->ops.get(), topt);
         ParGridFunction &src = isA ? *S.a_gf : *S.b_gf;
         ParGridFunction &dst = isA ? *D.a_gf : *D.b_gf;
         src.SetFromTrueDofs(isA ? a_t : b_t);
         plan.Apply(src, dst);
         dst.GetTrueDofs(isA ? a_t : b_t);
         if (coulomb)
         {
            gs = D.Gauge().Apply(a_t, 1e-12, true);
            t_gauge = gs.time;
         }
         if (isA) { ComposeB(*D.mc->ops, b0_t, a_t, b_t); }
         const TransferStats &st = plan.Stats();
         nnf_total += st.nnotfound;
         cur = 1 - cur;
         const double t_rez = WallTime() - tr0;
         t_cum_transfer += t_rez;
         res.t_transfer_total = t_cum_transfer;
         auto mr = measure(k, 1, t, s);
         record("step", mr.first, mr.second, 1, t_rez, st.t_setup, st.t_apply, t_gauge, st.iters,
                st.inner_iters, (coulomb || coulomb_pre) ? &gs : nullptr);
         if (rank == 0 && cfg.verbose)
         {
            std::cout << "  rezone at k=" << k << " (t=" << t << "): " << t_rez << " s (setup " << st.t_setup
                      << ", apply " << st.t_apply << ", gauge " << t_gauge << "), errB " << mr.first.errB
                      << ", div " << mr.first.div_scaled << ", label Newton resid " << lres << " (mean its "
                      << mi << ")" << std::endl;
         }
      }
   }

   // final comparison on the uniform mesh at T
   {
      Lag &m = *M[cur];
      const StepRec &r = res.rec.back();
      res.aT = a_t;
      res.bT = b_t;
      Vector d = b_t;
      d -= b_init;
      ParGridFunction g(m.mc->sp->RT.get());
      g.SetFromTrueDofs(d);
      const int ord = qorder;
      const double nd = std::sqrt(InnerGF(g, g, ord));
      g.SetFromTrueDofs(b_init);
      const double n0 = std::sqrt(InnerGF(g, g, ord));
      res.final_dB = nd / n0;
      res.final_E_ratio = r.E_B / E0;
      res.final_H_ratio = r.H / H0;
      res.final_flux_err = r.flux_err;
      res.final_div = r.div_scaled;
      res.final_errB = r.errB;
      res.final_errA = r.errA;
      res.nnotfound = nnf_total;
      {
         Vector da = a_t; da -= a_init;
         Vector db = b_t; db -= b_init;
         res.a_change_inf = da.Normlinf();
         res.b_change_inf = db.Normlinf();
      }
      StepRec fr = r;
      if (!cfg.csv.empty())
      {
         CsvRow row = make_row("final", K, 0, fr, 0, 0, 0, 0, 0, 0, nullptr, nullptr);
         row.Set("final_dB", res.final_dB);
         row.Set("final_E_ratio", res.final_E_ratio);
         row.Set("final_H_ratio", res.final_H_ratio);
         row.Set("final_flux_err", res.final_flux_err);
         row.Set("final_div", res.final_div);
         row.Set("final_errB", res.final_errB);
         row.Set("a_change_inf", res.a_change_inf);
         row.Set("t_transfer_total", res.t_transfer_total);
         CsvAppend(comm, cfg.csv, row);
      }
   }
   (void)kCols;
}


/// Projection floors at the rezone times (independent of the transfer operator): at every rezone time
/// t_r the relative error of the integrated-dof projection Pi_RT^int B_pf(t_r) of the exact push-forward
/// onto (a) the deformed source mesh just before the rezone and (b) the new uniform mesh.  This is the
/// best any transfer could give on that mesh ("floor" of E3).  One CSV row per rezone.
void RunFloors(const Cfg &cfg0, MPI_Comm comm)
{
   Cfg cfg = cfg0;
   if (cfg.q < 0) { cfg.q = cfg.p; }
   const int p = cfg.p, q = cfg.q, N = cfg.N, K = cfg.K;
   Problem pr(cfg.field, cfg.b0, 0.0, 2.0);
   PushForward pf;
   pf.F = pr.F; pf.B0 = pr.B0eff; pf.variant = cfg.variant;
   const int qorder = 2 * p + q + 1;
   Lag M0(comm, N, p, q), M1(comm, N, p, q);
   Lag *M[2] = {&M0, &M1};
   int cur = 0;
   auto floor_on = [&](Lag &m, real_t c)
   {
      pf.c = c;
      PushForwardCoef cf(pf, false);
      ParGridFunction b(m.mc->sp->RT.get());
      ProjectB(b, cf, ProjMode::Integrated, kExactNq);
      Vector bt;
      b.GetTrueDofs(bt);
      const Meas ms = Sweep(m, pf, nullptr, bt, qorder);
      return ms.errB / ms.normB;
   };
   for (int r = 1; r <= cfg.nrez; r++)
   {
      const int k = (int)std::llround((double)r * K / cfg.nrez);
      const double t = (k == K) ? cfg.T : k * (cfg.T / K);
      const double s = (k == K) ? 0.0 : SProfile(t, cfg.T);
      Lag &S = *M[cur];
      Lag &D = *M[1 - cur];
      MoveNodes(S.pm(), S.lab, cfg.eps * s, cfg.variant);
      const double fpre = floor_on(S, cfg.eps * s);
      SetNodes(D.pm(), D.X0);
      ComputeLabels(D.pm(), cfg.variant, cfg.eps * s, D.X0, D.lab);
      const double fpost = floor_on(D, cfg.eps * s);
      cur = 1 - cur;
      CsvRow row;
      row.Set("exp", "exp5_floor"); row.Set("p", p); row.Set("q", q); row.Set("N", N);
      row.Set("eps", (double)cfg.eps); row.Set("nrez", cfg.nrez); row.Set("rezone_index", r);
      row.Set("step", k); row.Set("t", t); row.Set("s", s);
      row.Set("minDetJ_src", JacobianStats(S.pm()).min_det);
      row.Set("floor_pre", fpre); row.Set("floor_post", fpost);
      if (!cfg.csv.empty()) { CsvAppend(comm, cfg.csv, row); }
      if (Mpi::WorldRank() == 0)
      {
         std::cout << "rezone " << r << " t=" << t << " floor on deformed source " << fpre << ", on uniform target " << fpost << std::endl;
      }
   }
}

} // namespace

// ---------------------------------------------------------------------------------------------

int SelfTest(MPI_Comm comm, int N, int p, int q, int nsteps)
{
   int fails = 0;
   const int rank = Mpi::WorldRank();
   auto check = [&](bool ok, const std::string &what)
   {
      if (rank == 0) { std::cout << (ok ? "[PASS] " : "[FAIL] ") << what << std::endl; }
      if (!ok) { fails++; }
   };
   // (i) gradient of f by central finite differences
   {
      std::mt19937 gen(5);
      std::uniform_real_distribution<real_t> U(0.0, 1.0);
      double worst = 0.0;
      for (int variant = 1; variant <= 3; variant++)
      {
         for (int n = 0; n < 200; n++)
         {
            real_t x[3] = {U(gen), U(gen), U(gen)}, G[9];
            GradF(variant, x, G);
            const real_t h = 1e-6;
            for (int j = 0; j < 3; j++)
            {
               real_t xp[3] = {x[0], x[1], x[2]}, xm[3] = {x[0], x[1], x[2]}, fp[3], fm[3];
               xp[j] += h; xm[j] -= h;
               DisplacementField(variant, xp, fp);
               DisplacementField(variant, xm, fm);
               for (int i = 0; i < 3; i++)
               {
                  worst = std::max(worst, (double)std::fabs((fp[i] - fm[i]) / (2 * h) - G[3 * i + j]));
               }
            }
         }
      }
      if (rank == 0) { std::cout << "  (i) max |FD - analytic grad f| = " << worst << " (variants 1,2,3)" << std::endl; }
      check(worst < 1e-8, "(i) grad f finite-difference check");
   }
   // (ii) Newton inverse, eps = 0.6, s = 1 (and 0.3)
   {
      std::mt19937 gen(6);
      std::uniform_real_distribution<real_t> U(-0.1, 1.1);
      for (real_t eps : {0.3, 0.6})
      {
         double worst = 0.0, mean_it = 0, worst_it = 0;
         const int n = 20000;
         for (int i = 0; i < n; i++)
         {
            real_t y[3] = {U(gen), U(gen), U(gen)}, X[3], res;
            const int it = InverseMap(1, eps * 1.0, y, X, 1e-14, &res);
            real_t r[3];
            worst = std::max(worst, (double)PhiResidual(1, eps, X, y, r));
            mean_it += it; worst_it = std::max(worst_it, (double)it);
         }
         if (rank == 0)
         {
            std::cout << "  (ii) eps=" << eps << " s=1: max residual " << worst << ", mean Newton its " << mean_it / n
                      << ", max " << worst_it << std::endl;
         }
         check(worst < 1e-13, "(ii) Newton inverse residual < 1e-13, eps=" + std::to_string(eps));
      }
   }
   // (iii) K=0 Lagrangian run
   for (int b0 = 0; b0 < 2; b0++)
   {
      Cfg cfg;
      cfg.N = N; cfg.p = p; cfg.q = q; cfg.K = nsteps; cfg.nrez = 0; cfg.eps = 0.6; cfg.op = "A_int";
      cfg.b0 = (b0 == 1);
      cfg.verbose = false;
      Result res;
      Run(cfg, res, comm);
      const StepRec &r0 = res.rec.front(), &rT = res.rec.back();
      double dH = 0, dFlux = 0, maxdiv = 0, maxrelerrE = 0, minDet = 1e300, peakErrB = 0, peakStep = 0, peakE = 0;
      const double relE0 = std::fabs(r0.E_B - r0.E_ex) / r0.E_ex;
      for (const StepRec &r : res.rec)
      {
         if (!cfg.b0) { dH = std::max(dH, std::fabs(r.H / r0.H - 1.0)); }
         dFlux = std::max(dFlux, r.flux_err);
         maxdiv = std::max(maxdiv, r.div_scaled);
         peakE = std::max(peakE, std::fabs(r.E_B - r.E_ex) / r.E_ex);
         maxrelerrE = std::max(maxrelerrE, std::fabs(r.E_B - r.E_ex) / r.E_ex / relE0);
         minDet = std::min(minDet, r.minDetJ);
         if (r.errB > peakErrB) { peakErrB = r.errB; peakStep = r.step; }
      }
      const std::string lab = std::string("K=0 eps=0.6 p=") + std::to_string(p) + (cfg.b0 ? " b0 on" : " b0 off");
      if (rank == 0)
      {
         std::cout << "  (iii) " << lab << ": ||a_T-a_0||_inf=" << res.a_change_inf << " ||b_T-b_0||_inf=" << res.b_change_inf
                   << "\n        errB(0)=" << r0.errB << " errB(T)=" << rT.errB << " peak errB=" << peakErrB << " at step "
                   << peakStep << " (ratio " << peakErrB / r0.errB << ")"
                   << "\n        max|H/H0-1|=" << dH << " max flux err=" << dFlux << " max scaled div=" << maxdiv
                   << " min detJ=" << minDet
                   << "\n        E err rel: t=0 " << relE0 << ", worst over t " << peakE << " (x" << maxrelerrE << " of t=0)"
                   << std::endl;
      }
      check(res.a_change_inf == 0.0 && res.b_change_inf == 0.0, "(iii) dofs constant (inf-norm change == 0), " + lab);
      check(std::fabs(rT.errB - r0.errB) <= 1e-14 * r0.errB + 1e-300, "(iii) final B error equals initial to 1e-14, " + lab);
      check(std::fabs(rT.E_B - r0.E_B) <= 1e-14 * r0.E_B, "(iii) final energy equals initial, " + lab);
      if (!cfg.b0) { check(dH < 1e-12, "(iii) helicity constant to roundoff, " + lab); }
      check(dFlux < 1e-12, "(iii) slice fluxes constant to roundoff, " + lab);
      check(maxdiv < 1e-12, "(iii) scaled div at roundoff at every step, " + lab);
      // energy tracking: informational (printed above)
   }
   // (iv) geometry independence of C_h, D_h, G_h after moving; D_h C_h = 0; label interpolation vs Newton
   {
      Lag M0(comm, N, p, q);
      MoveNodes(M0.pm(), M0.lab, 0.6 * 1.0, 1);
      Operators ops2(*M0.mc->sp);
      auto diff = [&](const HypreParMatrix &A, const HypreParMatrix &B)
      {
         std::unique_ptr<HypreParMatrix> Bm(new HypreParMatrix(B));
         *Bm *= -1.0;
         std::unique_ptr<HypreParMatrix> S(ParAdd(&A, Bm.get()));
         return MaxAbsEntry(*S);
      };
      const double dC = diff(*M0.mc->ops->C, *ops2.C), dG = diff(*M0.mc->ops->G, *ops2.G);
      // D_h maps to L2 point values of div (VALUE map): D' = diag(r) D with r = detJ_0/detJ_current(x_i)
      // (row scaling).  Check row proportionality with two random RT vectors, and D' C = 0.
      std::mt19937 gen(11);
      std::uniform_real_distribution<real_t> U(-1.0, 1.0);
      Vector b1(ops2.D->Width()), b2(ops2.D->Width());
      for (int i = 0; i < b1.Size(); i++) { b1(i) = U(gen); b2(i) = U(gen); }
      Vector y1(ops2.D->Height()), y2(y1), z1(y1), z2(y1);
      ops2.D->Mult(b1, y1); ops2.D->Mult(b2, y2);
      M0.mc->ops->D->Mult(b1, z1); M0.mc->ops->D->Mult(b2, z2);
      double prop = 0.0, rmin = 1e300, rmax = -1e300;
      for (int i = 0; i < y1.Size(); i++)
      {
         const double r1 = y1(i) / z1(i), r2 = y2(i) / z2(i);
         prop = std::max(prop, std::fabs(r1 - r2) / std::fabs(r1));
         rmin = std::min(rmin, r1); rmax = std::max(rmax, r1);
      }
      const double dc = ops2.MaxAbsDC(), cg = ops2.MaxAbsCG(), dc0 = M0.mc->ops->MaxAbsDC();
      if (rank == 0)
      {
         std::cout << "  (iv) after moving to eps=0.6,s=1: max|C-C'|=" << dC << " max|G-G'|=" << dG
                   << " ; D' = diag(r) D with r in [" << rmin << ", " << rmax << "] (row-proportionality defect " << prop
                   << ") ; max|D' C|=" << dc << " (D C stale: " << dc0 << ", entries of D are O(N^3 p^2)) max|C G|=" << cg << std::endl;
      }
      check(dC == 0.0 && dG == 0.0, "(iv) C_h and G_h independent of the node positions (bitwise)");
      check(prop < 1e-9 && dc < 1e-9 * rmax && dc0 < 1e-9, "(iv) D_h = row scaling of the reference D_h, D_h C_h = 0 to roundoff");
      // label interpolation vs Newton at quadrature points
      GridFunction labgf(const_cast<FiniteElementSpace *>(M0.pm().GetNodalFESpace()));
      labgf = M0.lab;
      const IntegrationRule &ir = IntRules.Get(Geometry::CUBE, 2 * p + 2 * q + 2);
      double worst = 0.0;
      Vector xv(3), Xl(3);
      for (int e = 0; e < M0.pm().GetNE(); e++)
      {
         ElementTransformation *T = M0.pm().GetElementTransformation(e);
         for (int k = 0; k < ir.GetNPoints(); k++)
         {
            const IntegrationPoint &ip = ir.IntPoint(k);
            T->SetIntPoint(&ip);
            T->Transform(ip, xv);
            labgf.GetVectorValue(e, ip, Xl);
            real_t X[3], xx[3] = {xv(0), xv(1), xv(2)};
            InverseMap(1, 0.6, xx, X);
            for (int d = 0; d < 3; d++) { worst = std::max(worst, (double)std::fabs(X[d] - Xl(d))); }
         }
      }
      if (rank == 0)
      {
         std::cout << "  (iv) info: max |label interpolated from the label dofs - Phi^{-1}(x_h)| at quadrature points (s=1, eps=0.6, q="
                   << q << ", N=" << N << ") = " << worst << std::endl;
      }
   }
   if (rank == 0) { std::cout << (fails ? "SELFTEST FAILED" : "SELFTEST PASSED") << " (" << fails << " failures)" << std::endl; }
   return fails;
}

int main(int argc, char *argv[])
{
   Mpi::Init(argc, argv);
   Hypre::Init();
   MPI_Comm comm = MPI_COMM_WORLD;
   const int rank = Mpi::WorldRank();

   Cfg cfg;
   bool selftest = false, floors = false;
   double eps = 0.3;
   OptionsParser args(argc, argv);
   args.AddOption(&cfg.N, "-N", "--N", "Cells per direction.");
   args.AddOption(&cfg.p, "-p", "--order", "ND order p.");
   args.AddOption(&cfg.q, "-q", "--geom-order", "Geometric order (default p).");
   args.AddOption(&eps, "-eps", "--eps", "Displacement amplitude.");
   args.AddOption(&cfg.K, "-K", "--steps", "Number of time steps over [0,T].");
   args.AddOption(&cfg.nrez, "-rezones", "--rezones", "Number of rezones (equally spaced in (0,T], last at T).");
   args.AddOption(&cfg.op, "-op", "--op", "A_int|A_l2|A_pt|B_int|B_l2|B_l2c.");
   args.AddOption(&cfg.gauge, "-gauge", "--gauge", "none | coulomb (A-route; Coulomb fix of the transferred a on the new mesh) | coulomb_pre (fix on the deformed source mesh before the transfer).");
   args.AddOption(&cfg.field, "-field", "--field", "abc | mod | mod2.");
   args.AddOption(&cfg.b0, "-b0", "--b0", "-no-b0", "--no-b0", "Mean field B0 = (0.3,-0.2,0.5).");
   args.AddOption(&cfg.variant, "-fv", "--fv", "Displacement variant (1,2,3).");
   args.AddOption(&cfg.every, "-every", "--every", "Log every k-th step (rezone steps and the last step always).");
   args.AddOption(&cfg.csv, "-csv", "--csv", "Append rows to this CSV file.");
   args.AddOption(&cfg.tag, "-tag", "--tag", "Series label.");
   args.AddOption(&floors, "-floors", "--floors", "-no-floors", "--no-floors", "Only compute the projection floors at the rezone times (no transfer).");
   args.AddOption(&selftest, "-selftest", "--selftest", "-no-selftest", "--no-selftest", "Run the self tests and exit.");
   args.Parse();
   if (!args.Good()) { if (rank == 0) { args.PrintUsage(std::cout); } return 1; }
   if (rank == 0) { args.PrintOptions(std::cout); }
   cfg.eps = eps;
   if (selftest)
   {
      const int q = cfg.q < 0 ? (cfg.p >= 2 ? 2 : cfg.p) : cfg.q;
      return SelfTest(comm, cfg.N, cfg.p, q, 20) == 0 ? 0 : 1;
   }
   if (floors) { RunFloors(cfg, comm); return 0; }
   Result res;
   const double t0 = WallTime();
   Run(cfg, res, comm);
   if (rank == 0)
   {
      std::cout << std::scientific << std::setprecision(4)
                << "FINAL op=" << cfg.op << " gauge=" << cfg.gauge << " p=" << cfg.p << " eps=" << cfg.eps
                << " rezones=" << cfg.nrez << " b0=" << cfg.b0 << ": ||b_T-b_0||/||b_0||=" << res.final_dB
                << " errB(T)=" << res.final_errB << " E(T)/E(0)=" << res.final_E_ratio << " H(T)/H(0)=" << res.final_H_ratio
                << " flux_err=" << res.final_flux_err << " div=" << res.final_div << " transfer_time=" << res.t_transfer_total
                << " s, notfound=" << res.nnotfound << ", wall " << WallTime() - t0 << " s" << std::defaultfloat << std::endl;
   }
   return 0;
}
