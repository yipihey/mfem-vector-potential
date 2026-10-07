// gauge_check.cpp -- unit checks of vp::CoulombGauge and the "B_l2c == curl of least-squares A" identity.
//
//   mpirun -np 1 bin/gauge_check -N 8 -p 2 -eps 0.3 [-lsq]
#include "exp_transfer_common.hpp"
#include "vp_gauge.hpp"
#include <random>

using namespace vp;


// ---------------------------------------------------------------------------------------------
// B_l2c (M-orthogonal projection of the L2-transferred b onto ker D_h) vs the curl of the
// least-squares vector potential  argmin_a ||C a - b_l2||_M  (Coulomb gauge, curl-curl + CG + AMS).
// Source mesh = f1 eps, target mesh = uniform; field abc, no B0 (so the harmonic part is zero).
// ---------------------------------------------------------------------------------------------
static void RunLsqCheck(MPI_Comm comm, int N, int p, real_t eps)
{
   const int rank = Mpi::WorldRank();
   const int q = p;
   MeshSpec s1; s1.N = N; s1.eps = eps; s1.fv = 1;
   MeshSpec s2; s2.N = N; s2.eps = 0.0;
   Problem pr("abc", false, 0.0, 2.0);
   MeshCase M1(comm, s1, p, q), M2(comm, s2, p, q);
   MeshState S1(M1, pr), S2(M2, pr);
   Vector a1, b1, aA, bA;
   S1.ProjectExactA(pr, false, a1);
   S1.Compose(a1, b1);
   RemoteEvaluator ev(*M1.pm);
   TransferOptions topt;
   ParGridFunction a1gf(M1.sp->ND.get()), b1gf(M1.sp->RT.get());
   a1gf.SetFromTrueDofs(a1);
   b1gf.SetFromTrueDofs(b1);
   Transfer Tl2(TransferKind::B_l2, ev, *M2.sp->RT, M2.ops.get(), topt);
   Transfer Tl2c(TransferKind::B_l2c, ev, *M2.sp->RT, M2.ops.get(), topt);
   Transfer TAl2(TransferKind::A_l2, ev, *M2.sp->ND, M2.ops.get(), topt);
   ParGridFunction g_l2(M2.sp->RT.get()), g_l2c(M2.sp->RT.get()), g_Al2(M2.sp->ND.get());
   Tl2.Apply(b1gf, g_l2);
   Tl2c.Apply(b1gf, g_l2c);
   TAl2.Apply(a1gf, g_Al2);
   Vector b_l2, b_l2c, a_Al2, b_Al2;
   g_l2.GetTrueDofs(b_l2);
   g_l2c.GetTrueDofs(b_l2c);
   g_Al2.GetTrueDofs(a_Al2);
   S2.Compose(a_Al2, b_Al2);

   const int qq = M2.pm->GetNodalFESpace()->GetMaxElementOrder();
   const IntegrationRule &ir = IntRules.Get(Geometry::CUBE, 2 * p + 2 * qq + 2);
   ParBilinearForm MRT(M2.sp->RT.get());
   auto *vi = new VectorFEMassIntegrator;
   vi->SetIntRule(&ir);
   MRT.AddDomainIntegrator(vi);
   MRT.Assemble();
   MRT.Finalize();
   std::unique_ptr<HypreParMatrix> Mrt(MRT.ParallelAssemble());
   ParBilinearForm MND(M2.sp->ND.get());
   auto *vi2 = new VectorFEMassIntegrator;
   vi2->SetIntRule(&ir);
   MND.AddDomainIntegrator(vi2);
   MND.Assemble();
   MND.Finalize();
   std::unique_ptr<HypreParMatrix> Mnd(MND.ParallelAssemble());
   auto Mnorm = [&](const Vector &x)
   {
      Vector y(x.Size());
      Mrt->Mult(x, y);
      return std::sqrt(std::max(InnerProduct(comm, x, y), (real_t)0.0));
   };
   // (1) D_h b_l2c = 0, D_h b_l2 != 0
   const DivDiag dc = DivergenceDiag(*M2.ops, b_l2c), dl = DivergenceDiag(*M2.ops, b_l2);
   // (2) orthogonality of (b_l2c - b_l2) to ker D_h = range C (+ harmonic): random a
   Vector diff(b_l2c);
   diff -= b_l2;
   Vector Mdiff(diff.Size());
   Mrt->Mult(diff, Mdiff);
   std::mt19937 gen(777);
   std::uniform_real_distribution<double> U(-1.0, 1.0);
   real_t max_orth = 0.0;
   for (int t = 0; t < 5; t++)
   {
      Vector arand(M2.sp->ND->GetTrueVSize()), Ca(M2.sp->RT->GetTrueVSize());
      for (int i = 0; i < arand.Size(); i++) { arand(i) = U(gen); }
      M2.ops->C->Mult(arand, Ca);
      const real_t ip = InnerProduct(comm, Mdiff, Ca);
      max_orth = std::max(max_orth, std::fabs(ip) / (Mnorm(diff) * Mnorm(Ca)));
   }
   // also the "unit-vector" version: (b_l2c - b_l2) is orthogonal to the full ker D: any v = b_l2c itself
   const real_t orth_self = std::fabs(InnerProduct(comm, Mdiff, b_l2c)) / (Mnorm(diff) * Mnorm(b_l2c));
   // (3) least-squares A: (C^T M C + delta M_ND) a = C^T M b_l2  (CG + AMS), then b_ls = C a
   std::unique_ptr<HypreParMatrix> CMC(RAP(Mrt.get(), M2.ops->C.get()));
   const real_t delta = 1e-9;
   std::unique_ptr<HypreParMatrix> Asys(Add(1.0, *CMC, delta, *Mnd));
   Vector rhs(M2.sp->ND->GetTrueVSize()), Mb(b_l2.Size()), a_ls(rhs.Size());
   Mrt->Mult(b_l2, Mb);
   M2.ops->C->MultTranspose(Mb, rhs);
   a_ls = 0.0;
   const double t0 = WallTime();
   HypreAMS ams(*Asys, M2.sp->ND.get());
   ams.SetPrintLevel(0);
   CGSolver cgs(comm);
   cgs.SetOperator(*Asys);
   cgs.SetPreconditioner(ams);
   cgs.SetRelTol(1e-13);
   cgs.SetAbsTol(0.0);
   cgs.SetMaxIter(500);
   cgs.SetPrintLevel(-1);
   cgs.Mult(rhs, a_ls);
   const double t_ams = WallTime() - t0;
   Vector b_ls(b_l2.Size());
   M2.ops->C->Mult(a_ls, b_ls);
   Vector d1(b_ls), d2(b_Al2), d3(b_l2);
   d1 -= b_l2c;
   d2 -= b_l2c;
   d3 -= b_l2c;
   const real_t nrm = Mnorm(b_l2c);
   // gauge fix of a_ls (should already be Coulomb since the shifted system kills the gradient part)
   CoulombGauge cg2(*M2.sp, *M2.ops);
   // floors: b_l2c error to the exact field
   const StateDiag dref = S2.diag->Evaluate(nullptr, b_l2c, false);
   if (rank == 0)
   {
      std::cout << std::scientific << std::setprecision(3)
                << "\nLSQ check: p=" << p << " N=" << N << " source f1 eps=" << eps << " -> uniform, abc, no B0\n"
                << "  div_h(b_l2c) max = " << dc.maxabs << "   div_h(b_l2) max = " << dl.maxabs << "\n"
                << "  M-orthogonality of (b_l2c - b_l2) to C a_rand (5 random a): max normalised = " << max_orth
                << ",  to b_l2c itself: " << orth_self << "\n"
                << "  LSQ solve (CG+AMS, shift " << delta << "): iters " << cgs.GetNumIterations() << " converged "
                << cgs.GetConverged() << " time " << t_ams << " s; coulomb residual of a_ls "
                << cg2.CoulombResidual(a_ls) << ", grad fraction of a_ls " << cg2.GradFraction(a_ls) << "\n"
                << "  ||C a_ls - b_l2c||_M / ||b_l2c||_M = " << Mnorm(d1) / nrm << "\n"
                << "  ||b_l2 - b_l2c||_M / ||b_l2c||_M   = " << Mnorm(d3) / nrm << "  (the divergence-cleaning change)\n"
                << "  ||C a_A_l2 - b_l2c||_M / ||b_l2c||_M = " << Mnorm(d2) / nrm
                << "  (A_l2 then B = C a: different operator)\n"
                << "  B error of b_l2c vs exact: " << dref.errB_rel << std::defaultfloat << std::endl;
   }
}

int main(int argc, char *argv[])
{
   Mpi::Init(argc, argv);
   Hypre::Init();
   MPI_Comm comm = MPI_COMM_WORLD;
   const int rank = Mpi::WorldRank();
   int N = 8, p = 2;
   real_t eps = 0.3;
   bool lsq = false;
   OptionsParser args(argc, argv);
   args.AddOption(&N, "-N", "--N", "cells");
   args.AddOption(&p, "-p", "--order", "order");
   args.AddOption(&eps, "-eps", "--eps", "deformation of the mesh");
   args.AddOption(&lsq, "-lsq", "--lsq", "-no-lsq", "--no-lsq", "B_l2c vs least-squares-A check");
   args.Parse();
   if (!args.Good()) { if (rank == 0) { args.PrintUsage(std::cout); } return 1; }
   const int q = p;

   MeshSpec sB; sB.N = N; sB.eps = eps; sB.fv = 1;
   MeshSpec sA; sA.N = N; sA.eps = 0.0;
   Problem pr("abc", false, 0.0, 2.0);
   MeshCase MB(comm, sB, p, q);
   MeshState SB(MB, pr);
   CoulombGauge cg(*MB.sp, *MB.ops);
   if (rank == 0)
   {
      std::cout << "p=" << p << " N=" << N << " eps=" << eps << " ND=" << MB.sp->NdofsND() << " H1="
                << MB.sp->NdofsH1() << " gauge setup " << cg.SetupTime() << " s" << std::endl;
   }
   const real_t mism = cg.StiffnessMismatch();
   const real_t cgmax = MB.ops->MaxAbsCG();
   if (rank == 0)
   {
      std::cout << "max|L - K|/max|K| (L = G^T M_ND G, K = DiffusionIntegrator) = " << mism << "\n"
                << "max|C G| = " << cgmax << std::endl;
   }
   for (real_t gamp : {0.0, 0.5})
   {
      for (real_t gk : {2.0, 6.0})
      {
         if (gamp == 0.0 && gk != 2.0) { continue; }
         Problem prg("abc", false, gamp, gk);
         Vector a, b, a2, b2;
         SB.ProjectExactA(prg, gamp != 0.0, a);
         SB.Compose(a, b);
         const StateDiag d0 = SB.diag->Evaluate(&a, b, true);
         a2 = a;
         const GaugeFixStats st = cg.Apply(a2, 1e-12);
         SB.Compose(a2, b2);
         const StateDiag d1 = SB.diag->Evaluate(&a2, b2, true);
         Vector zero(b.Size());
         zero = 0.0;
         const real_t dB = SB.diag->RTDistance(b, b2) / SB.diag->RTDistance(b, zero);
         if (rank == 0)
         {
            std::cout << std::scientific << std::setprecision(3) << "gauge g=" << gamp << " k=" << gk << "pi: iters " << st.iters
                      << " time " << st.time << " resid " << st.resid_before << " -> " << st.resid_after
                      << " gradfrac " << st.grad_fraction_before << " -> " << st.grad_fraction_after
                      << " (diag: " << d0.grad_frac << " -> " << d1.grad_frac << ")\n   |dB|/|B| = " << dB
                      << " H " << d0.helicity << " -> " << d1.helicity << " (dH/hel_scale = "
                      << (d1.helicity - d0.helicity) / d0.hel_scale << ")  |a| " << d0.a_L2 << " -> " << d1.a_L2
                      << std::defaultfloat << std::endl;
         }
         // Jacobi variant
         for (int k : {5, 20})
         {
            Vector a3 = a;
            const GaugeFixStats sj = cg.JacobiGaugeSmooth(a3, k);
            if (rank == 0)
            {
               std::cout << std::scientific << std::setprecision(3) << "   jacobi:" << k << " gradfrac " << sj.grad_fraction_before
                         << " -> " << sj.grad_fraction_after << " time " << sj.time << std::defaultfloat << std::endl;
            }
         }
      }
   }
   if (lsq) { RunLsqCheck(comm, N, p, eps); }
   return 0;
}
