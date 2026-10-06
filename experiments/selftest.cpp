// selftest.cpp -- consistency checks of the vp library.  Prints PASS/FAIL per
// check and an overall verdict; exit code 0 iff everything passed.
//  1. finite-difference check of the analytic curls (abc, mod, mod2)
//  2. mesh: fully periodic (3 N^3 faces), slice families have N^2 faces
//  3. D_h C_h = 0 and C_h G_h = 0 on an eps=0.25 deformed mesh, p=1..4
//  4. constant B0: Pi_RT(B0) exact (L2 error, div, slice fluxes, quadrature
//     cross-check of the flux routine)
//  5. integrated projection commutes: C_h Pi_ND^int A = Pi_RT^int curl A,
//     G_h Pi_H1 chi = Pi_ND^int grad chi, and the integrated ND dofs are
//     edge-continuous (same function after GetTrueDofs round trip)
//  6. helicity of the ABC field and matrix energy vs matrix-free energy
#include "vp_core.hpp"
#include <iostream>
#include <iomanip>

using namespace vp;

static int nfail = 0, ncheck = 0;
static int myrank = 0;
static void Check(bool ok, const std::string &what, double val, double tol)
{
   ncheck++;
   if (!ok) { nfail++; }
   if (myrank == 0)
   {
      std::cout << (ok ? "PASS " : "FAIL ") << std::left << std::setw(70) << what
                << " value=" << std::scientific << std::setprecision(3) << val
                << " tol=" << tol << std::endl;
   }
}

int main(int argc, char *argv[])
{
   Mpi::Init(argc, argv);
   Hypre::Init();
   myrank = Mpi::WorldRank();
   MPI_Comm comm = MPI_COMM_WORLD;
   FieldParams fp;

   // 1. curl finite-difference test
   {
      bool pass = false;
      real_t w = SelfTestCurl(fp, 1e-6, &pass);
      Check(pass, "analytic curl vs central differences (abc,mod,mod2)", w, 1e-6);
   }

   const int N = 4;
   // 2./3. deformed meshes
   for (int p = 1; p <= 4; p++)
   {
      const int q = std::max(2, std::min(p, 3));
      std::unique_ptr<ParMesh> pm(MakePeriodicBox(comm, N, q));
      auto fams = RecordSlices(*pm, N, {0, N / 2});
      bool ok_f = true;
      for (auto &f : fams) { ok_f = ok_f && (f.nfaces_global >= N * N); }
      DeformMesh(*pm, 0.25, 1);
      auto js = JacobianStats(*pm);
      std::string tag = " p=" + std::to_string(p) + " q=" + std::to_string(q);
      Check(js.n_neg == 0, "eps=0.25 mesh valid (no detJ<=0)" + tag, js.min_det, 0);
      Spaces sp(*pm, p);
      Operators ops(sp);
      const real_t dc = ops.MaxAbsDC(), cg = ops.MaxAbsCG();
      Check(dc < 1e-12, "max|D_h C_h| deformed" + tag, dc, 1e-12);
      Check(cg < 1e-12, "max|C_h G_h| deformed" + tag, cg, 1e-12);
      if (p == 1)
      {
         Check(ok_f, "slice families have >= N^2 faces", (double)fams[0].nfaces_global, N * N);
      }

      // random-vector action D_h C_h a, C_h G_h c
      Vector a(sp.ND->GetTrueVSize()), b(sp.RT->GetTrueVSize());
      a.Randomize(3);
      ops.ApplyCurl(a, b);
      auto dd = DivergenceDiag(ops, b);
      Check(dd.maxabs < 1e-11 * std::max(1.0, b.Normlinf()), "max|D_h C_h a| random a" + tag, dd.maxabs, 1e-11);
      Vector c(sp.H1->GetTrueVSize()), ga(sp.ND->GetTrueVSize()), cga(sp.RT->GetTrueVSize());
      c.Randomize(5);
      ops.G->Mult(c, ga);
      ops.C->Mult(ga, cga);
      real_t m = cga.Normlinf();
      MPI_Allreduce(MPI_IN_PLACE, &m, 1, MPITypeMap<real_t>::mpi_type, MPI_MAX, comm);
      Check(m < 1e-11, "max|C_h G_h c| random c" + tag, m, 1e-11);
   }

   // 4. constant B0 exactness, slices (deformed mesh too)
   for (int p = 1; p <= 4; p++)
   {
      for (real_t eps : {0.0, 0.25})
      {
         const int q = 2;
         std::unique_ptr<ParMesh> pm(MakePeriodicBox(comm, N, q));
         auto fams = RecordSlices(*pm, N, {0, N / 2});
         DeformMesh(*pm, eps, 1);
         Spaces sp(*pm, p);
         Operators ops(sp);
         Vector B0v(3); B0v(0) = 0.3; B0v(1) = -0.2; B0v(2) = 0.5;
         VectorConstantCoefficient B0c(B0v);
         ParGridFunction b0(sp.RT.get());
         // On a non-affine mesh a constant field is NOT in the mapped RT space
         // (contravariant Piola); only for eps=0 (affine) is the projection
         // exact.  The integrated (flux) projection is divergence-free exactly
         // on any mesh and has exact slice fluxes.
         const ProjMode mode = (eps == 0.0) ? ProjMode::Pointwise : ProjMode::Integrated;
         ProjectB(b0, B0c, mode);
         const real_t e = L2ErrorVec(b0, B0c, 2 * p + 2 * q + 2);
         Vector b0t; b0.GetTrueDofs(b0t);
         auto dd = DivergenceDiag(ops, b0t);
         std::string tag = " p=" + std::to_string(p) + " eps=" + std::to_string(eps);
         if (eps == 0.0)
         {
            Check(e < 1e-13, "||Pi_RT(B0) - B0||_L2 (affine mesh)" + tag, e, 1e-13);
         }
         else if (myrank == 0)
         {
            std::cout << "INFO  non-affine mesh: ||Pi_RT^int(B0) - B0||_L2 = " << e << tag << std::endl;
         }
         Check(dd.maxabs < 1e-12, "max|D_h Pi_RT(B0)|" + tag, dd.maxabs, 1e-12);
         real_t worst = 0, worstq = 0;
         for (auto &f : fams)
         {
            const real_t fl = SliceFlux(b0, f);
            worst = std::max(worst, std::fabs(fl - B0v(f.dir)));
            if (Mpi::WorldSize() == 1)
            {
               const real_t fq = SliceFluxQuadrature(b0, f, 2 * p + 4);
               worstq = std::max(worstq, std::fabs(fq - fl));
            }
         }
         Check(worst < 1e-12, "slice flux of Pi_RT(B0) == B0_i" + tag, worst, 1e-12);
         if (Mpi::WorldSize() == 1)
         {
            Check(worstq < 1e-12, "SliceFlux (dofs) == SliceFluxQuadrature" + tag, worstq, 1e-12);
         }
         // flux of random curl field must not change the flux
         Vector a(sp.ND->GetTrueVSize()), bb(sp.RT->GetTrueVSize());
         a.Randomize(11);
         ops.ApplyCurl(a, bb);
         ParGridFunction bg(sp.RT.get());
         bg.SetFromTrueDofs(bb);
         real_t wf = 0;
         for (auto &f : fams) { wf = std::max(wf, std::fabs(SliceFlux(bg, f))); }
         Check(wf < 1e-11, "slice flux of random C_h a vanishes" + tag, wf, 1e-11);
         if (Mpi::WorldSize() == 1)
         {
            real_t wq = 0;
            for (auto &f : fams)
            { wq = std::max(wq, std::fabs(SliceFlux(bg, f) - SliceFluxQuadrature(bg, f, 2 * p + 4))); }
            Check(wq < 1e-10, "SliceFlux == quadrature for random curl field" + tag, wq, 1e-10);
         }
      }
   }

   // 5. integrated projections commute
   for (int p = 1; p <= 4; p++)
   {
      const int q = 2;
      std::unique_ptr<ParMesh> pm(MakePeriodicBox(comm, N, q));
      DeformMesh(*pm, 0.2, 3);
      Spaces sp(*pm, p);
      Operators ops(sp);
      Field F = MakeField("abc", fp);
      VectorFunctionCoefficient Ac(3, F.A), Bc(3, F.curlA);
      ParGridFunction a(sp.ND.get()), b(sp.RT.get());
      ProjectA(a, Ac, ProjMode::Integrated);
      ProjectB(b, Bc, ProjMode::Integrated);
      Vector at, bt, cb(sp.RT->GetTrueVSize());
      a.GetTrueDofs(at); b.GetTrueDofs(bt);
      ops.ApplyCurl(at, cb);
      cb -= bt;
      real_t m = cb.Normlinf();
      MPI_Allreduce(MPI_IN_PLACE, &m, 1, MPITypeMap<real_t>::mpi_type, MPI_MAX, comm);
      std::string tag = " p=" + std::to_string(p);
      Check(m < 1e-11, "C_h Pi_ND^int A == Pi_RT^int curl A" + tag, m, 1e-11);
      // pointwise version commutes only up to O(h^p)
      ParGridFunction ap(sp.ND.get());
      ProjectA(ap, Ac, ProjMode::Pointwise);
      Vector apt; ap.GetTrueDofs(apt);
      ops.ApplyCurl(apt, cb);
      cb -= bt;
      real_t mp = cb.Normlinf();
      MPI_Allreduce(MPI_IN_PLACE, &mp, 1, MPITypeMap<real_t>::mpi_type, MPI_MAX, comm);
      if (myrank == 0) { std::cout << "INFO  pointwise ND: max|C_h Pi^pt A - Pi_RT^int B| = " << mp << tag << std::endl; }

      // gradient commutation
      Gauge gg(0.3, 2.0);
      FunctionCoefficient chic([&](const Vector &x) { return gg.chi(x); });
      VectorFunctionCoefficient gradc(3, [&](const Vector &x, Vector &v) { gg.grad(x, v); });
      ParGridFunction chi(sp.H1.get()), ga(sp.ND.get());
      chi.ProjectCoefficient(chic);
      ProjectA(ga, gradc, ProjMode::Integrated);
      Vector ct, gt, g2(sp.ND->GetTrueVSize());
      chi.GetTrueDofs(ct); ga.GetTrueDofs(gt);
      ops.G->Mult(ct, g2);
      g2 -= gt;
      real_t mg = g2.Normlinf();
      MPI_Allreduce(MPI_IN_PLACE, &mg, 1, MPITypeMap<real_t>::mpi_type, MPI_MAX, comm);
      Check(mg < 1e-11, "G_h Pi_H1 chi == Pi_ND^int grad chi" + tag, mg, 1e-11);
      ParGridFunction gp(sp.ND.get());
      ProjectA(gp, gradc, ProjMode::Pointwise);
      Vector gpt; gp.GetTrueDofs(gpt);
      gpt -= g2; // g2 = G chi - Pi^int grad chi ; just report pointwise defect
      gp.GetTrueDofs(gpt);
      ops.G->Mult(ct, g2);
      gpt -= g2;
      real_t mgp = gpt.Normlinf();
      MPI_Allreduce(MPI_IN_PLACE, &mgp, 1, MPITypeMap<real_t>::mpi_type, MPI_MAX, comm);
      if (myrank == 0) { std::cout << "INFO  pointwise ND: max|Pi^pt grad chi - G_h Pi chi| = " << mgp << tag << std::endl; }
   }

   // 6. helicity and energy via matrices vs. matrix-free; exact values
   for (int p = 1; p <= 3; p++)
   {
      const int q = 1;
      std::unique_ptr<ParMesh> pm(MakePeriodicBox(comm, 8, q));
      Spaces sp(*pm, p);
      Operators ops(sp);
      Field F = MakeField("abc", fp);
      VectorFunctionCoefficient Ac(3, F.A), Bc(3, F.curlA);
      ParGridFunction a(sp.ND.get()), b(sp.RT.get());
      ProjectA(a, Ac, ProjMode::Integrated);
      Vector at, bt(sp.RT->GetTrueVSize());
      a.GetTrueDofs(at);
      ops.ApplyCurl(at, bt);
      b.SetFromTrueDofs(bt);
      const real_t Hm = HelicityMass(sp, at, bt);
      const real_t Hg = InnerGF(a, b, 2 * p + 2 * q + 2);
      const real_t Hex = ExactHelicity("abc", fp);
      const real_t Em = EnergyMass(sp, bt);
      const real_t Eg = 0.5 * InnerGF(b, b, 2 * p + 2 * q + 2);
      const real_t Eq = 0.5 * InnerCoef(*pm, Bc, Bc, 2 * p + 2 * q + 6);
      const real_t Eex = ExactEnergy("abc", fp);
      std::string tag = " p=" + std::to_string(p) + " N=8";
      Check(std::fabs(Hm - Hg) < 1e-12 * std::fabs(Hex), "helicity: mixed-mass == matrix-free" + tag, Hm - Hg, 1e-12);
      Check(std::fabs(Em - Eg) < 1e-12 * Eex, "energy: mass matrix == matrix-free" + tag, Em - Eg, 1e-12);
      Check(std::fabs(Eq - Eex) < 1e-10 * Eex, "exact energy: quadrature == closed form" + tag, Eq - Eex, 1e-10);
      const real_t tolH = 0.5 * std::pow(1.0 / 8, p) * 20;
      Check(std::fabs(Hm - Hex) < tolH * Hex, "helicity(abc) ~ k(a^2+b^2+c^2)" + tag, (Hm - Hex) / Hex, tolH);
      if (myrank == 0)
      {
         std::cout << "INFO  H_h=" << std::setprecision(12) << Hm << " H_exact=" << Hex
                   << " E_h=" << Em << " E_exact=" << Eex << tag << std::endl;
      }
   }
   // exact energies of mod / mod2 by quadrature (undeformed mesh)
   {
      std::unique_ptr<ParMesh> pm(MakePeriodicBox(comm, 4, 1));
      for (const char *nm : {"abc", "mod", "mod2"})
      {
         Field F = MakeField(nm, fp);
         VectorFunctionCoefficient Bc(3, F.curlA);
         const real_t Eq = 0.5 * InnerCoef(*pm, Bc, Bc, 14);
         const real_t Eex = ExactEnergy(nm, fp);
         Check(std::fabs(Eq - Eex) < 1e-10 * Eex, std::string("closed-form energy ") + nm, Eq - Eex, 1e-10);
         VectorFunctionCoefficient Ac(3, F.A);
         const real_t Hq = InnerCoef(*pm, Ac, Bc, 14);
         const real_t Hex = ExactHelicity(nm, fp);
         Check(std::fabs(Hq - Hex) < 1e-10, std::string("closed-form helicity ") + nm, Hq - Hex, 1e-10);
      }
   }

   if (myrank == 0)
   {
      std::cout << "\n" << (nfail == 0 ? "ALL PASS" : "SOME FAIL") << " (" << ncheck - nfail
                << "/" << ncheck << " checks passed)" << std::endl;
   }
   return nfail == 0 ? 0 : 1;
}
