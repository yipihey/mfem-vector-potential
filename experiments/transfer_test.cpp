// transfer_test.cpp -- sanity checks of vp_transfer: (1) transfer onto the SAME mesh must be the
// identity (all 7 ops), (2) transfer deformed -> uniform: points found, periodic images, div.
#include "vp_transfer.hpp"
#include <iostream>
#include <iomanip>
using namespace vp;
static int nfail = 0;
static void Check(bool ok, const std::string &w, double v, double tol)
{
   if (!ok) { nfail++; }
   if (Mpi::Root())
      std::cout << (ok ? "PASS " : "FAIL ") << std::left << std::setw(60) << w << std::scientific
                << std::setprecision(3) << " value=" << v << " tol=" << tol << std::endl;
}
int main(int argc, char *argv[])
{
   Mpi::Init(argc, argv);
   Hypre::Init();
   MPI_Comm comm = MPI_COMM_WORLD;
   int N = 4, p = 2, q = 2;
   real_t eps = 0.3;
   OptionsParser args(argc, argv);
   args.AddOption(&N, "-N", "", "");
   args.AddOption(&p, "-p", "", "");
   args.AddOption(&q, "-q", "", "");
   args.AddOption(&eps, "-eps", "", "");
   args.Parse();
   FieldParams fp;
   Field F = MakeField("abc", fp);
   VectorFunctionCoefficient Ac(3, F.A), Bc(3, F.curlA);
   for (int deformed = 0; deformed < 2; deformed++)
   {
      std::unique_ptr<ParMesh> pm(MakePeriodicBox(comm, N, q));
      if (deformed) { DeformMesh(*pm, eps, 1); }
      Spaces sp(*pm, p);
      Operators ops(sp);
      ParGridFunction a(sp.ND.get()), b(sp.RT.get());
      ProjectA(a, Ac, ProjMode::Integrated);
      ProjectB(b, Bc, ProjMode::Integrated);
      RemoteEvaluator ev(*pm);
      for (TransferKind k : ParseKindList("all"))
      {
         const bool isA = IsARoute(k);
         ParGridFunction dst(isA ? sp.ND.get() : sp.RT.get());
         const ParGridFunction &src = isA ? a : b;
         Transfer T(k, ev, *dst.ParFESpace(), &ops);
         T.Apply(src, dst);
         Vector d(dst), s(src);
         d -= s;
         real_t err = d.Normlinf(), sc = s.Normlinf();
         const auto &st = T.Stats();
         std::string nm = std::string(deformed ? "deformed " : "uniform  ") + KindName(k);
         Check(st.nnotfound == 0, nm + " nnotfound", (double)st.nnotfound, 0);
         Check(err < 2e-8 * sc, nm + " identity max|dst-src|/max|src| (shifted=" + std::to_string(st.nshifted) + ", iters=" + std::to_string(st.iters) + ", inner=" + std::to_string(st.inner_iters) + ")", err / sc, 2e-8);
         if (Mpi::Root()) std::cout << "      npts=" << st.npoints << " t_setup=" << st.t_setup << " (locate " << st.t_locate << ") t_apply=" << st.t_apply << " eval=" << st.t_eval << " dofs=" << st.t_dofs << " solve=" << st.t_solve << std::endl;
      }
   }
   if (Mpi::Root()) std::cout << (nfail ? "SOME FAILED" : "ALL PASSED") << std::endl;
   return nfail != 0;
}
