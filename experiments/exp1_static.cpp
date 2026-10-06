// exp1_static.cpp -- E1: static representation B_h = Pi_RT(B0) + C_h a on the
// undeformed periodic box: errors, divergence, energy, helicity, slice flux,
// dofs, timings, memory.  One CSV row per run.  See exp_common.hpp.
#include "exp_common.hpp"
int main(int argc, char *argv[])
{
   return vp::RunStatic(argc, argv, /*distort=*/false, "exp1_static");
}
