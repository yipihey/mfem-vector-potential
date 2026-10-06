// exp2_distort.cpp -- E2: same as exp1 on meshes deformed by
// x' = x + eps f_variant(x) (periodic, fixed topology); also reports Jacobian
// statistics.  Meshes with non-positive detJ are recorded as valid=0.
#include "exp_common.hpp"
int main(int argc, char *argv[])
{
   return vp::RunStatic(argc, argv, /*distort=*/true, "exp2_distort");
}
