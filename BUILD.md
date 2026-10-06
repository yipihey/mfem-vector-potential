# Build and reproduction notes

All paths are under `/home/user/mfem-vector-potential` (`$REPO`); `EXT=$REPO/external`.
Everything below was actually run on a 4-core / 15 GB Ubuntu 24.04.5 container.
`scripts/build_all.sh` mirrors these steps from a clean `external/` (not re-executed after writing).
Raw logs are in `external/logs/`.

## Pinned sources

| Component | URL | Commit |
|---|---|---|
| mfem (upstream master, v4.10.1 dev) | https://github.com/mfem/mfem.git | `ed7e19c2acfa0cae13d8f8c453237764d4ae7374` |
| mfem-develop (fork, MFEM 4.8.1 based) | https://github.com/ruijie-xi/mfem-develop.git | `842fac6a3b038eaee848921f096e80ceb5f895fa` |
| MHD-ALE | https://github.com/ruijie-xi/MHD-ALE.git | `0e0564ed92acff6ad84d3f86f6609bb7ad916ce8` |
| gslib (master, `v1.0.9-1-g375eda2`) | https://github.com/Nek5000/gslib.git | `375eda25b04205d1f749ab0e00a0bb3f22225b22` |

Note: `external/MHD-ALE/mfem-develop/` is an untracked extra copy of the fork (same commit 842fac6); it is not used by any build step.

## Toolchain

- gcc / g++ / mpicxx: 13.3.0 (Ubuntu 13.3.0-6ubuntu2~24.04.1)
- mpirun: Open MPI 4.1.6
- cmake 3.28.3, ninja 1.11.1, GNU make
- apt: `libhypre-2.28.0` / `libhypre-dev` 2.28.0-8build2, `libmetis5` / `libmetis-dev` 5.1.0.dfsg-7build3,
  `openmpi-bin` / `libopenmpi3t64` / `libopenmpi-dev` 4.1.6-7ubuntu2,
  `liblapack-dev` / `libblas-dev` 3.12.0-3build1.1

## Install layout

| What | Path |
|---|---|
| GSLIB (static, MPI) | `external/gslib/build/{lib/libgs.a,include/gslib.h,include/gslib/}` (no separate install) |
| hypre / metis shims | `external/deps/hypre`, `external/deps/metis` (symlinks to apt files, see below) |
| Upstream MFEM 4.10.1 | `external/install/mfem` (`include/mfem/mfem.hpp`, `lib/libmfem.a`, `share/mfem/config.mk`, `lib/cmake/mfem/MFEMConfig.cmake`) |
| Fork MFEM 4.8.1 (for MHD-ALE) | `external/install/mfem-develop` (same layout) |
| Upstream build tree (ex3/ex3p/ex4/ex4p, gslib miniapps) | `external/build-mfem` |
| MHD-ALE executable | `external/MHD-ALE/MHD` |

## Commands

### 1. GSLIB
```
cd $EXT && git clone https://github.com/Nek5000/gslib.git gslib      # 375eda25...
cd gslib && make CC=mpicc MPI=1 -j4 > ../logs/gslib-build.log 2>&1   # 6 s
```
Produces `build/lib/libgs.a` and `build/include/gslib.h` (+ `gslib/config.h`). Compiled with `-DGSLIB_USE_MPI`.

### Workaround: hypre / metis shim directories
MFEM's `FindHYPRE.cmake` / `FindMETIS.cmake` expect `<DIR>/include` and `<DIR>/lib`. With `-DHYPRE_DIR=/usr` the header
was found (version 22800) but the library was not (Debian multiarch: `/usr/lib/x86_64-linux-gnu`, headers in `/usr/include/hypre`),
so shims were created:
```
mkdir -p $EXT/deps/hypre/{include,lib} $EXT/deps/metis/{include,lib}
ln -sf /usr/include/hypre/*                  $EXT/deps/hypre/include/
ln -sf /usr/lib/x86_64-linux-gnu/libHYPRE.so $EXT/deps/hypre/lib/libHYPRE.so
ln -sf /usr/include/metis.h                  $EXT/deps/metis/include/metis.h
ln -sf /usr/lib/x86_64-linux-gnu/libmetis.so $EXT/deps/metis/lib/libmetis.so
```

### 2. Upstream MFEM (parallel, Release, GSLIB)
```
cd $EXT && mkdir build-mfem && cd build-mfem
cmake -G Ninja ../mfem -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX=$EXT/install/mfem \
  -DMFEM_USE_MPI=YES -DMFEM_USE_METIS=YES -DMFEM_USE_METIS_5=YES -DMFEM_USE_GSLIB=YES -DMFEM_USE_LAPACK=YES \
  -DHYPRE_DIR=$EXT/deps/hypre -DMETIS_DIR=$EXT/deps/metis -DGSLIB_DIR=$EXT/gslib/build
cmake . -DMFEM_ENABLE_EXAMPLES=YES -DMFEM_ENABLE_MINIAPPS=YES     # needed to get ex3p etc. as targets
cmake --build . -j4 --target mfem ex3 ex3p ex4 ex4p field-interp findpts pfindpts   # 7m13s wall
# install workaround, see below:
cmake . -DMFEM_ENABLE_EXAMPLES=NO -DMFEM_ENABLE_MINIAPPS=NO
cmake --build . -j4 && cmake --install .
cmake . -DMFEM_ENABLE_EXAMPLES=YES -DMFEM_ENABLE_MINIAPPS=YES     # restore
```
Workaround: with examples enabled, `cmake --install` fails ("cannot find .../examples/ex0") because it wants every example
built. Toggling examples/miniapps off for the install avoids building ~150 targets; the already-built binaries stay in
`build-mfem/examples` and `build-mfem/miniapps/gslib`. The configured MFEM has MPI, hypre 2.28, METIS 5, LAPACK, GSLIB; the version string is `MFEM v4.10.1 (development)`.

### 4. Fork MFEM (for MHD-ALE) -- cmake install did produce config.mk, so no GNU-make fallback was needed
```
cd $EXT && mkdir build-mfem-develop && cd build-mfem-develop
cmake -G Ninja ../mfem-develop -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX=$EXT/install/mfem-develop \
  -DMFEM_USE_MPI=YES -DMFEM_USE_METIS=YES -DMFEM_USE_METIS_5=YES -DMFEM_USE_GSLIB=YES -DMFEM_USE_LAPACK=YES \
  -DHYPRE_DIR=$EXT/deps/hypre -DMETIS_DIR=$EXT/deps/metis -DGSLIB_DIR=$EXT/gslib/build
cmake --build . -j4          # 4m05s wall (library only; examples off by default)
cmake --install .
```
`install/mfem-develop/share/mfem/config.mk` reports `MFEM_VERSION_STRING = 4.8.1`, `MFEM_USE_MPI/GSLIB = YES`, and
`include/mfem/mfem.hpp` exists, matching the MHD-ALE makefile (`-I$(MFEM_INSTALL_DIR)/include/mfem`). The fork sources were not modified.

### 5. MHD-ALE
```
cd $EXT/MHD-ALE
make MFEM_INSTALL_DIR=$EXT/install/mfem-develop -j4        # 19 s, no errors (warnings only); makefile untouched
```
Compiled with `-std=c++11 -O3` (taken from the fork's config.mk). `./MHD -h` lists the options.

## Reproduction of upstream examples (step 3)
OpenMPI refuses to run as root unless `OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1` are exported (done); `--oversubscribe` used for 2 ranks.
```
cd $EXT/build-mfem/examples
mpirun -np 2 --oversubscribe ./ex3p -m $EXT/mfem/data/beam-hex.mesh -o 2 -no-vis > $EXT/logs/ex3p.log
mpirun -np 2 --oversubscribe ./ex4p -m $EXT/mfem/data/beam-hex.mesh -o 2 -no-vis > $EXT/logs/ex4p.log
./ex3 -m $EXT/mfem/data/beam-hex.mesh -o 2 -no-vis > $EXT/logs/ex3.log
cd ../miniapps/gslib
mpirun -np 2 --oversubscribe ./pfindpts -no-vis > $EXT/logs/pfindpts.log
./field-interp -no-vis > $EXT/logs/field-interp.log        # also ./field-interp -no-vis -fts 3 -ft 3 (ND source/target)
```
Key output (all exit code 0):
- `ex3p` (H(curl), 2 ranks, ND order 2): `Number of finite element unknowns: 821568`, `PCG Iterations = 15`, `Final PCG Relative Residual Norm = 8.91258e-13`, `|| E_h - E ||_{L^2} = 0.000150664` (PCG wall time 25.7 s on the oversubscribed container).
- `ex4p` (H(div), 2 ranks, RT order 2): `Number of finite element unknowns: 803840`, 14 iterations, `Average reduction factor = 0.129289`, `|| F_h - F ||_{L^2} = 0.00287333`.
- `ex3` (serial): 821568 unknowns, 291+ iterations, `|| E_h - E ||_{L^2} = 0.000253425` (the serial example uses a different preconditioner than the parallel AMS-based run, hence the iteration count and slightly different error).
- `pfindpts`: `Searched total points: 200 ... Points not found: 0`, `Max interp error: 0.0055`, `Max dist^2 (of found): 2.7e-27` (GSLIB with MPI works).
- `field-interp` runs for H1 and ND (H(curl)) fields.

Check of the install for user programs (`/tmp`-style test, `#include "mfem.hpp"`, `Mpi::Init`): compiled and printed `MFEM v4.10.1 (development)`.

## MHD-ALE smoke tests (step 6)
Run from `$EXT/MHD-ALE` with the same OMPI env vars. Commands are in `external/logs/mhd-2d-cmd.txt` and `mhd-3d-cmd.txt`; outputs in `external/MHD-ALE/output/` (git-ignored) with copies of `MHD.out` in `external/logs/mhd-*.out`.

2-D Taylor-Green (test_TaylorGreen.sh settings, refine 0, 1 proc, `-tf 0.1`, `-fsri 10`):
```
mpirun -np 1 ./MHD -p 0 -rs 0 -ptf plot_time/plot_time_TaylorGreen_test.dat -s 3 -m ./mesh/Disk-4x4-quad.mesh -pv \
  -od output/TaylorGreen-Dim2-Refine0-tf0.1 -pvb TaylorGreen -ale -tf 0.1 -mst 1 -ok 2 -ot 1 -or 1 -oa 3 -cfl 0.5 -dim 2 \
  -g1 0.85 -g2 1.02 -fsr -fsri 10 -rmv 1 -rme 1 -rma 2 -rmr 3 -cgt 1e-12 -bpt 1
```
Runs in < 1 s (4 time steps, 16 zones; no remap triggered since fsri 10 > steps; `-rma 2` prints "Helicity preserving remap is only for 3D problems. Switch to DG convection."). Tail of MHD.out:
```
L_inf  error of rho: 0.0203825   L_2 error of rho: 0.00615624
L_inf  error of v: 0.00834239    L_2 error of v: 0.00435944
L_inf  error of e: 0.094995      L_2 error of e: 0.0278165
L_inf  error of B: 0.0182107     L_2 error of B: 0.0075586
```
Files: `B_divergence_error.dat` (~1.5e-14), `helicity.dat`, `ale.dat`, `rho_min.dat`, `time_step.dat`, `gfdata/Cycle_*`, `TaylorGreen/` (ParaView). MHD.err empty.
A second 2-D run (`-rs 1 -fsri 2`, dir `...Dim2-Refine1-tf0.1-fsri2`) exercised the remap (ale.dat shows remap at steps 2 and 4): finished, L_2 errors rho 3.9e-3, v 2.0e-3, e 7.1e-3, B 2.3e-3.

3-D (`-dim 3`, ND vector potential `-oa 2`, RT B, helicity-preserving remap `-rma 2`, 2 procs). Problem 0 (Taylor-Green) is used; the suggested option set was accepted unchanged apart from adding `-p 0 -s 3 -ale -mst 1 -rmv 1 -rme 1 -rmr 3 -cgt 1e-12 -bpt 1`:
```
mpirun -np 2 ./MHD -p 0 -dim 3 -m ./mesh/Cube-4x4x4-hex.mesh -rs 1 -tf 0.2 -oa 2 -ok 1 -ot 0 -or 0 -s 3 -ale -mst 1 \
  -fsr -fsri 2 -rmv 1 -rme 1 -rma 2 -rmr 3 -cgt 1e-12 -bpt 1 -pv -pvb TaylorGreen \
  -od output/TaylorGreen-Dim3-Refine1-tf0.2-fsri2 -ptf plot_time/plot_time_TaylorGreen.dat
```
Wall time 20.6 s, MHD.err empty, 7 steps with remap at steps 2, 4, 6 (`Remapping rho/velocity/internal energy ...` in MHD.out; the 3-D helicity-preserving path for A ran). Tail:
```
L_inf  error of rho: 0.0299795   L_2: 0.0131142
L_inf  error of v: 0.0720900     L_2: 0.0308485
L_inf  error of e: 0.239437      L_2: 0.0851272
L_inf  error of B: 0.0441488     L_2: 0.0184916
```
`B_divergence_error.dat` stays ~1e-14-level; `helicity.dat` values ~1e-19 (the Taylor-Green field has zero helicity, so this checks drift only). A smaller 3-D run (`-rs 0 -tf 0.05 -fsri 5`) is a single time step and does not hit the remap, hence the larger run above.
In the repo's `plot_time/*.dat` files nothing was modified (`git status` of MHD-ALE shows only the untracked `mfem-develop/`).

## Compiling a program against the upstream install
GNU make (as MHD-ALE does):
```
include /home/user/mfem-vector-potential/external/install/mfem/share/mfem/config.mk
MFEM_FLAGS += -I$(MFEM_INSTALL_DIR)/include/mfem      # needed: mfem.hpp lives in include/mfem/
prog: prog.cpp
	$(MFEM_CXX) $(MFEM_FLAGS) $< -o $@ $(MFEM_LIBS)
```
`MFEM_CXX=/usr/bin/c++`, `MFEM_CXXFLAGS=-std=c++17 -O3 -DNDEBUG`, static `libmfem.a`; the MPI include flags are in `MFEM_TPLFLAGS`, MPI/hypre/metis/lapack/blas/libgs.a in `MFEM_EXT_LIBS` (rpaths included), so `c++` works without `mpicxx`.
CMake: `find_package(mfem REQUIRED PATHS /home/user/mfem-vector-potential/external/install/mfem/lib/cmake/mfem)` then `target_link_libraries(prog mfem)`.
