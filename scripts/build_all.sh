#!/usr/bin/env bash
# Reproduces the full external/ build: GSLIB, upstream MFEM (parallel), the
# MFEM fork required by MHD-ALE, MHD-ALE itself, and the smoke tests.
# Mirrors the commands actually used (see BUILD.md). Assumes Ubuntu 24.04 with
# apt packages: build-essential cmake ninja-build libopenmpi-dev openmpi-bin
# libhypre-dev libmetis-dev liblapack-dev libblas-dev git.
set -euo pipefail

REPO=${REPO:-/home/user/mfem-vector-potential}
EXT=$REPO/external
JOBS=${JOBS:-4}

MFEM_COMMIT=ed7e19c2acfa0cae13d8f8c453237764d4ae7374          # mfem/mfem master (v4.10.x dev)
MFEMDEV_COMMIT=842fac6a3b038eaee848921f096e80ceb5f895fa       # ruijie-xi/mfem-develop (MFEM 4.8.1 based)
MHDALE_COMMIT=0e0564ed92acff6ad84d3f86f6609bb7ad916ce8        # ruijie-xi/MHD-ALE
GSLIB_COMMIT=375eda25b04205d1f749ab0e00a0bb3f22225b22         # Nek5000/gslib master (v1.0.9-1)

# Running as root inside the container: let OpenMPI do so.
export OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1

mkdir -p "$EXT/logs" "$EXT/install"
cd "$EXT"

clone_at() { # url dir commit
  [ -d "$2/.git" ] || git clone "$1" "$2"
  git -C "$2" checkout -q "$3"
}

# ---------------------------------------------------------------- sources
clone_at https://github.com/mfem/mfem.git            mfem         $MFEM_COMMIT
clone_at https://github.com/ruijie-xi/mfem-develop.git mfem-develop $MFEMDEV_COMMIT
clone_at https://github.com/ruijie-xi/MHD-ALE.git    MHD-ALE      $MHDALE_COMMIT
clone_at https://github.com/Nek5000/gslib.git        gslib        $GSLIB_COMMIT

# ---------------------------------------------------------------- 1. GSLIB (MPI=1 is the default)
( cd gslib && make CC=mpicc MPI=1 -j"$JOBS" > "$EXT/logs/gslib-build.log" 2>&1 )
# -> gslib/build/lib/libgs.a, gslib/build/include/gslib.h

# ---------------------------------------------------------------- shims for apt hypre / metis
# MFEM's FindHYPRE/FindMETIS look for <DIR>/include and <DIR>/lib; Debian puts the
# libraries in /usr/lib/x86_64-linux-gnu and the hypre headers in /usr/include/hypre.
mkdir -p deps/hypre/include deps/hypre/lib deps/metis/include deps/metis/lib
ln -sf /usr/include/hypre/*                       deps/hypre/include/
ln -sf /usr/lib/x86_64-linux-gnu/libHYPRE.so      deps/hypre/lib/libHYPRE.so
ln -sf /usr/include/metis.h                       deps/metis/include/metis.h
ln -sf /usr/lib/x86_64-linux-gnu/libmetis.so      deps/metis/lib/libmetis.so

COMMON_CMAKE=(
  -G Ninja -DCMAKE_BUILD_TYPE=Release
  -DMFEM_USE_MPI=YES -DMFEM_USE_METIS=YES -DMFEM_USE_METIS_5=YES
  -DMFEM_USE_GSLIB=YES -DMFEM_USE_LAPACK=YES
  -DHYPRE_DIR=$EXT/deps/hypre -DMETIS_DIR=$EXT/deps/metis -DGSLIB_DIR=$EXT/gslib/build
)

# ---------------------------------------------------------------- 2. upstream MFEM
mkdir -p build-mfem
cmake -S mfem -B build-mfem "${COMMON_CMAKE[@]}" \
  -DCMAKE_INSTALL_PREFIX=$EXT/install/mfem \
  -DMFEM_ENABLE_EXAMPLES=YES -DMFEM_ENABLE_MINIAPPS=YES \
  > logs/mfem-cmake-configure.log 2>&1
cmake --build build-mfem -j"$JOBS" --target mfem ex3 ex3p ex4 ex4p field-interp findpts pfindpts \
  > logs/mfem-build.log 2>&1
# 'cmake --install' with examples enabled requires *all* examples to be built, so
# toggle examples/miniapps off for the install, then back on (build tree stays usable).
cmake build-mfem -DMFEM_ENABLE_EXAMPLES=NO -DMFEM_ENABLE_MINIAPPS=NO >> logs/mfem-cmake-configure.log 2>&1
cmake --build build-mfem -j"$JOBS" > /dev/null
cmake --install build-mfem > logs/mfem-install.log 2>&1
cmake build-mfem -DMFEM_ENABLE_EXAMPLES=YES -DMFEM_ENABLE_MINIAPPS=YES >> logs/mfem-cmake-configure.log 2>&1

# ---------------------------------------------------------------- 3. reproduce upstream examples
M=$EXT/mfem/data
( cd build-mfem/examples
  mpirun -np 2 --oversubscribe ./ex3p -m $M/beam-hex.mesh -o 2 -no-vis > $EXT/logs/ex3p.log 2>&1
  mpirun -np 2 --oversubscribe ./ex4p -m $M/beam-hex.mesh -o 2 -no-vis > $EXT/logs/ex4p.log 2>&1
  ./ex3 -m $M/beam-hex.mesh -o 2 -no-vis > $EXT/logs/ex3.log 2>&1 )
( cd build-mfem/miniapps/gslib
  mpirun -np 2 --oversubscribe ./pfindpts -no-vis > $EXT/logs/pfindpts.log 2>&1
  ./field-interp -no-vis > $EXT/logs/field-interp.log 2>&1 )

# ---------------------------------------------------------------- 4. MFEM fork (for MHD-ALE)
mkdir -p build-mfem-develop
cmake -S mfem-develop -B build-mfem-develop "${COMMON_CMAKE[@]}" \
  -DCMAKE_INSTALL_PREFIX=$EXT/install/mfem-develop \
  > logs/mfem-develop-cmake-configure.log 2>&1
cmake --build build-mfem-develop -j"$JOBS" > logs/mfem-develop-build.log 2>&1
cmake --install build-mfem-develop > logs/mfem-develop-install.log 2>&1
# -> install/mfem-develop/share/mfem/config.mk and include/mfem/mfem.hpp

# ---------------------------------------------------------------- 5. MHD-ALE
( cd MHD-ALE
  make MFEM_INSTALL_DIR=$EXT/install/mfem-develop -j"$JOBS" > "$EXT/logs/mhd-ale-build.log" 2>&1 )

# ---------------------------------------------------------------- 6. MHD-ALE smoke tests
cd MHD-ALE
run_mhd() { # outdir nproc args...
  local od=$1 np=$2; shift 2; mkdir -p "$od"
  mpirun -np "$np" ./MHD "$@" -pv -pvb TaylorGreen -od "$od" > "$od/MHD.out" 2> "$od/MHD.err"
}
# 2-D Taylor-Green, refine 0, 1 proc (test_TaylorGreen.sh settings, -tf 0.1)
run_mhd output/TaylorGreen-Dim2-Refine0-tf0.1 1 \
  -p 0 -rs 0 -ptf plot_time/plot_time_TaylorGreen.dat -s 3 -m ./mesh/Disk-4x4-quad.mesh -ale \
  -tf 0.1 -mst 1 -ok 2 -ot 1 -or 1 -oa 3 -cfl 0.5 -dim 2 -g1 0.85 -g2 1.02 -fsr -fsri 10 \
  -rmv 1 -rme 1 -rma 2 -rmr 3 -cgt 1e-12 -bpt 1
# 2-D, exercising remap every 2 steps
run_mhd output/TaylorGreen-Dim2-Refine1-tf0.1-fsri2 1 \
  -p 0 -rs 1 -ptf plot_time/plot_time_TaylorGreen.dat -s 3 -m ./mesh/Disk-4x4-quad.mesh -ale \
  -tf 0.1 -mst 1 -ok 2 -ot 1 -or 1 -oa 3 -cfl 0.5 -dim 2 -g1 0.85 -g2 1.02 -fsr -fsri 2 \
  -rmv 1 -rme 1 -rma 2 -rmr 3 -cgt 1e-12 -bpt 1
# 3-D Taylor-Green: ND vector potential (oa 2), RT B, helicity-preserving remap (-rma 2), 2 procs
run_mhd output/TaylorGreen-Dim3-Refine1-tf0.2-fsri2 2 \
  -p 0 -dim 3 -m ./mesh/Cube-4x4x4-hex.mesh -rs 1 -tf 0.2 -oa 2 -ok 1 -ot 0 -or 0 -s 3 -ale -mst 1 \
  -fsr -fsri 2 -rmv 1 -rme 1 -rma 2 -rmr 3 -cgt 1e-12 -bpt 1 -ptf plot_time/plot_time_TaylorGreen.dat
echo "build_all.sh finished"
