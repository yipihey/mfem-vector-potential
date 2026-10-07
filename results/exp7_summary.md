# E7: performance and scalability of one remap

One remap = M1 (f1, eps 0.3, geometric order q = p) -> M2 (uniform), abc field, b0 on, periodic unit box, N^3 hexahedra, ND order p (RT_{p-1}, L2_{p-1}, H1_p).
Driver `experiments/exp7_perf.cpp`, matrix `scripts/run_exp7.sh` (stages 1, 2, 3), analysis `scripts/plot_exp7.py`.
Raw data: `results/exp7_perf.csv` (223 rows, one per run and repeat), `results/exp7_baseline.csv` (mesh/space/operator construction only, 32 rows), logs `results/exp7_logs/`.
All tables (7.1-7.10) are in `results/exp7_tables.md`; figures `results/exp7_time_vs_nd.png`, `exp7_fractions.png`, `exp7_curl_memory.png`, `exp7_strong_scaling.png`, `exp7_cost_vs_p.png`.

## Method

* Every stage is timed with `MPI_Wtime` between barriers (max over ranks): mesh generation + deformation, Spaces, Operators (G, C, D), `RemoteEvaluator` (GSLIB setup on M1), the transfer plan (collect target points, FindPoints incl. periodic-image retries, mass operator / AMG setup), `Apply` (source evaluation at the located points, dof functionals / rhs, CG solves), the gauge fix, 10 x `C_h a`, and `b = b0 + C a`.
* "one-shot" = RemoteEvaluator setup + plan + apply + gauge apply + compose (mesh/space/operator build and the gauge object setup excluded, see below); "steady" = plan reused (eval + dofs + solves + gauge apply + compose; this is the E4 convention).
  One untimed warm-up apply absorbs lazy AMG setup; 2 timed applies; the whole process was run twice (tags r1, r2; r2 skipped when r1 took > 300 s, which happened for only two configurations: B_l2c p=2 N=32 and A_int+Coulomb p=4 N=12) and the minimum is reported. Repeat-to-repeat spread: median 3.8 %, max 27 % (small runs); at 4 ranks median 2.4 %.
* Memory: VmRSS / VmHWM from `/proc/self/status` (max and sum over ranks). **Peak RSS is the peak of the whole process** and therefore includes both meshes, both sets of spaces and operators, the exact source projection (initial data) and the GSLIB structures; `baseline` rows give the same without any transfer so the increment can be formed. Diagnostics (rel. error of B vs the exact field, max |D_h b|) are computed after the last memory snapshot and are outside all timings.
* Scheduling: `exp5_ale` was idle at the start of every stage (2 consecutive checks 60 s apart; stage 1/2 and stage 3 each waited 120 s). Nothing else was run by me during the matrix. I cannot exclude short activity of the other agent during a run; the small repeat spreads (above) and the 4-rank efficiencies (75-95 %) suggest the runs were mostly clean.
* All points were found (0 not found in every run); 3-9 % of the target points needed a periodic-image retry (Table 7.8). Sanity: rel. error of B after the remap 2.9e-2 / 7.3e-3 / 3.3e-3 / 1.9e-3 for p=2, N = 8/16/24/32 (A_int, order 2), div_h b at roundoff for A ops and B_l2c (1e-12 .. 1e-11), 0.03 - 26 (raw max, not scaled) for B_l2.
* Bounds of the study: one mesh pair, one field, one machine (4 cores, 15 GB, 1 rank = 1 core, OMP_NUM_THREADS = 1), MFEM + hypre + GSLIB as built here; absolute numbers are for this build (-O2).

## Compact results (1 rank; seconds: one-shot / steady; peak RSS of the process)

| p | N | ND dofs | A_pt | A_int | A_int+Coulomb | A_int+Jacobi:20 | A_l2+Coulomb | B_l2 | B_l2c | peak RSS A_int / B_l2c |
|---|---|---|---|---|---|---|---|---|---|---|
| 2 | 16 | 98 304 | 1.01 / 0.28 | 3.99 / 1.09 | 4.50 / 1.42 | 4.29 / 1.16 | 5.12 / 2.17 | 4.42 / 1.95 | 12.8 / 10.1 | 282 / 403 MB |
| 2 | 24 | 331 776 | 3.76 / 1.00 | 15.0 / 3.72 | 15.0 / 4.77 | 14.7 / 3.92 | 17.2 / 7.69 | 16.0 / 6.49 | 63.0 / 53.5 | 892 / 1300 MB |
| 2 | 32 | 786 432 | 9.02 / 2.41 | 37.0 / 9.28 | 40.4 / 13.7 | 36.8 / 10.0 | 48.3 / 20.8 | 40.2 / 13.5 | 191 / 164 | 2060 / 3020 MB |
| 3 | 16 | 331 776 | 4.40 / 1.44 | 24.5 / 8.15 | 26.7 / 11.3 | 25.2 / 8.50 | 21.5 / 11.1 | 16.4 / 7.13 | 58.2 / 46.9 | 871 / 1040 MB |
| 4 | 12 | 331 776 | 7.23 / 2.32 | 43.8 / 15.7 | 51.2 / 22.8 | 43.8 / 16.3 | 32.1 / 18.4 | 24.8 / 10.7 | 58.3 / 44.9 | 966 / 966 MB |
| 1 | 32 | 98 304 | 1.49 / 0.27 | 3.49 / 0.62 | 3.52 / 0.78 | 3.52 / 0.70 | 8.35 / 3.84 | 7.70 / 2.52 | 43.1 / 37.8 | 395 / 805 MB |

Mesh + spaces + operators build (both meshes) is small: 0.35 s (p=2 N=16), 3.6 s (N=32), 3.1 s (p=4 N=12). Peak RSS of that baseline alone: 120 MB (p=2 N=16), 797 MB (p=2 N=32), 401 MB (p=4 N=12) = about 1.0 kB per ND dof for two meshes (p=2).

## Scaling slopes (log-log of time vs ND dofs, 1 rank, N varied at fixed p; Table 7.1, 7.9)

| op | slope (p=2, 4 sizes, 8 .. 32) | range over p=1..4 (2-point fits for p != 2) |
|---|---|---|
| A_pt, A_int, A_int+Jacobi, B_l2 | 0.98 - 1.03 | 0.93 - 1.09 |
| A_int+Coulomb, A_l2+Coulomb | 1.03 - 1.07 (steady 1.07 / 1.06) | 0.97 - 1.17 |
| B_l2c | **1.24** (steady 1.30) | 1.24 - 1.47 |
| curl apply C_h a | time per ND dof 3.2 / 3.9 / 3.9 / 4.5 ns for N = 8/16/24/32 at p=2 (about 1.0 - 1.1) | 2.0 ns (p=1) .. 7.5 ns (p=4) per dof |
| memory (peak RSS) | linear: 2.6 kB per ND dof (A_int), 3.8-3.9 (B_l2c); of which 1.0 is the baseline (two meshes + spaces + operators); p=2 | |

Cost per ND dof vs p at fixed N (Table 7.7, µs per dof, one-shot, N=8, p=2/3/4): A_pt 12.6/14.7/20.2, A_int 41.7/67.1/138, B_l2 44.5/48.4/80.7, A_l2+Coulomb 48.6/57.1/89.6, B_l2c 88.8/97.5/130. At N=16, p=1/2/3: A_int 37.5/40.6/73.9, B_l2 77.3/45.0/49.3, B_l2c 191/131/175. Per point, FindPoints costs 2.0-2.3 µs (p=1), 2.8-3.4 µs (p=2), 4.4 µs (p=3), 7.4-8.0 µs (p=4) (Newton on higher-order elements) and the number of points per ND dof is roughly constant (9 for A_int at p=2..4, 12 at p=1; the quadrature of the integrated dofs, nq = p+2, is the multiplier).

## Where the time goes (Table 7.2, one-shot, 1 rank)

* A-route with integrated dofs (A_int): FindPoints 53-65 % (+ evaluator setup), source evaluation (collect + interpolate) 25-39 %, dof functionals 5-8 %, gauge 0 - 11 %, `b = b0 + C a` below 0.1 %. These are the **mesh-to-mesh costs shared by both routes** and they dominate everything except B_l2c.
* L2 variants add the mass solve (1-26 % of the remap: 7-10 % at p=2, 26 % for A_l2 at p=1 N=32, 1-4 % at p=3..4) and the rhs assembly (13-18 %): A_l2/B_l2 cost about the same as A_int (4.4-5.1 s vs 4.0 s at p=2 N=16).
* B_l2c: the div-clean solve is 69 % (p=2 N=16), 77-79 % (N=24/32), 70 % (p=3 N=16) of its time; it costs 3.2x (N=16) to 5.2x (N=32) of A_int one-shot at p=2 (steady, plan reused: 164 s vs 9.3 s, i.e. 18x at N=32).
* Gauge fix, apply only: Coulomb 4-16 % of the one-shot time (p=2: 0.225 s N=16, 3.7 s N=32; p=3 N=16 3.0 s; p=4 N=12 6.8 s; 13-23 AMG-CG iterations), Jacobi:20 1-2 % (0.06 s .. 0.76 s at p=2; p=4 N=12 0.77 s), i.e. 4-9x cheaper than Coulomb.
* **Caveat (Table 7.10)**: the gauge object is built per mesh (ND mass assembly, L = G^T M G triple product, AMG): 0.14 s .. 45 s (p=2 N=16..32), 55 s (p=3 N=16), 165 s (p=4 N=12), i.e. 2.6-10x the whole steady remap including the gauge apply; it is not in the totals above. On a moving mesh it would be paid every remap (Coulomb and, as implemented, Jacobi, which uses the same L), and it dominates the peak memory (A_int+Coulomb p=2: 7.1 kB per dof, p=3: 11 kB, p=4: 18 kB, versus 2.6, 2.6, 2.9 kB without gauge). This is an inefficiency of this implementation (L equals the H1 stiffness matrix; the Jacobi variant could be matrix-free), not an intrinsic cost, but it was not measured otherwise.

## Iteration counts (1 rank, Table 7.5)

| quantity | p=1 | p=2 | p=3 | p=4 |
|---|---|---|---|---|
| L2 projection, ND mass (A_l2), Jacobi-CG to 1e-12 | 29 (N=16), 30 (32) | 18, 19, 19, 18 (N=8,16,24,32) | 15, 15 (N=8,16) | 13, 12 (N=8,12) |
| L2 projection, RT mass (B_l2) | 9, 17 (N=16,32) | 6, 10, 12, 12 | 6, 10 | 6, 8 |
| B_l2c outer (Schur, AMG-PCG) / inner (M_RT CG total) | 28/232, 32/528 | 29/180, 32/330, 32/462, 34/630 | 39/240, 58/590 | 53/324, 69/560 |
| Coulomb AMG-CG | 13, 16 | 14, 16, 16, 16 | 21, 22 | 22, 23 |

The mass solves are bounded (saturating, independent of N; they get smaller with p because Jacobi-CG on the uniform target mesh sees a narrower spectrum, a property of this test, not a statement about distorted meshes). The Coulomb solve saturates at 16 (p=2), 22 (p=3), 23 (p=4) iterations. The Schur outer count is flat in N at p=1,2 (28-34) but grows with N at p=3 (39 -> 58) and p=4 (53 -> 69) and with p; the inner count per outer iteration grows with N (6-12 for the RT mass solve at 1e-14), so total B_l2c work scales as N^3 x N^0.3 ... N^0.5 (slopes 1.24-1.47).

## Strong scaling (Table 7.6, 4 cores, ranks oversubscribed together with the harness; indicative)

One-shot speedup 1 -> 2 -> 4 ranks: p=2 N=16: 1.6-2.0 -> 3.0-3.7 (A_pt 1.62 / 2.96, A_int 1.77 / 3.09, A_int+Coulomb 2.01 / 3.42, A_l2+Coulomb 1.83 / 3.29, B_l2 1.87 / 3.68, B_l2c 1.80 / 3.29); p=3 N=16: 1.7-1.9 -> 3.1-3.6 (B_l2c 1.82 / 3.56); p=2 N=24: 1.7-1.9 -> 3.2-3.6 (B_l2c 1.93 / 3.48). FindPoints scales as well as the rest (S(4) = 3.0-4.0), the Coulomb apply 2.3-3.3 (p=2) and 4.1-5.0 (p=3, superlinear, cache effects), the curl apply 2.5-3.4, the solves 2.7-4.4. Peak memory per rank drops 403 -> 164 MB (p=2 N=16, B_l2c) from 1 to 4 ranks, summed over ranks it grows 403 -> 636 MB (ghost layers, per-rank copies of the GSLIB hash / setup). 4 ranks use all 4 cores plus the launcher, so efficiencies of 75-95 % are an upper bound for what a real distributed run would give; communication latency (shared memory here) is not probed.

## Narrative

The A-route's own operations are local: the dof functionals (4-8 % of a remap), the source interpolation, the curl apply (2-8 ns per dof, 1.3 ms for 3.3e5 dofs at p=2), the composition b = b0 + C a (below 0.1 %) and the Jacobi gauge smoothing (20 mat-vec sweeps, 1-2 % of a remap, but it only damps the rough part of the gauge, E4) need no global communication beyond a halo exchange. What dominates every operator except B_l2c is the mesh-to-mesh machinery shared by both routes, GSLIB FindPoints (53-65 % for A_int, 2-8 us per point, 9-12 points per ND dof, 3-9 % of the points needing a periodic-image retry) plus evaluating the source at the points (25-39 %); it is linear in the number of dofs (slope 0.98-1.03), scales with ranks like everything else (S(4) = 3.0-4.0 here) and grows with p (A_int cost per dof x3.3 from p=2 to p=4) through more points per element and costlier Newton searches. The global solves are: the L2 projections (Jacobi-PCG on a well-conditioned mass matrix, bounded and N-independent iteration counts 6-30, 1-26 % of the remap), the Coulomb gauge (Poisson with AMG-CG, 13-23 iterations independent of N, 4-16 % of the remap, but a per-mesh setup of 2.6-10 steady remaps in this implementation and 7-18 kB per dof of memory), and the exact div-clean of B_l2c (AMG-preconditioned outer CG on the Schur complement with tight inner M_RT solves: 28-69 outer and 180-630 inner iterations, 70-80 % of its time). B_l2c is 3.2x (N=16) to 5.2x (N=32) the cost of A_int one-shot and 18x at N=32 with plans reused, and it is the only operator with a super-linear slope (1.24-1.47), because the outer count grows with N at p >= 3 (39 -> 58, 53 -> 69) and the inner count per outer step grows with N. The discrete operators are cheap in memory (nnz C_h = 2(p+1) per ND dof, nnz D_h = p+1 per RT dof, independent of N), and a whole remap peaks at 2.6-3.0 kB per ND dof for A_int (of which 1.0 kB is the two meshes + spaces + operators) and 3.8-3.9 kB for B_l2c. Extrapolating linearly to 1e9 dofs at p=2, A_int needs about 47 us per dof one-shot (4.7e4 core-seconds, about 0.5 s on 1e5 cores at the 75-95 % efficiency seen on 4 cores, which says nothing yet about a distributed GSLIB search) and 2.6 TB; B_l2c needs 243 us per dof already at 8e5 dofs, probably 5-6x more per dof at 1e9 if the N^1.25 slope persists, 3.8 TB, and hundreds of dependent global inner iterations. A cosmological-scale code could therefore accept the local operations, the integrated-dof transfer, the L2 mass projections and Jacobi gauge smoothing, but would hesitate at the Coulomb gauge (unless L is assembled cheaply or applied matrix-free) and above all at the Schur-complement div-clean at every remap. Caveats on reliability: no other job was seen during the matrix, the repeat-to-repeat spread is 2-4 % (max 27 % on sub-second runs), the 4-rank runs use all four cores and shared-memory MPI (no communication latency), p=3/p=4 slopes are 2-point fits, and the B_l2c slope was measured only up to 8e5 dofs.

## Deviations and not done

* Peak RSS includes initial-data projection and both meshes; the "increment" columns in Table 7.3 subtract the baseline run (separate process).
* Gauge setup is reported separately (Table 7.10) and not added to the remap totals.
* `A_pt` was run for reference only; B_pt, B_int, A_l2 without gauge were not run (A_l2 shares everything with A_l2+Coulomb apart from the gauge).
* The N=32 (p=2) and N=12 (p=4) series have one repeat only where the first run took > 300 s (B_l2c, A_l2+Coulomb, A_int+Coulomb/Jacobi); in those cases "min of 2" reduces to the 2 timed applies inside one process.
* p=3 and p=4 slopes are 2-point fits.
* `src/` and `external/` untouched; no commit.
