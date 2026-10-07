# E6: MHD-ALE 3-D Taylor-Green, A-route remap (-rma 0 / -rma 2) vs frequent forced remaps

All numbers below are produced by `mhdale/analyze.py` from the instrumented MHD-ALE runs
(`external/MHD-ALE/output/E6*/<run>/{step_diag,remap_diag}.dat`, exact commands in `mhdale/all_run_commands.txt`).
Instrumentation: `mhdale/instrumentation.patch` (does not change the evolution: final errors, helicity.dat,
B_divergence_error.dat, time_step.dat, rho_min.dat byte-identical to the uninstrumented binary for `-rma 0` and `-rma 2`, 2 ranks).

## What the test problem is
* `exact_B` of `testcase_TaylorGreen.cpp` is **time independent**: B = sqrt(4 pi) beta (sin pi x cos pi y, -cos pi x sin pi y, 0), beta = 0.5, mu = 4 pi.
  The TG state (v, e, B) is a *steady* smooth equilibrium (energy source term + Lorentz force balance the pressure), so the "exact" solution at time t equals the initial condition.
  The error vs exact therefore measures drift of the numerical solution away from the steady state caused by mesh motion/remap/time stepping; there is no physical decay (ideal MHD, no resistivity).
* In 3-D, A = (0, 0, A_z), B in the xy plane, so **helicity int A.B is identically 0** (|H| < 1e-13 in every run). The helicity diagnostic is therefore uninformative for this test; no helicity drift can be inferred.
* Mesh smoother in 3-D: `-mst 1` = `MeshSmoothType::INITIAL`: the target mesh is always the *initial* Cartesian mesh. The Lagrangian mesh follows the flow (|v| up to 1), so the nodes
  move a lot between remaps: max|newnodes - nodes| / h_min = **0.31 per 2 steps** and **1.65 per 10 steps** (i.e. more than one cell width at -fsri 10). Mesh quality criteria
  (detJ in [0.4,5], ratio < 4) never triggered a remap before the forced one in any run (remap counts equal floor(steps/interval)).
* Without ALE (`-no-ale`) the Lagrangian mesh distorts, dt shrinks from 0.017 to below 0.008 and the run needs 66 steps for tf = 0.5.

## Runs (3-D, `-p 0 -dim 3 -m mesh/Cube-4x4x4-hex.mesh -rs 1`, `-ok 2 -ot 1 -or 1 -oa 3`, `-s 3 -cfl 0.5 -mst 1 -rmv 1 -rme 1 -rmr 3 -cgt 1e-12 -bpt 1`, 2 MPI ranks)
(A: 45000 ND dofs, B: 43200 RT dofs, 512 hexes.)  -rs 2 in 3-D was not run: one `-rma 2` remap costs ~100 s at -rs 1 on the shared machine.
The fsri-2 `-rma 2` run was cut to tf = 0.25 (a tf = 0.5 attempt was killed after 4 remaps, ~100 s each, partial output in `external/MHD-ALE/output/E6_killed`); the fsri-2 `-rma 0` and no-ALE runs exist for both tf = 0.25 and 0.5.
Definitions: dE_B = relative change of E_B = 0.5 int B^2/mu over the run; "remap sum (A)" = sum of per-remap jumps (E_after-E_before)/E_B(0) of the route actually used (A -> B = curl A);
"shadow sum" = same for the shadow direct-B transfer (RT field point-interpolated by GSLIB `InterpolateRemap`, not used in the evolution); the "mean abs E jump" columns are mean |E_after-E_before|/E_before per remap.
"mean rel. B-err change" = mean over remaps of (L2err(B) after - before)/before.  Max divB columns: L2 norm over the whole run (A-route, every step) and over remap events (shadow); "div dof" = max |dof| of DivergenceInterpolator(B).
Wall-time columns: remap overhead = (all remap work: rho, v, e, A, mesh smoothing) / (ODE steps + remaps), instrumentation excluded; A-remap frac = only the A remap.  The machine was shared with other jobs, so wall times are noisy (an ODE step took ~2 s when idle, 3.5-6.4 s during these runs).

| run | steps | #remaps | final L2 err B | dE_B total | remap sum (A) | shadow sum | dKE total | max divB L2 (A-route) | max divB L2 (shadow) | max abs div dof (A / shadow) | mean abs E jump A / shadow | mean rel. B-err change at remap A / shadow | remap overhead | A-remap frac | total wall [s] |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| fsri2_rma0_rs1_tf0.25 | 15 | 7 | 4.037e-03 | 5.423e-05 | 2.364e-05 | -2.562e-04 | -1.265e-03 | 4.646e-14 | 1.938e-02 | 5.2e-13 / 1.9e-01 | 3.4e-06 / 4.2e-05 | -0.05 / +0.06 | 0.670 | 0.011 | 341 |
| fsri2_rma2_rs1_tf0.25 | 15 | 7 | 3.990e-03 | 2.981e-05 | 3.990e-07 | 1.808e-05 | -1.275e-03 | 4.605e-14 | 3.034e-02 | 5.1e-13 / 3.1e-01 | 1.2e-07 / 9.4e-05 | -0.01 / +0.06 | 0.883 | 0.650 | 874 |
| noale_rs1_tf0.25 | 21 | 0 | 7.473e-03 | -4.312e-05 | n/a | n/a | 4.353e-05 | 4.114e-14 | n/a | 4.7e-13 / - | - | - | 0 | 0 | 191 |
| fsri10_rma0_rs1_tf0.5 | 37 | 3 | 7.792e-03 | 2.248e-04 | -8.580e-06 | 4.357e-04 | -9.202e-04 | 4.598e-14 | 6.889e-02 | 4.7e-13 / 5.4e-01 | 2.9e-06 / 1.5e-04 | -0.08 / -0.00 | 0.714 | 0.002 | 496 |
| fsri10_rma2_rs1_tf0.5 | 37 | 3 | 8.148e-03 | 2.326e-04 | 3.127e-06 | 4.960e-04 | -9.296e-04 | 4.572e-14 | 9.685e-02 | 4.7e-13 / 6.9e-01 | 1.0e-06 / 1.7e-04 | +0.03 / -0.00 | 0.862 | 0.720 | 1025 |
| fsri2_rma0_rs1_tf0.5 | 30 | 14 | 6.572e-03 | 6.864e-04 | 4.760e-05 | 4.369e-04 | -3.026e-03 | 4.669e-14 | 1.938e-02 | 5.2e-13 / 1.9e-01 | 3.4e-06 / 7.0e-05 | -0.03 / +0.07 | 0.694 | 0.012 | 613 |
| noale_rs1_tf0.5 | 66 | 0 | 2.006e-02 | 2.120e-04 | n/a | n/a | -4.706e-04 | 4.117e-14 | n/a | 4.7e-13 / - | - | - | 0 | 0 | 311 |

## Interpretation (3-D)
1. **Remap frequency.** Remapping matters a lot relative to the pure Lagrangian run (B error 2.0e-2 without ALE vs 6.6e-3 / 7.8e-3 with remaps every 2 / 10 steps at t = 0.5; at t = 0.25: 7.5e-3 vs 4.0e-3), because the Lagrangian mesh distorts. Between the two remap frequencies the effect on the B error is modest (more frequent is ~15 % better at t = 0.5).
   E_B itself changes only at the 1e-4 level in all runs (2e-4 no-ALE and fsri 10, 7e-4 at fsri 2 over t = 0.5) and the part due to the A-remaps is small: sum of jumps 4.8e-5 (rma 0, fsri 2) or <1e-5 (others) of E_B(0); most of the E_B change comes from the Lagrangian steps. The visible cost of frequent remapping is in the *other* fields: kinetic energy changes by -3.0e-3 (fsri 2) vs -9.2e-4 (fsri 10) vs -4.7e-4 (no ALE), i.e. the DG v/e remaps (not the A remap) dissipate KE.
2. **-rma 0 vs -rma 2.** Final B errors agree within a few percent (4.04e-3 vs 3.99e-3 at fsri 2/t = 0.25; 7.79e-3 vs 8.15e-3 at fsri 10/t = 0.5), i.e. no measurable accuracy difference at this resolution. HPRemap does give a much smaller E_B jump per remap (mean |jump| 1.2e-7 vs 3.4e-6 at fsri 2, 1.0e-6 vs 2.9e-6 at fsri 10; cumulative 4e-7 vs 2.4e-5 at t = 0.25) but this is small compared with the E_B drift from the Lagrangian steps (3e-5 to 5e-5 at t = 0.25), so it is not visible in E_B(t) or in the errors. Cost: the A remap takes 0.3-0.4 s with -rma 0 and 76 s (fsri 2) / 235 s (fsri 10, mesh motion 1.65 h) with -rma 2, i.e. ~200-700x more; the whole remap goes 27 s -> 104 s (fsri 2); total wall 341 s -> 874 s. Both keep div_h B at roundoff (L2 4.6e-14, max |dof| 5e-13) since B = curl_h A.
   Note that even with -rma 0 the remap overhead is 67-71 % of the wall time, but only ~1 % of it is the A remap: the non-A parts (rho L2 projection, DG velocity/energy remaps, mesh smoothing; not separately timed) dominate, ~25 s per remap on 2 ranks on the loaded machine.
3. **Shadow direct-B transfer vs A-route.** Interpolating the RT field directly (same GSLIB machinery, no evolution) gives div B of L2 norm 2e-2 to 1e-1 (max |dof| 0.19-0.69) instead of 5e-14 / 5e-13, i.e. ~11 orders of magnitude worse in absolute terms, but note that its size is comparable to the discretization error scale (|B| ~ 1.25, error 4e-3-8e-3; div is a derivative so its norm is larger). The energy jumps of the shadow route are 12-50x larger than the A-route ones (mean |jump| 4.2e-5 and 1.5e-4 vs 3.4e-6 and 2.9e-6 at fsri 2 and 10) and at fsri 10 they all have the same sign (+8e-5, +1.7e-4, +1.8e-4), so repeated use would give a secular E_B gain (sum 4.4e-4, twice the total E_B change of the actual run), whereas at fsri 2 the jumps alternate in sign (net -2.6e-4 over 7 remaps). The shadow B differs from the A-route B by 5e-4-2.7e-3 relative L2, the same order as the B discretization error, and its error vs exact right after the remap is ~6 % worse at fsri 2 (ratio shadow/before 1.06) while the A-route is 3-5 % better with -rma 0 (ratio 0.95-0.97; 1 % with -rma 2); remapping onto the initial mesh also removes mesh distortion error.
   Caveats: the shadow transfer is a one-shot comparison starting from the A-route's (div-free) state and is *not* compounded; it uses plain point interpolation of RT dof functionals (not flux-conserving, no div cleaning), so it represents the naive direct-B transfer, not the best possible one. Its cost is small (0.2-0.3 s) but it is also comparable with the rma-0 A transfer (shadow/A-remap time 0.6-0.7).
4. **Limits of this experiment.** One resolution (-rs 1), smooth steady problem with planar B (helicity trivially 0), at most 14 remaps, time horizon t = 0.5 (about 30-37 steps), and machine-shared noisy timings. Conclusions about secular drift over many remaps are suggestive only (e.g. fsri 10 shadow jumps are same-signed), not demonstrated.

## 2-D side note (scalar H1 potential, `-dim 2 -m mesh/Disk-4x4-quad.mesh -rs 2 -ok 2 -ot 1 -or 1 -oa 3`, same other options, 2 ranks, tf = 0.3 and 1.0; `-rma 0` = interpolation, `-rma 1` = DG convection)
| run | steps | #remaps | final L2 err B | dE_B total | remap sum (A) | shadow sum | dKE total | max divB L2 (A-route) | max divB L2 (shadow) | max abs div dof (A / shadow) | mean abs E jump A / shadow | mean rel. B-err change at remap A / shadow | remap overhead | A-remap frac | total wall [s] |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| fsri10_rma0_rs2_tf0.3 | 39 | 3 | 1.090e-03 | -6.444e-06 | -4.233e-07 | 2.444e-05 | -2.526e-05 | 2.323e-13 | 1.737e-02 | 1.6e-12 / 1.9e-01 | 1.4e-07 / 8.1e-06 | -0.10 / -0.00 | 0.660 | 0.001 | 21 |
| fsri10_rma1_rs2_tf0.3 | 39 | 3 | 9.686e-04 | -9.474e-06 | -2.376e-06 | 7.259e-06 | -2.477e-05 | 3.320e-13 | 1.648e-02 | 1.8e-12 / 1.7e-01 | 1.0e-06 / 8.5e-06 | -0.12 / -0.00 | 0.881 | 0.663 | 57 |
| fsri2_rma0_rs2_tf0.3 | 35 | 17 | 1.463e-03 | 1.016e-04 | 7.201e-06 | -2.231e-04 | -8.807e-04 | 3.549e-13 | 6.029e-03 | 2.5e-12 / 8.2e-02 | 4.2e-07 / 1.4e-05 | -0.03 / +0.04 | 0.715 | 0.006 | 24 |
| fsri2_rma1_rs2_tf0.3 | 35 | 17 | 1.471e-03 | 9.619e-05 | -3.106e-07 | -2.244e-04 | -8.821e-04 | 3.467e-13 | 7.154e-03 | 2.0e-12 / 9.3e-02 | 3.7e-07 / 1.4e-05 | -0.01 / +0.04 | 0.907 | 0.667 | 62 |
| noale_rs2_tf0.3 | 54 | 0 | 2.710e-03 | -2.166e-06 | n/a | n/a | -1.014e-06 | 2.212e-13 | n/a | 1.4e-12 / - | - | - | 0 | 0 | 13 |
| fsri10_rma0_rs2_tf1.0 | 130 | 12 | 2.963e-03 | -2.007e-05 | -6.365e-06 | 1.623e-04 | -7.277e-05 | 2.899e-13 | 2.647e-02 | 1.6e-12 / 3.0e-01 | 5.3e-07 / 1.4e-05 | -0.06 / +0.01 | 0.697 | 0.001 | 75 |
| fsri10_rma1_rs2_tf1.0 | 130 | 12 | 2.535e-03 | -1.441e-04 | -8.825e-05 | 5.107e-06 | -7.215e-05 | 3.320e-13 | 2.405e-02 | 1.8e-12 / 3.1e-01 | 7.4e-06 / 9.8e-06 | -0.07 / +0.00 | 0.894 | 0.655 | 88 |
| fsri2_rma0_rs2_tf1.0 | 116 | 57 | 3.993e-03 | 3.611e-04 | 2.412e-05 | 5.166e-04 | -2.772e-03 | 3.549e-13 | 6.029e-03 | 2.5e-12 / 8.2e-02 | 4.2e-07 / 1.7e-05 | -0.01 / +0.02 | 0.716 | 0.005 | 26 |
| fsri2_rma1_rs2_tf1.0 | 116 | 57 | 3.665e-03 | 2.251e-04 | -6.026e-05 | -7.875e-05 | -2.758e-03 | 3.467e-13 | 7.154e-03 | 2.0e-12 / 9.3e-02 | 1.2e-06 / 6.8e-06 | +0.00 / +0.02 | 0.907 | 0.671 | 127 |
| noale_rs2_tf1.0 | 788 | 0 | 5.448e-02 | 3.934e-04 | n/a | n/a | -5.459e-03 | 2.215e-13 | n/a | 1.4e-12 / - | - | - | 0 | 0 | 164 |

* Same qualitative picture. B error at t = 1: no ALE 5.4e-2 (788 steps), remap/10 steps 3.0e-3 (rma 0) / 2.5e-3 (rma 1), remap/2 steps 4.0e-3 / 3.7e-3, i.e. in 2-D more frequent remapping is *not* better (3-D showed a modest gain from more frequent remapping) and the DG remap is 8-14 % more accurate than interpolation at up to 5x the total cost (A remap ~66 % of the wall time, total 127 s vs 26 s at fsri 2; 88 s vs 75 s at fsri 10).
  The remap-induced E_B jump of the DG remap is *larger* than that of interpolation (mean |jump| 1.2e-6 vs 4.2e-7 at fsri 2; 7.4e-6 vs 5.3e-7 at fsri 10).
  The A-route keeps max |div dof| at 1e-12 while the shadow RT interpolation gives L2 div norm 6e-3-3e-2 and E_B jumps 10-30x larger than interpolating A with -rma 0 (1.4e-5 vs 4.2e-7 at fsri 2).
  (The `-no-ale` 2-D run has 788 steps because dt collapses as the Lagrangian mesh distorts.)

## Appendix: full auto-generated tables
### 3d runs (output root `external/MHD-ALE/output/E6`)

E_B = 0.5 int B^2/mu (mu = 4 pi). dE_B total = (E_B(end)-E_B(0))/E_B(0). 'remap sum' = sum over remaps of (E_after-E_before)/E_B(0) for the A-route (what the code does) and 'shadow sum' the same for the direct-B (RT GSLIB interpolation) shadow transfer, which was NOT used for evolution. max divB columns: L2 norm of div B over the whole run for the A-route (every step, from step_diag.dat), and the max over remap events of the L2 norm / max|dof| of div B of the shadow B. 'remap overhead' = sum of (A + rho + v + e + mesh) remap wall time / (sum ODE-step + remap wall time), instrumentation excluded; 'A-remap' = only the A-remap part.

| run | steps | #remaps | final L2 err B | final L2 err v | final L2 err e | dE_B total | remap sum (A) | shadow sum | dKE total | max abs helicity | max divB L2 (A-route) | max divB L2 (shadow) | max abs div dof (A / shadow) | mean L2diff(shadow,A)/normB | mean abs E jump A / shadow | mean rel. B-err change at remap A / shadow | remap overhead | A-remap frac | shadow cost / A-remap | total wall [s] |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| fsri2_rma0_rs1_tf0.25 | 15 | 7 | 4.037e-03 | 2.397e-03 | 8.568e-03 | 5.423e-05 | 2.364e-05 | -2.562e-04 | -1.265e-03 | 2.347e-15 | 4.646e-14 | 1.938e-02 | 5.2e-13 / 1.9e-01 | 5.082e-04 | 3.4e-06 / 4.2e-05 | -0.05 / +0.06 | 0.670 | 0.011 | 0.726 | 341 |
| fsri2_rma2_rs1_tf0.25 | 15 | 7 | 3.990e-03 | 2.424e-03 | 8.567e-03 | 2.981e-05 | 3.990e-07 | 1.808e-05 | -1.275e-03 | 1.335e-17 | 4.605e-14 | 3.034e-02 | 5.1e-13 / 3.1e-01 | 9.158e-04 | 1.2e-07 / 9.4e-05 | -0.01 / +0.06 | 0.883 | 0.650 | 0.004 | 874 |
| noale_rs1_tf0.25 | 21 | 0 | 7.473e-03 | 2.447e-03 | 1.273e-02 | -4.312e-05 | n/a | n/a | 4.353e-05 | 1.335e-17 | 4.114e-14 | n/a | 4.7e-13 / - | n/a | - | - | 0 | 0 | n/a | 191 |
| fsri10_rma0_rs1_tf0.5 | 37 | 3 | 7.792e-03 | 3.707e-03 | 1.172e-02 | 2.248e-04 | -8.580e-06 | 4.357e-04 | -9.202e-04 | 3.911e-15 | 4.598e-14 | 6.889e-02 | 4.7e-13 / 5.4e-01 | 1.273e-03 | 2.9e-06 / 1.5e-04 | -0.08 / -0.00 | 0.714 | 0.002 | 0.578 | 496 |
| fsri10_rma2_rs1_tf0.5 | 37 | 3 | 8.148e-03 | 3.855e-03 | 1.175e-02 | 2.326e-04 | 3.127e-06 | 4.960e-04 | -9.296e-04 | 1.489e-17 | 4.572e-14 | 9.685e-02 | 4.7e-13 / 6.9e-01 | 2.702e-03 | 1.0e-06 / 1.7e-04 | +0.03 / -0.00 | 0.862 | 0.720 | 0.001 | 1025 |
| fsri2_rma0_rs1_tf0.5 | 30 | 14 | 6.572e-03 | 3.472e-03 | 9.058e-03 | 6.864e-04 | 4.760e-05 | 4.369e-04 | -3.026e-03 | 1.634e-13 | 4.669e-14 | 1.938e-02 | 5.2e-13 / 1.9e-01 | 6.757e-04 | 3.4e-06 / 7.0e-05 | -0.03 / +0.07 | 0.694 | 0.012 | 0.638 | 613 |
| noale_rs1_tf0.5 | 66 | 0 | 2.006e-02 | 9.253e-03 | 1.729e-02 | 2.120e-04 | n/a | n/a | -4.706e-04 | 1.335e-17 | 4.117e-14 | n/a | 4.7e-13 / - | n/a | - | - | 0 | 0 | n/a | 311 |


### 2d runs (output root `external/MHD-ALE/output/E6_2d`)

E_B = 0.5 int B^2/mu (mu = 4 pi). dE_B total = (E_B(end)-E_B(0))/E_B(0). 'remap sum' = sum over remaps of (E_after-E_before)/E_B(0) for the A-route (what the code does) and 'shadow sum' the same for the direct-B (RT GSLIB interpolation) shadow transfer, which was NOT used for evolution. max divB columns: L2 norm of div B over the whole run for the A-route (every step, from step_diag.dat), and the max over remap events of the L2 norm / max|dof| of div B of the shadow B. 'remap overhead' = sum of (A + rho + v + e + mesh) remap wall time / (sum ODE-step + remap wall time), instrumentation excluded; 'A-remap' = only the A-remap part.

| run | steps | #remaps | final L2 err B | final L2 err v | final L2 err e | dE_B total | remap sum (A) | shadow sum | dKE total | max abs helicity | max divB L2 (A-route) | max divB L2 (shadow) | max abs div dof (A / shadow) | mean L2diff(shadow,A)/normB | mean abs E jump A / shadow | mean rel. B-err change at remap A / shadow | remap overhead | A-remap frac | shadow cost / A-remap | total wall [s] |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| fsri10_rma0_rs2_tf0.3 | 39 | 3 | 1.090e-03 | 5.438e-04 | 2.168e-03 | -6.444e-06 | -4.233e-07 | 2.444e-05 | -2.526e-05 | 0.000e+00 | 2.323e-13 | 1.737e-02 | 1.6e-12 / 1.9e-01 | 1.751e-04 | 1.4e-07 / 8.1e-06 | -0.10 / -0.00 | 0.660 | 0.001 | 1.642 | 21 |
| fsri10_rma1_rs2_tf0.3 | 39 | 3 | 9.686e-04 | 5.514e-04 | 2.168e-03 | -9.474e-06 | -2.376e-06 | 7.259e-06 | -2.477e-05 | 0.000e+00 | 3.320e-13 | 1.648e-02 | 1.8e-12 / 1.7e-01 | 2.501e-04 | 1.0e-06 / 8.5e-06 | -0.12 / -0.00 | 0.881 | 0.663 | 0.001 | 57 |
| fsri2_rma0_rs2_tf0.3 | 35 | 17 | 1.463e-03 | 8.260e-04 | 2.221e-03 | 1.016e-04 | 7.201e-06 | -2.231e-04 | -8.807e-04 | 0.000e+00 | 3.549e-13 | 6.029e-03 | 2.5e-12 / 8.2e-02 | 1.038e-04 | 4.2e-07 / 1.4e-05 | -0.03 / +0.04 | 0.715 | 0.006 | 1.824 | 24 |
| fsri2_rma1_rs2_tf0.3 | 35 | 17 | 1.471e-03 | 8.361e-04 | 2.222e-03 | 9.619e-05 | -3.106e-07 | -2.244e-04 | -8.821e-04 | 0.000e+00 | 3.467e-13 | 7.154e-03 | 2.0e-12 / 9.3e-02 | 1.132e-04 | 3.7e-07 / 1.4e-05 | -0.01 / +0.04 | 0.907 | 0.667 | 0.006 | 62 |
| noale_rs2_tf0.3 | 54 | 0 | 2.710e-03 | 7.138e-04 | 4.113e-03 | -2.166e-06 | n/a | n/a | -1.014e-06 | 0.000e+00 | 2.212e-13 | n/a | 1.4e-12 / - | n/a | - | - | 0 | 0 | n/a | 13 |
| fsri10_rma0_rs2_tf1.0 | 130 | 12 | 2.963e-03 | 1.038e-03 | 2.773e-03 | -2.007e-05 | -6.365e-06 | 1.623e-04 | -7.277e-05 | 0.000e+00 | 2.899e-13 | 2.647e-02 | 1.6e-12 / 3.0e-01 | 2.791e-04 | 5.3e-07 / 1.4e-05 | -0.06 / +0.01 | 0.697 | 0.001 | 1.953 | 75 |
| fsri10_rma1_rs2_tf1.0 | 130 | 12 | 2.535e-03 | 1.010e-03 | 2.776e-03 | -1.441e-04 | -8.825e-05 | 5.107e-06 | -7.215e-05 | 0.000e+00 | 3.320e-13 | 2.405e-02 | 1.8e-12 / 3.1e-01 | 3.486e-04 | 7.4e-06 / 9.8e-06 | -0.07 / +0.00 | 0.894 | 0.655 | 0.001 | 88 |
| fsri2_rma0_rs2_tf1.0 | 116 | 57 | 3.993e-03 | 1.094e-03 | 2.821e-03 | 3.611e-04 | 2.412e-05 | 5.166e-04 | -2.772e-03 | 0.000e+00 | 3.549e-13 | 6.029e-03 | 2.5e-12 / 8.2e-02 | 1.339e-04 | 4.2e-07 / 1.7e-05 | -0.01 / +0.02 | 0.716 | 0.005 | 2.243 | 26 |
| fsri2_rma1_rs2_tf1.0 | 116 | 57 | 3.665e-03 | 1.107e-03 | 2.824e-03 | 2.251e-04 | -6.026e-05 | -7.875e-05 | -2.758e-03 | 0.000e+00 | 3.467e-13 | 7.154e-03 | 2.0e-12 / 9.3e-02 | 1.339e-04 | 1.2e-06 / 6.8e-06 | +0.00 / +0.02 | 0.907 | 0.671 | 0.006 | 127 |
| noale_rs2_tf1.0 | 788 | 0 | 5.448e-02 | 4.305e-02 | 5.801e-02 | 3.934e-04 | n/a | n/a | -5.459e-03 | 0.000e+00 | 2.215e-13 | n/a | 1.4e-12 / - | n/a | - | - | 0 | 0 | n/a | 164 |

