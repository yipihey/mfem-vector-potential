# E5 summary: smooth ALE cycle with rezones (frozen-in field, material motion x -> x + eps sin^2(pi t/T) f1(x))

Data: `results/exp5_ale.csv` (one row per logged step; at a rezone step two rows: the state on the deformed source mesh *before* the transfer (`rezone=0`) and on the new uniform mesh *after* it (`rezone=1`); `kind=final` rows carry the summary quantities),
`results/exp5_floors.csv` (operator-independent projection floors at the rezone times), per-run parts/logs `results/exp5_parts/`, `results/exp5_logs/`.
Driver `experiments/exp5_ale.cpp` (+ `experiments/exp5_common.hpp`), runs `scripts/run_exp5.sh`, `scripts/run_exp5_floors.sh`, figures `results/exp5_*.png` (`scripts/plot_exp5.py`), full tables `results/exp5_tables.md` (`scripts/make_tables5.py`).
Self tests: `bin/exp5_ale -selftest -N 8 -p 2 -q 2` (logs `results/exp5_logs/selftest_p*_q*.log`): all PASS.
Timing caveat as in E3/E4: 1 MPI rank per job, 2 jobs at once on a 4-core VM shared with another agent (load average 3-4); wall times are inflated by an unknown factor, only ratios mean something.

## Set-up (what exactly was run)

* Periodic unit box, N = 8, field `abc` (Beltrami, helicity != 0) unless `B0` is on, f = `DisplacementField` variant 1, T = 1, K = 100 uniform steps, eps = 0.3 and 0.6 (min detJ of the Lagrangian mesh 0.90 and 0.60), p = 2 and 3 with **geometric order q = 2 for both** (q = 3 control for p = 3, see below).
* State: A-route `a = Pi_ND^int(A_exact)` (b = b0 + C_h a recomputed after every transfer), B-route `b = b0 + C_h a`, i.e. the same discrete field. During the Lagrangian phase only the nodes move: nodes = X + eps s(t) f(X) per node dof from the material labels X (label vectors kept per mesh);
  `C_h`, `G_h` are bitwise independent of the node positions and `D_h C_h = 0` to roundoff (selftest iv; note `D_h` itself is *not* geometry independent: it maps to the L2 *point values* of div, so it is a row scaling `diag(detJ_0/detJ(x_i))` of the reference matrix with the same kernel. I therefore evaluate the divergence at quadrature points as `div_ref/detJ_current`, verified against `GetDivergence` to 1e-12.)
* Rezone r of R: at t_r = r T/R (the last one at t = T) a fresh uniform mesh with labels `X_u = Phi_{t_r}^{-1}(y)` (Newton, 1e-14, mean 3-4 iterations, residual < 1e-14 verified at every rezone) is built, the state is moved to it with the chosen `vp::Transfer` plan (GSLIB point location on the current deformed mesh, rebuilt at every rezone and included in the reported transfer time), nodes continue from the labels `X_u`.
* Exact field: Cauchy push-forward `B = F B0(X)/det F`, `A = F^{-T} A0(X)` with `X = Phi_t^{-1}(x)` obtained by Newton **at the physical quadrature point** (deviation from the task text, which suggested interpolating a label grid function: with q = 2 the interpolated labels differ from `Phi^{-1}(x_h)` by up to 9e-4 at s = 1, eps = 0.6, which is the size of the p = 3 error, selftest iv). Relative errors are `||b_h - B_pf||/||B_pf||` (L2, quadrature order 2p+q+1; orders 7..10 change errB by < 1e-3 relative).
* Energy, helicity and the divergence are quadrature sweeps on the moving mesh (`E = 1/2 int |b_h|^2`, `H = int a_h . b_h`) instead of assembled mass matrices (same rule for the exact and discrete energies). Slice fluxes: `vp::SliceFlux` of all 6 slice families.
* "dB" below = `||b_T - b_0||/||b_0||`, both on the uniform mesh at T (dof vectors of the same topology), b_0 = initial discrete state (not the exact field).
* B0 = (0.3,-0.2,0.5) series: B-route only (see "A-route and a mean field" below).

Run matrix (thinned for the ~2 h budget, all done): R = 0 (Lagrangian reference; all operators are identical, run once with A_int), R = 2, 5, 20 for p = 2, 3 x eps = 0.3, 0.6 and the operators A_int (no gauge), A_int + Coulomb, A_l2 + Coulomb, B_int, B_l2, B_l2c; **R = 1** only for p = 2, eps = 0.3 (a rezone at t = T from the unmoved uniform mesh is an exact identity transfer: dB = 1e-15..2e-14, as it must be; I added **R = 2** (t = T/2 and T) as the smallest case with genuine transfers);
**B_int at p = 3, eps = 0.6, R = 20 was not run** (17 min per run); A_pt only p = 2, eps = 0.3, R = 1, 5, 20; A_l2 without gauge fixing only p = 2, eps = 0.3 (R = 1, 5, 20). p = 3 runs are logged every 2nd step (every 5th at R = 20; all rezone steps pre/post and t = T always), p = 2 every step.
Extras beyond the task: `-gauge coulomb_pre` (Coulomb fix on the *deformed source* mesh before the transfer; A_int/A_l2, p = 2 all, p = 3 eps = 0.6 R = 5, 20), a q = 3 control (p = 3, eps = 0.6, R = 5; and Lagrangian p = 3 q = 3), and the projection-floor runs. In r = 5 the rezone at t = 0.6 is an exact identity (s(0.6) = s(0.4), the mesh is uniform again), so 4 transfers are genuine; r = 2 and r = 20 have none such.

## Headline

**Lagrangian motion alone (R = 0).** a and b are constant: ||a_T - a_0||_inf = ||b_T - b_0||_inf = 0 (by construction, nothing is touched), and the invariants computed on the *moving* geometry stay at roundoff: max |H/H0 - 1| = 3e-16, max slice-flux error 7e-16 (4e-16 p = 3), max scaled div 1.2e-15 (3.2e-15 p = 3), for eps up to 0.6 (min detJ 0.60) -- i.e. the discrete geometric conservation holds for the compatible dofs. The error at T equals the error at t = 0 bit for bit (selftest iii). What does change is the representation error of the (frozen) field on the distorted mesh:

| p, q, eps | err B(0) = err B(T) | peak err B (at s = 1) | peak/initial | energy error (E-E_ex)/E_ex: t = 0 / worst over t |
|---|---|---|---|---|
| 2, 2, 0.3 | 2.28e-2 | 2.35e-2 | 1.03 | -2.6e-4 / 3.1e-4 |
| 2, 2, 0.6 | 2.28e-2 | 2.78e-2 | 1.22 | -2.6e-4 / 4.4e-4 (x1.7) |
| 3, 2, 0.3 | 1.52e-3 | 5.47e-3 | 3.6 | -7.7e-7 / 5.6e-5 |
| 3, 2, 0.6 | 1.52e-3 | 1.32e-2 | 8.7 | -7.7e-7 / 2.6e-4 (x340) |
| 3, 3, 0.3 | 1.52e-3 | 1.55e-3 | 1.03 | -7.7e-7 / 9.7e-7 |
| 3, 3, 0.6 | 1.52e-3 | 1.73e-3 | 1.14 | -7.7e-7 / 1.6e-6 (x2.1) |

At p = 2 the peak is mild (FE error dominated). At p = 3 with q = 2 the peak is 3.6x / 8.7x the initial error, but with q = 3 it is only 1.03x / 1.14x: **the p = 3, q = 2 peak is the geometry error of representing the flow map f1 by Q2 nodes, not the ND/RT space.** (The frozen state is the exact discrete Lagrangian state; what deviates is the continuous push-forward it is compared with.) This matters for every rezone result at p = 3, q = 2 below. The tracking of E_exact(t) by E_B(t) within a factor 2 of the t = 0 error holds for p = 2 (x1.7) and p = 3 with q = 3 (x2.1), not for p = 3, q = 2 (x340, geometry).

**Final error after R rezones** (error vs the exact field at T on the uniform mesh, `errB(T)`, shown as R=5 / R=20; the t = 0 error, i.e. the floor of the uniform mesh, is 2.28e-2 at p = 2 and 1.52e-3 at p = 3):

| operator | p=2, eps=0.3 | p=2, eps=0.6 | p=3, eps=0.3 | p=3, eps=0.6 |
|---|---|---|---|---|
| A_pt | 3.27e-02 / 2.87e-02 | -- | -- | -- |
| A_int | 2.78e-02 / 2.81e-02 | 4.17e-02 / 4.89e-02 | 4.09e-03 / 3.30e-03 | 2.20e-02 / 1.22e-02 |
| A_int + Coulomb | 2.70e-02 / 2.71e-02 | 3.46e-02 / 3.75e-02 | 3.36e-03 / 2.79e-03 | 1.74e-02 / 7.39e-03 |
| A_int + Coulomb (pre) | 2.75e-02 / 2.73e-02 | 3.57e-02 / 3.75e-02 | -- | 1.50e-02 / 6.37e-03 |
| A_l2 | 2.68e-02 / 2.57e-02 | -- | -- | -- |
| A_l2 + Coulomb | 2.56e-02 / 2.38e-02 | 3.76e-02 / 3.44e-02 | 5.84e-03 / 3.09e-03 | 2.68e-02 / 1.08e-02 |
| A_l2 + Coulomb (pre) | 2.63e-02 / 2.43e-02 | 3.50e-02 / 3.45e-02 | -- | 1.80e-02 / 7.09e-03 |
| B_int | 2.75e-02 / 2.67e-02 | 3.37e-02 / 3.75e-02 | 2.91e-03 / 2.67e-03 | 1.50e-02 / -- |
| B_l2 | 2.66e-02 / 2.64e-02 | 3.10e-02 / 3.40e-02 | 2.87e-03 / 2.66e-03 | 1.36e-02 / 6.50e-03 |
| B_l2c | 2.66e-02 / 2.63e-02 | 3.07e-02 / 3.31e-02 | 2.82e-03 / 2.64e-03 | 1.34e-02 / 5.97e-03 |

(Full tables for R = 2, 5, 20, dB, energy, helicity, flux, divergence, time: `results/exp5_tables.md`, Tables 5.2-5.8.)

Conservation / cost at R = 5 (dB / E(T)/E(0)-1 / H(T)/H(0)-1 / max slice-flux error / scaled div / total transfer wall time in s, q = 2):

| operator | p=2, eps=0.6: dB / E(T)/E(0)-1 / H(T)/H(0)-1 / flux err / scaled div / transfer s | p=3, eps=0.6: same |
|---|---|---|
| A_pt | -- | -- |
| A_int | 3.49e-02 / -2.81e-03 / -4.04e-03 / 4.44e-16 / 1.20e-15 / 10 | 2.19e-02 / 6.50e-04 / 1.73e-04 / 6.11e-16 / 3.17e-15 / 37 |
| A_int + Coulomb | 2.60e-02 / -2.97e-03 / -3.66e-03 / 3.89e-16 / 1.20e-15 / 16 | 1.73e-02 / 4.34e-04 / 1.36e-04 / 1.15e-15 / 3.12e-15 / 87 |
| A_int + Coulomb (pre) | 2.75e-02 / -3.04e-03 / -3.81e-03 / 7.08e-16 / 1.23e-15 / 22 | 1.49e-02 / 2.16e-04 / -5.46e-06 / 5.13e-16 / 3.17e-15 / 150 |
| A_l2 | -- | -- |
| A_l2 + Coulomb | 3.00e-02 / 7.82e-04 / -1.05e-04 / 5.55e-16 / 1.20e-15 / 17 | 2.68e-02 / 7.37e-04 / 2.26e-05 / 1.26e-15 / 3.14e-15 / 70 |
| A_l2 + Coulomb (pre) | 2.67e-02 / 5.64e-04 / -1.40e-04 / 4.72e-16 / 1.23e-15 / 24 | 1.79e-02 / 3.22e-04 / 3.94e-06 / 2.78e-16 / 3.15e-15 / 94 |
| B_int | 2.48e-02 / -3.37e-03 / -- / 3.22e-03 / 1.88e-02 / 26 | 1.49e-02 / 2.90e-04 / -- / 1.13e-04 / 7.99e-03 / 139 |
| B_l2 | 2.10e-02 / -1.99e-03 / -- / 5.51e-03 / 2.82e-02 / 11 | 1.35e-02 / 5.14e-05 / -- / 8.07e-03 / 2.09e-02 / 20 |
| B_l2c | 2.05e-02 / -1.93e-03 / -- / 1.58e-14 / 1.12e-14 / 20 | 1.33e-02 / 5.11e-05 / -- / 1.16e-14 / 8.03e-15 / 51 |

* **Divergence and flux.** All A-route variants: scaled div 1e-15 and slice flux 1e-15 at every step, whatever the operator (even A_pt) -- no cleaning solve. B-route without cleaning: div 2e-3..5e-2 (scaled), slice-flux error 1e-4..2e-2 (B_l2 at p = 3, eps = 0.6: 8e-3). B_l2c: div 1e-15..1e-14, flux 1e-14 without mean field; with B0 (p = 2, eps = 0.3, R = 5) B_l2c keeps the flux to 1e-6, B_l2 1.6e-4, B_int 8.8e-4 (the A-route Lagrangian: 9e-16), see the B0 table below.
* **Energy.** |E(T)/E(0) - 1| <= 1e-2 for everything (p = 3: <= 4e-4 at R = 20, up to 2e-3 (A_int) and 7e-3 (A_l2 + Coulomb) at R = 2, eps = 0.6); no operator is energy conserving. At p = 2 the rezoned B-route and A_int lose 0.2-1% (eps = 0.6: B_int -1.0%, B_l2c -0.9%, A_int -0.4%), A_l2 + Coulomb stays at +0.03%.
* **Helicity** (A-route): |H(T)/H(0) - 1| = 2e-3..6e-3 for A_int at p = 2 (R = 20), 1.6e-4..2.8e-4 at p = 3, <= 4e-4 for A_l2 + Coulomb. Not conserved by any transfer, but small. (Not defined for the B-route here.)

## Interpretation (honest, one set of 10 sentences)

1. **Mesh motion alone costs nothing in the dofs** (R = 0: dofs, helicity, fluxes and div are invariant to roundoff) -- this demonstrates the discrete geometric conservation for the compatible dofs; the only motion-related error is the representation error of the field on the distorted mesh (peak at s = 1: +3%/+22% at p = 2, but a factor 3.6/8.7 at p = 3 with q = 2, which turns out to be a geometric-order effect: with q = 3 it is +3%/+14%).
2. The dominant error in the rezoned runs is **not accumulation but the projection floor of the new mesh at the moment of the rezone**. At a rezone at large deformation the field on the uniform mesh is more demanding (compressed, sheared): the floor `||Pi_RT B_pf - B_pf||` of the exact push-forward on the uniform mesh is 2.6e-2 / 5.5e-2 (p = 2, eps = 0.3 / 0.6, t = 0.5), 2.3e-3 / 1.4e-2 (p = 3, q = 2) and 8.7e-3 (p = 3, q = 3, eps = 0.6, t = 0.4) against 2.3e-2 / 1.5e-3 at t = 0 (Table 5.15). Right after the rezone at maximal deformation the achieved error is, in units of this floor: B_l2c / B_l2 1.01 (p = 2, eps = 0.6, t = 0.5), 1.04 (p = 2, eps = 0.6, t = 0.4), 1.00 (p = 3, q = 3, eps = 0.6, t = 0.4); A_int + Coulomb 1.04, 1.08, 1.02; B_int 1.04, 1.08, --; A_l2 + Coulomb 1.46, 1.28, 1.40 (p = 3, q = 2 at t = 0.5: 1.18, 1.30, 1.25, 2.72 for B_l2c, A_int + Coulomb, B_int, A_l2 + Coulomb, where the pre-rezone error already contains the geometry error).
3. **Information lost is not recovered when the mesh moves back**: the later rezones towards the original configuration have a floor of 1.8e-3 (p = 3) but the state stays at 4.7e-3 (q = 3, R = 5, t = 0.8) or at 3.7e-3 at T, 2.4-2.6x the floor of the uniform mesh at T -- so the final error is set by the *worst* rezone (the one at max compression/shear), then reduced a little by the later ones.
4. **Each rezone adds a small, deformation-proportional error and the total does not grow linearly in the number of rezones.** The mean jump in err B across a transfer (Table 5.9) at R = 20, p = 2, eps = 0.3 is +1.7e-4 (B_l2c, 0.7% of the error), +1.8e-4 (B_int), +2.6e-4 (A_int), +2.1e-4 (A_int + Coulomb), +4.3e-5 (A_l2 + Coulomb), +2.9e-4 (A_pt); it falls roughly like 1/R (R = 2: 1.8e-3, R = 5: 7.1e-4, R = 20: 1.7e-4 for B_l2c), so R x jump is about constant (3.4-3.6e-3) -- the cumulative cost of rezoning a given motion is set by the total mesh displacement that has to be absorbed, not by the number of transfers; accordingly `errB(T)` is flat or decreasing in R (B_l2c, p = 2, eps = 0.3: 2.68, 2.66, 2.63e-2 for R = 2, 5, 20; p = 3, eps = 0.6, q = 2: 3.3e-2, 1.3e-2, 6.0e-3, because the floor and the Q2-geometry error that get baked in at a rezone are smaller the smaller the deformation between rezones). (At p = 3, q = 2 the per-rezone jumps are negative because rezoning removes the Q2-geometry error, so they cannot be read as transfer errors there.)
5. The A-route with Coulomb gauge: **A_int + Coulomb is within 1.5-3% of B_l2c at p = 2, eps = 0.3 (R = 5, 20), 13% at eps = 0.6, 6-30% at p = 3 (R = 20: eps = 0.3: 2.79e-3 vs 2.64e-3, eps = 0.6: 7.4e-3 vs 6.0e-3; R = 5: 3.4e-3 vs 2.8e-3 and 1.7e-2 vs 1.3e-2), but 9-46% worse at R = 2** (where the one big rezone dominates: p = 2, eps = 0.6: 4.9e-2 vs 3.5e-2; p = 3, eps = 0.3: 7.4e-3 vs 5.8e-3, eps = 0.6: 4.8e-2 vs 3.3e-2). It buys exact divergence and exact flux with no cleaning solve, which B_l2c obtains with a global Schur solve (B_l2 and B_l2c have the same B error to 1-2%, so the cleaning costs no accuracy, only 2-3x time). Costs per rezone are comparable (p = 2: A_int + Coulomb 2.8 s vs B_l2c 4.5 s; p = 3: 19 s vs 18 s, contended); the Coulomb solve itself is 14-23 CG iterations, 0.11 s (p = 2) / 0.7-1.4 s (p = 3) per rezone, i.e. 4-12% of the rezone.
6. **Gauge fixing helps the A-route here, and where it is applied matters.** The pushed-forward potential `F^{-T}A_0` is not divergence free in the deformed geometry, so the gradient part of the frozen `a` grows with the deformation since the last gauge fix (fraction `||P_grad a||/||a||` before the fix: 0.19 / 0.36 at R = 2 for eps = 0.3 / 0.6, 0.07 / 0.135 (R = 5), 0.019 / 0.038 (R = 20); Table 5.11), exactly the pollution mechanism of E3, and fixing the gauge after each transfer reduces `errB(T)` of A_int by 3-39% vs no fix (p = 3, eps = 0.6: R = 5 2.2e-2 -> 1.7e-2, R = 20 1.2e-2 -> 7.4e-3; p = 2, eps = 0.6: 4.2e-2 -> 3.5e-2, 4.9e-2 -> 3.8e-2; p = 2, eps = 0.3: -3%) and that of A_l2 by 7% (p = 2, eps = 0.3, R = 20: 2.6e-2 -> 2.4e-2, dB 1.2e-2 -> 6.7e-3). Applying it *before* the transfer, on the deformed source mesh (`coulomb_pre`, extra), is better at p = 3 (eps = 0.6: A_int 1.5e-2 / 6.4e-3 and A_l2 1.8e-2 / 7.1e-3 at R = 5 / 20 vs 1.7e-2 / 7.4e-3 and 2.7e-2 / 1.1e-2, i.e. A_int + Coulomb(pre) is within 12% / 7% of B_l2c) and neutral at p = 2, but costs a new AMG setup per rezone (1.7-2.9x the total transfer time at p = 3).
7. Ranking of the cheap options: A_l2 + Coulomb is the best at p = 2, eps = 0.3 (R = 20: errB(T) 2.38e-2 = 1.04x the floor, B_l2c 1.15x) but the worst of the gauge-fixed variants in the large-deformation cases (R = 2: 1.06-2.5x B_l2c; R = 5, 20 at p = 3: 1.2-2.1x); A_int is the most robust A-route operator. A_pt, which E4 showed to blow up after ~15 ping-pong remaps, is stable for 20 smooth rezones (p = 2, eps = 0.3: errB(T) 3.3e-2 / 2.9e-2 vs 2.8e-2 / 2.8e-2 for A_int at R = 5 / 20, i.e. +18% / +2%) because consecutive meshes differ little. B_int is the most expensive operator (p = 3: about 7x B_l2 per rezone) and is never better than B_l2.
8. **Mean field.** With B0 on (p = 2, eps = 0.3, R = 5, B-route) the errors are those of the B0-off runs to < 1% (dB 1.35e-2 for B_l2/B_l2c, 1.53e-2 B_int); the net slice flux, exact for the A-route (9e-16 at all times in the Lagrangian reference), drifts for every B-route operator: 1.0e-6 (B_l2c), 1.6e-4 (B_l2), 8.8e-4 (B_int) at T (maximum during the cycle 1.1e-6, 1.0e-3, 1.2e-3), while B_l2 and B_int also have div 1e-2. (Table 5.12, `results/exp5_b0.png`.)
9. **A-route and a mean field is not solved by a plain transfer (finding, not implemented).** Under a frozen-in deformation the mean field is not constant: `B_pf = F B0/det F` has physical mean exactly B0 (net flux through every homologous slice is invariant), but a fluctuating part with L2 norm 22% (eps = 0.3) / 48% (eps = 0.6) of |B0| at s = 1. In the A-route state `(a, b0_t)` the dofs `b0_t` carry this deformed harmonic part during the Lagrangian phase, so the motion is fine; but a rezone transfers only `a`, and resetting `b0_u = Pi_RT^int(B0)` (constant) throws the deformed fluctuation away, while carrying `b0_t` over as dofs is meaningless on a different mesh. A correct A-route rezone needs the fluctuation expressed through a potential, `C_h a_B0 = Pi_RT(B_pf - B0)` on the new mesh, i.e. a curl inversion (a global Coulomb-gauge solve with an AMS-type preconditioner), or a direct B transfer of the harmonic part. The driver therefore aborts for `-b0` with an A-route operator and R > 0 (the Lagrangian R = 0 case with B0 runs: error, flux, div identical to the B0-off case).
10. Caveats: one deformation family (f1, a standing-wave shear that returns to the identity), one field family (Beltrami `abc`, k = 2 pi on N = 8: 8 cells per wavelength), N = 8 only, 1 MPI rank, q = 2 (the p = 3 numbers at q = 2 are dominated by geometry error as shown above; the q = 3 control covers only p = 3, eps = 0.6, R = 5), and the final-state comparisons are L2 relative errors against an analytic push-forward whose geometry is exact while the mesh is a Q2/Q3 interpolation of the map. Nothing here says anything about long cycles with many return trips (E4 covers ping-pong, where the A-route drifts); a 20-rezone cycle is the longest tested.

## Selected tables (all tables in `results/exp5_tables.md`)

### Table 5.1: Lagrangian reference (no rezone), N=8, abc, K=100 steps

| tag | p | q | eps | B0 | err B(t=0) | peak err B (t) | peak/initial | err B(T) | err A(0) | peak err A | (E-E_ex)/E_ex at 0 | worst over t | min detJ | max &#124;H/H0-1&#124; | max flux err | max scaled div | &#124;&#124;a_T-a_0&#124;&#124;_inf |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| lag | 2 | 2 | 0.3 | 0 | 2.28e-02 | 2.35e-02 (0.50) | 1.03 | 2.28e-02 | 2.76e-03 | 5.59e-03 | -2.61e-04 | 3.09e-04 | 0.901 | 3.33e-16 | 6.80e-16 | 1.18e-15 | 0.00e+00 |
| lag | 2 | 2 | 0.6 | 0 | 2.28e-02 | 2.78e-02 (0.50) | 1.22 | 2.28e-02 | 2.76e-03 | 1.15e-02 | -2.61e-04 | 4.44e-04 | 0.604 | 3.33e-16 | 6.80e-16 | 1.18e-15 | 0.00e+00 |
| lag | 3 | 2 | 0.3 | 0 | 1.52e-03 | 5.47e-03 (0.50) | 3.61 | 1.52e-03 | 1.26e-04 | 4.87e-03 | -7.67e-07 | 5.56e-05 | 0.901 | 3.33e-16 | 4.30e-16 | 3.12e-15 | 0.00e+00 |
| lag | 3 | 2 | 0.6 | 0 | 1.52e-03 | 1.32e-02 (0.50) | 8.74 | 1.52e-03 | 1.26e-04 | 1.11e-02 | -7.67e-07 | 2.60e-04 | 0.604 | 3.33e-16 | 4.30e-16 | 3.12e-15 | 0.00e+00 |
| lag_q3 | 3 | 3 | 0.3 | 0 | 1.52e-03 | 1.55e-03 (0.50) | 1.03 | 1.52e-03 | 1.25e-04 | 3.48e-04 | -7.67e-07 | 9.67e-07 | 0.910 | 4.44e-16 | 7.49e-16 | 3.19e-15 | 0.00e+00 |
| lag_q3 | 3 | 3 | 0.6 | 0 | 1.52e-03 | 1.73e-03 (0.50) | 1.14 | 1.52e-03 | 1.25e-04 | 7.20e-04 | -7.67e-07 | 1.62e-06 | 0.639 | 4.44e-16 | 7.49e-16 | 3.19e-15 | 0.00e+00 |
| lag_b0 | 2 | 2 | 0.3 | 1 | 2.27e-02 | 2.35e-02 (0.50) | 1.03 | 2.27e-02 | 2.76e-03 | 5.59e-03 | -2.59e-04 | 3.08e-04 | 0.901 | -- | 8.88e-16 | 1.19e-15 | 0.00e+00 |

### Table 5.9: error added by one rezone (jump of the relative B error across the transfer; mean over the rezones of a run)

Rows: operator; columns: (p, eps) at r = 5 and 20 rezones.  Entry: mean(err_post - err_pre) with err = ||b_h - B_pf||/||B_pf|| measured on the source mesh right before and on the uniform target right after the transfer; in parentheses: mean err_post/err_pre; in brackets: sum of the jumps over all rezones of the run (note that for r = 5 the rezone at t=0.6 is an exact identity because s(0.6)=s(0.4): 4 genuine transfers; r=20 has none, r=2 none).

| operator | p=2 e=0.3 r=2 | p=2 e=0.3 r=5 | p=2 e=0.3 r=20 | p=2 e=0.6 r=2 | p=2 e=0.6 r=5 | p=2 e=0.6 r=20 | p=3 e=0.3 r=2 | p=3 e=0.3 r=5 | p=3 e=0.3 r=20 | p=3 e=0.6 r=2 | p=3 e=0.6 r=5 | p=3 e=0.6 r=20 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| A_pt | -- | +1.9e-03 (1.08) [sum +9.6e-03] | +2.9e-04 (1.01) [sum +5.8e-03] | -- | -- | -- | -- | -- | -- | -- | -- | -- |
| A_int | +3.6e-03 (1.15) [sum +7.1e-03] | +9.6e-04 (1.05) [sum +4.8e-03] | +2.6e-04 (1.01) [sum +5.2e-03] | +1.3e-02 (1.50) [sum +2.5e-02] | +3.7e-03 (1.14) [sum +1.8e-02] | +1.3e-03 (1.04) [sum +2.6e-02] | +4.5e-04 (1.06) [sum +9.1e-04] | +2.5e-05 (1.02) [sum +1.3e-04] | +2.1e-05 (1.01) [sum +4.3e-04] | +4.0e-03 (1.21) [sum +8.0e-03] | +5.2e-04 (1.06) [sum +2.6e-03] | +1.4e-04 (1.02) [sum +2.8e-03] |
| A_int+Coulomb | +3.0e-03 (1.13) [sum +6.0e-03] | +7.9e-04 (1.04) [sum +3.9e-03] | +2.1e-04 (1.01) [sum +4.3e-03] | +1.1e-02 (1.46) [sum +2.1e-02] | +2.3e-03 (1.11) [sum +1.2e-02] | +7.5e-04 (1.03) [sum +1.5e-02] | -3.3e-04 (0.96) [sum -6.7e-04] | -1.0e-04 (0.98) [sum -5.1e-04] | +3.4e-06 (1.00) [sum +6.7e-05] | +1.5e-03 (1.15) [sum +3.0e-03] | -4.0e-04 (1.01) [sum -2.0e-03] | -8.4e-05 (1.00) [sum -1.7e-03] |
| A_int+Coulomb(pre) | +3.8e-03 (1.17) [sum +7.6e-03] | +8.9e-04 (1.04) [sum +4.4e-03] | +2.2e-04 (1.01) [sum +4.4e-03] | +1.4e-02 (1.54) [sum +2.7e-02] | +2.5e-03 (1.11) [sum +1.2e-02] | +7.5e-04 (1.03) [sum +1.5e-02] | -- | -- | -- | -- | -8.9e-04 (0.97) [sum -4.4e-03] | -1.3e-04 (0.99) [sum -2.5e-03] |
| A_l2 | -- | +7.3e-04 (1.04) [sum +3.6e-03] | +1.4e-04 (1.01) [sum +2.8e-03] | -- | -- | -- | -- | -- | -- | -- | -- | -- |
| A_l2+Coulomb | +2.5e-03 (1.14) [sum +5.0e-03] | +5.1e-04 (1.04) [sum +2.5e-03] | +4.3e-05 (1.00) [sum +8.6e-04] | +1.7e-02 (1.84) [sum +3.4e-02] | +2.7e-03 (1.15) [sum +1.3e-02] | +5.7e-04 (1.03) [sum +1.1e-02] | +2.5e-03 (1.35) [sum +4.9e-03] | +3.5e-04 (1.08) [sum +1.8e-03] | +4.2e-05 (1.01) [sum +8.4e-04] | +2.2e-02 (2.07) [sum +4.5e-02] | +1.9e-03 (1.20) [sum +9.3e-03] | +1.8e-04 (1.03) [sum +3.6e-03] |
| A_l2+Coulomb(pre) | +2.9e-03 (1.15) [sum +5.8e-03] | +6.4e-04 (1.04) [sum +3.2e-03] | +6.8e-05 (1.00) [sum +1.4e-03] | +1.2e-02 (1.54) [sum +2.4e-02] | +2.1e-03 (1.12) [sum +1.1e-02] | +5.7e-04 (1.03) [sum +1.1e-02] | -- | -- | -- | -- | -1.9e-04 (1.05) [sum -9.4e-04] | +1.8e-05 (1.01) [sum +3.6e-04] |
| B_int | +2.7e-03 (1.12) [sum +5.5e-03] | +9.0e-04 (1.04) [sum +4.5e-03] | +1.8e-04 (1.01) [sum +3.7e-03] | +7.1e-03 (1.40) [sum +1.4e-02] | +2.1e-03 (1.10) [sum +1.1e-02] | +7.4e-04 (1.03) [sum +1.5e-02] | -9.4e-04 (0.87) [sum -1.9e-03] | -1.8e-04 (0.96) [sum -9.0e-04] | -9.5e-07 (1.00) [sum -1.9e-05] | -3.1e-03 (1.04) [sum -6.3e-03] | -9.2e-04 (0.97) [sum -4.6e-03] | -- |
| B_l2 | +1.8e-03 (1.08) [sum +3.7e-03] | +7.3e-04 (1.03) [sum +3.6e-03] | +1.8e-04 (1.01) [sum +3.5e-03] | +3.8e-03 (1.32) [sum +7.5e-03] | +1.5e-03 (1.09) [sum +7.7e-03] | +5.4e-04 (1.02) [sum +1.1e-02] | -1.2e-03 (0.83) [sum -2.4e-03] | -2.6e-04 (0.94) [sum -1.3e-03] | -1.3e-05 (1.00) [sum -2.5e-04] | -5.0e-03 (0.97) [sum -1.0e-02] | -1.2e-03 (0.94) [sum -6.0e-03] | -2.0e-04 (0.97) [sum -4.0e-03] |
| B_l2c | +1.8e-03 (1.08) [sum +3.6e-03] | +7.1e-04 (1.03) [sum +3.6e-03] | +1.7e-04 (1.01) [sum +3.3e-03] | +3.6e-03 (1.32) [sum +7.2e-03] | +1.5e-03 (1.08) [sum +7.5e-03] | +4.8e-04 (1.02) [sum +9.7e-03] | -1.3e-03 (0.82) [sum -2.5e-03] | -2.7e-04 (0.93) [sum -1.4e-03] | -1.4e-05 (1.00) [sum -2.9e-04] | -5.2e-03 (0.96) [sum -1.0e-02] | -1.2e-03 (0.93) [sum -6.2e-03] | -2.3e-04 (0.97) [sum -4.5e-03] |

### Table 5.12: mean field B0 = (0.3,-0.2,0.5) on (p=2, eps=0.3, 5 rezones; B-route; A_int Lagrangian reference with B0 for comparison)

| operator | rezones | ||b_T-b_0||/||b_0|| | err B(T) | E(T)/E(0)-1 | slice-flux error at T | scaled div | max flux err over t | transfer time (s) |
|---|---|---|---|---|---|---|---|---|
| B_int | 5 | 1.53e-02 | 2.75e-02 | -3.32e-04 | 8.81e-04 | 9.82e-03 | 1.15e-03 | 26.0 |
| B_l2 | 5 | 1.35e-02 | 2.65e-02 | -3.55e-04 | 1.64e-04 | 1.29e-02 | 1.03e-03 | 9.7 |
| B_l2c | 5 | 1.35e-02 | 2.65e-02 | -3.22e-04 | 9.55e-07 | 2.84e-15 | 1.14e-06 | 19.9 |
| A_int | 0 | 0.00e+00 | 2.27e-02 | 0.00e+00 | 8.88e-16 | 1.19e-15 | 8.88e-16 | 0.0 |

### Table 5.15: error right after each rezone versus the projection floor of the exact push-forward on the new uniform mesh

floor = ||Pi_RT^int B_pf(t_r) - B_pf(t_r)||/||B_pf|| on the new uniform mesh (best any transfer can do; operator independent, `-floors` mode), floor(src) = same on the deformed source mesh just before the rezone. Columns: relative B error right after the transfer (and in brackets error/floor). r=5: the rezone at t=0.6 is an identity (s(0.6)=s(0.4)).

**p=2, eps=0.3, r=2 rezones (q=2)**

| rezone | t | min detJ (source) | floor(src) | floor(new) | B_l2c | B_l2 | B_int | A_int+Coulomb | A_int | A_l2+Coulomb |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 0.50 | 0.901 | 2.35e-02 | 2.62e-02 | 2.90e-02 (1.11) | 2.90e-02 (1.11) | 3.07e-02 (1.17) | 3.07e-02 (1.17) | 3.07e-02 (1.17) | 3.34e-02 (1.28) |
| 2 | 1.00 | 0.917 | 2.58e-02 | 2.28e-02 | 2.68e-02 (1.18) | 2.68e-02 (1.18) | 2.87e-02 (1.26) | 2.93e-02 (1.29) | 3.04e-02 (1.33) | 2.84e-02 (1.25) |

**p=2, eps=0.3, r=5 rezones (q=2)**

| rezone | t | min detJ (source) | floor(src) | floor(new) | B_l2c | B_l2 | B_int | A_int+Coulomb | A_int | A_l2+Coulomb |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 0.20 | 0.988 | 2.29e-02 | 2.31e-02 | 2.45e-02 (1.06) | 2.45e-02 (1.06) | 2.49e-02 (1.08) | 2.49e-02 (1.08) | 2.49e-02 (1.08) | 2.50e-02 (1.08) |
| 2 | 0.40 | 0.932 | 2.35e-02 | 2.54e-02 | 2.91e-02 (1.15) | 2.91e-02 (1.15) | 3.06e-02 (1.21) | 3.07e-02 (1.21) | 3.09e-02 (1.22) | 3.28e-02 (1.29) |
| 3 | 0.60 | 1.000 | 2.54e-02 | 2.54e-02 | 2.91e-02 (1.15) | 2.91e-02 (1.15) | 3.06e-02 (1.21) | 3.07e-02 (1.21) | 3.09e-02 (1.22) | 3.28e-02 (1.29) |
| 4 | 0.80 | 0.940 | 2.51e-02 | 2.31e-02 | 2.52e-02 (1.09) | 2.52e-02 (1.09) | 2.54e-02 (1.10) | 2.49e-02 (1.08) | 2.62e-02 (1.13) | 2.41e-02 (1.04) |
| 5 | 1.00 | 0.989 | 2.31e-02 | 2.28e-02 | 2.66e-02 (1.17) | 2.66e-02 (1.17) | 2.75e-02 (1.21) | 2.70e-02 (1.18) | 2.78e-02 (1.22) | 2.56e-02 (1.12) |

**p=2, eps=0.6, r=2 rezones (q=2)**

| rezone | t | min detJ (source) | floor(src) | floor(new) | B_l2c | B_l2 | B_int | A_int+Coulomb | A_int | A_l2+Coulomb |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 0.50 | 0.604 | 2.78e-02 | 5.55e-02 | 5.63e-02 (1.01) | 5.61e-02 (1.01) | 5.78e-02 (1.04) | 5.80e-02 (1.04) | 5.80e-02 (1.04) | 8.10e-02 (1.46) |
| 2 | 1.00 | 0.731 | 5.36e-02 | 2.28e-02 | 3.48e-02 (1.53) | 3.54e-02 (1.55) | 4.15e-02 (1.82) | 4.85e-02 (2.13) | 5.25e-02 (2.30) | 6.34e-02 (2.78) |

**p=2, eps=0.6, r=5 rezones (q=2)**

| rezone | t | min detJ (source) | floor(src) | floor(new) | B_l2c | B_l2 | B_int | A_int+Coulomb | A_int | A_l2+Coulomb |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 0.20 | 0.953 | 2.31e-02 | 2.42e-02 | 2.70e-02 (1.12) | 2.70e-02 (1.12) | 2.81e-02 (1.16) | 2.82e-02 (1.17) | 2.82e-02 (1.17) | 2.94e-02 (1.22) |
| 2 | 0.40 | 0.716 | 2.71e-02 | 4.44e-02 | 4.63e-02 (1.04) | 4.63e-02 (1.04) | 4.80e-02 (1.08) | 4.80e-02 (1.08) | 4.88e-02 (1.10) | 5.69e-02 (1.28) |
| 3 | 0.60 | 1.000 | 4.44e-02 | 4.44e-02 | 4.63e-02 (1.04) | 4.63e-02 (1.04) | 4.80e-02 (1.08) | 4.80e-02 (1.08) | 4.88e-02 (1.10) | 5.69e-02 (1.28) |
| 4 | 0.80 | 0.800 | 4.10e-02 | 2.42e-02 | 3.14e-02 (1.30) | 3.17e-02 (1.31) | 3.46e-02 (1.43) | 3.65e-02 (1.51) | 4.55e-02 (1.88) | 3.66e-02 (1.51) |
| 5 | 1.00 | 0.959 | 2.40e-02 | 2.28e-02 | 3.07e-02 (1.35) | 3.10e-02 (1.36) | 3.37e-02 (1.48) | 3.46e-02 (1.52) | 4.17e-02 (1.83) | 3.76e-02 (1.65) |

## Bugs / deviations / unfinished

* Deviations are listed in "Set-up" (Newton evaluation at physical points instead of label GF interpolation; quadrature-based energy/helicity/div; thinned matrix; p = 3 logging stride; R = 1 only at p = 2 / eps = 0.3; R = 2 added; B_int p = 3, eps = 0.6, R = 20 skipped; the `D_h` geometry dependence). The A-route with B0 and R > 0 is not implemented (point 9).
* While developing, the first lag runs used quadrature order 2p+2q and the GetDivergence-based div; they were rerun with the final code (all numbers here are from the final binary). The two `rc=1` entries for `lagq3_*` in `results/exp5_logs/_summary.log` are the first attempts that I killed to change the logging stride; they were rerun.
* Not done: more than 1 rank, N != 8, other fields/deformations, p = 1, time-step (K) sensitivity (K = 100 only; the Lagrangian phase has no time discretisation error by construction), error attribution by tracking the same field through a *single* long Lagrangian segment with q = 3 for the whole matrix (only the control above), an A-route curl-inversion for the mean field.
