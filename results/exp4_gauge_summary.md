# E4-gauge summary: Coulomb (discrete minimum-norm) gauge fixing of the A-route under repeated remaps

Code: `src/vp_gauge.{hpp,cpp}` (`vp::CoulombGauge`), wired into `experiments/exp4_repeat.cpp` (`-gauge none|coulomb|jacobi:k`) and
`experiments/exp3_remap.cpp` (`-gfix none|coulomb|jacobi:k`; exp3 already has `-gauge <amplitude>` for the gauge-perturbation test, so the new option is `-gfix`).
Checks / least-squares test: `experiments/gauge_check.cpp` (`bin/gauge_check -N 8 -p P [-lsq]`, outputs `results/exp4_gauge_check_p{2,3}.txt`).
Runs: `scripts/run_exp4_gauge.sh` (parts g_p2 g_p3 g_long g_jac g_b0 g_cost e3). Data: `results/exp4_gauge.csv` (333 rows, all old exp4 columns kept + `gauge, t_transfer_last, gauge_iters_*, gauge_time_*, gauge_resid_*, gauge_gf_*, gauge_dB_rel, gauge_dH_scaled`),
`results/exp3_gfix.csv` (one remap, 9 rows; extra columns `gfix_*`, `pre_grad_frac`, `pre_helicity`). Plots: `scripts/plot_exp4_gauge.py` -> `results/exp4_gauge_p2.png`, `_p3.png`, `_long.png`, `_misc.png`. Old files untouched.
Set-up as in exp4: N=8, abc, uniform <-> f1 eps 0.3, q=p, state on mesh A (even n), "no gauge" and B_l2c numbers are from `results/exp4_repeat.csv` (same code path, same set-up).
The gauge fix is applied after every A-route transfer (the B route is untouched); diagnostics (grad fraction before/after, helicity before/after, ||B change||) are evaluated outside the timers and only at logged steps (every 5, or 20 in the long run).
Machine: 1 rank per job, 3-4 jobs on 4 cores (the other agent's exp5 jobs were running), so ALL timings are inflated/noisy; only same-run ratios mean anything.

## What the gauge fix is

`CoulombGauge`: L = G_h^T M_ND G_h built by a Hypre RAP from the fully assembled ND mass matrix (same quadrature order 2p+2q+2 as the diagnostics); a <- a - G_h chi, L chi = G_h^T M_ND a; AMG(BoomerAMG)-preconditioned CG with the preconditioner wrapped by projections onto mean-free vectors (the constant null space; without the wrapper the CG stalled at 500 iterations for p=1 because the AMG coarse solve of the singular matrix injects constants; found by the first exp3 p=1 run, fixed, all results below are with the fix). rtol 1e-12.
* **L equals the H1 DiffusionIntegrator stiffness to roundoff**: max|L-K|/max|K| = 2.1e-15 (p=2, both meshes), 2.5e-15 / 2.9e-15 (p=3); `CoulombGauge::StiffnessMismatch()`. (As expected: the covariant Piola map of ND makes grad(H1) mass integrand identical to the Laplacian integrand.) max|C_h G_h| = 0 (the sparse product has no nonzero entry).
* After Apply: ||G^T M a|| ~ 1e-15 .. 1e-13 (p=1,2,3), gradient fraction ~1e-11 (p=2: 1e-11 measured by the independent diag solve, 0 at p=3 where the diag CG stops at its initial residual).
* **B invariance and helicity (exact-sequence test).** Over all gauge solves logged (p=2: 100/400 remaps, p=3: 100 remaps, with/without B0, jacobi variants, and exp3 one-shot) `||C_h a_fixed - C_h a||_M/||B||` <= 1.0e-15 and the helicity change |H_after - H_before|/(||a||_L2 ||B||_L2) <= 4.0e-14 (typically 1e-15..1e-14), i.e. roundoff: int grad(chi_h).B_h = 0 for div_h-free B_h holds to machine precision (also in the unit check: gauge perturbation g=0.5, 2pi and 6pi: dB 6e-16..5e-15, dH 3e-16..3e-14, ||a|| 2.4 -> 1.414). Therefore the gauge fix never changes the one-shot accuracy: exp3 (`results/exp3_gfix.csv`, p=1,2,3, eps 0.3 -> uniform, abc+B0, A_pt/A_int/A_l2) reproduces the old B errors in `exp3_remap.csv` to all printed digits (e.g. p=2: A_pt 3.4709e-2, A_int 2.9143e-2, A_l2 3.0212e-2; p=3: 4.351e-3, 2.432e-3, 3.105e-3), energies and helicity identical; only the gradient fraction drops (p=2: 8.4e-3 -> 3e-11, p=3: 1e-3 -> ~0) and a(.) changes. Its only effect is on the NEXT transfer.

## Headline: does Coulomb gauge fixing rescue the A-route?

**Partly. It removes the gauge growth completely, but it does not remove the instability of the pointwise operator, only delays/slows the A_int instability, and it makes A_l2 at p=3 *worse*.** Numbers (state on mesh A):

p=2 (B-error floor 2.28e-2):

| op | n | E/E_ref | B err (x floor) | grad frac | H(n)/H(0) |
|---|---|---|---|---|---|
| A_pt + Coulomb | 20 / 50 / 100 | 1.024 / 21 / 1.1e8 | 0.155 / 4.5 / 1.1e4 | 1e-11 | 1.0005 / 0.978 / -3.7e4 |
| A_pt, no gauge | 20 / 50 | 1.65 / 2.3e4 | 0.80 / 151 | 0.2 / 0.77 | |
| A_int + Coulomb | 100 / 200 / 300 / 400 | 1.009 / 1.035 / 1.143 / 3.31 | 0.107 (4.7x) / 0.203 (8.9x) / 0.418 (18x) / 1.58 (69x) | 4e-11 | 0.9983 / 0.9949 / 0.957 / 0.705 |
| A_int, no gauge | 100 / 200 / 300 / 400 | 1.013 / 1.23 / 18.8 / 1765 | 0.136 (5.9x) / 0.49 / 4.2 / 42 | 0.078 .. 0.95 | 0.9951 .. -0.57 |
| A_l2 + Coulomb | 100 / 200 / 300 / 400 | 1.0039 / 1.0063 / 1.0109 / 1.0188 | 0.1085 (4.8x) / 0.145 (6.3x) / 0.178 (7.8x) / 0.214 (9.4x) | 1e-11 | 0.9926 / 0.9859 / 0.9797 / 0.9740 |
| A_l2, no gauge | 100 / 400 | 1.010 / 1.58 | 0.120 (5.3x) / 0.77 | 0.016 .. 0.07 | 0.9961 / ~0.99 |
| B_l2c (exp4) | 100 / 200 / 300 / 400 | 0.940 / 0.906 / 0.890 / 0.8875 | 0.0873 (3.8x) / 0.143 (6.3x) / 0.190 (8.3x) / 0.233 (10.2x) | -- | n/a (B-route state) |

p=3 (floor 1.51e-3), n=100: A_int + Coulomb: E/E_ref 1.0037, B err 3.88e-2 (25.6x) [no gauge: 1.023, 0.137 (90x)]; A_l2 + Coulomb: 1.0007, 2.83e-2 (18.7x) [no gauge: 1.0003, 1.87e-2 (12.3x)]; A_pt + Coulomb: blows up (E/E_ref 1.04 at n=50, 1.8e3 at n=100, B err 0.20 -> 42); B_l2c: 0.9997, 3.94e-3 (2.6x).

* **A_pt: not rescued.** The pointwise-interpolation transfer of the *curl part* is unstable by itself (B err growth x1.1-1.15 per remap with the fix vs ~x1.2 without, at p=2; the fix delays the B err > 3 point by roughly 20-30 remaps at p=2 and p=3 but the field is garbage by n=50-100 either way). div_h B stays at roundoff throughout (exactly div-free but wrong).
* **A_int: gauge growth is gone (grad fraction constant at 1e-11), the error and energy drift are reduced (p=2 n=100: B err 0.107 vs 0.136; p=3: 0.039 vs 0.137), the n=400 blow-up is strongly reduced (E/E_ref 3.3 vs 1765, B err 1.6 vs 42) but not removed:** the error still grows geometrically after n~250 (x1.014 per remap, vs x1.02 without gauge), and the helicity collapses (-29% at n=400). So the A_int instability observed in exp4 was *aggravated* by the gauge component (feedback) but has a gauge-independent part: the transfer is not contractive on the curl part either. I did not find the mechanism (e.g. spectral radius of the A_int ping-pong map restricted to the Coulomb subspace is not computed here).
* **A_l2: drift becomes slow and saturating** at p=2, n=400: E/E_ref 1.019, B err 0.214 (9.4x floor) - *better than B_l2c in both energy (+1.9% vs -11.3%) and error (0.214 vs 0.233)* at n=400 (B_l2c is better below n~200: 0.087 vs 0.109 at n=100); helicity drifts -2.6% over 400 remaps (no gauge fix: -1.2%, i.e. *worse* with the fix). A_l2 + Coulomb n=400 error growth slows (0.178 -> 0.214 per 100 remaps).
* **p=3, A_l2 + Coulomb is worse than A_l2 without gauge fix** (n=100: 2.83e-2 vs 1.87e-2; n=20: 1.30e-2 vs 9.5e-3; at n=10: 9.2e-3 vs 7.3e-3, i.e. already after 9 transfers). Since B is invariant under the fix, this is a genuine effect of the removed gradient content on the next L2 projection (removing O(h^p) gradient content of the interpolant changes what the L2-projection of a does to its curl). I did not investigate further; both the p=3 A_l2 drift and its worsening are reproducible (deterministic code). Conclusion: gauge fixing is not a universal improvement; for A_l2 at p=3 the un-fixed state is more accurate over 100 remaps (its grad fraction was only 2e-3 there).
* **Energy and helicity drift per 100 remaps (p=2, first 100)**: E: A_int +0.90% (no gauge +1.28%), A_l2 +0.39% (+1.0%), A_pt blow-up, B_l2c -6.0%. H/H(0)-1: A_int -0.17% (-0.49%), A_l2 -0.74% (-0.39%), p=3: A_int +0.22% (+0.39%), A_l2 -0.012% (-0.0045%). B_l2c has no A, hence no helicity. After 400 remaps: A_l2 + Coulomb +1.9% energy, -2.6% helicity; B_l2c -11.3% energy.
* **B0 on (p=2, A_int + Coulomb, n=100, `g_b0`)**: identical to the B0-off run (B err 0.1068 vs 0.1070; E/E_ref 1.0090; flux error 3e-16..1e-15 for all 100 remaps: mean flux exact, B_l2c: 6.4e-4); the fix does not interact with B0 (it acts on A only, B0 enters b only through Pi_RT B0).
* **Local variant (EXPERIMENTAL, `jacobi:k` = k steps of Chebyshev-accelerated Jacobi on [lmax/30, 1.05 lmax] of D^-1 L, zero initial guess, only mat-vecs)**, A_int p=2, n=100: jacobi:5: grad fraction stays flat at 1.7e-3 (vs 7.8e-2 no gauge, 4e-11 Coulomb), B err 0.1104, E/E_ref 1.0069; jacobi:20: grad fraction 8e-6, B err 0.10702 (= Coulomb 0.10702 to 6 digits), E/E_ref 1.0090. Cost 0.012 s (5 sweeps) / 0.028 s (20 sweeps) vs 0.1 s (Coulomb, same contended machine). So a *local* smoother suffices to stop the gauge feedback at n=100: the gradient content injected per remap is rough (high-frequency), and a few sweeps remove it. Not tested at n=400 or p=3, and the long-time remaining instability of A_int (see above) would obviously be unaffected.
* **AMG-CG iterations stay bounded**: 14-17 (p=2, all 400 remaps), 20-22 (p=3), 10 (p=1), independent of n and of the amplitude of the gauge component (0 .. 1). One-time setup of the gauge object (both meshes): 3.2 s (p=2), 23 s (p=3; dominated by assembling the full ND mass matrix and the RAP, contended machine; reduce with a direct H1-stiffness build, since L == K) plus the first-solve AMG setup.

## Cost per remap (`g_cost`: n=20, every 4, p=2/3, run serially in one script with B_l2c; median over the steady-state logged steps 2..20; contended machine)

| p | op | transfer (s) | gauge solve (s, its) | total (s) | vs B_l2c |
|---|---|---|---|---|---|
| 2 | B_l2c | 2.33 (29 outer + 180 inner CG) | | 2.33 | 1 |
| 2 | A_int | 0.49 | | 0.49 | 0.21 |
| 2 | A_int + Coulomb | 0.53 | 0.083 (14) | 0.61 | 0.26 |
| 2 | A_l2 | 0.81 (18 CG) | | 0.81 | 0.35 |
| 2 | A_l2 + Coulomb | 0.84 | 0.066 (14) | 0.91 | 0.39 |
| 3 | B_l2c | 4.49 (40 + 246) | | 4.49 | 1 |
| 3 | A_int | 1.80 | | 1.80 | 0.40 |
| 3 | A_int + Coulomb | 2.24 | 0.42 (21) | 2.62 | 0.58 |
| 3 | A_l2 | 2.29 (16) | | 2.29 | 0.51 |
| 3 | A_l2 + Coulomb | 2.10 | 0.45 (21) | 2.57 | 0.57 |

The gauge solve costs 8-16% (p=2) / 19-21% (p=3) of the transfer and 3-4% (p=2) / 9-10% (p=3) of B_l2c; the first transfer of B_l2c costs 3.5-4x its steady state (lazy AMG setup), the gauge object's AMG setup is also lazy. Timing noise: the same gauge solve took 0.1 s (p=2) / 1.0 s (p=3) during the contended g_p2/g_p3 runs, so +-2x. The previous exp4 cost numbers (B_l2c 3.3 s p=2) are of the same order.

## B_l2c vs "A_l2 then B = C a" vs the constrained least-squares A (`gauge_check -lsq`, f1 0.3 -> uniform, abc, no B0)

* B_l2c = M_RT-orthogonal projection of the L2-transferred b_l2 onto ker D_h (the construction: b <- b - M^-1 D^T M_L2 phi gives M(b_l2 - b_l2c) = D^T(...), orthogonal to ker D_h). On the torus ker D_h = range C_h + 3-dim harmonic space (constant fluxes). The constrained projection argmin ||C a - b_l2||_M gives C a = projection onto range C_h only, i.e. equals B_l2c up to the harmonic (mean-flux) part, which is zero here (zero mean flux source).
* Verified: p=2 / p=3: max|D_h b_l2c| = 2.0e-12 / 5.3e-13 (b_l2: 3.3 / 0.95); (b_l2c - b_l2) is M-orthogonal to C_h a_rand (5 random a): 2.8e-16 / 2.2e-15 (normalised), and to b_l2c itself 4.6e-16 / 1.2e-15; least-squares A from the curl-curl system (C^T M C + 1e-9 M_ND) a = C^T M b_l2 (shift kills the gradient part, CG + HypreAMS, 33 / 67 iterations, tol 1e-13): ||C a_ls - b_l2c||_M/||b_l2c||_M = **2.5e-11 / 2.5e-11** (shift-limited). So the constrained A-route *is* B_l2c; its curl-curl solve (AMS, 33-67 CG its, 3 s at p=2, 19 s at p=3 here, contended) is far costlier than the transfer itself, so nobody would implement it that way.
* "A_l2 then C a" is a different operator: ||C a_{A_l2} - b_l2c||/||b_l2c|| = 1.15e-2 (p=2) / 2.5e-3 (p=3), comparable to the discretisation error (2.6e-2 / 1.8e-3). A_l2 minimises the L2 distance of a, not of C a.

## Honest assessment: is the Coulomb-fixed A-route competitive with B_l2c?

* **Accuracy.** Only A_l2 (+Coulomb) is: p=2 it is 1.25x worse than B_l2c at n=100 (0.109 vs 0.087), equal at n~200 (0.145 vs 0.143) and better beyond (0.178 vs 0.190 at n=300, 0.214 vs 0.233 at n=400), with a smaller energy drift (+1.9% vs -11% at n=400). At p=3 it is 7x worse (2.8e-2 vs 3.9e-3 at n=100; 18.7x vs 2.6x floor) and the gauge fix makes it *worse* than plain A_l2 (1.87e-2). A_int + Coulomb: 1.2x worse than B_l2c at p=2 / n=100, 10x worse at p=3, and unstable at n>250 (p=2). A_pt: unusable. Single-remap accuracy is unchanged by the fix and A_int is as good as B_l2c in one shot (exp3), so the difference appears only over many remaps.
* **Stability.** Gauge fixing removes the gauge-growth mechanism entirely and is cheap, but only A_l2 becomes unconditionally well-behaved over the tested 400 remaps (p=2); A_int and A_pt retain instabilities of the curl part. B_l2c is stable but loses energy (-11% after 400 at p=2, -0.03%/100 at p=3). (One mesh pair, one field; thresholds depend on these.)
* **Cost.** A_int/A_l2 + Coulomb is 0.26-0.58 of B_l2c per remap in this measurement (gauge solve 3-10% of B_l2c), i.e. the A-route stays 2-4x cheaper per remap and has exact mean flux, helicity-gauge control and exact div_h B without any divergence-cleaning solve; but the gauge fix adds a *global* (AMG-CG, 14-22 iterations) solve per remap, so the A-route loses its "local/cheap" selling point; the local `jacobi:k` variant (EXPERIMENTAL) gave identical results at n=100 for k=20 and almost identical for k=5, at 1/4-1/10 of the Coulomb cost, and is the cheaper route to keep A bounded.
* Bottom line: with gauge fixing, A_l2 is competitive with B_l2c at p=2 for the tested 400 remaps (similar error, better energy, exact mean flux, cheaper) but not at p=3, where B_l2c is 7x more accurate; A_int/A_pt are not competitive. The brief's question "is the A-route worth it for AREPO (B-based)" does not become yes: the gain is bounded by B_l2c, which is just the exact constrained least-squares and needs no gauge.

## Deviations, bugs, not done

* exp3 option is `-gfix` (not `-gauge`, which was taken by the gauge-perturbation amplitude).
* Bug found and fixed during the work: AMG-CG on the singular L stalled (500 its) at p=1; fixed with the mean-free preconditioner; all reported runs were redone with the fixed solver (first-pass runs with the old solver were discarded; they showed the same physics, only the iteration counts, 8-9 instead of 14-17, were different).
* `GaugeFixStats::time` excludes the diagnostics (before/after fractions, requiring 1-2 extra solves) and uses rtol 1e-12; grad_fraction_after is evaluated with the 1e-10 rtol solve (at p=3 it returns 0 since the residual is below the CG initial-residual threshold: do not read it as exact zero).
* The p=3 CSV `g_p3` contains A_pt/A_int/A_l2 n=100 only (no n=400, p=3 long run, no p=1 repeated-remap run); jacobi variants only at p=2, n=100; no `-gauge` run on exp5 by me (the other agent uses the interface).
* The A_int late instability with gauge (n>250) is not explained. The p=3 A_l2 worsening by gauge fixing is not explained.
* Timings are from a machine with 3-4 concurrent jobs; ratios within `g_cost` only.
