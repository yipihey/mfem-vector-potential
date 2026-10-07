# E4 summary: repeated remaps (ping-pong between a uniform mesh A and a deformed mesh B)

Data: `results/exp4_repeat.csv` (782 rows, one per operator and logged step; per-part files `results/exp4_parts/`, logs `results/exp4_logs/`), driver `experiments/exp4_repeat.cpp`,
runs `scripts/run_exp4.sh`, figures `results/exp4_s1.png` (p = 1,2,3), `results/exp4_s2.png` (with mean field B0), `results/exp4_long.png` (n = 400), plots by `scripts/plot_exp4.py`, tables by `scripts/make_tables4.py`
(also `results/exp4_tables.md`).

**Set-up.** M_A = uniform N=8, M_B = f1 eps 0.3 (N=8), q = p, field `abc` (Beltrami, helicity != 0, the discrete helicity of the initial state is the reference), state at n = 0 = projection of the exact field onto M_A.
Remap n = 1 goes A -> B, n = 2 goes B -> A, etc.; "n remaps" means n transfers. The state is a (A-route; b = B0 + C_h a is recomputed after every transfer) or the total B (B-route). All tables/plots show the state **on mesh A (even n)**; the
quantities are normalised with the same-mesh reference ("ref" = projection of the exact field onto that mesh, so E/E_ref = 1 and B err/floor = 1 for a perfect remap).
Transfer plans (GSLIB point location, mass operators, AMG hierarchy) are built once per direction and reused (setup time in Table 4.6). Quadrature: nq = p+2 for A_int/B_int. 
Series: s1 = p = 1,2,3 without B0, n = 100, all seven operators (diagnostics every 5 steps); s2 = p=2 with B0 = (0.3,-0.2,0.5), n=100; s3 = p=2, n = 400 for A_pt, A_int, B_int, B_l2c (every 20 steps), plus s3b: A_l2, B_l2 at n=400 (added after
the s3 result below). The full matrix was run (n was not reduced).
Timings: 1 rank, 2 jobs concurrently on a 4-core machine that was partially occupied by another agent's MHD-ALE jobs; wall times are inflated and only ratios matter.

## Headline results

(All numbers: p = 2 unless stated. "Per 100 remaps" is taken from n = 100, Table 4.7.)

1. **Pointwise operators (A_pt, B_pt) are unstable at p >= 2; so is B_int at p=2.** The B error of A_pt/B_pt grows geometrically, ~x1.25 (A_pt) and ~x1.5-1.6 (B_pt) per remap,
   until the field is garbage within 20-30 remaps (the B error first exceeds 50% at n = 20 (A_pt, p=2), 10 (B_pt, p=2), 30 (A_pt, p=3; 15% at n=20), 20 (B_pt, p=3); logged every 5 steps); this does not depend on B0 (s2 = s1).
   A_pt's divergence stays at roundoff throughout (until the amplitude is 1e14, where the roundoff of the scaled number is 1e-9) -- exactly div-free but wrong. B_pt additionally shows div ~ O(1) after 10 remaps.
   **B_int (p=2)** is stable for ~20 remaps (B err 5.7e-2 at n=20), then blows up: 0.29 at n=50, 9.1 at n=100, geometric growth x1.07-1.08 per remap afterwards (n = 200: 1.9e4); the scaled divergence grows with it (0.37 at n=50, 13 at n=100).
   B_int at p=3 is still stable at n=100 (error 2.0e-2 = 13x floor, energy +0.3%, div 2.3e-2 and growing), at p=1 it is "stable" only in the sense of a bounded, very inaccurate state (see 5).
2. **A_int and A_l2 are not stable either, but much slower; the cause is the gauge.** The gradient fraction of a_n, ||P_grad a_n||/||a_n|| (P_grad = L2 (M_ND)-orthogonal projection onto range G_h, via a Poisson solve; initial value 8e-11 = 0),
   rises monotonically for every A-route operator, and so does the Coulomb residual ||G_h^T M a_n||: p=2: A_int 0.013 (n=10), 0.021 (20), 0.045 (50), 0.078 (100), 0.25 (200), 0.88 (300), 0.95 (400);
   A_l2 0.0067, 0.0077, 0.011, 0.016 (100), 0.025 (200), 0.041 (300), 0.070 (400); A_pt 0.047 (n=10), 0.215 (20), 0.77 (50). For A_int the long run (s3) shows the onset of an instability after ~150-200 remaps: ||a_n|| stays at 1.00 until
   n ~ 150, then 1.036 (200), 2.7 (300), 25 (400); E_B(n)/E_ref = 1.013 (100), 1.23 (200), 19 (300), 1765 (400); B error 0.14, 0.49, 4.2, 42, geometric growth x1.02 per remap; div_h b stays at roundoff (3.6e-14 at n=400) throughout.
   A_l2 (n=400): E/E_ref 1.010, 1.042, 1.15, 1.58 at n = 100,200,300,400, B error 0.12, 0.22, 0.40, 0.77: the same drift, 3-4 times slower, accelerating. So the answer to "does ||a_n||/the gauge component drift": the *total* ||a_n|| stays flat (p>=2: within 0.2% for 100 remaps)
   because the gauge part is small compared with the physical part, but the gauge part grows without bound (it is never damped, since the transfers do not map gradients to gradients and nothing projects it out), and at some point it feeds back
   into B. **No gauge fixing was applied (as specified); these results show that it would be needed for the A-route.** A Coulomb gauge projection every k remaps (a Poisson solve) is the obvious remedy but was not tested.
3. **B_l2c (L2 projection + div clean) is the only operator that is stable in all runs** (p=2: n=400; p=3: n=100): div_h at roundoff at every step (scaled 9e-15, 8.7e-15..9.9e-15 over the long run, 9e-16 at p=3), energy
   monotonically decreasing and saturating: E/E_ref = 0.9925 (n=10), 0.940 (100), 0.906 (200), 0.8896 (300), 0.8875 (400) at p=2, **-6.0% per 100 remaps initially, ~ -11% after 400**; at p=3 only -0.03% per 100 remaps (0.9997 at n=100). B error 8.7e-2 (n=100) -> 0.23 (400) at p=2 (3.8x -> 10x the floor),
   3.9e-3 at p=3 (n=100, 2.6x floor). The slice flux of the total field is *not* exactly preserved: with B0 (s2) the net flux error grows to 6.4e-4 per 100 remaps (without B0 it is 2e-13 -- there is nothing to lose when the mean is zero). B_l2 (no clean) at p=2 is stable to n=200 (E 0.968, B err 0.18, div 0.34) but at n=400 the unconstrained
   B starts to gain energy (E/E_ref = 1.117, div 0.81) -- the divergence error accumulates (0.2 at n=100) and eventually feeds back.
4. **Energy/helicity.** No operator shows *secular* gain or loss of energy that is small: per 100 remaps (p=2) A_int +1.3% and accelerating, A_l2 +1.0% accelerating (+58% at n=400), B_l2 -3.6%, B_l2c -6.0% (saturating at -11%), A_pt/B_pt/B_int blow up. At p=3: A_int +2.3%, A_l2 +0.03%, B_int +0.34%, B_l2 -0.02%, B_l2c -0.034%.
   At p=1 (under-resolved: 8 cells per wavelength, floor 22%) *every* operator destroys the field: after 100 remaps the energy is 0.002 (B_l2c), 0.012 (B_l2), 0.25 (A_int), 0.88 (A_l2), 2.1 (A_pt: gauge growth), 1.2 (B_int), 32 (B_pt) times the reference,
   and the B error is 0.94-5.4 (i.e. no information left) for all but A_l2 (0.30). A_l2 is the only operator that works at p=1 (energy -12%, B error 0.30 after 100 remaps).
   Helicity (A-route only, relative change (H_n - H_ref)/(||a||||B||)): p=2: A_int -0.49%, A_l2 -0.39% at n=100 (-1.2% at n=400 for A_l2); p=3: A_int +0.39%, A_l2 -0.0045%; p=1: A_pt -87%, A_int -96%, A_l2 -15%. Helicity is *not* preserved exactly by any operator, the
   sign and size of the drift depend on p; A_l2 at p=3 is the best (-5e-5 after 100 remaps).
5. **Accuracy after n remaps (state on mesh A, p=3, n=100; floor 1.5e-3):** B_l2c 3.9e-3 (2.6x), B_l2 4.2e-3 (2.7x), A_l2 1.9e-2 (12x), B_int 2.0e-2 (13x), A_int 0.137 (90x), A_pt/B_pt diverged. At p=2 (n=100, floor 2.3e-2): B_l2c 8.7e-2 (3.8x), B_l2 0.10 (4.5x),
   A_l2 0.12 (5.3x), A_int 0.135 (6.0x). The error of B_l2c/B_l2 grows sub-linearly in n (B_l2c at p=2: 0.087, 0.143, 0.19, 0.233 at n = 100, 200, 300, 400), that of A_int/A_l2 super-linearly (accelerating, see 2). **The B-route with L2 projection is 1.4-5x more accurate than the best A-route operator after 100 remaps (p=2: 0.087 vs 0.12; p=3: 3.9e-3 vs 1.9e-2).**
   dB_first (||B_n - B_first||/||B||, difference to the first state on the same mesh) is within ~4% of the error to the exact field for all non-diverged operators, i.e. the loss is real error, not oscillation between the two meshes.
6. **Mean field B0 (s2)** does not change anything for the fluctuating field (errors, energy, helicity, gradient fraction agree with s1 to better than 0.5%); it only changes the flux: A-route 1e-15 forever; B_l2c 6.4e-4 per 100 remaps (p=2), B_l2 4.4e-2, B_int 0.25 (blow-up), i.e. the A-route (which carries B0 separately and exactly)
   keeps the mean field exact over arbitrarily many remaps, while all B-route operators let it drift.
7. **Cost per remap (Table 4.6; plan setup excluded, contended machine).** p=2: A_pt 0.09 s, B_pt 0.05, A_int 0.31, B_l2 0.48, A_l2 0.63, B_int 0.82, B_l2c 3.3; p=3: A_pt 0.26, B_pt 0.22, A_int 1.2, A_l2 1.2, B_l2 2.1, B_int 7.4, B_l2c 12.9 (40 outer + 730 inner iterations
   per remap; the L2 projections need 15 (RT) / 25-28 (ND) CG iterations). So the only operator that is stable *and* exactly div-free over many remaps with the best accuracy (B_l2c) costs per remap 7x / 6x B_l2, 10x / 11x A_int and 5x / 10x A_l2 (p = 2 / p = 3),
   and its cost is dominated by the global Schur-complement solve (inner mass CGs); the cheapest stable-ish choices (A_l2, A_int) are 5-11x cheaper but drift into the gauge.

## Which operators are usable (as measured here, N=8, f1 0.3 <-> uniform, abc, p = 2,3)

| operator | n = 100 | n = 400 (p=2) | exact div_h b | remark |
|---|---|---|---|---|
| A_pt, B_pt | unstable (blow-up after 15-30 remaps) | -- | A_pt yes, B_pt no | do not use |
| B_int | p=2 unstable after ~50; p=3 degrading (13x floor) | unstable | no (2e-2 .. O(1)) | needs many Gauss points, expensive |
| A_int | p=2: 6x floor, p=3: 90x floor, gauge fraction 8% / 2% | unstable after ~200 | yes | gauge growth |
| A_l2 | p=2: 5x, p=3: 12x floor | slow drift (E +58%, B err 0.77) | yes | slowest A drift, best at p=1 |
| B_l2 | p=2: 4.5x, p=3: 2.7x | marginal (E +12%, div 0.8) | no (1e-2..0.2) | divergence accumulates |
| B_l2c | p=2: 3.8x, p=3: 2.6x | stable (E -11%, B err 0.23 = 10x floor) | yes (9e-15) | global solve, 6-10x costlier |

## Honest caveats

* One mesh pair (uniform <-> f1 0.3, N=8, q = p), one field family (abc Beltrami), 1 rank; the instability thresholds (number of remaps) depend on p, N, eps and the field -- the blow-up of the pointwise operators is geometric and fast, the A_int/A_l2 drift is slow and
  the A_int onset time (~150-200 at p=2) is not a prediction for other settings (at p=3 the A_int error is already 90x the floor at n=100 but gauge fraction only 2%; its blow-up time was not measured).
* The ping-pong mesh pair is one particular (smooth, moderate-distortion) pair; a real ALE rezoning sequence has meshes that differ little between steps, where every operator is much closer to the identity (that is E5, not done here).
* "p=1" is a deliberately bad case (8 cells per wavelength, 22% floor): conclusions about p=1 are about robustness for under-resolved fields.
* The quadrature rule of the integrated operators (nq = p+2) was not varied in E4; Exp. 3(e) shows A_int accuracy improves ~10-20% and B_int div improves ~2x with nq=8, which would not change the stability picture (not tested).
* GSLIB located every point in every step (0 not found, summed over all runs, in both directions).
* The long run for the L2 variants (s3b) was added after the A_int result; B_pt was not run to n=400 (it is garbage after ~15 remaps).

## Tables

### Table 4.1: p=1, N=8, abc, no B0, uniform <-> f1 0.3, state on mesh A (even n)

| op | n | E/E_ref | B err | B err / floor | div (scaled L2) | &#124;&#124;B_n-B_first&#124;&#124;/&#124;&#124;B&#124;&#124; | flux err | (H-H_ref)/(&#124;a&#124;&#124;B&#124;) | &#124;a&#124;/&#124;a_ref&#124; | grad frac | Coulomb res. | first n with B err > 0.5 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| A_pt | 10 | 0.5943 | 4.331e-01 | 1.93 | 1.2e-16 | 3.80e-01 | 2.6e-16 | -4.54e-01 | 0.7429 | 1.518e-01 | 7.84e-02 |  |
| A_pt | 20 | 0.4924 | 6.537e-01 | 2.91 | 1.1e-16 | 6.30e-01 | 5.2e-16 | -6.64e-01 | 0.5983 | 2.692e-01 | 1.11e-01 |  |
| A_pt | 50 | 0.7025 | 1.029e+00 | 4.59 | 1.2e-16 | 1.03e+00 | 9.2e-16 | -8.26e-01 | 0.524 | 5.594e-01 | 2.14e-01 |  |
| A_pt | 100 | 2.1 | 1.602e+00 | 7.14 | 2.0e-16 | 1.63e+00 | 3.9e-16 | -8.65e-01 | 1.063 | 8.699e-01 | 7.24e-01 | 15 |
| A_int | 10 | 0.5764 | 3.980e-01 | 1.77 | 1.2e-16 | 3.37e-01 | 5.4e-16 | -4.47e-01 | 0.7378 | 7.639e-02 | 3.79e-02 |  |
| A_int | 20 | 0.3879 | 5.593e-01 | 2.49 | 9.4e-17 | 5.26e-01 | 6.6e-16 | -6.67e-01 | 0.5693 | 1.095e-01 | 3.99e-02 |  |
| A_int | 50 | 0.2385 | 8.200e-01 | 3.65 | 7.9e-17 | 8.09e-01 | 4.9e-16 | -8.92e-01 | 0.3326 | 1.797e-01 | 4.23e-02 |  |
| A_int | 100 | 0.2475 | 9.586e-01 | 4.27 | 7.7e-17 | 9.56e-01 | 9.4e-16 | -9.57e-01 | 0.281 | 4.916e-01 | 1.10e-01 | 20 |
| A_l2 | 10 | 0.9847 | 2.392e-01 | 1.07 | 1.5e-16 | 8.50e-02 | 4.3e-16 | -2.14e-02 | 0.9896 | 3.179e-02 | 2.07e-02 |  |
| A_l2 | 20 | 0.9716 | 2.510e-01 | 1.12 | 1.5e-16 | 1.15e-01 | 1.9e-16 | -3.90e-02 | 0.981 | 4.302e-02 | 2.66e-02 |  |
| A_l2 | 50 | 0.9346 | 2.751e-01 | 1.23 | 1.5e-16 | 1.63e-01 | 1.0e-15 | -8.48e-02 | 0.958 | 5.590e-02 | 3.20e-02 |  |
| A_l2 | 100 | 0.8793 | 3.049e-01 | 1.36 | 1.6e-16 | 2.12e-01 | 3.9e-16 | -1.51e-01 | 0.9245 | 6.111e-02 | 3.32e-02 | >100 |
| B_pt | 10 | 1.199 | 4.688e-01 | 2.09 | 4.0e-01 | 4.22e-01 | 2.2e-02 | n/a | n/a | n/a | n/a |  |
| B_pt | 20 | 1.662 | 8.046e-01 | 3.59 | 7.9e-01 | 7.93e-01 | 3.1e-02 | n/a | n/a | n/a | n/a |  |
| B_pt | 50 | 4.7 | 1.862e+00 | 8.3 | 1.9e+00 | 1.90e+00 | 6.1e-02 | n/a | n/a | n/a | n/a |  |
| B_pt | 100 | 32.49 | 5.430e+00 | 24.2 | 5.7e+00 | 5.57e+00 | 1.9e-02 | n/a | n/a | n/a | n/a | 15 |
| B_int | 10 | 0.5614 | 4.067e-01 | 1.81 | 8.0e-02 | 3.48e-01 | 3.5e-01 | n/a | n/a | n/a | n/a |  |
| B_int | 20 | 0.3743 | 5.753e-01 | 2.56 | 7.6e-02 | 5.44e-01 | 4.3e-01 | n/a | n/a | n/a | n/a |  |
| B_int | 50 | 0.3146 | 8.649e-01 | 3.85 | 5.4e-02 | 8.57e-01 | 3.4e-01 | n/a | n/a | n/a | n/a |  |
| B_int | 100 | 1.212 | 1.275e+00 | 5.68 | 6.4e-02 | 1.29e+00 | 1.5e-01 | n/a | n/a | n/a | n/a | 15 |
| B_l2 | 10 | 0.5264 | 3.705e-01 | 1.65 | 1.7e-01 | 3.02e-01 | 5.7e-02 | n/a | n/a | n/a | n/a |  |
| B_l2 | 20 | 0.2975 | 5.333e-01 | 2.38 | 2.5e-01 | 4.96e-01 | 2.1e-03 | n/a | n/a | n/a | n/a |  |
| B_l2 | 50 | 0.07226 | 7.987e-01 | 3.56 | 2.5e-01 | 7.87e-01 | 1.5e-01 | n/a | n/a | n/a | n/a |  |
| B_l2 | 100 | 0.01238 | 9.376e-01 | 4.18 | 1.5e-01 | 9.34e-01 | 2.8e-01 | n/a | n/a | n/a | n/a | 20 |
| B_l2c | 10 | 0.5114 | 3.643e-01 | 1.62 | 1.4e-14 | 2.94e-01 | 2.2e-14 | n/a | n/a | n/a | n/a |  |
| B_l2c | 20 | 0.2659 | 5.307e-01 | 2.37 | 9.9e-15 | 4.94e-01 | 1.6e-14 | n/a | n/a | n/a | n/a |  |
| B_l2c | 50 | 0.03858 | 8.200e-01 | 3.65 | 3.6e-15 | 8.09e-01 | 7.0e-15 | n/a | n/a | n/a | n/a |  |
| B_l2c | 100 | 0.001711 | 9.638e-01 | 4.29 | 7.9e-16 | 9.62e-01 | 4.4e-15 | n/a | n/a | n/a | n/a | 20 |

### Table 4.2: p=2, N=8, abc, no B0, uniform <-> f1 0.3, state on mesh A (even n)

| op | n | E/E_ref | B err | B err / floor | div (scaled L2) | &#124;&#124;B_n-B_first&#124;&#124;/&#124;&#124;B&#124;&#124; | flux err | (H-H_ref)/(&#124;a&#124;&#124;B&#124;) | &#124;a&#124;/&#124;a_ref&#124; | grad frac | Coulomb res. | first n with B err > 0.5 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| A_pt | 10 | 1.024 | 1.419e-01 | 6.23 | 1.2e-15 | 1.40e-01 | 6.1e-16 | +4.59e-03 | 1.004 | 4.703e-02 | 3.14e-02 |  |
| A_pt | 20 | 1.653 | 7.981e-01 | 35 | 1.3e-15 | 7.98e-01 | 3.1e-16 | +1.55e-02 | 1.048 | 2.150e-01 | 1.26e-01 |  |
| A_pt | 50 | 2.272e+04 | 1.507e+02 | 6.61e+03 | 8.0e-14 | 1.51e+02 | 1.1e-15 | +1.50e+01 | 42.34 | 7.703e-01 | 2.13e+01 |  |
| A_pt | 100 | 1.177e+14 | 1.085e+07 | 4.76e+08 | 5.7e-09 | 1.08e+07 | 8.0e-11 | +4.93e+10 | 3.139e+06 | 7.751e-01 | 1.68e+06 | 20 |
| A_int | 10 | 0.9988 | 3.440e-02 | 1.51 | 1.2e-15 | 2.56e-02 | 4.5e-16 | -1.93e-03 | 0.9992 | 1.293e-02 | 9.81e-03 |  |
| A_int | 20 | 0.9989 | 4.143e-02 | 1.82 | 1.2e-15 | 3.45e-02 | 4.0e-16 | -2.36e-03 | 0.9991 | 2.084e-02 | 1.31e-02 |  |
| A_int | 50 | 1.001 | 7.166e-02 | 3.14 | 1.2e-15 | 6.79e-02 | 5.3e-16 | -3.47e-03 | 0.9995 | 4.453e-02 | 2.37e-02 |  |
| A_int | 100 | 1.013 | 1.355e-01 | 5.95 | 1.2e-15 | 1.34e-01 | 9.2e-16 | -4.88e-03 | 1.001 | 7.812e-02 | 4.12e-02 | >100 |
| A_l2 | 10 | 1.001 | 4.163e-02 | 1.83 | 1.2e-15 | 3.48e-02 | 3.4e-16 | -5.69e-04 | 0.9998 | 6.733e-03 | 5.48e-03 |  |
| A_l2 | 20 | 1.001 | 5.234e-02 | 2.3 | 1.2e-15 | 4.71e-02 | 6.4e-16 | -9.92e-04 | 0.9996 | 7.747e-03 | 5.82e-03 |  |
| A_l2 | 50 | 1.004 | 7.881e-02 | 3.46 | 1.2e-15 | 7.55e-02 | 3.7e-16 | -2.17e-03 | 0.9991 | 1.096e-02 | 7.05e-03 |  |
| A_l2 | 100 | 1.01 | 1.197e-01 | 5.25 | 1.2e-15 | 1.17e-01 | 5.0e-16 | -3.90e-03 | 0.9986 | 1.556e-02 | 9.38e-03 | >100 |
| B_pt | 10 | 1.612 | 5.252e-01 | 23 | 9.2e-01 | 5.25e-01 | 3.6e-02 | n/a | n/a | n/a | n/a |  |
| B_pt | 20 | 1840 | 4.281e+01 | 1.88e+03 | 7.1e+01 | 4.28e+01 | 4.9e-01 | n/a | n/a | n/a | n/a |  |
| B_pt | 50 | 1.081e+16 | 1.040e+08 | 4.56e+09 | 2.1e+08 | 1.04e+08 | 2.9e+06 | n/a | n/a | n/a | n/a |  |
| B_pt | 100 | 2.58e+37 | 5.079e+18 | 2.23e+20 | 1.1e+19 | 5.08e+18 | 1.8e+17 | n/a | n/a | n/a | n/a | 10 |
| B_int | 10 | 1.002 | 3.713e-02 | 1.63 | 3.5e-02 | 2.92e-02 | 5.8e-03 | n/a | n/a | n/a | n/a |  |
| B_int | 20 | 1.01 | 5.735e-02 | 2.52 | 7.0e-02 | 5.26e-02 | 1.1e-02 | n/a | n/a | n/a | n/a |  |
| B_int | 50 | 1.14 | 2.907e-01 | 12.8 | 3.7e-01 | 2.90e-01 | 1.4e-02 | n/a | n/a | n/a | n/a |  |
| B_int | 100 | 83.81 | 9.054e+00 | 397 | 1.3e+01 | 9.05e+00 | 2.5e-01 | n/a | n/a | n/a | n/a | 60 |
| B_l2 | 10 | 0.9937 | 3.104e-02 | 1.36 | 4.6e-02 | 2.09e-02 | 2.6e-03 | n/a | n/a | n/a | n/a |  |
| B_l2 | 20 | 0.9881 | 3.788e-02 | 1.66 | 6.9e-02 | 3.02e-02 | 5.0e-04 | n/a | n/a | n/a | n/a |  |
| B_l2 | 50 | 0.9753 | 6.281e-02 | 2.76 | 1.2e-01 | 5.85e-02 | 1.5e-02 | n/a | n/a | n/a | n/a |  |
| B_l2 | 100 | 0.9636 | 1.027e-01 | 4.51 | 2.0e-01 | 1.00e-01 | 4.3e-02 | n/a | n/a | n/a | n/a | >100 |
| B_l2c | 10 | 0.9925 | 2.988e-02 | 1.31 | 8.9e-15 | 1.92e-02 | 2.2e-14 | n/a | n/a | n/a | n/a |  |
| B_l2c | 20 | 0.9852 | 3.460e-02 | 1.52 | 8.9e-15 | 2.59e-02 | 4.2e-14 | n/a | n/a | n/a | n/a |  |
| B_l2c | 50 | 0.9659 | 5.427e-02 | 2.38 | 8.8e-15 | 4.92e-02 | 9.9e-14 | n/a | n/a | n/a | n/a |  |
| B_l2c | 100 | 0.9403 | 8.729e-02 | 3.83 | 8.7e-15 | 8.43e-02 | 2.0e-13 | n/a | n/a | n/a | n/a | >100 |

### Table 4.3: p=3, N=8, abc, no B0, uniform <-> f1 0.3, state on mesh A (even n)

| op | n | E/E_ref | B err | B err / floor | div (scaled L2) | &#124;&#124;B_n-B_first&#124;&#124;/&#124;&#124;B&#124;&#124; | flux err | (H-H_ref)/(&#124;a&#124;&#124;B&#124;) | &#124;a&#124;/&#124;a_ref&#124; | grad frac | Coulomb res. | first n with B err > 0.5 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| A_pt | 10 | 1.001 | 2.228e-02 | 14.7 | 3.2e-15 | 2.22e-02 | 7.1e-16 | +7.62e-05 | 1 | 3.297e-03 | 2.20e-03 |  |
| A_pt | 20 | 1.022 | 1.472e-01 | 97.1 | 3.2e-15 | 1.47e-01 | 5.7e-16 | +5.58e-05 | 1 | 1.450e-02 | 8.79e-03 |  |
| A_pt | 50 | 1.724e+04 | 1.313e+02 | 8.67e+04 | 1.2e-13 | 1.31e+02 | 1.7e-15 | -4.04e+00 | 15.61 | 6.492e-01 | 6.08e+00 |  |
| A_pt | 100 | 2.644e+14 | 1.626e+07 | 1.07e+10 | 1.4e-08 | 1.63e+07 | 3.1e-10 | -8.14e+10 | 1.953e+06 | 6.549e-01 | 7.72e+05 | 30 |
| A_int | 10 | 1 | 5.422e-03 | 3.58 | 3.1e-15 | 5.21e-03 | 4.0e-16 | +4.16e-05 | 1 | 1.856e-03 | 1.23e-03 |  |
| A_int | 20 | 1 | 8.056e-03 | 5.32 | 3.1e-15 | 7.91e-03 | 5.5e-16 | +2.23e-04 | 1 | 3.091e-03 | 1.96e-03 |  |
| A_int | 50 | 1.002 | 2.119e-02 | 14 | 3.2e-15 | 2.11e-02 | 1.0e-15 | +1.09e-03 | 1.001 | 7.338e-03 | 4.48e-03 |  |
| A_int | 100 | 1.023 | 1.370e-01 | 90.4 | 3.2e-15 | 1.37e-01 | 7.2e-16 | +3.92e-03 | 1.002 | 2.238e-02 | 1.45e-02 | >100 |
| A_l2 | 10 | 1 | 7.291e-03 | 4.81 | 3.2e-15 | 7.14e-03 | 5.7e-16 | -8.28e-06 | 1 | 1.092e-03 | 7.05e-04 |  |
| A_l2 | 20 | 1 | 9.530e-03 | 6.29 | 3.2e-15 | 9.42e-03 | 5.0e-16 | -1.62e-05 | 1 | 1.391e-03 | 8.50e-04 |  |
| A_l2 | 50 | 1 | 1.360e-02 | 8.98 | 3.2e-15 | 1.35e-02 | 2.5e-16 | -3.31e-05 | 1 | 1.721e-03 | 1.00e-03 |  |
| A_l2 | 100 | 1 | 1.867e-02 | 12.3 | 3.2e-15 | 1.86e-02 | 3.4e-16 | -4.47e-05 | 1 | 1.960e-03 | 1.09e-03 | >100 |
| B_pt | 10 | 1.005 | 1.731e-02 | 11.4 | 6.0e-02 | 1.72e-02 | 2.7e-03 | n/a | n/a | n/a | n/a |  |
| B_pt | 20 | 1.586 | 6.866e-01 | 453 | 2.7e+00 | 6.87e-01 | 1.1e-01 | n/a | n/a | n/a | n/a |  |
| B_pt | 50 | 3.481e+10 | 1.866e+05 | 1.23e+08 | 7.2e+05 | 1.87e+05 | 4.6e+04 | n/a | n/a | n/a | n/a |  |
| B_pt | 100 | 1.019e+29 | 3.192e+14 | 2.11e+17 | 1.2e+15 | 3.19e+14 | 8.1e+13 | n/a | n/a | n/a | n/a | 20 |
| B_int | 10 | 1 | 3.330e-03 | 2.2 | 2.2e-03 | 2.97e-03 | 3.1e-04 | n/a | n/a | n/a | n/a |  |
| B_int | 20 | 1 | 4.525e-03 | 2.99 | 3.2e-03 | 4.27e-03 | 3.8e-04 | n/a | n/a | n/a | n/a |  |
| B_int | 50 | 1.001 | 9.086e-03 | 6 | 7.1e-03 | 8.96e-03 | 2.6e-04 | n/a | n/a | n/a | n/a |  |
| B_int | 100 | 1.003 | 2.002e-02 | 13.2 | 2.3e-02 | 2.00e-02 | 4.3e-03 | n/a | n/a | n/a | n/a | >100 |
| B_l2 | 10 | 1 | 2.915e-03 | 1.92 | 4.6e-03 | 2.49e-03 | 1.2e-03 | n/a | n/a | n/a | n/a |  |
| B_l2 | 20 | 1 | 3.226e-03 | 2.13 | 5.6e-03 | 2.85e-03 | 1.4e-03 | n/a | n/a | n/a | n/a |  |
| B_l2 | 50 | 0.9999 | 3.625e-03 | 2.39 | 7.8e-03 | 3.29e-03 | 1.8e-03 | n/a | n/a | n/a | n/a |  |
| B_l2 | 100 | 0.9998 | 4.156e-03 | 2.74 | 1.1e-02 | 3.87e-03 | 2.5e-03 | n/a | n/a | n/a | n/a | >100 |
| B_l2c | 10 | 1 | 2.922e-03 | 1.93 | 8.8e-16 | 2.50e-03 | 1.3e-14 | n/a | n/a | n/a | n/a |  |
| B_l2c | 20 | 1 | 3.228e-03 | 2.13 | 8.9e-16 | 2.85e-03 | 2.6e-14 | n/a | n/a | n/a | n/a |  |
| B_l2c | 50 | 0.9999 | 3.561e-03 | 2.35 | 9.0e-16 | 3.22e-03 | 6.5e-14 | n/a | n/a | n/a | n/a |  |
| B_l2c | 100 | 0.9997 | 3.942e-03 | 2.6 | 9.0e-16 | 3.64e-03 | 1.3e-13 | n/a | n/a | n/a | n/a | >100 |

### Table 4.4: p=2, N=8, abc WITH mean field B0 = (0.3,-0.2,0.5)

| op | n | E/E_ref | B err | B err / floor | div (scaled L2) | &#124;&#124;B_n-B_first&#124;&#124;/&#124;&#124;B&#124;&#124; | flux err | (H-H_ref)/(&#124;a&#124;&#124;B&#124;) | &#124;a&#124;/&#124;a_ref&#124; | grad frac | Coulomb res. | first n with B err > 0.5 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| A_pt | 10 | 1.024 | 1.416e-01 | 6.23 | 1.2e-15 | 1.40e-01 | 5.0e-16 | +4.59e-03 | 1.004 | 4.703e-02 | 3.14e-02 |  |
| A_pt | 20 | 1.65 | 7.961e-01 | 35 | 1.3e-15 | 7.96e-01 | 6.7e-16 | +1.55e-02 | 1.048 | 2.150e-01 | 1.26e-01 |  |
| A_pt | 50 | 2.261e+04 | 1.503e+02 | 6.61e+03 | 8.2e-14 | 1.50e+02 | 3.9e-15 | +1.50e+01 | 42.34 | 7.703e-01 | 2.13e+01 |  |
| A_pt | 100 | 1.171e+14 | 1.082e+07 | 4.76e+08 | 5.7e-09 | 1.08e+07 | 8.0e-11 | +4.93e+10 | 3.139e+06 | 7.751e-01 | 1.68e+06 | 20 |
| A_int | 10 | 0.9988 | 3.431e-02 | 1.51 | 1.2e-15 | 2.56e-02 | 8.9e-16 | -1.93e-03 | 0.9992 | 1.293e-02 | 9.81e-03 |  |
| A_int | 20 | 0.9989 | 4.133e-02 | 1.82 | 1.3e-15 | 3.44e-02 | 4.4e-16 | -2.36e-03 | 0.9991 | 2.084e-02 | 1.31e-02 |  |
| A_int | 50 | 1.001 | 7.149e-02 | 3.14 | 1.3e-15 | 6.77e-02 | 7.2e-16 | -3.47e-03 | 0.9995 | 4.453e-02 | 2.37e-02 |  |
| A_int | 100 | 1.013 | 1.352e-01 | 5.95 | 1.2e-15 | 1.33e-01 | 1.3e-15 | -4.88e-03 | 1.001 | 7.812e-02 | 4.12e-02 | >100 |
| A_l2 | 10 | 1.001 | 4.153e-02 | 1.83 | 1.2e-15 | 3.47e-02 | 4.4e-16 | -5.69e-04 | 0.9998 | 6.733e-03 | 5.48e-03 |  |
| A_l2 | 20 | 1.001 | 5.221e-02 | 2.3 | 1.2e-15 | 4.70e-02 | 8.0e-16 | -9.92e-04 | 0.9996 | 7.747e-03 | 5.82e-03 |  |
| A_l2 | 50 | 1.004 | 7.862e-02 | 3.46 | 1.3e-15 | 7.53e-02 | 5.6e-16 | -2.17e-03 | 0.9991 | 1.096e-02 | 7.05e-03 |  |
| A_l2 | 100 | 1.01 | 1.194e-01 | 5.25 | 1.2e-15 | 1.17e-01 | 3.3e-16 | -3.90e-03 | 0.9986 | 1.556e-02 | 9.38e-03 | >100 |
| B_pt | 10 | 1.609 | 5.240e-01 | 23 | 9.2e-01 | 5.24e-01 | 3.6e-02 | n/a | n/a | n/a | n/a |  |
| B_pt | 20 | 1832 | 4.271e+01 | 1.88e+03 | 7.1e+01 | 4.27e+01 | 5.0e-01 | n/a | n/a | n/a | n/a |  |
| B_pt | 50 | 1.076e+16 | 1.037e+08 | 4.56e+09 | 2.1e+08 | 1.04e+08 | 2.9e+06 | n/a | n/a | n/a | n/a |  |
| B_pt | 100 | 2.568e+37 | 5.067e+18 | 2.23e+20 | 1.1e+19 | 5.07e+18 | 1.9e+17 | n/a | n/a | n/a | n/a | 10 |
| B_int | 10 | 1.002 | 3.704e-02 | 1.63 | 3.5e-02 | 2.91e-02 | 5.8e-03 | n/a | n/a | n/a | n/a |  |
| B_int | 20 | 1.01 | 5.722e-02 | 2.52 | 7.0e-02 | 5.24e-02 | 1.1e-02 | n/a | n/a | n/a | n/a |  |
| B_int | 50 | 1.139 | 2.900e-01 | 12.8 | 3.7e-01 | 2.89e-01 | 1.4e-02 | n/a | n/a | n/a | n/a |  |
| B_int | 100 | 83.42 | 9.032e+00 | 397 | 1.3e+01 | 9.03e+00 | 2.5e-01 | n/a | n/a | n/a | n/a | 60 |
| B_l2 | 10 | 0.9937 | 3.097e-02 | 1.36 | 4.5e-02 | 2.09e-02 | 2.7e-03 | n/a | n/a | n/a | n/a |  |
| B_l2 | 20 | 0.9882 | 3.779e-02 | 1.66 | 6.9e-02 | 3.01e-02 | 7.8e-04 | n/a | n/a | n/a | n/a |  |
| B_l2 | 50 | 0.9754 | 6.266e-02 | 2.76 | 1.2e-01 | 5.84e-02 | 1.5e-02 | n/a | n/a | n/a | n/a |  |
| B_l2 | 100 | 0.9637 | 1.024e-01 | 4.51 | 2.0e-01 | 9.99e-02 | 4.4e-02 | n/a | n/a | n/a | n/a | >100 |
| B_l2c | 10 | 0.9926 | 2.981e-02 | 1.31 | 8.9e-15 | 1.91e-02 | 6.4e-05 | n/a | n/a | n/a | n/a |  |
| B_l2c | 20 | 0.9853 | 3.452e-02 | 1.52 | 8.8e-15 | 2.59e-02 | 1.3e-04 | n/a | n/a | n/a | n/a |  |
| B_l2c | 50 | 0.9661 | 5.414e-02 | 2.38 | 8.8e-15 | 4.91e-02 | 3.2e-04 | n/a | n/a | n/a | n/a |  |
| B_l2c | 100 | 0.9406 | 8.708e-02 | 3.83 | 8.7e-15 | 8.41e-02 | 6.4e-04 | n/a | n/a | n/a | n/a | >100 |

### Table 4.5: long run, p=2, N=8, no B0

| op | n | E/E_ref | B err | B err / floor | div (scaled L2) | &#124;&#124;B_n-B_first&#124;&#124;/&#124;&#124;B&#124;&#124; | flux err | (H-H_ref)/(&#124;a&#124;&#124;B&#124;) | &#124;a&#124;/&#124;a_ref&#124; | grad frac | Coulomb res. | first n with B err > 0.5 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| A_pt | 100 | 1.177e+14 | 1.085e+07 | 4.76e+08 | 5.7e-09 | 1.08e+07 | 8.0e-11 | +4.93e+10 | 3.139e+06 | 7.751e-01 | 1.68e+06 |  |
| A_pt | 200 | 6.351e+34 | 2.520e+17 | 1.11e+19 | 1.3e+02 | 2.52e+17 | 1.8e+00 | +5.33e+31 | 6.461e+16 | 7.168e-01 | 3.30e+16 |  |
| A_pt | 300 | 3.648e+55 | 6.039e+27 | 2.65e+29 | 3.1e+12 | 6.04e+27 | 3.9e+10 | +3.87e+52 | 1.526e+27 | 7.086e-01 | 7.72e+26 |  |
| A_pt | 400 | 2.101e+76 | 1.449e+38 | 6.36e+39 | 7.4e+22 | 1.45e+38 | 1.0e+21 | +2.28e+73 | 3.652e+37 | 7.068e-01 | 1.84e+37 | 20 |
| A_int | 100 | 1.013 | 1.355e-01 | 5.95 | 1.2e-15 | 1.34e-01 | 9.2e-16 | -4.88e-03 | 1.001 | 7.812e-02 | 4.12e-02 |  |
| A_int | 200 | 1.228 | 4.888e-01 | 21.4 | 1.3e-15 | 4.88e-01 | 4.9e-16 | -8.61e-03 | 1.036 | 2.549e-01 | 1.70e-01 |  |
| A_int | 300 | 18.79 | 4.225e+00 | 185 | 3.8e-15 | 4.23e+00 | 1.0e-15 | -5.47e-02 | 2.709 | 8.844e-01 | 1.61e+00 |  |
| A_int | 400 | 1765 | 4.200e+01 | 1.84e+03 | 3.6e-14 | 4.20e+01 | 3.6e-15 | -1.57e+00 | 25.2 | 9.470e-01 | 1.60e+01 | 220 |
| A_l2 | 100 | 1.01 | 1.197e-01 | 5.25 | 1.2e-15 | 1.17e-01 | 5.0e-16 | -3.90e-03 | 0.9986 | 1.556e-02 | 9.38e-03 |  |
| A_l2 | 200 | 1.042 | 2.206e-01 | 9.68 | 1.2e-15 | 2.19e-01 | 5.1e-16 | -6.65e-03 | 0.9984 | 2.535e-02 | 1.56e-02 |  |
| A_l2 | 300 | 1.15 | 3.976e-01 | 17.4 | 1.2e-15 | 3.97e-01 | 4.6e-16 | -8.82e-03 | 1.001 | 4.113e-02 | 2.59e-02 |  |
| A_l2 | 400 | 1.582 | 7.680e-01 | 33.7 | 1.3e-15 | 7.68e-01 | 5.0e-16 | -1.17e-02 | 1.012 | 6.975e-02 | 4.40e-02 | 340 |
| B_int | 100 | 83.81 | 9.054e+00 | 397 | 1.3e+01 | 9.05e+00 | 2.5e-01 | n/a | n/a | n/a | n/a |  |
| B_int | 200 | 3.612e+08 | 1.900e+04 | 8.34e+05 | 2.8e+04 | 1.90e+04 | 8.0e+02 | n/a | n/a | n/a | n/a |  |
| B_int | 300 | 1.698e+15 | 4.120e+07 | 1.81e+09 | 6.2e+07 | 4.12e+07 | 1.7e+06 | n/a | n/a | n/a | n/a |  |
| B_int | 400 | 8.046e+21 | 8.969e+10 | 3.93e+12 | 1.4e+11 | 8.97e+10 | 3.6e+09 | n/a | n/a | n/a | n/a | 60 |
| B_l2 | 100 | 0.9636 | 1.027e-01 | 4.51 | 2.0e-01 | 1.00e-01 | 4.3e-02 | n/a | n/a | n/a | n/a |  |
| B_l2 | 200 | 0.9677 | 1.759e-01 | 7.72 | 3.4e-01 | 1.74e-01 | 1.0e-01 | n/a | n/a | n/a | n/a |  |
| B_l2 | 300 | 1.011 | 2.584e-01 | 11.3 | 5.3e-01 | 2.57e-01 | 1.6e-01 | n/a | n/a | n/a | n/a |  |
| B_l2 | 400 | 1.117 | 3.747e-01 | 16.4 | 8.1e-01 | 3.74e-01 | 2.0e-01 | n/a | n/a | n/a | n/a | >400 |
| B_l2c | 100 | 0.9403 | 8.729e-02 | 3.83 | 8.7e-15 | 8.43e-02 | 2.0e-13 | n/a | n/a | n/a | n/a |  |
| B_l2c | 200 | 0.9064 | 1.426e-01 | 6.26 | 8.9e-15 | 1.41e-01 | 3.9e-13 | n/a | n/a | n/a | n/a |  |
| B_l2c | 300 | 0.8896 | 1.895e-01 | 8.31 | 9.3e-15 | 1.88e-01 | 5.8e-13 | n/a | n/a | n/a | n/a |  |
| B_l2c | 400 | 0.8875 | 2.332e-01 | 10.2 | 9.9e-15 | 2.32e-01 | 7.8e-13 | n/a | n/a | n/a | n/a | >400 |

### Table 4.6: cost per remap (transfer + composition, without diagnostics; contended machine) and CG iterations per remap

| series | op | s/remap | outer or L2 iterations / remap | inner iterations / remap | plan setup [s] |
|---|---|---|---|---|---|
| s1_p1 | A_pt | 0.009 | - | - | 0.06 |
| s1_p1 | A_int | 0.024 | - | - | 0.13 |
| s1_p1 | A_l2 | 0.101 | 25.0 | - | 0.18 |
| s1_p1 | B_pt | 0.004 | - | - | 0.02 |
| s1_p1 | B_int | 0.036 | - | - | 0.18 |
| s1_p1 | B_l2 | 0.076 | 14.5 | - | 0.19 |
| s1_p1 | B_l2c | 0.815 | 24.5 | 426 | 0.15 |
| s1_p2 | A_pt | 0.087 | - | - | 0.30 |
| s1_p2 | A_int | 0.311 | - | - | 1.20 |
| s1_p2 | A_l2 | 0.628 | 27.4 | - | 0.99 |
| s1_p2 | B_pt | 0.052 | - | - | 0.20 |
| s1_p2 | B_int | 0.821 | - | - | 2.97 |
| s1_p2 | B_l2 | 0.482 | 15.0 | - | 1.00 |
| s1_p2 | B_l2c | 3.266 | 29.9 | 562 | 0.95 |
| s1_p3 | A_pt | 0.259 | - | - | 0.89 |
| s1_p3 | A_int | 1.176 | - | - | 3.05 |
| s1_p3 | A_l2 | 1.240 | 25.0 | - | 1.74 |
| s1_p3 | B_pt | 0.224 | - | - | 0.44 |
| s1_p3 | B_int | 7.390 | - | - | 23.04 |
| s1_p3 | B_l2 | 2.081 | 15.0 | - | 3.34 |
| s1_p3 | B_l2c | 12.946 | 40.4 | 730 | 3.55 |
| s3 | A_pt | 0.081 | - | - | 0.30 |
| s3 | A_int | 0.332 | - | - | 1.05 |
| s3 | A_l2 | 0.612 | 27.7 | - | 1.06 |
| s3 | B_int | 0.809 | - | - | 3.04 |
| s3 | B_l2 | 0.493 | 15.0 | - | 1.00 |
| s3 | B_l2c | 3.432 | 29.8 | 559 | 0.85 |

### Table 4.7: drift per 100 remaps (state on mesh A), from the long run (n=400, p=2) and the n=100 runs

| series | op | (E_n/E_0-1) per 100 | B err growth per 100 | grad frac at end | &#124;a_n&#124;/&#124;a_0&#124; at end |
|---|---|---|---|---|---|
| s1_p1 | A_pt | +1.10e+00 | +1.38e+00 | 8.699e-01 | 1.063 |
| s1_p1 | A_int | -7.53e-01 | +7.34e-01 | 4.916e-01 | 0.281 |
| s1_p1 | A_l2 | -1.21e-01 | +8.05e-02 | 6.111e-02 | 0.9245 |
| s1_p1 | B_pt | diverged (E/E_0 = 32, B err = 5.4 at n=100) | | | |
| s1_p1 | B_int | +2.12e-01 | +1.05e+00 | n/a | n/a |
| s1_p1 | B_l2 | -9.88e-01 | +7.13e-01 | n/a | n/a |
| s1_p1 | B_l2c | -9.98e-01 | +7.39e-01 | n/a | n/a |
| s1_p2 | A_pt | diverged (E/E_0 = 1.2e+14, B err = 1.1e+07 at n=100) | | | |
| s1_p2 | A_int | +1.28e-02 | +1.13e-01 | 7.812e-02 | 1.001 |
| s1_p2 | A_l2 | +9.96e-03 | +9.69e-02 | 1.556e-02 | 0.9986 |
| s1_p2 | B_pt | diverged (E/E_0 = 2.6e+37, B err = 5.1e+18 at n=100) | | | |
| s1_p2 | B_int | diverged (E/E_0 = 84, B err = 9.1 at n=100) | | | |
| s1_p2 | B_l2 | -3.64e-02 | +7.99e-02 | n/a | n/a |
| s1_p2 | B_l2c | -5.97e-02 | +6.45e-02 | n/a | n/a |
| s1_p3 | A_pt | diverged (E/E_0 = 2.6e+14, B err = 1.6e+07 at n=100) | | | |
| s1_p3 | A_int | +2.26e-02 | +1.35e-01 | 2.238e-02 | 1.002 |
| s1_p3 | A_l2 | +2.99e-04 | +1.71e-02 | 1.960e-03 | 1 |
| s1_p3 | B_pt | diverged (E/E_0 = 1e+29, B err = 3.2e+14 at n=100) | | | |
| s1_p3 | B_int | +3.44e-03 | +1.85e-02 | n/a | n/a |
| s1_p3 | B_l2 | -2.00e-04 | +2.64e-03 | n/a | n/a |
| s1_p3 | B_l2c | -3.38e-04 | +2.43e-03 | n/a | n/a |
| s2 | A_pt | diverged (E/E_0 = 1.2e+14, B err = 1.1e+07 at n=100) | | | |
| s2 | A_int | +1.27e-02 | +1.12e-01 | 7.812e-02 | 1.001 |
| s2 | A_l2 | +9.92e-03 | +9.66e-02 | 1.556e-02 | 0.9986 |
| s2 | B_pt | diverged (E/E_0 = 2.6e+37, B err = 5.1e+18 at n=100) | | | |
| s2 | B_int | diverged (E/E_0 = 83, B err = 9 at n=100) | | | |
| s2 | B_l2 | -3.63e-02 | +7.97e-02 | n/a | n/a |
| s2 | B_l2c | -5.94e-02 | +6.43e-02 | n/a | n/a |
| s3 | A_pt | diverged (E/E_0 = 2.1e+76, B err = 1.4e+38 at n=400) | | | |
| s3 | A_int | diverged (E/E_0 = 1.8e+03, B err = 42 at n=400) | | | |
| s3 | A_l2 | +1.46e-01 | +1.86e-01 | 6.975e-02 | 1.012 |
| s3 | B_int | diverged (E/E_0 = 8e+21, B err = 9e+10 at n=400) | | | |
| s3 | B_l2 | +2.93e-02 | +8.80e-02 | n/a | n/a |
| s3 | B_l2c | -2.81e-02 | +5.26e-02 | n/a | n/a |
