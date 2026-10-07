# E3 summary: one remap M1 -> M2, A-route versus B-route

Data: `results/exp3_remap.csv` (258 rows = one per transfer operator and configuration; per-part files in `results/exp3_parts/`, logs in `results/exp3_logs/`),
driver `experiments/exp3_remap.cpp`, runs `scripts/run_exp3.sh`, figures `results/exp3_*.png` (`scripts/plot_exp3.py`), tables below by `scripts/make_tables3.py`
(also in `results/exp3_tables.md`). Library: `src/vp_transfer.{hpp,cpp}`; sanity test `bin/transfer_test` (see "Verification").

**Timing caveat.** All runs were 1 MPI rank, two at a time, on a 4-core VM that was shared with another agent's 2-rank MHD-ALE jobs (load average 2.5-4), so every
wall time below is inflated by an unknown factor (roughly 1.5-3x judging from user-time vs wall-time) and only the ratios between operators are meaningful.
Peak RSS of any run: 3.8 GB (p=3, N=16).

## Set-up (what exactly is compared)

* M1 = periodic unit box, N^3 hexes, geometric order q = p, deformed x' = x + eps f1(x). The A-route state is a1 = Pi_ND^int(A_exact) (integrated dofs, 6 Gauss points per
  sub-edge), b1 = Pi_RT^int(B0) + C_h a1. The B-route state is the same b1 (total B including the mean field). So both routes start from the *identical* discrete field.
* M2 in {uniform (a), f2 eps 0.15 (b), uniform N2 = 12, 16 (d)}. After the transfer the A-route composes b2 = Pi_RT^int(B0) + C_h a2.
* "floor" = ||B_exact - (B0 + C_h Pi_ND^int A_exact)||_L2/||B|| on M2 = error of the projection of the exact field onto M2 (for integrated dofs this equals the RT^int interpolation error,
  by the commuting diagram), i.e. the best a transfer could hope for. "err/floor" is therefore the cost of the remap in units of the discretisation error.
* Operators (all in `vp::Transfer`, evaluation of the M1 field at the required points by GSLIB `FindPointsGSLIB` + `GetVectorValue`):
  A_pt / B_pt: point values of the tangential / normal component at the dof nodes of M2 (MFEM `Project`); A_int / B_int: the integrated functionals of `vp::ProjectA/B`
  (sub-edge circulations, sub-face fluxes) with nq = p+2 Gauss points per direction (default; see (e)); A_l2 / B_l2: global L2 projection (rhs by quadrature of order 2p+2q+2 at points of M2,
  partial-assembly mass with the *same* rule, Jacobi-CG rel. tol 1e-12; the same-mesh transfer is the identity to 1e-12);
  B_l2c: B_l2 followed by the M-orthogonal projection onto ker D_h (see "Cost").
* Periodic point location: query points are wrapped into [0,1)^3; on deformed meshes the unit-cell copy of the mesh sticks out of [0,1)^3 so some wrapped points are in no element
  and are searched again at the 26 periodic images (up to 7e5 points per run needed a shift). Not-found points: **0 in all 258 runs**. No inward perturbation of 1e-12 was needed.
  Note: GSLIB's default tolerance for "point on the border" is |dx| <= 1e-4 (dist^2 <= 1e-8); it is tightened to 1e-9 in `RemoteEvaluator`, otherwise points slightly outside an element are
  silently clamped to it instead of being searched at the image. GSLIB's general (non-H1) interpolation path handles ND and RT (it calls `GetVectorValue`); no changes were needed.

## Findings

1. **A-route divergence is at roundoff for every operator**: max scaled ||D_h b2||/(||B||/h) = 1.3e-14, max |D_h b2| entry = 7e-12 (p = 4), net slice flux of the transferred field equals B0 to <= 4.1e-15 --
   exactly as for a field that was never transferred. This is by construction (b2 = B0 + C_h a2) and it holds for A_pt too, whose B error is the worst. 
2. **B-route divergence is not at roundoff unless it is cleaned.** Scaled L2 divergence (p=2, N=8, eps 0.3): B_pt 2.6e-2, B_l2 9.0e-3, B_int 6.1e-3 (nq = 4) -- of the same size as the B error itself
   (2.6e-2 .. 3.5e-2), decreasing roughly like h^(p-1): p = 1: 0.10, 0.045, 0.056; p = 4: 3.7e-4, 1.8e-4, 5.6e-5. **B_l2c** brings it to 6e-15 (max over all runs 6e-14 scaled, 1.9e-13 in
   max|D_h b| relative to ||B||/h) with a global solve. The B-route also does **not** preserve the net slice flux exactly: errors 3e-7 .. 0.4 (B_pt), 5e-7 .. 0.19 (B_int), 1.3e-6 .. 4e-2 (B_l2),
   and after cleaning 9e-12 .. 1.3e-2 (the clean projection fixes div but not the mean; the mean error is largest at p=1). The A-route keeps the mean to 1e-15.
3. **T_B_int divergence vs quadrature (Table 3.6, (e)).** "Small but not roundoff" is confirmed, but "small" is relative: the scaled L2 divergence of B_int decreases only algebraically with the
   number of Gauss points per direction nq, roughly like nq^-1 (p=2: 1.2e-2, 9.2e-3, 6.1e-3, 3.4e-3, 1.9e-3 for nq = 2,3,4,6,8; p=3: 3.2e-3, 1.9e-3, 1.4e-3, 8.4e-4, 6.2e-4; p=1 is noisier),
   while the cost grows like nq^2 (3.5M points at p=3, nq=8). The cause is the kink of the source B1_h's normal component at M1 faces cutting M2 sub-faces; no practical nq approaches 1e-14.
   The B error is almost independent of nq for B_int (p=2: 2.80e-2 .. 2.81e-2), whereas A_int's B error does depend on nq (p=3: 4.2e-3 at nq=2, 3.0e-3 (3), 2.6e-3 (4), 2.4e-3 at the default nq=5, 2.3e-3 (6), 2.15e-3 at nq=8; floor 1.5e-3), so the default nq = p+2 leaves a
   few tens of percent on the table for A_int at p>=3.
4. **B error: the A-route is not better, and is somewhat worse than the B-route, increasingly so with p and eps.** Over all (a) configurations the ratio (A-route error)/(B-route error) of the
   same family has median 1.10 (pt), 1.03 (int), 1.12 (l2), ranges 0.87-4.4 (pt), 0.98-2.0 (int), 0.98-4.6 (l2); in (b) (f1 -> f2 mod2) medians 0.94, 1.01, 1.15; in (d) 1.25, 1.05, 1.29 (maxima 3.6, 1.9, 2.8).
   In units of the floor (eps 0.3, N=8): p=2: A_pt 1.53, A_int 1.28, A_l2 1.33, B_pt 1.47, B_int 1.24, B_l2 1.16, B_l2c 1.15; p=3: 2.88, 1.61, 2.05, 1.55, 1.27, 1.21, 1.20;
   p=4: 8.9, 2.95, 5.05, 2.17, 1.68, 1.42, 1.42. So the B-route with L2 projection is within 15-40% of the floor, A_int within 1.3-3x, A_pt up to 9x (p=4) and 8x (p=3, eps 0.5).
   At p = 1 everything is within 1-4% of the (large, 22%) floor and the differences are irrelevant. The difference grows with eps (p=3: A_int 1.14, 1.61, 3.08 x floor for eps = 0.1, 0.3, 0.5; B_l2c 1.09, 1.20, 1.40).
   Plausible cause (not separately verified): the curl of the interpolation error of A is controlled by the interpolation error of the *derivative* of E1 A1_h, which has kinks at M1 faces, whereas B is interpolated
   directly; B_int and B_l2 interpolate/ project a piecewise polynomial B1_h that is already in the target space type.
   The best A-route operator is A_int, then A_l2, then A_pt, the same order as in the gauge test. All operators lose accuracy relative to the floor when the target is finer than the source (d): err(M2)/err(M1) = 0.75-0.8 (B-route and A, p=1) and 0.89-1.04 (B-route, p=2; A_int 0.99-1.01, A_pt 1.26-1.28, A_l2 1.1-1.2), but 1.35-1.6 (A_int, p=3), 1.8-2.2 (A_l2), 3.1-3.6 (A_pt) versus 0.76-1.04 for the B-route at p = 3 (the
   field cannot be recovered beyond what M1 stored, so about 1 is the best possible) -- the A-route *amplifies* the source error when refining at p = 3.
5. **Energy.** The remap changes the energy by an amount comparable with the discretisation error: relative to the energy of the projection of the exact field on M2, (E2 - Eref)/Eref (eps 0.3, N = 8, Table 3.1):
   p=2: +1.4e-4 (A_pt), -4.2e-4 (A_int), +3.3e-4 (A_l2), +6.3e-3 (B_pt), -4.8e-4, -3.8e-4, -3.8e-4 (B_int, B_l2, B_l2c); p=3: 1.3e-5, -3.1e-6, 7.4e-6, -3.7e-5, -6.8e-6, -4.1e-6, -4.1e-6; p=4: 6.7e-7, -2.0e-7, 1.2e-7, -9.3e-6, -8.8e-8, 5.4e-7, 5.4e-7.
   No operator is systematically energy-conserving; B_pt has the largest energy error at p>=2 (0.6% at p=2) and the A-route energy errors are of the same size as the B_int/B_l2 ones (A_pt p=3: 1.3e-5 versus 4e-6 for B_l2). At p = 1 (where the field is under-resolved, floor 22%) the
   energy is *lost* (-6.5% w.r.t. Eref) by A_pt, A_int, B_int, B_l2, B_l2c, while A_l2 (-0.35%) and B_pt (+0.5%) keep it. The L2 projections of B (B_l2, B_l2c) are contractions in the energy norm, so they can only lose energy up to the quadrature error of the kinked source (at p = 4 E2 exceeds E1 by 5e-7 relative for that reason; at p<=3 E2 < E1).
   Helicity (b, mod2, exact H = 0, scaled by ||a||·||B||): |H2 - H1|/(||a||||B||) (N=8) = 7e-5 / 1.4e-3 / 1.8e-3 (p=1: A_pt / A_int / A_l2), 4.2e-4 / 1.8e-4 / 4.2e-5 (p=2), 3.3e-6 / 1.7e-5 / 9.0e-6 (p=3), versus a discrete helicity on M1 itself of 2e-5, 7e-6, 1e-7 (p = 1,2,3); i.e. the remap perturbs the helicity at 1-2 orders above the discretisation level of H on M1 for p = 1,3, with no operator consistently best.
6. **Gauge pollution (c, Table 3.4).** A1 = Pi_ND^int(A + grad chi) is transferred with a smooth (k = 2 pi) or a rough (k = 6 pi) gauge g = 0.5 (so that the gradient part carries ~80% of ||A||_L2 on M1: columns "grad frac"),
   B unchanged. ||b2(A+grad chi) - b2(A)||/||b2||: smooth gauge: p=1: A_pt 0.35, A_int 0.14, A_l2 0.17; p=2: 0.127, 0.041, 0.051; p=3: 0.023, 0.0079, 0.013 (the pollution falls by roughly 3-5x per order, i.e. like the interpolation error of
   grad chi_h along the new edges; the floor for the B error is 1.5e-3 at p=3, so 0.8% gauge pollution is 5x larger). Rough gauge: p=1: 0.93, 0.54, 0.33; p=2: 0.95, 0.75, 0.80; p=3: 0.66, 0.15, 0.50 -- i.e. order one
   (a rough gauge function on N=8 cells is not resolved by M1, so the discrete gradient is a poor representation of it and nothing commutes). A_int is the best at p>=2, A_pt the worst everywhere except p=1/rough where A_l2 wins;
   **none of the transfers is gauge invariant**. The result with the gauge of the main run is the "B err (with gauge)" column. The integrated functionals commute with the gradient only up to the quadrature error of the kink of grad chi_h along an
   M2 edge crossing M1 faces, which is the same mechanism as the divergence error of B_int in finding 3. (Check: ||b1(gauge) - b1(no gauge)|| on M1 itself is <= 1.5e-8, i.e. quadrature level, so the pollution comes entirely from the transfer.)
7. **Cost (Tables 3.1, 3.7; contended machine).** Apply time (without point location) at p=3, N=8 (55k RT dofs): B_pt 0.4 s, A_pt 0.7, B_l2 2.8, A_l2 3.7, A_int 3.4, B_l2c 10.4, B_int 10.8; setup (collect + FindPoints) is 1.7x-2.3x the apply time for the
   local ones and 1.3x for the L2 variants. The number of evaluation points drives everything: A_pt 74k, B_pt 55k, A_l2 = B_l2 = 262k, A_int 369k, B_int 1.4M (nq^2 points per face dof against nq per edge dof). One-shot totals (setup + apply): p=2: A_pt 0.5 s, B_pt 0.3, A_int 1.8, A_l2 2.1, B_l2 1.7,
   B_l2c 3.9, B_int 4.0; p=4: A_pt 4.8, B_pt 3.5, B_l2 15.8, A_l2 16.6, B_l2c 28.4, A_int 28.8, B_int 137.6.
   Iterations of the Jacobi-CG mass solve of the L2 projections: RT 5-21, ND 11-35 (the larger values for the deformed targets of (b) and for N=16).
   **B_l2c** is a global solve: AMG-preconditioned CG on S = D M^-1 D^T (matrix-free, M^-1 by an inner Jacobi-CG to 1e-14), 23-53 outer iterations and 96-974 inner iterations in total (N=8; 28-54 and 232-1375 at N=16; the largest numbers are for the deformed target of (b)), i.e. 1.8-2.3x the one-shot cost of B_l2, and 3.6 / 2.1 / 1.5 / 1.0 x that of A_int for p = 1,2,3,4 (3x apply-only at p=3: 10.4 s vs 3.4 s); cleaning changes the B error by < 1% (it lowers it slightly) and the energy (relative to Eref) by < 0.05 percentage points. Implementation note: the *plain* CG of the specification (Jacobi, no preconditioner) needed 86-355 outer iterations at p=2-4
   (about 10x more); with the AMG preconditioner the preconditioned-residual criterion only reduces ||D b||_2 to 1e-9 .. 1e-11, so the implementation repeats up to a few defect-correction rounds until ||D b'||_2 <= 1e-12 ||D b||_2 (a few rounds). The achieved
   scaled divergence is 1e-15 .. 6e-14.
8. **Honest summary of the one-remap comparison.** On smooth fields (abc/mod2 and mean field) with these meshes, going through A buys exactly one thing: b2 is divergence-free (and carries the exact mean flux) for any transfer operator, even the cheapest (A_pt),
   without a global solve. It does not buy accuracy: at equal point-location cost the B-route is 0-4x more accurate; the cheap pointwise variants (A_pt especially) lose a lot of accuracy at high p and for distorted source meshes; A is gauge dependent and every A transfer pollutes B with the gauge part
   (0.8% .. 35% for a O(1) smooth gauge, O(1) for a rough one). If divergence-free output is required of the B-route one pays a global solve (B_l2c) of 1-3.6x the one-shot cost of A_int (3x at p=3, apply only), but it is more accurate than every A-route operator at p >= 2 (1.15-1.4 x floor versus 1.3-3 x for A_int).
   The only comparison where the A-route is competitive: A_int vs B_int (error ratio median 1.03; the energy error is lower for A_int), at a 3-5x lower cost than B_int (fewer evaluation points) and with exact div.
   Limitations: one deformation family, two field families, uniform/low-contrast meshes, periodic torus, q = p, 1 rank; the p>=3 A-route/B-route accuracy gap depends on the default nq of A_int (finding 3), so A_int at nq = 8 would be closer to the B-route.

## Verification (`bin/transfer_test`)

For all seven operators, a transfer from a mesh onto itself (uniform N = 4, p = 2 and deformed eps 0.3 with q = 2, and p = 3, q = 3 with 2 ranks) reproduces the dofs: max |dst - src|/max|src| <= 1e-11 (L2 variants 5e-12, others 1e-14), with
points found at periodic images (up to 8064 shifted points) and zero not found. Also checked: 1-rank and 2-rank results coincide in that test.

## Tables

### Table 3.1 (a): M1 = f1 eps -> uniform, N=8, abc + B0

| eps | p | op | B err | floor | err/floor | div (scaled L2) | max&#124;D b&#124; abs | dE/E1 | (E2-Eref)/Eref | flux err | t one-shot [s] | iters (outer/inner) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0.1 | 1 | A_pt | 2.25e-01 | 2.24e-01 | 1.00 | 1.64e-16 | 4.97e-14 | -2.51e-02 | -2.51e-02 | 5.55e-16 | 0.07 | - |
| 0.1 | 1 | A_int | 2.25e-01 | 2.24e-01 | 1.00 | 1.49e-16 | 4.26e-14 | -2.50e-02 | -2.50e-02 | 7.22e-16 | 0.18 | - |
| 0.1 | 1 | A_l2 | 2.24e-01 | 2.24e-01 | 1.00 | 1.66e-16 | 4.97e-14 | +3.52e-04 | +3.73e-04 | 1.11e-15 | 0.34 | 15/0 |
| 0.1 | 1 | B_pt | 2.28e-01 | 2.24e-01 | 1.02 | 3.34e-02 | 5.60e+00 | +1.14e-04 | +1.35e-04 | 1.14e-02 | 0.03 | - |
| 0.1 | 1 | B_int | 2.27e-01 | 2.24e-01 | 1.01 | 3.56e-02 | 5.81e+00 | -5.18e-03 | -5.16e-03 | 1.21e-02 | 0.24 | - |
| 0.1 | 1 | B_l2 | 2.26e-01 | 2.24e-01 | 1.01 | 4.35e-02 | 1.14e+01 | -2.56e-02 | -2.56e-02 | 2.33e-02 | 0.29 | 5/0 |
| 0.1 | 1 | B_l2c | 2.26e-01 | 2.24e-01 | 1.01 | 5.13e-14 | 1.06e-11 | -2.60e-02 | -2.59e-02 | 8.03e-07 | 0.61 | 23/96 |
| 0.1 | 2 | A_pt | 2.45e-02 | 2.27e-02 | 1.08 | 1.24e-15 | 4.97e-13 | -2.59e-04 | -2.53e-04 | 7.22e-16 | 0.74 | - |
| 0.1 | 2 | A_int | 2.45e-02 | 2.27e-02 | 1.08 | 1.23e-15 | 4.21e-13 | -2.62e-04 | -2.56e-04 | 6.11e-16 | 2.85 | - |
| 0.1 | 2 | A_l2 | 2.46e-02 | 2.27e-02 | 1.08 | 1.24e-15 | 3.85e-13 | +6.65e-05 | +7.23e-05 | 9.44e-16 | 3.22 | 16/0 |
| 0.1 | 2 | B_pt | 2.55e-02 | 2.27e-02 | 1.12 | 3.37e-03 | 1.07e+00 | +2.43e-03 | +2.44e-03 | 6.10e-04 | 0.46 | - |
| 0.1 | 2 | B_int | 2.44e-02 | 2.27e-02 | 1.07 | 2.98e-03 | 1.34e+00 | -4.13e-04 | -4.08e-04 | 6.11e-04 | 6.74 | - |
| 0.1 | 2 | B_l2 | 2.41e-02 | 2.27e-02 | 1.06 | 4.43e-03 | 1.64e+00 | -2.65e-04 | -2.59e-04 | 7.89e-04 | 2.52 | 6/0 |
| 0.1 | 2 | B_l2c | 2.41e-02 | 2.27e-02 | 1.06 | 2.13e-15 | 6.11e-13 | -2.65e-04 | -2.59e-04 | 5.99e-10 | 5.68 | 29/180 |
| 0.1 | 3 | A_pt | 2.08e-03 | 1.51e-03 | 1.37 | 3.18e-15 | 1.54e-12 | +3.21e-05 | +3.23e-05 | 1.30e-15 | 2.07 | - |
| 0.1 | 3 | A_int | 1.72e-03 | 1.51e-03 | 1.14 | 3.17e-15 | 1.51e-12 | +3.00e-05 | +3.02e-05 | 8.88e-16 | 10.03 | - |
| 0.1 | 3 | A_l2 | 1.75e-03 | 1.51e-03 | 1.16 | 3.19e-15 | 1.48e-12 | -4.68e-07 | -3.30e-07 | 1.11e-15 | 7.68 | 14/0 |
| 0.1 | 3 | B_pt | 1.99e-03 | 1.51e-03 | 1.32 | 2.94e-03 | 1.93e+00 | +1.30e-04 | +1.30e-04 | 1.52e-04 | 1.41 | - |
| 0.1 | 3 | B_int | 1.68e-03 | 1.51e-03 | 1.11 | 1.01e-03 | 8.62e-01 | +3.00e-05 | +3.02e-05 | 1.52e-04 | 36.96 | - |
| 0.1 | 3 | B_l2 | 1.65e-03 | 1.51e-03 | 1.09 | 1.59e-03 | 6.97e-01 | +3.59e-06 | +3.73e-06 | 1.66e-04 | 6.07 | 6/0 |
| 0.1 | 3 | B_l2c | 1.64e-03 | 1.51e-03 | 1.09 | 1.19e-15 | 4.19e-13 | +3.56e-06 | +3.70e-06 | 9.66e-11 | 14.61 | 39/240 |
| 0.1 | 4 | A_pt | 1.72e-04 | 7.49e-05 | 2.30 | 6.65e-15 | 4.41e-12 | +1.99e-07 | +2.00e-07 | 3.44e-15 | 7.02 | - |
| 0.1 | 4 | A_int | 1.02e-04 | 7.49e-05 | 1.36 | 6.64e-15 | 4.01e-12 | +2.35e-07 | +2.37e-07 | 4.11e-15 | 42.79 | - |
| 0.1 | 4 | A_l2 | 1.31e-04 | 7.49e-05 | 1.74 | 6.62e-15 | 4.21e-12 | +2.29e-08 | +2.45e-08 | 3.61e-15 | 23.39 | 11/0 |
| 0.1 | 4 | B_pt | 1.06e-04 | 7.49e-05 | 1.41 | 1.07e-04 | 1.72e-01 | -8.14e-06 | -8.14e-06 | 2.91e-06 | 5.46 | - |
| 0.1 | 4 | B_int | 9.28e-05 | 7.49e-05 | 1.24 | 2.60e-05 | 2.82e-02 | +2.61e-07 | +2.63e-07 | 7.67e-07 | 192.80 | - |
| 0.1 | 4 | B_l2 | 8.50e-05 | 7.49e-05 | 1.14 | 5.82e-05 | 6.43e-02 | +6.27e-07 | +6.28e-07 | 2.76e-06 | 23.15 | 6/0 |
| 0.1 | 4 | B_l2c | 8.49e-05 | 7.49e-05 | 1.13 | 9.09e-16 | 5.13e-13 | +6.27e-07 | +6.28e-07 | 8.89e-12 | 42.56 | 53/324 |
| 0.3 | 1 | A_pt | 2.33e-01 | 2.24e-01 | 1.04 | 1.60e-16 | 4.26e-14 | -6.60e-02 | -6.56e-02 | 7.77e-16 | 0.06 | - |
| 0.3 | 1 | A_int | 2.30e-01 | 2.24e-01 | 1.03 | 1.62e-16 | 4.97e-14 | -6.56e-02 | -6.52e-02 | 2.78e-16 | 0.18 | - |
| 0.3 | 1 | A_l2 | 2.26e-01 | 2.24e-01 | 1.01 | 1.65e-16 | 4.26e-14 | -3.90e-03 | -3.51e-03 | 6.38e-16 | 0.38 | 15/0 |
| 0.3 | 1 | B_pt | 2.59e-01 | 2.24e-01 | 1.16 | 1.02e-01 | 1.79e+01 | +4.43e-03 | +4.83e-03 | 1.02e-01 | 0.02 | - |
| 0.3 | 1 | B_int | 2.31e-01 | 2.24e-01 | 1.03 | 5.57e-02 | 1.21e+01 | -6.74e-02 | -6.70e-02 | 1.09e-01 | 0.28 | - |
| 0.3 | 1 | B_l2 | 2.28e-01 | 2.24e-01 | 1.02 | 4.50e-02 | 1.05e+01 | -6.47e-02 | -6.43e-02 | 2.56e-02 | 0.32 | 5/0 |
| 0.3 | 1 | B_l2c | 2.28e-01 | 2.24e-01 | 1.02 | 1.71e-14 | 4.07e-12 | -6.50e-02 | -6.47e-02 | 5.55e-06 | 0.65 | 24/100 |
| 0.3 | 2 | A_pt | 3.47e-02 | 2.27e-02 | 1.53 | 1.24e-15 | 3.84e-13 | +7.46e-05 | +1.38e-04 | 8.33e-16 | 0.48 | - |
| 0.3 | 2 | A_int | 2.91e-02 | 2.27e-02 | 1.28 | 1.23e-15 | 4.69e-13 | -4.84e-04 | -4.21e-04 | 1.05e-15 | 1.81 | - |
| 0.3 | 2 | A_l2 | 3.02e-02 | 2.27e-02 | 1.33 | 1.21e-15 | 4.26e-13 | +2.68e-04 | +3.31e-04 | 8.33e-16 | 2.06 | 18/0 |
| 0.3 | 2 | B_pt | 3.35e-02 | 2.27e-02 | 1.47 | 2.60e-02 | 1.24e+01 | +6.27e-03 | +6.33e-03 | 5.92e-03 | 0.32 | - |
| 0.3 | 2 | B_int | 2.81e-02 | 2.27e-02 | 1.24 | 6.06e-03 | 1.99e+00 | -5.40e-04 | -4.77e-04 | 3.90e-04 | 4.00 | - |
| 0.3 | 2 | B_l2 | 2.63e-02 | 2.27e-02 | 1.16 | 8.99e-03 | 3.39e+00 | -4.39e-04 | -3.76e-04 | 8.56e-04 | 1.66 | 6/0 |
| 0.3 | 2 | B_l2c | 2.62e-02 | 2.27e-02 | 1.15 | 6.23e-15 | 1.97e-12 | -4.41e-04 | -3.78e-04 | 9.36e-09 | 3.85 | 29/180 |
| 0.3 | 3 | A_pt | 4.35e-03 | 1.51e-03 | 2.88 | 3.19e-15 | 1.40e-12 | +1.15e-05 | +1.29e-05 | 1.61e-15 | 1.82 | - |
| 0.3 | 3 | A_int | 2.43e-03 | 1.51e-03 | 1.61 | 3.20e-15 | 1.73e-12 | -4.57e-06 | -3.13e-06 | 1.28e-15 | 9.78 | - |
| 0.3 | 3 | A_l2 | 3.10e-03 | 1.51e-03 | 2.05 | 3.22e-15 | 1.39e-12 | +5.99e-06 | +7.43e-06 | 1.55e-15 | 7.97 | 15/0 |
| 0.3 | 3 | B_pt | 2.35e-03 | 1.51e-03 | 1.55 | 5.00e-03 | 3.21e+00 | -3.84e-05 | -3.70e-05 | 1.17e-03 | 1.29 | - |
| 0.3 | 3 | B_int | 1.92e-03 | 1.51e-03 | 1.27 | 1.05e-03 | 7.32e-01 | -8.18e-06 | -6.75e-06 | 5.04e-05 | 35.32 | - |
| 0.3 | 3 | B_l2 | 1.82e-03 | 1.51e-03 | 1.21 | 1.97e-03 | 9.51e-01 | -5.54e-06 | -4.10e-06 | 1.89e-05 | 6.40 | 6/0 |
| 0.3 | 3 | B_l2c | 1.81e-03 | 1.51e-03 | 1.20 | 1.58e-15 | 5.48e-13 | -5.57e-06 | -4.14e-06 | 9.79e-09 | 14.56 | 39/240 |
| 0.3 | 4 | A_pt | 6.65e-04 | 7.49e-05 | 8.88 | 6.63e-15 | 4.52e-12 | +6.58e-07 | +6.74e-07 | 3.77e-15 | 4.83 | - |
| 0.3 | 4 | A_int | 2.21e-04 | 7.49e-05 | 2.95 | 6.57e-15 | 4.35e-12 | -2.20e-07 | -2.04e-07 | 3.50e-15 | 28.78 | - |
| 0.3 | 4 | A_l2 | 3.78e-04 | 7.49e-05 | 5.05 | 6.70e-15 | 5.79e-12 | +1.00e-07 | +1.16e-07 | 3.61e-15 | 16.55 | 13/0 |
| 0.3 | 4 | B_pt | 1.62e-04 | 7.49e-05 | 2.17 | 3.70e-04 | 2.89e-01 | -9.30e-06 | -9.28e-06 | 2.80e-07 | 3.46 | - |
| 0.3 | 4 | B_int | 1.26e-04 | 7.49e-05 | 1.68 | 5.61e-05 | 3.37e-02 | -1.05e-07 | -8.82e-08 | 5.66e-07 | 137.60 | - |
| 0.3 | 4 | B_l2 | 1.07e-04 | 7.49e-05 | 1.42 | 1.77e-04 | 1.25e-01 | +5.19e-07 | +5.35e-07 | 1.33e-06 | 15.82 | 6/0 |
| 0.3 | 4 | B_l2c | 1.06e-04 | 7.49e-05 | 1.42 | 9.11e-16 | 5.12e-13 | +5.19e-07 | +5.35e-07 | 8.15e-10 | 28.36 | 53/324 |
| 0.5 | 1 | A_pt | 2.64e-01 | 2.24e-01 | 1.18 | 1.62e-16 | 4.26e-14 | -8.45e-02 | -8.23e-02 | 1.67e-16 | 0.04 | - |
| 0.5 | 1 | A_int | 2.39e-01 | 2.24e-01 | 1.07 | 1.60e-16 | 4.26e-14 | -9.73e-02 | -9.52e-02 | 5.00e-16 | 0.11 | - |
| 0.5 | 1 | A_l2 | 2.29e-01 | 2.24e-01 | 1.02 | 1.65e-16 | 4.26e-14 | -1.20e-02 | -9.71e-03 | 8.33e-16 | 0.22 | 15/0 |
| 0.5 | 1 | B_pt | 3.05e-01 | 2.24e-01 | 1.36 | 3.12e-01 | 7.31e+01 | -1.34e-02 | -1.11e-02 | 2.82e-01 | 0.02 | - |
| 0.5 | 1 | B_int | 2.40e-01 | 2.24e-01 | 1.07 | 8.48e-02 | 2.02e+01 | -9.95e-02 | -9.74e-02 | 7.59e-02 | 0.18 | - |
| 0.5 | 1 | B_l2 | 2.33e-01 | 2.24e-01 | 1.04 | 5.39e-02 | 1.18e+01 | -1.00e-01 | -9.80e-02 | 4.51e-03 | 0.19 | 5/0 |
| 0.5 | 1 | B_l2c | 2.32e-01 | 2.24e-01 | 1.04 | 2.13e-14 | 3.84e-12 | -1.00e-01 | -9.83e-02 | 3.15e-06 | 0.44 | 24/100 |
| 0.5 | 2 | A_pt | 4.51e-02 | 2.27e-02 | 1.98 | 1.22e-15 | 3.98e-13 | +2.04e-04 | +4.39e-04 | 6.66e-16 | 0.35 | - |
| 0.5 | 2 | A_int | 3.05e-02 | 2.27e-02 | 1.34 | 1.20e-15 | 4.41e-13 | -7.36e-04 | -5.00e-04 | 5.00e-16 | 1.34 | - |
| 0.5 | 2 | A_l2 | 3.50e-02 | 2.27e-02 | 1.54 | 1.20e-15 | 3.69e-13 | +3.11e-04 | +5.46e-04 | 8.88e-16 | 1.50 | 19/0 |
| 0.5 | 2 | B_pt | 3.39e-02 | 2.27e-02 | 1.49 | 3.64e-02 | 1.21e+01 | -2.56e-03 | -2.32e-03 | 7.99e-03 | 0.23 | - |
| 0.5 | 2 | B_int | 2.86e-02 | 2.27e-02 | 1.26 | 5.45e-03 | 1.71e+00 | -9.46e-04 | -7.11e-04 | 2.79e-04 | 3.15 | - |
| 0.5 | 2 | B_l2 | 2.68e-02 | 2.27e-02 | 1.18 | 1.41e-02 | 5.64e+00 | -9.82e-04 | -7.47e-04 | 6.23e-06 | 1.27 | 6/0 |
| 0.5 | 2 | B_l2c | 2.67e-02 | 2.27e-02 | 1.17 | 7.66e-15 | 2.02e-12 | -9.86e-04 | -7.50e-04 | 1.98e-07 | 2.72 | 29/180 |
| 0.5 | 3 | A_pt | 1.22e-02 | 1.51e-03 | 8.08 | 3.15e-15 | 1.66e-12 | +1.39e-04 | +1.44e-04 | 1.22e-15 | 1.45 | - |
| 0.5 | 3 | A_int | 4.66e-03 | 1.51e-03 | 3.08 | 3.19e-15 | 1.72e-12 | +1.10e-05 | +1.63e-05 | 1.28e-15 | 6.71 | - |
| 0.5 | 3 | A_l2 | 5.99e-03 | 1.51e-03 | 3.97 | 3.22e-15 | 1.45e-12 | +1.98e-05 | +2.51e-05 | 8.88e-16 | 5.45 | 16/0 |
| 0.5 | 3 | B_pt | 3.13e-03 | 1.51e-03 | 2.07 | 8.23e-03 | 3.82e+00 | +1.47e-04 | +1.52e-04 | 4.04e-05 | 1.10 | - |
| 0.5 | 3 | B_int | 2.38e-03 | 1.51e-03 | 1.58 | 1.78e-03 | 1.47e+00 | +9.36e-06 | +1.47e-05 | 2.62e-05 | 24.73 | - |
| 0.5 | 3 | B_l2 | 2.13e-03 | 1.51e-03 | 1.41 | 3.21e-03 | 1.35e+00 | +3.08e-06 | +8.41e-06 | 5.76e-05 | 5.37 | 6/0 |
| 0.5 | 3 | B_l2c | 2.11e-03 | 1.51e-03 | 1.40 | 1.96e-15 | 6.96e-13 | +2.99e-06 | +8.33e-06 | 4.73e-08 | 9.37 | 39/240 |
| 0.5 | 4 | A_pt | 1.18e-03 | 7.49e-05 | 15.73 | 6.61e-15 | 4.27e-12 | +2.39e-06 | +2.46e-06 | 3.50e-15 | 4.92 | - |
| 0.5 | 4 | A_int | 3.94e-04 | 7.49e-05 | 5.26 | 6.62e-15 | 4.36e-12 | +2.37e-07 | +3.01e-07 | 3.39e-15 | 30.07 | - |
| 0.5 | 4 | A_l2 | 7.45e-04 | 7.49e-05 | 9.95 | 6.62e-15 | 4.64e-12 | +6.99e-07 | +7.63e-07 | 3.77e-15 | 25.78 | 14/0 |
| 0.5 | 4 | B_pt | 2.67e-04 | 7.49e-05 | 3.57 | 8.27e-04 | 6.44e-01 | +1.43e-06 | +1.49e-06 | 1.55e-06 | 5.86 | - |
| 0.5 | 4 | B_int | 1.95e-04 | 7.49e-05 | 2.61 | 1.20e-04 | 9.17e-02 | +5.74e-08 | +1.21e-07 | 5.07e-07 | 208.07 | - |
| 0.5 | 4 | B_l2 | 1.61e-04 | 7.49e-05 | 2.15 | 3.87e-04 | 1.94e-01 | -1.05e-06 | -9.91e-07 | 1.47e-06 | 25.16 | 6/0 |
| 0.5 | 4 | B_l2c | 1.59e-04 | 7.49e-05 | 2.12 | 9.44e-16 | 5.40e-13 | -1.06e-06 | -9.91e-07 | 1.47e-09 | 43.74 | 52/318 |

### Table 3.2 (a): N=16 (p<=2), M1 = f1 eps -> uniform N=16

| eps | p | op | B err | floor | err/floor | div (scaled L2) | dE/E1 | t one-shot [s] | iters |
|---|---|---|---|---|---|---|---|---|---|
| 0.1 | 1 | A_pt | 1.13e-01 | 1.13e-01 | 1.01 | 2.60e-16 | -1.30e-02 | 0.42 | - |
| 0.1 | 1 | A_int | 1.14e-01 | 1.13e-01 | 1.01 | 2.59e-16 | -1.28e-02 | 1.45 | - |
| 0.1 | 1 | A_l2 | 1.13e-01 | 1.13e-01 | 1.00 | 2.69e-16 | -3.11e-04 | 3.05 | 26/0 |
| 0.1 | 1 | B_pt | 1.23e-01 | 1.13e-01 | 1.09 | 1.88e-02 | +2.78e-04 | 0.22 | - |
| 0.1 | 1 | B_int | 1.16e-01 | 1.13e-01 | 1.03 | 3.32e-02 | -1.28e-02 | 2.10 | - |
| 0.1 | 1 | B_l2 | 1.14e-01 | 1.13e-01 | 1.01 | 2.04e-02 | -1.31e-02 | 2.49 | 9/0 |
| 0.1 | 1 | B_l2c | 1.14e-01 | 1.13e-01 | 1.01 | 1.66e-14 | -1.32e-02 | 6.94 | 28/232 |
| 0.1 | 2 | A_pt | 6.93e-03 | 5.72e-03 | 1.21 | 2.32e-15 | -1.36e-05 | 3.36 | - |
| 0.1 | 2 | A_int | 6.79e-03 | 5.72e-03 | 1.19 | 2.32e-15 | -2.23e-05 | 14.01 | - |
| 0.1 | 2 | A_l2 | 7.09e-03 | 5.72e-03 | 1.24 | 2.34e-15 | +1.70e-05 | 14.65 | 17/0 |
| 0.1 | 2 | B_pt | 7.98e-03 | 5.72e-03 | 1.39 | 1.74e-03 | +2.11e-03 | 2.34 | - |
| 0.1 | 2 | B_int | 6.77e-03 | 5.72e-03 | 1.18 | 6.19e-04 | -2.14e-05 | 37.05 | - |
| 0.1 | 2 | B_l2 | 6.49e-03 | 5.72e-03 | 1.13 | 9.29e-04 | -2.68e-06 | 14.44 | 10/0 |
| 0.1 | 2 | B_l2c | 6.49e-03 | 5.72e-03 | 1.13 | 9.39e-16 | -2.71e-06 | 46.15 | 32/330 |
| 0.3 | 1 | A_pt | 1.47e-01 | 1.13e-01 | 1.30 | 2.66e-16 | -1.60e-02 | 0.52 | - |
| 0.3 | 1 | A_int | 1.21e-01 | 1.13e-01 | 1.07 | 2.57e-16 | -2.42e-02 | 1.46 | - |
| 0.3 | 1 | A_l2 | 1.17e-01 | 1.13e-01 | 1.04 | 2.69e-16 | -9.85e-05 | 3.43 | 29/0 |
| 0.3 | 1 | B_pt | 1.53e-01 | 1.13e-01 | 1.36 | 1.33e-01 | -2.90e-03 | 0.25 | - |
| 0.3 | 1 | B_int | 1.18e-01 | 1.13e-01 | 1.04 | 5.05e-02 | -2.45e-02 | 2.02 | - |
| 0.3 | 1 | B_l2 | 1.15e-01 | 1.13e-01 | 1.02 | 3.72e-02 | -2.56e-02 | 2.48 | 9/0 |
| 0.3 | 1 | B_l2c | 1.14e-01 | 1.13e-01 | 1.01 | 2.97e-14 | -2.57e-02 | 9.04 | 28/232 |
| 0.3 | 2 | A_pt | 9.48e-03 | 5.72e-03 | 1.66 | 2.33e-15 | -2.56e-05 | 3.59 | - |
| 0.3 | 2 | A_int | 7.34e-03 | 5.72e-03 | 1.28 | 2.32e-15 | -3.42e-05 | 13.77 | - |
| 0.3 | 2 | A_l2 | 8.31e-03 | 5.72e-03 | 1.45 | 2.31e-15 | +2.90e-05 | 16.44 | 19/0 |
| 0.3 | 2 | B_pt | 8.24e-03 | 5.72e-03 | 1.44 | 4.65e-03 | -2.98e-04 | 2.23 | - |
| 0.3 | 2 | B_int | 7.11e-03 | 5.72e-03 | 1.24 | 8.71e-04 | -4.23e-05 | 37.16 | - |
| 0.3 | 2 | B_l2 | 6.56e-03 | 5.72e-03 | 1.15 | 2.23e-03 | -6.98e-05 | 12.33 | 10/0 |
| 0.3 | 2 | B_l2c | 6.55e-03 | 5.72e-03 | 1.14 | 9.93e-16 | -6.99e-05 | 43.56 | 32/330 |
| 0.5 | 1 | A_pt | 1.91e-01 | 1.13e-01 | 1.69 | 2.61e-16 | -1.70e-03 | 0.57 | - |
| 0.5 | 1 | A_int | 1.28e-01 | 1.13e-01 | 1.13 | 2.62e-16 | -2.20e-02 | 1.56 | - |
| 0.5 | 1 | A_l2 | 1.17e-01 | 1.13e-01 | 1.03 | 2.61e-16 | -2.95e-03 | 3.02 | 29/0 |
| 0.5 | 1 | B_pt | 1.52e-01 | 1.13e-01 | 1.35 | 1.53e-01 | -2.54e-03 | 0.22 | - |
| 0.5 | 1 | B_int | 1.19e-01 | 1.13e-01 | 1.05 | 3.59e-02 | -2.52e-02 | 2.07 | - |
| 0.5 | 1 | B_l2 | 1.15e-01 | 1.13e-01 | 1.02 | 2.39e-02 | -2.63e-02 | 2.99 | 9/0 |
| 0.5 | 1 | B_l2c | 1.14e-01 | 1.13e-01 | 1.01 | 1.78e-14 | -2.64e-02 | 8.86 | 28/232 |
| 0.5 | 2 | A_pt | 1.76e-02 | 5.72e-03 | 3.08 | 2.33e-15 | +2.54e-04 | 3.87 | - |
| 0.5 | 2 | A_int | 8.23e-03 | 5.72e-03 | 1.44 | 2.33e-15 | -3.33e-05 | 14.71 | - |
| 0.5 | 2 | A_l2 | 9.78e-03 | 5.72e-03 | 1.71 | 2.34e-15 | +3.93e-05 | 17.13 | 19/0 |
| 0.5 | 2 | B_pt | 9.17e-03 | 5.72e-03 | 1.60 | 8.17e-03 | +5.35e-04 | 2.56 | - |
| 0.5 | 2 | B_int | 7.55e-03 | 5.72e-03 | 1.32 | 8.77e-04 | -5.35e-05 | 39.84 | - |
| 0.5 | 2 | B_l2 | 6.91e-03 | 5.72e-03 | 1.21 | 3.12e-03 | -7.00e-05 | 15.25 | 10/0 |
| 0.5 | 2 | B_l2c | 6.90e-03 | 5.72e-03 | 1.21 | 2.36e-15 | -7.01e-05 | 45.51 | 31/320 |

### Table 3.3 (b): f1 0.3 -> f2 0.15 (mod2, no B0)

| N | p | op | B err | floor | err/floor | div (scaled L2) | dE/E1 | (E2-Eref)/Eref | H2 (exact 0) | (H2-H1)/(&#124;a&#124;&#124;B&#124;) | &#124;a2&#124;/&#124;a1&#124; | t [s] |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 8 | 1 | A_pt | 3.65e-01 | 3.19e-01 | 1.14 | 1.40e-16 | -1.81e-01 | -1.70e-01 | -5.98e-04 | -7.37e-05 | 0.927 | 0.06 |
| 8 | 1 | A_int | 3.64e-01 | 3.19e-01 | 1.14 | 1.28e-16 | -1.83e-01 | -1.72e-01 | +8.69e-03 | +1.38e-03 | 0.926 | 0.17 |
| 8 | 1 | A_l2 | 3.25e-01 | 3.19e-01 | 1.02 | 1.53e-16 | -1.49e-02 | -1.91e-03 | -1.17e-02 | -1.81e-03 | 0.997 | 0.42 |
| 8 | 1 | B_pt | 3.78e-01 | 3.19e-01 | 1.18 | 2.05e-01 | -1.02e-01 | -9.02e-02 | n/a | n/a | n/a | 0.03 |
| 8 | 1 | B_int | 3.65e-01 | 3.19e-01 | 1.14 | 6.43e-02 | -1.76e-01 | -1.65e-01 | n/a | n/a | n/a | 0.26 |
| 8 | 1 | B_l2 | 3.42e-01 | 3.19e-01 | 1.07 | 9.87e-02 | -1.02e-01 | -9.04e-02 | n/a | n/a | n/a | 0.35 |
| 8 | 1 | B_l2c | 3.41e-01 | 3.19e-01 | 1.07 | 6.12e-14 | -1.03e-01 | -9.17e-02 | n/a | n/a | n/a | 1.27 |
| 8 | 2 | A_pt | 7.40e-02 | 5.74e-02 | 1.29 | 8.30e-16 | -4.48e-03 | -4.23e-03 | -3.36e-03 | -4.23e-04 | 0.999 | 0.64 |
| 8 | 2 | A_int | 6.94e-02 | 5.74e-02 | 1.21 | 8.25e-16 | -5.54e-03 | -5.29e-03 | +1.48e-03 | +1.76e-04 | 0.998 | 3.05 |
| 8 | 2 | A_l2 | 7.01e-02 | 5.74e-02 | 1.22 | 8.24e-16 | +6.14e-04 | +8.66e-04 | +3.97e-04 | +4.17e-05 | 1.000 | 3.62 |
| 8 | 2 | B_pt | 8.21e-02 | 5.74e-02 | 1.43 | 1.07e-01 | +5.31e-03 | +5.57e-03 | n/a | n/a | n/a | 0.40 |
| 8 | 2 | B_int | 6.93e-02 | 5.74e-02 | 1.21 | 2.79e-02 | -5.31e-03 | -5.06e-03 | n/a | n/a | n/a | 6.68 |
| 8 | 2 | B_l2 | 6.49e-02 | 5.74e-02 | 1.13 | 5.09e-02 | -3.84e-03 | -3.59e-03 | n/a | n/a | n/a | 2.78 |
| 8 | 2 | B_l2c | 6.50e-02 | 5.74e-02 | 1.13 | 2.82e-14 | -3.90e-03 | -3.65e-03 | n/a | n/a | n/a | 14.71 |
| 8 | 3 | A_pt | 1.18e-02 | 8.63e-03 | 1.37 | 2.14e-15 | -1.40e-04 | -1.35e-04 | +2.73e-05 | +3.25e-06 | 1.000 | 2.11 |
| 8 | 3 | A_int | 1.07e-02 | 8.63e-03 | 1.23 | 2.12e-15 | -1.46e-04 | -1.42e-04 | +1.38e-04 | +1.69e-05 | 1.000 | 10.68 |
| 8 | 3 | A_l2 | 1.19e-02 | 8.63e-03 | 1.38 | 2.13e-15 | +1.34e-05 | +1.78e-05 | -7.24e-05 | -9.04e-06 | 1.000 | 8.26 |
| 8 | 3 | B_pt | 1.25e-02 | 8.63e-03 | 1.45 | 2.35e-02 | -3.13e-04 | -3.09e-04 | n/a | n/a | n/a | 1.51 |
| 8 | 3 | B_int | 1.05e-02 | 8.63e-03 | 1.22 | 4.09e-03 | -1.50e-04 | -1.45e-04 | n/a | n/a | n/a | 38.64 |
| 8 | 3 | B_l2 | 9.77e-03 | 8.63e-03 | 1.13 | 1.41e-02 | -1.33e-04 | -1.29e-04 | n/a | n/a | n/a | 6.58 |
| 8 | 3 | B_l2c | 9.75e-03 | 8.63e-03 | 1.13 | 8.72e-15 | -1.35e-04 | -1.31e-04 | n/a | n/a | n/a | 32.00 |
| 16 | 1 | A_pt | 1.95e-01 | 1.59e-01 | 1.23 | 1.95e-16 | -7.54e-02 | -7.23e-02 | +2.12e-03 | +2.98e-04 | 0.970 | 0.43 |
| 16 | 1 | A_int | 1.79e-01 | 1.59e-01 | 1.13 | 1.91e-16 | -8.15e-02 | -7.84e-02 | +6.09e-03 | +8.20e-04 | 0.970 | 1.41 |
| 16 | 1 | A_l2 | 1.62e-01 | 1.59e-01 | 1.02 | 1.98e-16 | -2.92e-03 | +4.85e-04 | -3.17e-03 | -3.97e-04 | 0.999 | 3.50 |
| 16 | 1 | B_pt | 2.19e-01 | 1.59e-01 | 1.38 | 1.85e-01 | -4.10e-02 | -3.77e-02 | n/a | n/a | n/a | 0.23 |
| 16 | 1 | B_int | 1.79e-01 | 1.59e-01 | 1.13 | 5.91e-02 | -7.79e-02 | -7.48e-02 | n/a | n/a | n/a | 2.13 |
| 16 | 1 | B_l2 | 1.68e-01 | 1.59e-01 | 1.06 | 4.88e-02 | -4.52e-02 | -4.20e-02 | n/a | n/a | n/a | 2.98 |
| 16 | 1 | B_l2c | 1.67e-01 | 1.59e-01 | 1.05 | 2.23e-14 | -4.55e-02 | -4.22e-02 | n/a | n/a | n/a | 21.33 |
| 16 | 2 | A_pt | 1.94e-02 | 1.45e-02 | 1.34 | 1.53e-15 | -3.95e-04 | -3.77e-04 | +3.21e-04 | +3.91e-05 | 1.000 | 3.56 |
| 16 | 2 | A_int | 1.80e-02 | 1.45e-02 | 1.24 | 1.54e-15 | -4.21e-04 | -4.02e-04 | +4.13e-05 | +4.64e-06 | 1.000 | 14.35 |
| 16 | 2 | A_l2 | 2.03e-02 | 1.45e-02 | 1.40 | 1.53e-15 | +1.21e-04 | +1.40e-04 | -4.65e-05 | -6.19e-06 | 1.000 | 15.93 |
| 16 | 2 | B_pt | 2.08e-02 | 1.45e-02 | 1.43 | 1.61e-02 | -8.32e-04 | -8.14e-04 | n/a | n/a | n/a | 2.30 |
| 16 | 2 | B_int | 1.79e-02 | 1.45e-02 | 1.23 | 3.81e-03 | -4.39e-04 | -4.20e-04 | n/a | n/a | n/a | 39.53 |
| 16 | 2 | B_l2 | 1.64e-02 | 1.45e-02 | 1.13 | 1.45e-02 | -4.18e-04 | -4.00e-04 | n/a | n/a | n/a | 15.49 |
| 16 | 2 | B_l2c | 1.64e-02 | 1.45e-02 | 1.13 | 1.24e-14 | -4.22e-04 | -4.04e-04 | n/a | n/a | n/a | 80.65 |
| 16 | 3 | A_pt | 2.46e-03 | 1.19e-03 | 2.07 | 4.01e-15 | +6.66e-07 | +5.55e-07 | -2.13e-05 | -2.63e-06 | 1.000 | 16.00 |
| 16 | 3 | A_int | 1.60e-03 | 1.19e-03 | 1.35 | 4.01e-15 | -2.58e-06 | -2.70e-06 | -1.30e-06 | -1.64e-07 | 1.000 | 78.94 |
| 16 | 3 | A_l2 | 2.10e-03 | 1.19e-03 | 1.77 | 4.00e-15 | +1.97e-06 | +1.86e-06 | -8.75e-06 | -1.08e-06 | 1.000 | 65.65 |
| 16 | 3 | B_pt | 1.68e-03 | 1.19e-03 | 1.41 | 4.16e-03 | -1.69e-05 | -1.71e-05 | n/a | n/a | n/a | 7.14 |
| 16 | 3 | B_int | 1.43e-03 | 1.19e-03 | 1.20 | 1.13e-03 | -4.36e-06 | -4.47e-06 | n/a | n/a | n/a | 203.78 |
| 16 | 3 | B_l2 | 1.39e-03 | 1.19e-03 | 1.17 | 2.13e-03 | +1.67e-06 | +1.56e-06 | n/a | n/a | n/a | 38.73 |
| 16 | 3 | B_l2c | 1.38e-03 | 1.19e-03 | 1.16 | 1.83e-15 | +1.62e-06 | +1.51e-06 | n/a | n/a | n/a | 217.00 |

### Table 3.4 (c): gauge pollution, g=0.5 (A1 += grad chi on M1), eps 0.3 -> uniform, N=8

| k | p | op | &#124;&#124;b2(g)-b2(0)&#124;&#124;/&#124;&#124;b2&#124;&#124; | &#124;&#124;a2(g)-a2(0)&#124;&#124;/&#124;&#124;a2(0)&#124;&#124; | B err (with gauge) | B err floor | grad frac a2 | grad frac a1 | div |
|---|---|---|---|---|---|---|---|---|---|
| 2pi | 1 | A_pt | 3.45e-01 | 1.152 | 4.17e-01 | 2.24e-01 | 0.749 | 0.786 | 2.04e-16 |
| 2pi | 1 | A_int | 1.36e-01 | 1.086 | 2.64e-01 | 2.24e-01 | 0.735 | 0.786 | 2.07e-16 |
| 2pi | 1 | A_l2 | 1.68e-01 | 1.200 | 2.81e-01 | 2.24e-01 | 0.768 | 0.786 | 2.40e-16 |
| 2pi | 2 | A_pt | 1.27e-01 | 1.361 | 1.33e-01 | 2.27e-02 | 0.806 | 0.806 | 2.05e-15 |
| 2pi | 2 | A_int | 4.05e-02 | 1.358 | 4.99e-02 | 2.27e-02 | 0.805 | 0.806 | 2.03e-15 |
| 2pi | 2 | A_l2 | 5.06e-02 | 1.359 | 5.90e-02 | 2.27e-02 | 0.805 | 0.806 | 1.97e-15 |
| 2pi | 3 | A_pt | 2.33e-02 | 1.360 | 2.37e-02 | 1.51e-03 | 0.806 | 0.806 | 5.21e-15 |
| 2pi | 3 | A_int | 7.90e-03 | 1.360 | 8.27e-03 | 1.51e-03 | 0.806 | 0.806 | 5.26e-15 |
| 2pi | 3 | A_l2 | 1.29e-02 | 1.360 | 1.33e-02 | 1.51e-03 | 0.806 | 0.806 | 5.26e-15 |
| 6pi | 1 | A_pt | 9.25e-01 | 0.891 | 2.31e+00 | 2.24e-01 | 0.569 | 0.859 | 5.12e-16 |
| 6pi | 1 | A_int | 5.35e-01 | 0.599 | 6.40e-01 | 2.24e-01 | 0.500 | 0.859 | 2.25e-16 |
| 6pi | 1 | A_l2 | 3.29e-01 | 1.013 | 4.08e-01 | 2.24e-01 | 0.710 | 0.859 | 4.21e-16 |
| 6pi | 2 | A_pt | 9.50e-01 | 3.720 | 3.03e+00 | 2.27e-02 | 0.961 | 0.968 | 4.48e-15 |
| 6pi | 2 | A_int | 7.48e-01 | 3.567 | 1.13e+00 | 2.27e-02 | 0.962 | 0.968 | 4.04e-15 |
| 6pi | 2 | A_l2 | 7.97e-01 | 3.734 | 1.32e+00 | 2.27e-02 | 0.965 | 0.968 | 4.11e-15 |
| 6pi | 3 | A_pt | 6.58e-01 | 4.057 | 8.74e-01 | 1.51e-03 | 0.971 | 0.971 | 1.31e-14 |
| 6pi | 3 | A_int | 1.54e-01 | 4.057 | 1.56e-01 | 1.51e-03 | 0.971 | 0.971 | 1.32e-14 |
| 6pi | 3 | A_l2 | 4.97e-01 | 4.059 | 5.73e-01 | 1.51e-03 | 0.971 | 0.971 | 1.33e-14 |

### Table 3.5 (d): resolution change M1 N=8 (f1 0.3) -> uniform N2

| N2 | p | op | B err (M2) | floor on M2 | B err on M1 | err2/err1 | div (scaled L2) | dE/E1 | t [s] |
|---|---|---|---|---|---|---|---|---|---|
| 12 | 1 | A_pt | 1.99e-01 | 1.50e-01 | 2.33e-01 | 0.85 | 2.08e-16 | -6.04e-02 | 0.11 |
| 12 | 1 | A_int | 1.86e-01 | 1.50e-01 | 2.33e-01 | 0.80 | 2.13e-16 | -6.49e-02 | 0.35 |
| 12 | 1 | A_l2 | 1.68e-01 | 1.50e-01 | 2.33e-01 | 0.72 | 2.23e-16 | -2.76e-02 | 0.71 |
| 12 | 1 | B_pt | 2.58e-01 | 1.50e-01 | 2.33e-01 | 1.11 | 2.02e-01 | -6.45e-03 | 0.06 |
| 12 | 1 | B_int | 1.85e-01 | 1.50e-01 | 2.33e-01 | 0.79 | 8.29e-02 | -6.49e-02 | 0.46 |
| 12 | 1 | B_l2 | 1.81e-01 | 1.50e-01 | 2.33e-01 | 0.78 | 5.47e-02 | -6.71e-02 | 0.75 |
| 12 | 1 | B_l2c | 1.81e-01 | 1.50e-01 | 2.33e-01 | 0.78 | 2.59e-14 | -6.74e-02 | 1.95 |
| 12 | 2 | A_pt | 3.15e-02 | 1.02e-02 | 2.49e-02 | 1.26 | 1.79e-15 | +2.85e-04 | 1.01 |
| 12 | 2 | A_int | 2.47e-02 | 1.02e-02 | 2.49e-02 | 0.99 | 1.78e-15 | -1.23e-04 | 4.01 |
| 12 | 2 | A_l2 | 3.03e-02 | 1.02e-02 | 2.49e-02 | 1.21 | 1.78e-15 | +2.84e-04 | 4.42 |
| 12 | 2 | B_pt | 2.58e-02 | 1.02e-02 | 2.49e-02 | 1.04 | 1.65e-02 | -8.99e-04 | 0.63 |
| 12 | 2 | B_int | 2.35e-02 | 1.02e-02 | 2.49e-02 | 0.94 | 4.36e-03 | -1.96e-04 | 10.79 |
| 12 | 2 | B_l2 | 2.22e-02 | 1.02e-02 | 2.49e-02 | 0.89 | 7.71e-03 | -3.61e-04 | 4.03 |
| 12 | 2 | B_l2c | 2.22e-02 | 1.02e-02 | 2.49e-02 | 0.89 | 4.65e-15 | -3.62e-04 | 11.03 |
| 12 | 3 | A_pt | 6.20e-03 | 4.50e-04 | 1.98e-03 | 3.13 | 4.63e-15 | +4.00e-05 | 6.67 |
| 12 | 3 | A_int | 2.68e-03 | 4.50e-04 | 1.98e-03 | 1.35 | 4.67e-15 | +8.71e-06 | 33.57 |
| 12 | 3 | A_l2 | 3.59e-03 | 4.50e-04 | 1.98e-03 | 1.81 | 4.67e-15 | +1.01e-05 | 24.69 |
| 12 | 3 | B_pt | 2.07e-03 | 4.50e-04 | 1.98e-03 | 1.04 | 4.70e-03 | +1.77e-05 | 4.76 |
| 12 | 3 | B_int | 1.62e-03 | 4.50e-04 | 1.98e-03 | 0.82 | 1.32e-03 | +3.41e-06 | 126.94 |
| 12 | 3 | B_l2 | 1.52e-03 | 4.50e-04 | 1.98e-03 | 0.76 | 2.10e-03 | -3.72e-06 | 23.99 |
| 12 | 3 | B_l2c | 1.50e-03 | 4.50e-04 | 1.98e-03 | 0.76 | 1.50e-15 | -3.77e-06 | 57.07 |
| 16 | 1 | A_pt | 2.05e-01 | 1.13e-01 | 2.33e-01 | 0.88 | 2.60e-16 | -3.61e-02 | 0.45 |
| 16 | 1 | A_int | 1.80e-01 | 1.13e-01 | 2.33e-01 | 0.77 | 2.53e-16 | -4.70e-02 | 1.43 |
| 16 | 1 | A_l2 | 1.98e-01 | 1.13e-01 | 2.33e-01 | 0.85 | 2.48e-16 | -1.62e-02 | 3.22 |
| 16 | 1 | B_pt | 2.44e-01 | 1.13e-01 | 2.33e-01 | 1.05 | 1.88e-01 | -6.94e-03 | 0.23 |
| 16 | 1 | B_int | 1.79e-01 | 1.13e-01 | 2.33e-01 | 0.77 | 6.38e-02 | -4.76e-02 | 2.04 |
| 16 | 1 | B_l2 | 1.76e-01 | 1.13e-01 | 2.33e-01 | 0.76 | 5.39e-02 | -4.86e-02 | 2.49 |
| 16 | 1 | B_l2c | 1.75e-01 | 1.13e-01 | 2.33e-01 | 0.75 | 1.84e-14 | -4.88e-02 | 7.35 |
| 16 | 2 | A_pt | 3.20e-02 | 5.72e-03 | 2.49e-02 | 1.28 | 2.32e-15 | +3.08e-04 | 3.59 |
| 16 | 2 | A_int | 2.51e-02 | 5.72e-03 | 2.49e-02 | 1.01 | 2.34e-15 | -3.06e-05 | 14.00 |
| 16 | 2 | A_l2 | 2.82e-02 | 5.72e-03 | 2.49e-02 | 1.13 | 2.32e-15 | +1.72e-04 | 16.98 |
| 16 | 2 | B_pt | 2.52e-02 | 5.72e-03 | 2.49e-02 | 1.01 | 1.30e-02 | -2.77e-04 | 2.33 |
| 16 | 2 | B_int | 2.37e-02 | 5.72e-03 | 2.49e-02 | 0.95 | 3.75e-03 | -1.02e-04 | 38.59 |
| 16 | 2 | B_l2 | 2.30e-02 | 5.72e-03 | 2.49e-02 | 0.92 | 5.77e-03 | -1.68e-04 | 15.09 |
| 16 | 2 | B_l2c | 2.29e-02 | 5.72e-03 | 2.49e-02 | 0.92 | 2.50e-15 | -1.69e-04 | 44.99 |
| 16 | 3 | A_pt | 7.16e-03 | 1.90e-04 | 1.98e-03 | 3.61 | 6.22e-15 | +4.63e-05 | 16.04 |
| 16 | 3 | A_int | 3.13e-03 | 1.90e-04 | 1.98e-03 | 1.58 | 6.21e-15 | +6.72e-06 | 78.25 |
| 16 | 3 | A_l2 | 4.40e-03 | 1.90e-04 | 1.98e-03 | 2.22 | 6.23e-15 | +1.65e-05 | 61.31 |
| 16 | 3 | B_pt | 1.99e-03 | 1.90e-04 | 1.98e-03 | 1.00 | 3.91e-03 | +2.13e-06 | 11.26 |
| 16 | 3 | B_int | 1.66e-03 | 1.90e-04 | 1.98e-03 | 0.84 | 1.07e-03 | -2.00e-06 | 289.47 |
| 16 | 3 | B_l2 | 1.59e-03 | 1.90e-04 | 1.98e-03 | 0.80 | 1.91e-03 | -3.19e-06 | 53.51 |
| 16 | 3 | B_l2c | 1.58e-03 | 1.90e-04 | 1.98e-03 | 0.80 | 1.54e-15 | -3.22e-06 | 172.77 |

### Table 3.6 (e): quadrature of the integrated transfers (N=8, eps 0.3 -> uniform, abc+B0)

| p | nq | op | B err | div (scaled L2) | max&#124;D b&#124; abs | dE/E1 | #points | t [s] |
|---|---|---|---|---|---|---|---|---|
| 1 | 2 | A_int | 2.33e-01 | 1.46e-16 | 4.26e-14 | -6.30e-02 | 12288 | 0.13 |
| 1 | 2 | B_int | 2.39e-01 | 1.38e-01 | 4.00e+01 | -5.18e-02 | 12288 | 0.12 |
| 1 | 3 | A_int | 2.30e-01 | 1.62e-16 | 4.97e-14 | -6.56e-02 | 18432 | 0.19 |
| 1 | 3 | B_int | 2.31e-01 | 5.57e-02 | 1.21e+01 | -6.74e-02 | 27648 | 0.27 |
| 1 | 4 | A_int | 2.30e-01 | 1.57e-16 | 3.55e-14 | -6.61e-02 | 24576 | 0.23 |
| 1 | 4 | B_int | 2.31e-01 | 4.22e-02 | 1.31e+01 | -6.17e-02 | 49152 | 0.45 |
| 1 | 6 | A_int | 2.30e-01 | 1.56e-16 | 5.68e-14 | -6.60e-02 | 36864 | 0.36 |
| 1 | 6 | B_int | 2.30e-01 | 3.43e-02 | 9.12e+00 | -6.48e-02 | 110592 | 1.07 |
| 1 | 8 | A_int | 2.30e-01 | 1.58e-16 | 4.09e-14 | -6.60e-02 | 49152 | 0.47 |
| 1 | 8 | B_int | 2.29e-01 | 1.30e-02 | 3.67e+00 | -6.64e-02 | 196608 | 1.88 |
| 2 | 2 | A_int | 2.97e-02 | 1.23e-15 | 5.22e-13 | -5.01e-04 | 55296 | 0.95 |
| 2 | 2 | B_int | 2.80e-02 | 1.24e-02 | 3.54e+00 | -8.17e-04 | 73728 | 1.14 |
| 2 | 3 | A_int | 2.98e-02 | 1.23e-15 | 4.69e-13 | -4.50e-04 | 82944 | 1.15 |
| 2 | 3 | B_int | 2.82e-02 | 9.18e-03 | 3.15e+00 | -5.26e-04 | 165888 | 2.75 |
| 2 | 4 | A_int | 2.91e-02 | 1.23e-15 | 4.69e-13 | -4.84e-04 | 110592 | 1.82 |
| 2 | 4 | B_int | 2.81e-02 | 6.06e-03 | 1.99e+00 | -5.40e-04 | 294912 | 4.87 |
| 2 | 6 | A_int | 2.86e-02 | 1.25e-15 | 3.69e-13 | -4.93e-04 | 165888 | 2.67 |
| 2 | 6 | B_int | 2.81e-02 | 3.38e-03 | 1.30e+00 | -5.41e-04 | 663552 | 9.15 |
| 2 | 8 | A_int | 2.84e-02 | 1.24e-15 | 3.98e-13 | -5.15e-04 | 221184 | 3.80 |
| 2 | 8 | B_int | 2.81e-02 | 1.93e-03 | 5.25e-01 | -5.44e-04 | 1179648 | 19.44 |
| 3 | 2 | A_int | 4.21e-03 | 3.16e-15 | 1.76e-12 | +7.06e-06 | 147456 | 2.36 |
| 3 | 2 | B_int | 2.00e-03 | 3.24e-03 | 2.75e+00 | +7.81e-06 | 221184 | 6.06 |
| 3 | 3 | A_int | 3.02e-03 | 3.22e-15 | 1.51e-12 | +1.64e-06 | 221184 | 6.12 |
| 3 | 3 | B_int | 1.96e-03 | 1.90e-03 | 9.89e-01 | +2.94e-06 | 497664 | 13.51 |
| 3 | 4 | A_int | 2.58e-03 | 3.17e-15 | 1.82e-12 | -5.80e-06 | 294912 | 8.07 |
| 3 | 4 | B_int | 1.94e-03 | 1.41e-03 | 8.26e-01 | -4.04e-06 | 884736 | 11.55 |
| 3 | 6 | A_int | 2.33e-03 | 3.19e-15 | 1.76e-12 | -5.08e-06 | 442368 | 5.96 |
| 3 | 6 | B_int | 1.93e-03 | 8.40e-04 | 5.92e-01 | -6.49e-06 | 1990656 | 26.83 |
| 3 | 8 | A_int | 2.15e-03 | 3.14e-15 | 1.48e-12 | -3.91e-06 | 589824 | 8.33 |
| 3 | 8 | B_int | 1.92e-03 | 6.16e-04 | 2.52e-01 | -8.88e-06 | 3538944 | 47.93 |

### Table 3.7: cost of the global variants (apply = t_apply, setup = point location etc.), N=8 eps 0.3 -> uniform

| p | op | #points | t_setup | (locate) | t_eval | t_dofs | t_solve | t_apply | CG iters | inner iters | max&#124;D b&#124; before cleaning |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | A_pt | 6144 | 0.046 | 0.037 | 0.011 | 0.001 | 0.000 | 0.011 | - | - | - |
| 1 | A_int | 18432 | 0.142 | 0.130 | 0.026 | 0.011 | 0.000 | 0.037 | - | - | - |
| 1 | A_l2 | 32768 | 0.240 | 0.229 | 0.050 | 0.035 | 0.050 | 0.135 | 15 | - | - |
| 1 | B_pt | 3072 | 0.020 | 0.020 | 0.002 | 0.001 | 0.000 | 0.002 | - | - | - |
| 1 | B_int | 27648 | 0.223 | 0.200 | 0.038 | 0.021 | 0.000 | 0.058 | - | - | - |
| 1 | B_l2 | 32768 | 0.228 | 0.216 | 0.049 | 0.025 | 0.022 | 0.095 | 5 | - | - |
| 1 | B_l2c | 32768 | 0.229 | 0.216 | 0.049 | 0.025 | 0.348 | 0.422 | 24 | 100 | 1.05e+01 |
| 2 | A_pt | 27648 | 0.349 | 0.311 | 0.106 | 0.025 | 0.000 | 0.131 | - | - | - |
| 2 | A_int | 110592 | 1.309 | 1.141 | 0.396 | 0.108 | 0.000 | 0.504 | - | - | - |
| 2 | A_l2 | 110592 | 1.209 | 1.153 | 0.396 | 0.276 | 0.182 | 0.854 | 18 | - | - |
| 2 | B_pt | 18432 | 0.237 | 0.203 | 0.061 | 0.023 | 0.000 | 0.083 | - | - | - |
| 2 | B_int | 294912 | 2.984 | 2.533 | 0.778 | 0.241 | 0.000 | 1.019 | - | - | - |
| 2 | B_l2 | 110592 | 0.968 | 0.921 | 0.358 | 0.263 | 0.071 | 0.692 | 6 | - | - |
| 2 | B_l2c | 110592 | 1.200 | 1.133 | 0.371 | 0.253 | 2.025 | 2.650 | 29 | 180 | 3.39e+00 |
| 3 | A_pt | 73728 | 1.140 | 0.982 | 0.506 | 0.179 | 0.000 | 0.684 | - | - | - |
| 3 | A_int | 368640 | 6.372 | 5.293 | 2.664 | 0.740 | 0.000 | 3.403 | - | - | - |
| 3 | A_l2 | 262144 | 4.278 | 4.037 | 1.944 | 1.360 | 0.392 | 3.696 | 15 | - | - |
| 3 | B_pt | 55296 | 0.889 | 0.745 | 0.309 | 0.095 | 0.000 | 0.404 | - | - | - |
| 3 | B_int | 1382400 | 24.474 | 20.809 | 8.534 | 2.314 | 0.000 | 10.848 | - | - | - |
| 3 | B_l2 | 262144 | 3.607 | 3.371 | 1.452 | 1.173 | 0.168 | 2.793 | 6 | - | - |
| 3 | B_l2c | 262144 | 4.180 | 3.939 | 1.666 | 1.233 | 7.479 | 10.379 | 39 | 240 | 9.51e-01 |
| 4 | A_pt | 153600 | 3.125 | 2.606 | 1.329 | 0.376 | 0.000 | 1.705 | - | - | - |
| 4 | A_int | 921600 | 18.759 | 15.850 | 7.922 | 2.102 | 0.000 | 10.024 | - | - | - |
| 4 | A_l2 | 512000 | 9.243 | 8.732 | 3.645 | 3.182 | 0.483 | 7.310 | 13 | - | - |
| 4 | B_pt | 122880 | 2.224 | 1.810 | 0.959 | 0.280 | 0.000 | 1.239 | - | - | - |
| 4 | B_int | 4423680 | 94.504 | 77.419 | 32.487 | 10.609 | 0.000 | 43.097 | - | - | - |
| 4 | B_l2 | 512000 | 8.713 | 8.178 | 3.996 | 2.899 | 0.214 | 7.108 | 6 | - | - |
| 4 | B_l2c | 512000 | 8.152 | 7.482 | 3.936 | 2.944 | 13.326 | 20.206 | 53 | 324 | 1.25e-01 |

