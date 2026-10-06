# E2 summary: the same representation on deformed periodic meshes

Data: `results/exp2_distort.csv` (467 rows, one per run, 1 MPI rank, `scripts/run_exp2.sh`), figures `results/exp2_distort.png`, `results/exp2_q_effect.png`
(`scripts/plot_exp2.py`), tables below by `scripts/make_tables.py exp2`.
Deformation x' = x + eps f(x) (variants f1, f2 (4pi), f3 = f1 + 0.5 f2), nodes in a discontinuous L2 space of order q, topology fixed; slice families recorded on the
undeformed mesh and reused. Fields: `abc` (a,b,c = 1,0.8,0.6) unless stated; A is projected on the DEFORMED mesh, B_exact is evaluated at physical points;
`int` = integrated (exact line integral) dofs, `pt` = pointwise dofs. min detJ is normalised so that the undeformed mesh has 1.
Provenance: the CSV was produced in three stages by `scripts/run_exp2.sh` (p = 1,2 main+extras; p = 3 main; p = 3 extras + p = 4) because of run time; between stages the only code change was the cap
of the integrated-dof quadrature at 8 Gauss points (`min(p+q+5, 8)`, accuracy <~ 1e-12, checked in `dbg` runs), and the p = 4 eps lists were thinned (see the script header). Runs overlapped on a 4-core machine; times are indicative.

## What was established

1. **div_h B_h = 0 is independent of geometry, as claimed.** Over all 400 valid runs that use a divergence-free B0 representation (p = 1..4, q = 1..4, f1/f2/f3, eps up to the point where min detJ drops to 0.09, N = 8, 16), ||D_h b||_L2/(||B||/h) <= 8.4e-15,
   max_i |(D_h b)_i| <= 1.8e-11 absolute (<= 2.5e-13 relative to ||B||/h), max|C_h G_h| <= 7.1e-15, max|D_h C_h| h_min^3 <= 2.9e-14, flux of the curl part through every slice <= 3.8e-15. No trend with eps, q or variant: the divergence is an algebraic property of the dofs, not of the geometry.
2. **B error grows only mildly with distortion; optimal order is kept.** Relative B error (f1, N=8, q=2): p=1 0.224 -> 0.239 (eps=0.3, min detJ 0.90) -> 0.289 (0.6, 0.60) -> 0.413 (0.9, 0.094);
   p=2 2.3e-2 -> 2.5e-2 -> 3.3e-2 -> 5.6e-2; p=3 1.5e-3 -> 1.8e-3 -> 3.0e-3 -> 6.1e-3; p=4 7.5e-5 -> 1.1e-4 (0.3) -> 2.2e-4 (0.6) -> 3.8e-4 (0.8). The degradation factor at fixed eps grows with p (a few percent at eps=0.3 for p=1, +50% for p=4),
   consistent with the loss of polynomial completeness of mapped RT_{p-1}/ND_p spaces on non-affine hexahedra. The h-convergence of B is O(h^p) up to eps=0.5 (min detJ 0.69-0.74): observed 0.97-0.99 (p=1), 1.97-1.99 (p=2) between N=8 and 16 (Table 2.8).
   **A converges worse than on the undeformed mesh:** its L2 order drops from p+1 to about p (p=1: 1.99 -> 1.09-1.4 at eps 0.3-0.5; p=2: 2.99 -> 2.0-2.2) -- A_h is not a better approximation of A than B_h is of B on distorted meshes, but B (the physical quantity) is unaffected in order.
3. **Energy and helicity**: energy error is the superconvergent 2p-type value at eps=0 and moves to the interpolation-error level as eps grows, changing sign near eps ~ 0.3-0.8 (e.g. p=2: -2.6e-4 -> +1.9e-3 at eps=0.9; p=4: 1.4e-9 -> 1.2e-7); for p=2 it stays below 1.9e-3 over the whole scan. The helicity error (abc) stays within a factor 2.4 of its eps=0 value (p=1: x1.3, p=2: 9.8e-3 -> 1.6e-2, p=4: 8.9e-8 -> 2.1e-7 at eps=0.8).
   For `mod2` (helicity exactly 0) the discrete helicity becomes non-zero on distorted meshes (p=2: 1.6e-6 at eps=0.1, 7.9e-4 at 0.6), i.e. the symmetry that made it vanish on the Cartesian mesh is broken (a small, convergent number, not a conservation failure).
4. **Mesh validity.** First invalid eps of the scan (non-positive detJ somewhere on the mesh): f1: q=1 0.8, q=2 1.0 (N=16 only; at N=8 the f1/q=2 series was stopped at eps=0.9 where min detJ = 0.094 < 0.15), q=3 1.0; f2: q=1 0.8, q=2 1.0 (N=16), q=3 1.0; f3: q=1 0.6 (p<=3; 0.8 for the thinned p=4 list), q=2 0.8, q=3 0.8, q=4 0.8 (Table 2.2). The normalised min detJ of f1 (q=2, N=8) is 0.99/0.90/0.60/0.30/0.09 at eps = 0.1/0.3/0.6/0.8/0.9: the prescribed range eps <= 0.3 is only a mild distortion (detJ >= 0.88, within-element min/max detJ ratio >= 0.87); a real stress test needs eps ~ 0.6-0.9 (this is where the plan's "min detJ ~ 0.1" is reached).
5. **Mean field and flux.** Net flux through all six slice families equals B0_i to <= 4.7e-15 for every run, deformed or not, for both B0 projections; the flux of C_h a_h through every slice is <= 3.8e-15. This is the discrete Stokes statement, and it holds exactly even for invalid-ish (min detJ 0.1) meshes.
   **But the plan's statement "Pi_RT(B0) is exact because constants lie in RT_0 ⊂ RT_{p-1}" is false on non-affine meshes** (contravariant Piola map: a constant physical field has reference field adj(J) B0, which is not constant). Table 2.5: ||Pi_RT(B0) - B0||_L2 = 1.4e-2 (p=1, eps=0.1), 9.4e-2 (eps=0.6), 1e-3 .. 8e-3 (p=2), 1.2e-5 .. 4.6e-4 (p=3), and only 2e-15 for p=4 with q=2 (then the mapped constant is in RT_3). With the **pointwise** projection the representation also **loses exact div-freeness**: max|D_h Pi^pt_RT(B0)| = 2.7e-3 (p=1, eps=0.1), 1.1e-1 (eps=0.6), 0.89 (f3, eps=0.6), 1.6e-4 .. 8.4e-3 (p=2); for p >= 3 (q=2) it remains at roundoff in these runs.
   With the **integrated (flux) projection** of B0, D_h b0 = 0 to <= 2.4e-13 for all p (commuting diagram) although b0 is not pointwise equal to B0; the total B error is not affected visibly (e.g. p=1, eps=0.3: 2.39e-1 in both cases). Recommendation: always project B0 with face fluxes (integrated dofs), never with `ProjectCoefficient`, on curved/non-affine meshes.
6. **Pointwise vs integrated A projection** (Table 2.4): for these smooth fields at N=8 the resulting B and A errors agree to <= 0.4% (and the energy errors agree except near their sign change); the pointwise projection's *commutation defect* is large only at dof level (selftest, f3 eps=0.2, q=2, N=4: max|C_h a^pt - b^int| = 1.7e-2 (p=1), 4.9e-3 (p=2), 5.7e-4 (p=3), 4.6e-5 (p=4) versus <= 1e-11 for the integrated one), which matters for gauge transfer (E3), not for the static accuracy.
7. **Cost** does not depend on eps: C_h application 1.0e-3 s (p=3, N=8, 41k dofs), 3.4-3.9e-3 s (p=4, 98k dofs), assembly 0.35 s / 1.5-1.8 s; RSS 62 MB / 116 MB.

## Open points / possible concerns

- Larger q gave *larger* B errors than q=1 for p >= 3 at equal (eps, f) (p=4, eps=0.3: 8.7e-5 (q=1), 1.1e-4 (q=2), 1.5e-4 (q=4); p=3: 1.67e-3 / 1.84e-3 / 1.99e-3). My diagnosis (not separately verified): with q=1 the f1/f2/f3 displacement is sampled only at the vertices, i.e. a smoother deformation with less intra-element variation of J; higher q resolves the true displacement and hence the stronger non-affinity. So q is not a pure "accuracy" knob here; geometry differs between q. q=1 meshes also become invalid earlier (f1 0.8 vs 1.0).
- Short-wavelength deformations hurt more: f2 (4pi, 4 cells per wavelength at N=8) at eps=0.3 gives 2.6e-2 (p=2, q=1) .. 3.4e-2 (q=2) vs 2.5e-2 for f1; for p=3, q=3: 4.9e-3 vs 2.0e-3.
- `max_DC` in the CSV scales with 1/detJ_min (use `max_DC_scaled`); on deformed meshes the L2 values of D_h b include 1/detJ so `div_max` grows with distortion at roundoff level only.
- All results for a single field family (abc, mod2) and a single mesh topology (Cartesian N^3 box with smooth periodic deformation); multi-rank runs were only exercised in `selftest` (1, 2, 4 ranks agree) and a spot check (np=1 vs 4, identical to 1e-12 except roundoff), not in the matrix.

## Tables

### Table 2.1 Distortion scan (f1, N=8, q=2, field abc, no B0, integrated projection)

min detJ is normalised to the undeformed value (quad points + corners); ratio = min over elements of (min detJ)/(max detJ) inside one element.

| p | eps | min detJ | min ratio | rel err B | rel err A | max&#124;D_h b&#124; | rel div L2 | energy rel err | helicity abs err | flux err | max&#124;D_hC_h&#124; scaled |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 0.00 | 1.000 | 1.000 | 2.244e-01 | 5.556e-02 | 4.3e-14 | 1.5e-16 | -5.04e-02 | 1.3e+00 | 8.9e-16 | 0.0e+00 |
| 1 | 0.05 | 0.997 | 0.996 | 2.248e-01 | 5.671e-02 | 3.6e-14 | 1.7e-16 | -5.03e-02 | 1.3e+00 | 2.1e-16 | 0.0e+00 |
| 1 | 0.10 | 0.989 | 0.985 | 2.260e-01 | 6.005e-02 | 5.2e-14 | 1.8e-16 | -5.01e-02 | 1.3e+00 | 4.3e-16 | 0.0e+00 |
| 1 | 0.15 | 0.975 | 0.966 | 2.280e-01 | 6.528e-02 | 4.8e-14 | 1.7e-16 | -4.98e-02 | 1.3e+00 | 7.8e-16 | 0.0e+00 |
| 1 | 0.20 | 0.956 | 0.941 | 2.308e-01 | 7.207e-02 | 5.0e-14 | 1.8e-16 | -4.93e-02 | 1.3e+00 | 7.1e-16 | 0.0e+00 |
| 1 | 0.25 | 0.931 | 0.909 | 2.345e-01 | 8.011e-02 | 6.4e-14 | 1.8e-16 | -4.86e-02 | 1.3e+00 | 3.1e-16 | 0.0e+00 |
| 1 | 0.30 | 0.901 | 0.872 | 2.390e-01 | 8.917e-02 | 5.0e-14 | 1.9e-16 | -4.77e-02 | 1.3e+00 | 4.3e-16 | 0.0e+00 |
| 1 | 0.40 | 0.824 | 0.776 | 2.511e-01 | 1.098e-01 | 5.7e-14 | 1.8e-16 | -4.50e-02 | 1.3e+00 | 3.1e-16 | 0.0e+00 |
| 1 | 0.50 | 0.725 | 0.657 | 2.675e-01 | 1.334e-01 | 4.3e-14 | 1.8e-16 | -4.07e-02 | 1.4e+00 | 7.2e-16 | 0.0e+00 |
| 1 | 0.60 | 0.604 | 0.521 | 2.892e-01 | 1.602e-01 | 6.4e-14 | 2.0e-16 | -3.36e-02 | 1.4e+00 | 5.4e-16 | 0.0e+00 |
| 1 | 0.70 | 0.461 | 0.375 | 3.176e-01 | 1.912e-01 | 5.7e-14 | 2.0e-16 | -2.23e-02 | 1.5e+00 | 9.2e-16 | 0.0e+00 |
| 1 | 0.80 | 0.296 | 0.224 | 3.559e-01 | 2.286e-01 | 8.5e-14 | 2.1e-16 | -3.35e-03 | 1.6e+00 | 1.1e-16 | 0.0e+00 |
| 1 | 0.90 | 0.094 | 0.076 | 4.133e-01 | 2.799e-01 | 1.1e-13 | 2.3e-16 | +3.31e-02 | 1.7e+00 | 1.0e-15 | 0.0e+00 |
| 2 | 0.00 | 1.000 | 1.000 | 2.279e-02 | 2.762e-03 | 4.1e-13 | 1.2e-15 | -2.61e-04 | 9.8e-03 | 3.3e-16 | 8.9e-16 |
| 2 | 0.05 | 0.997 | 0.996 | 2.285e-02 | 3.097e-03 | 3.5e-13 | 1.2e-15 | -2.59e-04 | 9.8e-03 | 8.7e-16 | 8.9e-16 |
| 2 | 0.10 | 0.989 | 0.985 | 2.302e-02 | 3.943e-03 | 4.2e-13 | 1.2e-15 | -2.55e-04 | 9.8e-03 | 6.4e-16 | 1.8e-15 |
| 2 | 0.15 | 0.975 | 0.966 | 2.331e-02 | 5.068e-03 | 4.2e-13 | 1.2e-15 | -2.47e-04 | 9.9e-03 | 4.0e-16 | 8.7e-16 |
| 2 | 0.20 | 0.956 | 0.941 | 2.373e-02 | 6.346e-03 | 4.4e-13 | 1.2e-15 | -2.35e-04 | 1.0e-02 | 8.3e-16 | 1.7e-15 |
| 2 | 0.25 | 0.931 | 0.909 | 2.429e-02 | 7.726e-03 | 4.0e-13 | 1.2e-15 | -2.19e-04 | 1.0e-02 | 5.0e-16 | 1.7e-15 |
| 2 | 0.30 | 0.901 | 0.872 | 2.500e-02 | 9.189e-03 | 4.5e-13 | 1.3e-15 | -1.97e-04 | 1.0e-02 | 6.5e-16 | 1.6e-15 |
| 2 | 0.40 | 0.824 | 0.776 | 2.693e-02 | 1.235e-02 | 4.5e-13 | 1.3e-15 | -1.32e-04 | 1.1e-02 | 4.3e-16 | 7.7e-16 |
| 2 | 0.50 | 0.725 | 0.657 | 2.968e-02 | 1.587e-02 | 5.6e-13 | 1.3e-15 | -2.44e-05 | 1.1e-02 | 3.5e-16 | 1.4e-15 |
| 2 | 0.60 | 0.604 | 0.521 | 3.343e-02 | 1.983e-02 | 5.8e-13 | 1.3e-15 | +1.50e-04 | 1.2e-02 | 5.7e-16 | 1.2e-15 |
| 2 | 0.70 | 0.461 | 0.375 | 3.852e-02 | 2.442e-02 | 6.8e-13 | 1.5e-15 | +4.35e-04 | 1.3e-02 | 2.6e-16 | 1.0e-15 |
| 2 | 0.80 | 0.296 | 0.224 | 4.555e-02 | 2.996e-02 | 1.1e-12 | 1.6e-15 | +9.25e-04 | 1.4e-02 | 4.6e-16 | 1.6e-15 |
| 2 | 0.90 | 0.094 | 0.076 | 5.617e-02 | 3.735e-02 | 1.9e-12 | 1.8e-15 | +1.88e-03 | 1.6e-02 | 4.9e-16 | 2.1e-15 |
| 3 | 0.00 | 1.000 | 1.000 | 1.515e-03 | 1.255e-04 | 1.9e-12 | 3.1e-15 | -7.67e-07 | 3.8e-05 | 6.7e-16 | 3.6e-15 |
| 3 | 0.05 | 0.997 | 0.996 | 1.524e-03 | 1.730e-04 | 1.9e-12 | 3.1e-15 | -7.41e-07 | 3.8e-05 | 5.4e-16 | 3.5e-15 |
| 3 | 0.10 | 0.989 | 0.985 | 1.550e-03 | 2.691e-04 | 1.6e-12 | 3.2e-15 | -6.61e-07 | 3.8e-05 | 8.8e-16 | 3.5e-15 |
| 3 | 0.15 | 0.975 | 0.966 | 1.594e-03 | 3.790e-04 | 1.5e-12 | 3.2e-15 | -5.24e-07 | 3.8e-05 | 6.8e-16 | 3.5e-15 |
| 3 | 0.20 | 0.956 | 0.941 | 1.656e-03 | 4.950e-04 | 1.4e-12 | 3.2e-15 | -3.24e-07 | 3.9e-05 | 1.1e-15 | 3.4e-15 |
| 3 | 0.25 | 0.931 | 0.909 | 1.738e-03 | 6.157e-04 | 1.6e-12 | 3.2e-15 | -5.26e-08 | 3.9e-05 | 6.2e-16 | 3.4e-15 |
| 3 | 0.30 | 0.901 | 0.872 | 1.840e-03 | 7.413e-04 | 2.4e-12 | 3.3e-15 | +3.02e-07 | 3.9e-05 | 4.6e-16 | 3.3e-15 |
| 3 | 0.40 | 0.824 | 0.776 | 2.113e-03 | 1.011e-03 | 1.7e-12 | 3.2e-15 | +1.33e-06 | 3.9e-05 | 1.0e-15 | 3.1e-15 |
| 3 | 0.50 | 0.725 | 0.657 | 2.493e-03 | 1.313e-03 | 2.0e-12 | 3.4e-15 | +2.98e-06 | 4.1e-05 | 7.9e-16 | 2.8e-15 |
| 3 | 0.60 | 0.604 | 0.521 | 3.009e-03 | 1.665e-03 | 2.8e-12 | 3.6e-15 | +5.63e-06 | 4.3e-05 | 6.4e-16 | 2.5e-15 |
| 3 | 0.70 | 0.461 | 0.375 | 3.708e-03 | 2.088e-03 | 2.4e-12 | 3.7e-15 | +1.00e-05 | 4.7e-05 | 5.0e-16 | 4.2e-15 |
| 3 | 0.80 | 0.296 | 0.224 | 4.676e-03 | 2.623e-03 | 3.4e-12 | 4.1e-15 | +1.76e-05 | 5.4e-05 | 5.1e-16 | 3.2e-15 |
| 3 | 0.90 | 0.094 | 0.076 | 6.121e-03 | 3.369e-03 | 6.8e-12 | 4.7e-15 | +3.24e-05 | 6.3e-05 | 8.5e-16 | 4.3e-15 |
| 4 | 0.00 | 1.000 | 1.000 | 7.505e-05 | 4.751e-06 | 4.1e-12 | 6.6e-15 | -1.41e-09 | 8.8e-08 | 3.3e-15 | 1.4e-14 |
| 4 | 0.10 | 0.989 | 0.985 | 7.970e-05 | 1.639e-05 | 4.5e-12 | 6.7e-15 | -6.93e-10 | 8.8e-08 | 3.7e-15 | 1.4e-14 |
| 4 | 0.20 | 0.956 | 0.941 | 9.257e-05 | 3.249e-05 | 4.5e-12 | 6.6e-15 | +1.49e-09 | 8.9e-08 | 3.5e-15 | 1.4e-14 |
| 4 | 0.30 | 0.901 | 0.872 | 1.120e-04 | 5.002e-05 | 4.7e-12 | 6.8e-15 | +5.35e-09 | 9.0e-08 | 3.3e-15 | 1.3e-14 |
| 4 | 0.40 | 0.824 | 0.776 | 1.378e-04 | 6.943e-05 | 4.7e-12 | 6.9e-15 | +1.14e-08 | 9.5e-08 | 3.3e-15 | 1.2e-14 |
| 4 | 0.60 | 0.604 | 0.521 | 2.179e-04 | 1.177e-04 | 8.1e-12 | 7.3e-15 | +3.76e-08 | 1.2e-07 | 3.4e-15 | 2.0e-14 |
| 4 | 0.80 | 0.296 | 0.224 | 3.770e-04 | 1.915e-04 | 1.2e-11 | 8.4e-15 | +1.25e-07 | 2.1e-07 | 3.4e-15 | 2.6e-14 |

### Table 2.2 Invalid / degenerate meshes (valid=0: some detJ <= 0)

| p | N | q | variant | eps | min detJ | #points detJ<=0 |
|---|---|---|---|---|---|---|
| 1 | 8 | 1 | f1 | 0.8 | -0.042 | 8 |
| 1 | 8 | 1 | f2 | 0.8 | -0.042 | 64 |
| 1 | 8 | 1 | f3 | 0.6 | -0.361 | 16 |
| 1 | 8 | 2 | f3 | 0.8 | -0.176 | 360 |
| 1 | 16 | 1 | f2 | 0.8 | -0.042 | 64 |
| 1 | 16 | 2 | f1 | 1.0 | -0.033 | 344 |
| 1 | 16 | 2 | f2 | 1.0 | -0.149 | 2752 |
| 1 | 16 | 2 | f3 | 0.8 | -0.112 | 2008 |
| 3 | 8 | 3 | f1 | 1.0 | -0.011 | 160 |
| 3 | 8 | 3 | f2 | 1.0 | -0.037 | 1536 |
| 3 | 8 | 3 | f3 | 0.8 | -0.073 | 368 |
| 4 | 8 | 1 | f3 | 0.8 | -1.643 | 112 |
| 4 | 8 | 4 | f3 | 0.8 | -0.076 | 680 |

### Table 2.3 Geometric order and variant (N=8, abc, int): relative B error at selected eps

| p | q | variant | eps | min detJ | rel err B | div rel L2 | energy rel err |
|---|---|---|---|---|---|---|---|
| 1 | 1 | f1 | 0.3 | 0.877 | 2.334e-01 | 1.9e-16 | -5.00e-02 |
| 1 | 1 | f2 | 0.3 | 0.877 | 2.400e-01 | 1.9e-16 | -4.55e-02 |
| 1 | 1 | f3 | 0.3 | 0.707 | 2.373e-01 | 1.7e-16 | -4.86e-02 |
| 1 | 1 | f1 | 0.6 | 0.451 | 2.644e-01 | 2.0e-16 | -4.58e-02 |
| 1 | 2 | f1 | 0.3 | 0.901 | 2.390e-01 | 1.9e-16 | -4.77e-02 |
| 1 | 2 | f2 | 0.3 | 0.878 | 2.594e-01 | 1.8e-16 | -3.69e-02 |
| 1 | 2 | f3 | 0.3 | 0.835 | 2.478e-01 | 1.8e-16 | -4.40e-02 |
| 1 | 2 | f1 | 0.6 | 0.604 | 2.892e-01 | 2.0e-16 | -3.36e-02 |
| 2 | 1 | f1 | 0.3 | 0.877 | 2.397e-02 | 1.2e-15 | -2.07e-04 |
| 2 | 1 | f2 | 0.3 | 0.877 | 2.557e-02 | 1.2e-15 | -1.27e-04 |
| 2 | 1 | f3 | 0.3 | 0.707 | 2.478e-02 | 1.2e-15 | -1.72e-04 |
| 2 | 1 | f1 | 0.6 | 0.451 | 2.805e-02 | 1.4e-15 | -1.43e-05 |
| 2 | 2 | f1 | 0.3 | 0.901 | 2.500e-02 | 1.3e-15 | -1.97e-04 |
| 2 | 2 | f2 | 0.3 | 0.878 | 3.425e-02 | 1.2e-15 | +2.44e-04 |
| 2 | 2 | f3 | 0.3 | 0.835 | 2.806e-02 | 1.2e-15 | -6.28e-05 |
| 2 | 2 | f1 | 0.6 | 0.604 | 3.343e-02 | 1.3e-15 | +1.50e-04 |
| 3 | 1 | f1 | 0.3 | 0.877 | 1.665e-03 | 3.2e-15 | -2.92e-07 |
| 3 | 1 | f2 | 0.3 | 0.877 | 1.891e-03 | 3.2e-15 | +5.13e-07 |
| 3 | 1 | f3 | 0.3 | 0.707 | 1.759e-03 | 3.2e-15 | +5.16e-08 |
| 3 | 1 | f1 | 0.6 | 0.451 | 2.113e-03 | 3.5e-15 | +1.38e-06 |
| 3 | 2 | f1 | 0.3 | 0.901 | 1.840e-03 | 3.3e-15 | +3.02e-07 |
| 3 | 2 | f2 | 0.3 | 0.878 | 3.570e-03 | 3.2e-15 | +9.55e-06 |
| 3 | 2 | f3 | 0.3 | 0.835 | 2.461e-03 | 3.3e-15 | +2.80e-06 |
| 3 | 2 | f1 | 0.6 | 0.604 | 3.009e-03 | 3.6e-15 | +5.63e-06 |
| 3 | 3 | f1 | 0.3 | 0.910 | 1.986e-03 | 3.2e-15 | +6.76e-07 |
| 3 | 3 | f2 | 0.3 | 0.907 | 4.873e-03 | 3.2e-15 | +1.79e-05 |
| 3 | 3 | f3 | 0.3 | 0.849 | 3.150e-03 | 3.2e-15 | +6.18e-06 |
| 3 | 3 | f1 | 0.6 | 0.639 | 3.606e-03 | 3.5e-15 | +8.70e-06 |
| 4 | 1 | f1 | 0.3 | 0.877 | 8.725e-05 | 6.7e-15 | +5.69e-10 |
| 4 | 1 | f1 | 0.6 | 0.451 | 1.197e-04 | 7.2e-15 | +7.27e-09 |
| 4 | 2 | f1 | 0.3 | 0.901 | 1.120e-04 | 6.8e-15 | +5.35e-09 |
| 4 | 2 | f1 | 0.6 | 0.604 | 2.179e-04 | 7.3e-15 | +3.76e-08 |
| 4 | 4 | f1 | 0.3 | 0.910 | 1.512e-04 | 6.8e-15 | +1.51e-08 |
| 4 | 4 | f1 | 0.6 | 0.640 | 3.638e-04 | 7.4e-15 | +1.18e-07 |

### Table 2.4 Pointwise vs integrated projection (N=8, q=2, f1, abc)

| p | eps | rel err B (pt) | rel err B (int) | rel err A (pt) | rel err A (int) | E err (pt) | E err (int) |
|---|---|---|---|---|---|---|---|
| 1 | 0.1 | 2.260e-01 | 2.260e-01 | 6.003e-02 | 6.005e-02 | -5.0e-02 | -5.0e-02 |
| 1 | 0.2 | 2.309e-01 | 2.308e-01 | 7.204e-02 | 7.207e-02 | -4.9e-02 | -4.9e-02 |
| 1 | 0.3 | 2.392e-01 | 2.390e-01 | 8.915e-02 | 8.917e-02 | -4.7e-02 | -4.8e-02 |
| 1 | 0.5 | 2.679e-01 | 2.675e-01 | 1.337e-01 | 1.334e-01 | -3.8e-02 | -4.1e-02 |
| 2 | 0.1 | 2.302e-02 | 2.302e-02 | 3.951e-03 | 3.943e-03 | -2.5e-04 | -2.5e-04 |
| 2 | 0.2 | 2.373e-02 | 2.373e-02 | 6.365e-03 | 6.346e-03 | -2.3e-04 | -2.4e-04 |
| 2 | 0.3 | 2.500e-02 | 2.500e-02 | 9.220e-03 | 9.189e-03 | -1.9e-04 | -2.0e-04 |
| 2 | 0.5 | 2.969e-02 | 2.968e-02 | 1.592e-02 | 1.587e-02 | +1.1e-05 | -2.4e-05 |
| 3 | 0.1 | 1.550e-03 | 1.550e-03 | 2.694e-04 | 2.691e-04 | -6.6e-07 | -6.6e-07 |
| 3 | 0.2 | 1.656e-03 | 1.656e-03 | 4.958e-04 | 4.950e-04 | -3.2e-07 | -3.2e-07 |
| 3 | 0.3 | 1.840e-03 | 1.840e-03 | 7.425e-04 | 7.413e-04 | +3.1e-07 | +3.0e-07 |
| 3 | 0.5 | 2.493e-03 | 2.493e-03 | 1.316e-03 | 1.313e-03 | +3.0e-06 | +3.0e-06 |
| 4 | 0.1 | 7.970e-05 | 7.970e-05 | 1.640e-05 | 1.639e-05 | -6.9e-10 | -6.9e-10 |
| 4 | 0.2 | 9.256e-05 | 9.257e-05 | 3.251e-05 | 3.249e-05 | +1.5e-09 | +1.5e-09 |
| 4 | 0.3 | 1.120e-04 | 1.120e-04 | 5.006e-05 | 5.002e-05 | +5.4e-09 | +5.3e-09 |

### Table 2.5 Mean field B0 = (0.3,-0.2,0.5), abc, N=8, q=2 (slice flux; B0 projection exactness)

| proj | p | variant | eps | b0-proj err L2 | b0-proj max&#124;D_h&#124; | flux err (max of 6 slices) | flux of curl part | max&#124;D_h b&#124; | rel err B |
|---|---|---|---|---|---|---|---|---|---|
| int | 1 | f1 | 0.0 | 1.7e-15 | 7.1e-15 | 5.6e-16 | 8.9e-16 | 4.3e-14 | 2.24e-01 |
| int | 1 | f1 | 0.1 | 1.4e-02 | 8.9e-15 | 6.4e-16 | 4.3e-16 | 5.0e-14 | 2.25e-01 |
| int | 1 | f1 | 0.2 | 2.8e-02 | 1.2e-14 | 1.2e-15 | 7.1e-16 | 5.6e-14 | 2.30e-01 |
| int | 1 | f1 | 0.3 | 4.2e-02 | 1.1e-14 | 3.6e-16 | 4.3e-16 | 5.7e-14 | 2.39e-01 |
| int | 1 | f1 | 0.4 | 5.8e-02 | 1.3e-14 | 8.3e-16 | 3.1e-16 | 5.7e-14 | 2.51e-01 |
| int | 1 | f1 | 0.5 | 7.5e-02 | 1.1e-14 | 6.1e-16 | 7.2e-16 | 5.7e-14 | 2.67e-01 |
| int | 1 | f1 | 0.6 | 9.4e-02 | 1.0e-14 | 8.3e-16 | 5.4e-16 | 7.1e-14 | 2.89e-01 |
| int | 1 | f3 | 0.2 | 3.7e-02 | 1.1e-14 | 8.3e-16 | 3.3e-16 | 5.7e-14 | 2.34e-01 |
| int | 1 | f3 | 0.4 | 8.0e-02 | 1.0e-14 | 6.1e-16 | 2.8e-16 | 5.0e-14 | 2.66e-01 |
| int | 1 | f3 | 0.6 | 1.3e-01 | 1.6e-14 | 1.1e-15 | 4.9e-16 | 6.4e-14 | 3.28e-01 |
| int | 2 | f1 | 0.0 | 2.3e-15 | 2.8e-14 | 9.4e-16 | 3.3e-16 | 3.9e-13 | 2.27e-02 |
| int | 2 | f1 | 0.1 | 1.0e-03 | 3.7e-14 | 6.1e-16 | 6.4e-16 | 4.1e-13 | 2.30e-02 |
| int | 2 | f1 | 0.2 | 2.0e-03 | 4.1e-14 | 1.0e-15 | 8.3e-16 | 4.7e-13 | 2.37e-02 |
| int | 2 | f1 | 0.3 | 3.2e-03 | 3.6e-14 | 7.2e-16 | 6.5e-16 | 4.4e-13 | 2.49e-02 |
| int | 2 | f1 | 0.4 | 4.5e-03 | 4.8e-14 | 3.9e-16 | 4.3e-16 | 4.2e-13 | 2.69e-02 |
| int | 2 | f1 | 0.5 | 5.9e-03 | 5.3e-14 | 8.9e-16 | 3.5e-16 | 5.8e-13 | 2.96e-02 |
| int | 2 | f1 | 0.6 | 7.7e-03 | 4.7e-14 | 6.1e-16 | 5.7e-16 | 5.9e-13 | 3.34e-02 |
| int | 2 | f3 | 0.2 | 4.4e-03 | 3.6e-14 | 5.0e-16 | 7.8e-16 | 4.8e-13 | 2.50e-02 |
| int | 2 | f3 | 0.4 | 1.0e-02 | 5.0e-14 | 1.0e-15 | 3.9e-16 | 4.4e-13 | 3.25e-02 |
| int | 2 | f3 | 0.6 | 1.8e-02 | 1.1e-13 | 7.2e-16 | 3.8e-16 | 1.8e-12 | 4.80e-02 |
| int | 3 | f1 | 0.0 | 2.0e-15 | 8.5e-14 | 8.9e-16 | 5.0e-16 | 1.3e-12 | 1.51e-03 |
| int | 3 | f1 | 0.1 | 1.2e-05 | 1.2e-13 | 5.0e-16 | 3.7e-16 | 1.5e-12 | 1.55e-03 |
| int | 3 | f1 | 0.2 | 4.8e-05 | 8.0e-14 | 4.4e-16 | 5.8e-16 | 1.5e-12 | 1.65e-03 |
| int | 3 | f1 | 0.3 | 1.1e-04 | 1.1e-13 | 1.1e-15 | 2.5e-16 | 1.5e-12 | 1.84e-03 |
| int | 3 | f1 | 0.4 | 2.0e-04 | 1.1e-13 | 7.2e-16 | 4.6e-16 | 1.8e-12 | 2.11e-03 |
| int | 3 | f1 | 0.5 | 3.1e-04 | 1.3e-13 | 8.3e-16 | 3.7e-16 | 1.8e-12 | 2.49e-03 |
| int | 3 | f1 | 0.6 | 4.6e-04 | 1.1e-13 | 3.9e-16 | 3.9e-16 | 2.2e-12 | 3.00e-03 |
| int | 3 | f3 | 0.2 | 1.2e-04 | 1.0e-13 | 3.9e-16 | 5.1e-16 | 1.3e-12 | 1.94e-03 |
| int | 3 | f3 | 0.4 | 5.2e-04 | 1.2e-13 | 6.7e-16 | 5.3e-16 | 2.1e-12 | 3.19e-03 |
| int | 3 | f3 | 0.6 | 1.3e-03 | 1.9e-13 | 8.3e-16 | 3.8e-16 | 4.6e-12 | 5.65e-03 |
| int | 4 | f1 | 0.0 | 2.2e-15 | 1.0e-13 | 4.7e-15 | 3.3e-15 | 4.3e-12 | 7.49e-05 |
| int | 4 | f1 | 0.2 | 2.3e-15 | 1.7e-13 | 3.6e-15 | 3.5e-15 | 4.6e-12 | 9.23e-05 |
| int | 4 | f1 | 0.4 | 2.3e-15 | 2.4e-13 | 4.1e-15 | 3.3e-15 | 4.7e-12 | 1.38e-04 |
| int | 4 | f1 | 0.6 | 2.5e-15 | 2.1e-13 | 3.9e-15 | 3.4e-15 | 8.3e-12 | 2.17e-04 |
| int | 4 | f3 | 0.4 | 2.3e-15 | 2.2e-13 | 3.7e-15 | 3.6e-15 | 6.6e-12 | 2.50e-04 |
| pt | 1 | f1 | 0.0 | 1.8e-15 | 0.0e+00 | 1.1e-15 | 1.9e-16 | 3.7e-14 | 2.24e-01 |
| pt | 1 | f1 | 0.1 | 1.4e-02 | 2.7e-03 | 1.4e-15 | 3.3e-16 | 2.7e-03 | 2.25e-01 |
| pt | 1 | f1 | 0.2 | 2.8e-02 | 1.1e-02 | 2.8e-16 | 2.9e-16 | 1.1e-02 | 2.30e-01 |
| pt | 1 | f1 | 0.3 | 4.2e-02 | 2.5e-02 | 1.1e-15 | 7.4e-16 | 2.5e-02 | 2.39e-01 |
| pt | 1 | f1 | 0.4 | 5.8e-02 | 4.6e-02 | 8.3e-16 | 3.6e-16 | 4.6e-02 | 2.51e-01 |
| pt | 1 | f1 | 0.5 | 7.5e-02 | 7.4e-02 | 7.2e-16 | 6.8e-16 | 7.4e-02 | 2.67e-01 |
| pt | 1 | f1 | 0.6 | 9.4e-02 | 1.1e-01 | 6.1e-16 | 7.6e-16 | 1.1e-01 | 2.89e-01 |
| pt | 1 | f3 | 0.2 | 3.8e-02 | 5.9e-02 | 3.3e-16 | 5.1e-16 | 5.9e-02 | 2.34e-01 |
| pt | 1 | f3 | 0.4 | 8.0e-02 | 2.7e-01 | 3.9e-16 | 6.4e-16 | 2.7e-01 | 2.68e-01 |
| pt | 1 | f3 | 0.6 | 1.4e-01 | 8.9e-01 | 7.8e-16 | 8.9e-16 | 8.9e-01 | 3.31e-01 |
| pt | 2 | f1 | 0.0 | 3.3e-15 | 1.6e-13 | 5.0e-16 | 1.1e-15 | 4.3e-13 | 2.27e-02 |
| pt | 2 | f1 | 0.1 | 1.0e-03 | 1.6e-04 | 5.0e-16 | 5.4e-16 | 1.6e-04 | 2.30e-02 |
| pt | 2 | f1 | 0.2 | 2.0e-03 | 6.7e-04 | 1.1e-15 | 4.7e-16 | 6.7e-04 | 2.37e-02 |
| pt | 2 | f1 | 0.3 | 3.2e-03 | 1.6e-03 | 7.5e-16 | 5.3e-16 | 1.6e-03 | 2.49e-02 |
| pt | 2 | f1 | 0.4 | 4.5e-03 | 2.9e-03 | 8.9e-16 | 1.2e-15 | 2.9e-03 | 2.69e-02 |
| pt | 2 | f1 | 0.5 | 6.0e-03 | 5.0e-03 | 5.0e-16 | 5.7e-16 | 5.0e-03 | 2.96e-02 |
| pt | 2 | f1 | 0.6 | 7.7e-03 | 8.4e-03 | 3.9e-16 | 5.1e-16 | 8.4e-03 | 3.34e-02 |
| pt | 2 | f3 | 0.2 | 4.5e-03 | 5.1e-03 | 7.2e-16 | 4.4e-16 | 5.1e-03 | 2.50e-02 |
| pt | 2 | f3 | 0.4 | 1.0e-02 | 2.4e-02 | 8.9e-16 | 3.8e-16 | 2.4e-02 | 3.25e-02 |
| pt | 2 | f3 | 0.6 | 1.8e-02 | 8.0e-02 | 5.0e-16 | 3.5e-16 | 8.0e-02 | 4.80e-02 |
| pt | 3 | f1 | 0.0 | 5.1e-15 | 3.1e-13 | 8.9e-16 | 3.1e-16 | 1.5e-12 | 1.51e-03 |
| pt | 3 | f1 | 0.1 | 1.2e-05 | 3.7e-13 | 3.3e-16 | 6.5e-16 | 1.5e-12 | 1.55e-03 |
| pt | 3 | f1 | 0.2 | 4.8e-05 | 4.4e-13 | 6.7e-16 | 7.8e-16 | 1.4e-12 | 1.65e-03 |
| pt | 3 | f1 | 0.3 | 1.1e-04 | 3.3e-13 | 3.9e-16 | 8.0e-16 | 1.5e-12 | 1.84e-03 |
| pt | 3 | f1 | 0.4 | 2.0e-04 | 4.1e-13 | 7.2e-16 | 6.8e-16 | 1.6e-12 | 2.11e-03 |
| pt | 3 | f1 | 0.5 | 3.1e-04 | 4.6e-13 | 8.9e-16 | 6.7e-16 | 1.7e-12 | 2.49e-03 |
| pt | 3 | f1 | 0.6 | 4.6e-04 | 5.0e-13 | 1.1e-15 | 8.6e-16 | 2.4e-12 | 3.00e-03 |
| pt | 3 | f3 | 0.2 | 1.2e-04 | 3.8e-13 | 4.4e-16 | 5.7e-16 | 1.6e-12 | 1.94e-03 |
| pt | 3 | f3 | 0.4 | 5.2e-04 | 3.8e-13 | 4.4e-16 | 6.0e-16 | 2.1e-12 | 3.19e-03 |
| pt | 3 | f3 | 0.6 | 1.3e-03 | 8.8e-13 | 6.1e-16 | 5.8e-16 | 3.1e-12 | 5.65e-03 |
| pt | 4 | f1 | 0.0 | 2.4e-15 | 6.4e-13 | 3.3e-15 | 3.3e-15 | 4.4e-12 | 7.49e-05 |
| pt | 4 | f1 | 0.2 | 2.6e-15 | 7.4e-13 | 3.2e-15 | 3.2e-15 | 6.8e-12 | 9.23e-05 |
| pt | 4 | f1 | 0.4 | 2.7e-15 | 8.7e-13 | 3.6e-15 | 3.3e-15 | 6.2e-12 | 1.37e-04 |
| pt | 4 | f1 | 0.6 | 2.9e-15 | 1.0e-12 | 3.1e-15 | 3.5e-15 | 8.2e-12 | 2.17e-04 |
| pt | 4 | f3 | 0.4 | 2.7e-15 | 8.5e-13 | 3.4e-15 | 3.1e-15 | 6.8e-12 | 2.50e-04 |

### Table 2.6 mod2 (multi-mode, non-Beltrami), N=8, q=2, f1, int

| p | eps | rel err B | rel div L2 | energy rel err | helicity abs err (exact 0) |
|---|---|---|---|---|---|
| 1 | 0.0 | 3.148e-01 | 1.6e-16 | -2.27e-01 | 7.4e-16 |
| 1 | 0.1 | 3.171e-01 | 1.6e-16 | -2.25e-01 | 1.4e-05 |
| 1 | 0.2 | 3.241e-01 | 1.4e-16 | -2.20e-01 | 1.2e-04 |
| 1 | 0.3 | 3.359e-01 | 1.5e-16 | -2.12e-01 | 4.2e-04 |
| 1 | 0.4 | 3.530e-01 | 1.6e-16 | -2.00e-01 | 1.1e-03 |
| 1 | 0.5 | 3.760e-01 | 1.6e-16 | -1.84e-01 | 2.4e-03 |
| 1 | 0.6 | 4.060e-01 | 1.6e-16 | -1.62e-01 | 4.7e-03 |
| 2 | 0.0 | 5.094e-02 | 7.7e-16 | -4.83e-03 | 1.2e-15 |
| 2 | 0.1 | 5.219e-02 | 8.5e-16 | -4.78e-03 | 1.6e-06 |
| 2 | 0.2 | 5.595e-02 | 8.4e-16 | -4.58e-03 | 1.5e-05 |
| 2 | 0.3 | 6.221e-02 | 8.8e-16 | -4.20e-03 | 6.0e-05 |
| 2 | 0.4 | 7.104e-02 | 8.7e-16 | -3.52e-03 | 1.7e-04 |
| 2 | 0.5 | 8.262e-02 | 8.9e-16 | -2.36e-03 | 4.0e-04 |
| 2 | 0.6 | 9.733e-02 | 9.7e-16 | -4.00e-04 | 7.9e-04 |
| 3 | 0.0 | 6.326e-03 | 2.1e-15 | -6.42e-05 | 5.1e-16 |
| 3 | 0.1 | 6.638e-03 | 2.1e-15 | -6.39e-05 | 1.4e-08 |
| 3 | 0.2 | 7.559e-03 | 2.1e-15 | -6.21e-05 | 9.8e-08 |
| 3 | 0.3 | 9.082e-03 | 2.2e-15 | -5.61e-05 | 2.6e-07 |
| 3 | 0.4 | 1.125e-02 | 2.2e-15 | -4.03e-05 | 5.2e-07 |
| 3 | 0.5 | 1.420e-02 | 2.2e-15 | -5.64e-06 | 1.0e-06 |
| 3 | 0.6 | 1.815e-02 | 2.4e-15 | +6.44e-05 | 2.6e-06 |

### Table 2.7 Cost vs eps (field abc, q=2, f1, int; N=8)

| p | eps | ND dofs | t assemble | t C_h apply | RSS MB | t total |
|---|---|---|---|---|---|---|
| 1 | 0.0 | 1536 | 0.014 | 3.96e-06 | 31 | 3.1 |
| 1 | 0.3 | 1536 | 0.002 | 1.30e-05 | 31 | 3.7 |
| 2 | 0.0 | 12288 | 0.039 | 8.68e-05 | 54 | 5.3 |
| 2 | 0.3 | 12288 | 0.034 | 8.20e-05 | 53 | 5.9 |
| 3 | 0.0 | 41472 | 0.339 | 1.08e-03 | 62 | 25.6 |
| 3 | 0.3 | 41472 | 0.368 | 9.57e-04 | 62 | 25.7 |
| 4 | 0.0 | 98304 | 1.825 | 3.86e-03 | 115 | 52.1 |
| 4 | 0.3 | 98304 | 1.528 | 3.36e-03 | 116 | 51.5 |

### Table 2.8 h-convergence on deformed meshes (f1, abc, int, no B0, N=8 -> 16; p <= 2)

| p | q | eps | min detJ (N=16) | err B N=8 | err B N=16 | order B | err A N=8 | err A N=16 | order A | max&#124;D_h b&#124; N=16 |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 1 | 0.0 | 1.000 | 1.994e+00 | 1.005e+00 | 0.99 | 7.857e-02 | 1.984e-02 | 1.99 | 1.7e-13 |
| 1 | 1 | 0.1 | 0.989 | 2.003e+00 | 1.010e+00 | 0.99 | 8.175e-02 | 2.282e-02 | 1.84 | 1.7e-13 |
| 1 | 1 | 0.2 | 0.956 | 2.029e+00 | 1.024e+00 | 0.99 | 9.059e-02 | 2.999e-02 | 1.59 | 1.9e-13 |
| 1 | 1 | 0.3 | 0.897 | 2.074e+00 | 1.050e+00 | 0.98 | 1.038e-01 | 3.927e-02 | 1.40 | 2.0e-13 |
| 1 | 1 | 0.5 | 0.694 | 2.230e+00 | 1.140e+00 | 0.97 | 1.400e-01 | 6.154e-02 | 1.19 | 1.9e-13 |
| 1 | 2 | 0.0 | 1.000 | 1.994e+00 | 1.005e+00 | 0.99 | 7.857e-02 | 1.984e-02 | 1.99 | 1.7e-13 |
| 1 | 2 | 0.1 | 0.990 | 2.008e+00 | 1.012e+00 | 0.99 | 8.492e-02 | 2.556e-02 | 1.73 | 2.1e-13 |
| 1 | 2 | 0.2 | 0.959 | 2.051e+00 | 1.035e+00 | 0.99 | 1.019e-01 | 3.807e-02 | 1.42 | 2.2e-13 |
| 1 | 2 | 0.3 | 0.908 | 2.124e+00 | 1.074e+00 | 0.98 | 1.261e-01 | 5.327e-02 | 1.24 | 1.9e-13 |
| 1 | 2 | 0.5 | 0.744 | 2.377e+00 | 1.209e+00 | 0.98 | 1.886e-01 | 8.870e-02 | 1.09 | 2.3e-13 |
| 2 | 1 | 0.0 | 1.000 | 2.025e-01 | 5.096e-02 | 1.99 | 3.906e-03 | 4.914e-04 | 2.99 | 1.8e-12 |
| 2 | 1 | 0.1 | 0.989 | 2.037e-01 | 5.127e-02 | 1.99 | 4.234e-03 | 6.440e-04 | 2.72 | 1.9e-12 |
| 2 | 1 | 0.2 | 0.956 | 2.071e-01 | 5.221e-02 | 1.99 | 5.161e-03 | 9.930e-04 | 2.38 | 1.8e-12 |
| 2 | 1 | 0.3 | 0.897 | 2.130e-01 | 5.382e-02 | 1.98 | 6.550e-03 | 1.439e-03 | 2.19 | 1.9e-12 |
| 2 | 1 | 0.5 | 0.694 | 2.333e-01 | 5.944e-02 | 1.97 | 1.044e-02 | 2.572e-03 | 2.02 | 2.8e-12 |
| 2 | 2 | 0.0 | 1.000 | 2.025e-01 | 5.096e-02 | 1.99 | 3.906e-03 | 4.914e-04 | 2.99 | 2.0e-12 |
| 2 | 2 | 0.1 | 0.990 | 2.045e-01 | 5.148e-02 | 1.99 | 5.576e-03 | 1.116e-03 | 2.32 | 1.8e-12 |
| 2 | 2 | 0.2 | 0.959 | 2.108e-01 | 5.310e-02 | 1.99 | 8.974e-03 | 2.093e-03 | 2.10 | 2.0e-12 |
| 2 | 2 | 0.3 | 0.908 | 2.221e-01 | 5.600e-02 | 1.99 | 1.299e-02 | 3.160e-03 | 2.04 | 2.0e-12 |
| 2 | 2 | 0.5 | 0.744 | 2.637e-01 | 6.671e-02 | 1.98 | 2.244e-02 | 5.589e-03 | 2.01 | 2.7e-12 |

Overall, all valid runs except the pointwise-B0 series (400 rows): max|C_h G_h| = 7.1e-15; max rel div L2 = 8.4e-15; max div_max/(|B|/h) = 2.5e-13; max max|D_h b| = 1.8e-11; max flux of curl part = 3.8e-15; max flux err = 4.7e-15; max (max_DC_scaled) = 2.9e-14
Pointwise-B0 series (35 rows): max rel div L2 = 1.6e-03, max max|D_h b| = 8.9e-01, max flux err = 3.6e-15
