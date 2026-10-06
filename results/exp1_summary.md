# E1 summary: static representation B_h = Pi_RT(B0) + C_h a_h on the undeformed periodic box

Data: `results/exp1_static.csv` (one row per run, 1 MPI rank, `scripts/run_exp1.sh`), figure `results/exp1_convergence.png`
(`scripts/plot_exp1.py`), tables below by `scripts/make_tables.py exp1`.
Setting: unit torus, N^3 hexahedra (N = 4,8,16, and 32 for p <= 2), ND_p / RT_{p-1} / L2_{p-1} / H1_p, fields `abc`
(Beltrami, a,b,c = 1, 0.8, 0.6) and `mod2` (non-Beltrami, 2pi and 4pi modes), geometric order q = 1 (the undeformed box is
affine; a check with q = 3 reproduces errB to 3e-11 and errA to 1e-7 relative, i.e. within roundoff of the tiny errors), error norms with quadrature order 2p+2q+2.
Timings were taken with up to two other single-rank jobs running on the 4-core machine, so they are indicative only.

## What was established

1. **Divergence is at roundoff for all p, N.** For b = C_h a_h: max over true dofs of |(D_h b)_i| <= 2.0e-11 (p=4, N=16), and
   ||D_h b||_L2 / (||B||_L2/h) <= 1.3e-14 (max-version <= 1.4e-13). The sparse products satisfy max|C_h G_h| <= 7.1e-15 (exactly 0
   for p <= 3) and max|D_h C_h| * h^3 <= 1.4e-14. The raw max|D_h C_h| grows like N^3 (5.8e-11 at p=4, N=16) only because
   `DivergenceInterpolator` returns L2 *values* of the divergence, i.e. contains the 1/detJ = N^3 factor; this is roundoff
   amplified by the geometry factor, not a defect (the scaled number is N-independent).
2. **Convergence orders.** B_h = C_h a_h converges as O(h^p) (observed 0.96-1.00, 1.96-2.00, 2.97-2.99, 3.97-3.99 for p = 1..4; finest pairs 1.00 / 2.00 / 2.99 / 3.99 for abc), i.e. the expected RT_{p-1} rate: no order is lost by going through A. A_h converges at the *better* rate O(h^{p+1})
   (finest pairs 2.00, 3.00, 3.99, 4.99; ND_p contains the full P_p), not O(h^p) as stated in the plan; the plan's h^p is only a lower bound.
   `mod2` is pre-asymptotic at N=4 (the 4pi mode has 2 cells per wavelength; e.g. p=1 A-order 1.65, p=2 2.70) and reaches the same rates by N=16.
3. **Energy and helicity converge at the superconvergent rate 2p** (abc: 1.99, 3.99, 5.99, 7.99 for the finest pair, both energy and
   helicity); at p=4, N=16 the relative energy error is 5.5e-12 and the helicity error 3.5e-10. For `mod`/`mod2` the helicity is exactly 0 by symmetry
   and the discrete value is <= 6e-15 (trivially satisfied, so these fields do not test helicity).
4. **Mean-field and flux.** On the affine undeformed mesh Pi_RT(B0) is exact (L2 error <= 1.4e-15, max|D_h Pi_RT(B0)| <= 4.8e-13 =
   N^3 x roundoff) and the net flux through all six slice families (x_i = 0 and x_i = 1/2) equals B0_i to <= 7e-15; the flux of C_h a_h
   through every slice is <= 7e-15, as the discrete Stokes argument predicts.
5. **Cost (1 rank).** One C_h application: 0.5 us (192 dofs) ... 6.5 ms (p=2, N=32, 786k ND dofs), 21 ms (p=4, N=16, 786k dofs); assembling G_h, C_h, D_h:
   2.0 s (p=2, N=32), 7.0 s (p=4, N=16). Jacobi-CG iterations on the ND mass matrix (rel. tol 1e-12, random rhs) grow slowly and saturate:
   6/15/34/40 (p=1), 10/21/28/29 (p=2), 10/20/24 (p=3), 9/19/22 (p=4) for increasing N, i.e. the ND mass matrix is well conditioned (iteration counts become independent of h; no condition numbers were measured directly). Peak
   RSS for 786k ND dofs is 0.7 GB (p=2/4) incl. the three sparse operators; total memory is far below the 10 GB budget.

## Caveats and surprising points

- `-proj pt` and `-proj int` give *identical* results on the undeformed mesh (to every printed digit). This is not a bug: for both test fields A_i does not depend on x_i
  (e.g. A_x = a sin kz + c cos ky), so on a Cartesian mesh every edge/interior dof functional (a line integral along a coordinate direction) equals the point value times the
  length. The two projections differ only on deformed meshes (exp2; commuting defect O(h^p) for the pointwise one) and for gauge-perturbed A. **E1 therefore cannot discriminate the two projections; use E2 for that.**
- MFEM's built-in `ND_HexahedronElement::ProjectIntegrated` / `RT_HexahedronElement::ProjectIntegrated` are `protected`, only reachable through elements built with the
  `IntegratedGLL` open basis, use a quadrature rule of order p only, and the curl/gradient interpolators assume the standard nodal representation. `vp::ProjectA/ProjectB(..., Integrated)` therefore
  evaluates the same functionals (sub-edge circulations, sub-face fluxes between Gauss-Lobatto points) with a (p+q+5)-point Gauss rule and converts them to the standard basis with the local matrix M^{-1}
  (verified in `selftest`: C_h Pi_ND^int A = Pi_RT^int curl A and G_h Pi_H1 chi = Pi_ND^int grad chi to <= 1e-11 on a deformed mesh).
- The design-note statement that constants are in RT_{p-1} is true only for affine elements (confirmed here, undeformed mesh); on deformed meshes see `exp2_summary.md`.
- `max_DC` absolute values in the CSV scale as N^3 (see item 1); use `max_DC_scaled`.
- Quadrature orders: errors 2p+2q+2, exact norms 2p+2q+6. Energy via mass matrix (full assembly, RT) agrees with the matrix-free energy to <= 1e-13 relative (checked for RT dofs <= 40000, `n/a` otherwise).

## Tables

### Table 1.1 Convergence (field, p, N): errors, observed order (between successive N)

pt and int rows are identical to all printed digits on the undeformed mesh (A_i is independent of x_i for both fields, so line integrals along coordinate edges equal point values times h); only `pt` is listed.

| field | p | N | h | ND dofs | err A (L2) | order A | err B (L2) | order B | rel err B |
|---|---|---|---|---|---|---|---|---|---|
| abc | 1 | 4 | 0.25 | 192 | 3.018e-01 | nan | 3.867e+00 | nan | 4.35e-01 |
| abc | 1 | 8 | 0.125 | 1536 | 7.857e-02 | 1.94 | 1.994e+00 | 0.96 | 2.24e-01 |
| abc | 1 | 16 | 0.0625 | 12288 | 1.984e-02 | 1.99 | 1.005e+00 | 0.99 | 1.13e-01 |
| abc | 1 | 32 | 0.03125 | 98304 | 4.973e-03 | 2.00 | 5.033e-01 | 1.00 | 5.66e-02 |
| abc | 2 | 4 | 0.25 | 1536 | 3.044e-02 | nan | 7.898e-01 | nan | 8.89e-02 |
| abc | 2 | 8 | 0.125 | 12288 | 3.906e-03 | 2.96 | 2.025e-01 | 1.96 | 2.28e-02 |
| abc | 2 | 16 | 0.0625 | 98304 | 4.914e-04 | 2.99 | 5.096e-02 | 1.99 | 5.73e-03 |
| abc | 2 | 32 | 0.03125 | 786432 | 6.153e-05 | 3.00 | 1.276e-02 | 2.00 | 1.44e-03 |
| abc | 3 | 4 | 0.25 | 5184 | 2.783e-03 | nan | 1.055e-01 | nan | 1.19e-02 |
| abc | 3 | 8 | 0.125 | 41472 | 1.775e-04 | 3.97 | 1.346e-02 | 2.97 | 1.52e-03 |
| abc | 3 | 16 | 0.0625 | 331776 | 1.115e-05 | 3.99 | 1.692e-03 | 2.99 | 1.90e-04 |
| abc | 4 | 4 | 0.25 | 12288 | 2.114e-04 | nan | 1.048e-02 | nan | 1.18e-03 |
| abc | 4 | 8 | 0.125 | 98304 | 6.719e-06 | 4.98 | 6.669e-04 | 3.97 | 7.51e-05 |
| abc | 4 | 16 | 0.0625 | 786432 | 2.109e-07 | 4.99 | 4.187e-05 | 3.99 | 4.71e-06 |
| mod2 | 1 | 4 | 0.25 | 192 | 4.066e-01 | nan | 5.892e+00 | nan | 6.57e-01 |
| mod2 | 1 | 8 | 0.125 | 1536 | 1.295e-01 | 1.65 | 2.825e+00 | 1.06 | 3.15e-01 |
| mod2 | 1 | 16 | 0.0625 | 12288 | 3.536e-02 | 1.87 | 1.380e+00 | 1.03 | 1.54e-01 |
| mod2 | 1 | 32 | 0.03125 | 98304 | 9.046e-03 | 1.97 | 6.837e-01 | 1.01 | 7.62e-02 |
| mod2 | 2 | 4 | 0.25 | 1536 | 5.750e-02 | nan | 1.744e+00 | nan | 1.94e-01 |
| mod2 | 2 | 8 | 0.125 | 12288 | 8.847e-03 | 2.70 | 4.571e-01 | 1.93 | 5.09e-02 |
| mod2 | 2 | 16 | 0.0625 | 98304 | 1.110e-03 | 2.99 | 1.149e-01 | 1.99 | 1.28e-02 |
| mod2 | 2 | 32 | 0.03125 | 786432 | 1.387e-04 | 3.00 | 2.874e-02 | 2.00 | 3.20e-03 |
| mod2 | 3 | 4 | 0.25 | 5184 | 1.049e-02 | nan | 4.325e-01 | nan | 4.82e-02 |
| mod2 | 3 | 8 | 0.125 | 41472 | 7.398e-04 | 3.83 | 5.677e-02 | 2.93 | 6.33e-03 |
| mod2 | 3 | 16 | 0.0625 | 331776 | 4.711e-05 | 3.97 | 7.172e-03 | 2.98 | 7.99e-04 |
| mod2 | 4 | 4 | 0.25 | 12288 | 1.644e-03 | nan | 8.418e-02 | nan | 9.38e-03 |
| mod2 | 4 | 8 | 0.125 | 98304 | 5.523e-05 | 4.90 | 5.521e-03 | 3.93 | 6.15e-04 |
| mod2 | 4 | 16 | 0.0625 | 786432 | 1.755e-06 | 4.98 | 3.491e-04 | 3.98 | 3.89e-05 |

### Table 1.2 Exactness of the discrete structure (field abc, proj pt)

max_DC_scaled = max|D_h C_h| * hmin^3 (D_h returns L2 *values* of div and so carries a 1/detJ = N^3 factor); div columns are for b = C_h a_h; flux error = max over the 6 slice families |flux - B0_i|.

| p | N | max|C_h G_h| | max|D_h C_h| | max_DC_scaled | max&#124;D_h b&#124; | &#124;D_h b&#124;_L2 | rel div (L2) | rel div (max) | flux err | B0-proj err (L2) | B0-proj div max |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 4 | 0.0e+00 | 0.0e+00 | 0.0e+00 | 1.4e-14 | 4.6e-15 | 1.3e-16 | 4.0e-16 | 6.1e-16 | 1.6e-16 | 0.0e+00 |
| 1 | 8 | 0.0e+00 | 0.0e+00 | 0.0e+00 | 2.8e-14 | 1.0e-14 | 1.4e-16 | 4.0e-16 | 1.9e-16 | 3.2e-16 | 0.0e+00 |
| 1 | 16 | 0.0e+00 | 0.0e+00 | 0.0e+00 | 1.4e-13 | 3.2e-14 | 2.3e-16 | 1.0e-15 | 9.4e-16 | 6.3e-16 | 0.0e+00 |
| 1 | 32 | 0.0e+00 | 0.0e+00 | 0.0e+00 | 6.8e-13 | 1.3e-13 | 4.5e-16 | 2.4e-15 | 5.9e-16 | 1.4e-15 | 0.0e+00 |
| 2 | 4 | 0.0e+00 | 5.7e-14 | 8.9e-16 | 9.5e-14 | 2.5e-14 | 7.2e-16 | 2.7e-15 | 3.9e-16 | 2.3e-16 | 3.6e-15 |
| 2 | 8 | 0.0e+00 | 4.5e-13 | 8.9e-16 | 4.1e-13 | 8.7e-14 | 1.2e-15 | 5.8e-15 | 7.2e-16 | 4.9e-16 | 1.9e-14 |
| 2 | 16 | 0.0e+00 | 3.6e-12 | 8.9e-16 | 1.8e-12 | 3.2e-13 | 2.3e-15 | 1.2e-14 | 8.9e-16 | 8.8e-16 | 8.7e-14 |
| 2 | 32 | 0.0e+00 | 2.9e-11 | 8.9e-16 | 8.2e-12 | 1.3e-12 | 4.6e-15 | 2.9e-14 | 1.7e-15 | 1.9e-15 | 4.3e-13 |
| 3 | 4 | 0.0e+00 | 1.1e-13 | 1.8e-15 | 3.4e-13 | 5.7e-14 | 1.6e-15 | 9.7e-15 | 3.2e-16 | 2.6e-16 | 1.3e-14 |
| 3 | 8 | 0.0e+00 | 1.8e-12 | 3.6e-15 | 1.6e-12 | 2.2e-13 | 3.1e-15 | 2.2e-14 | 5.8e-16 | 5.0e-16 | 6.9e-14 |
| 3 | 16 | 0.0e+00 | 1.5e-11 | 3.6e-15 | 7.2e-12 | 8.7e-13 | 6.1e-15 | 5.1e-14 | 6.3e-16 | 1.0e-15 | 3.2e-13 |
| 4 | 4 | 7.1e-15 | 9.1e-13 | 1.4e-14 | 9.4e-13 | 1.2e-13 | 3.5e-15 | 2.6e-14 | 1.7e-15 | 4.0e-16 | 3.0e-14 |
| 4 | 8 | 7.1e-15 | 7.3e-12 | 1.4e-14 | 4.7e-12 | 4.6e-13 | 6.5e-15 | 6.7e-14 | 4.2e-15 | 6.9e-16 | 1.2e-13 |
| 4 | 16 | 7.1e-15 | 5.8e-11 | 1.4e-14 | 2.0e-11 | 1.9e-12 | 1.3e-14 | 1.4e-13 | 7.0e-15 | 1.4e-15 | 4.8e-13 |

### Table 1.3 Energy and helicity (abc: E_exact = 0.5 k^2 (a^2+b^2+c^2) = 39.478418, H_exact = k(a^2+b^2+c^2) = 12.566371; mod2: E_exact = 40.267986, H_exact = 0)

| field | p | N | energy rel err | order | helicity abs err | order | matrix vs matrix-free energy diff |
|---|---|---|---|---|---|---|---|
| abc | 1 | 4 | -1.89e-01 | nan | 4.57e+00 | nan | -1.1e-14 |
| abc | 1 | 8 | -5.04e-02 | 1.91 | 1.25e+00 | 1.87 | 5.5e-14 |
| abc | 1 | 16 | -1.28e-02 | 1.98 | 3.21e-01 | 1.97 | 4.2e-14 |
| abc | 1 | 32 | -3.21e-03 | 1.99 | 8.06e-02 | 1.99 | -9.5e-14 |
| abc | 2 | 4 | -4.00e-03 | nan | 1.48e-01 | nan | -2.2e-15 |
| abc | 2 | 8 | -2.61e-04 | 3.94 | 9.78e-03 | 3.92 | -1.2e-14 |
| abc | 2 | 16 | -1.65e-05 | 3.99 | 6.20e-04 | 3.98 | 1.0e-13 |
| abc | 2 | 32 | -1.03e-06 | 4.00 | 3.89e-05 | 4.00 | n/a |
| abc | 3 | 4 | -4.74e-05 | nan | 2.36e-03 | nan | 2.3e-15 |
| abc | 3 | 8 | -7.67e-07 | 5.95 | 3.85e-05 | 5.94 | -2.6e-14 |
| abc | 3 | 16 | -1.21e-08 | 5.99 | 6.07e-07 | 5.98 | n/a |
| abc | 4 | 4 | -3.50e-07 | nan | 2.18e-05 | nan | 1.4e-14 |
| abc | 4 | 8 | -1.41e-09 | 7.96 | 8.85e-08 | 7.95 | n/a |
| abc | 4 | 16 | -5.54e-12 | 7.99 | 3.50e-10 | 7.98 | n/a |
| mod2 | 1 | 4 | -6.03e-01 | nan | 1.26e-17 | nan | 7.3e-15 |
| mod2 | 1 | 8 | -2.27e-01 | 1.41 | 2.23e-16 | -4.14 | -3.5e-14 |
| mod2 | 1 | 16 | -6.57e-02 | 1.79 | 2.25e-16 | -0.02 | 1.8e-15 |
| mod2 | 1 | 32 | -1.71e-02 | 1.94 | 6.09e-15 | -4.76 | -2.4e-14 |
| mod2 | 2 | 4 | -6.62e-02 | nan | 8.05e-17 | nan | -9.5e-15 |
| mod2 | 2 | 8 | -4.83e-03 | 3.78 | 5.66e-16 | -2.81 | 2.7e-14 |
| mod2 | 2 | 16 | -3.22e-04 | 3.91 | 4.81e-16 | 0.23 | -5.2e-14 |
| mod2 | 2 | 32 | -2.04e-05 | 3.98 | 3.21e-15 | -2.74 | n/a |
| mod2 | 3 | 4 | -3.43e-03 | nan | 2.96e-16 | nan | 1.4e-14 |
| mod2 | 3 | 8 | -6.42e-05 | 5.74 | 9.48e-16 | -1.68 | 4.9e-14 |
| mod2 | 3 | 16 | -1.05e-06 | 5.93 | 2.87e-15 | -1.60 | n/a |
| mod2 | 4 | 4 | -1.20e-04 | nan | 9.81e-17 | nan | -7.6e-15 |
| mod2 | 4 | 8 | -5.54e-07 | 7.76 | 5.41e-16 | -2.46 | n/a |
| mod2 | 4 | 16 | -2.26e-09 | 7.94 | 3.24e-15 | -2.58 | n/a |

### Table 1.4 Cost (field abc, proj pt, 1 rank; times in seconds, mean of repeated C_h applications)

| p | N | H1 dofs | ND dofs | RT dofs | L2 dofs | t assemble G,C,D | t C_h apply | CG its (Jacobi, ND mass) | RSS MB | t total |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 4 | 64 | 192 | 192 | 64 | 0.000 | 5.09e-07 | 6 | 26 | 0.2 |
| 1 | 8 | 512 | 1536 | 1536 | 512 | 0.002 | 3.49e-06 | 15 | 28 | 0.4 |
| 1 | 16 | 4096 | 12288 | 12288 | 4096 | 0.013 | 2.99e-05 | 34 | 49 | 2.0 |
| 1 | 32 | 32768 | 98304 | 98304 | 32768 | 0.099 | 2.88e-04 | 40 | 219 | 20.2 |
| 2 | 4 | 512 | 1536 | 1536 | 512 | 0.002 | 4.24e-06 | 10 | 29 | 0.3 |
| 2 | 8 | 4096 | 12288 | 12288 | 4096 | 0.018 | 4.20e-05 | 21 | 52 | 0.9 |
| 2 | 16 | 32768 | 98304 | 98304 | 32768 | 0.144 | 4.64e-04 | 28 | 241 | 5.1 |
| 2 | 32 | 262144 | 786432 | 786432 | 262144 | 1.959 | 6.48e-03 | 29 | 709 | 60.7 |
| 3 | 4 | 1728 | 5184 | 5184 | 1728 | 0.013 | 4.50e-05 | 10 | 52 | 1.1 |
| 3 | 8 | 13824 | 41472 | 41472 | 13824 | 0.308 | 4.67e-04 | 20 | 239 | 7.1 |
| 3 | 16 | 110592 | 331776 | 331776 | 110592 | 0.896 | 2.00e-03 | 24 | 272 | 12.6 |
| 4 | 4 | 4096 | 12288 | 12288 | 4096 | 0.098 | 1.79e-04 | 9 | 158 | 4.4 |
| 4 | 8 | 32768 | 98304 | 98304 | 32768 | 0.679 | 1.64e-03 | 19 | 114 | 7.8 |
| 4 | 16 | 262144 | 786432 | 786432 | 262144 | 7.047 | 2.09e-02 | 22 | 677 | 61.3 |

### Table 1.5 Mean field B0 = (0.3,-0.2,0.5), abc (flux error = max |flux - B0_i| over 6 slice families; curl part = max |flux of C_h a|)

| p | N | proj | flux err | flux of curl part | max|D_h b| | B0-proj L2 err |
|---|---|---|---|---|---|---|
| 1 | 4 | int | 5.0e-16 | 3.9e-16 | 1.4e-14 | 4.3e-16 |
| 1 | 4 | pt | 3.3e-16 | 6.1e-16 | 1.8e-14 | 1.6e-16 |
| 1 | 8 | int | 1.5e-15 | 6.7e-16 | 4.5e-14 | 5.0e-16 |
| 1 | 8 | pt | 1.1e-15 | 1.9e-16 | 3.7e-14 | 3.2e-16 |
| 2 | 4 | int | 3.9e-16 | 4.2e-16 | 8.5e-14 | 4.4e-16 |
| 2 | 4 | pt | 3.3e-16 | 3.9e-16 | 1.0e-13 | 2.3e-16 |
| 2 | 8 | int | 1.7e-15 | 5.0e-16 | 4.0e-13 | 5.5e-16 |
| 2 | 8 | pt | 9.4e-16 | 7.2e-16 | 4.4e-13 | 4.9e-16 |
| 3 | 4 | int | 6.4e-16 | 6.0e-16 | 3.3e-13 | 3.1e-16 |
| 3 | 4 | pt | 6.1e-16 | 3.2e-16 | 3.8e-13 | 2.6e-16 |
| 3 | 8 | int | 1.5e-15 | 3.8e-16 | 1.2e-12 | 4.5e-16 |
| 3 | 8 | pt | 8.3e-16 | 5.8e-16 | 1.6e-12 | 5.0e-16 |
| 4 | 4 | int | 2.0e-15 | 1.9e-15 | 1.1e-12 | 4.7e-16 |
| 4 | 4 | pt | 1.8e-15 | 1.7e-15 | 1.0e-12 | 4.0e-16 |
| 4 | 8 | int | 3.7e-15 | 3.5e-15 | 4.5e-12 | 6.9e-16 |
| 4 | 8 | pt | 4.8e-15 | 4.2e-15 | 4.7e-12 | 6.9e-16 |
