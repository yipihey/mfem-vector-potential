# Results and recommendation: vector-potential remapping for ALE MHD in MFEM

Companion to `PLAN.md` and `docs/design_note.md`. Per-experiment detail, tables
and figures: `results/exp1_summary.md` ... `results/exp5_summary.md`,
`results/exp4_gauge_summary.md`, `results/exp7_summary.md`,
`results/mhdale_summary.md`. Build provenance: `BUILD.md`.

## 0. The question and the short answer

Question: for an ALE code whose magnetic state is B (AREPO-like), is it worth
carrying the vector potential A in H(curl) so that rezoning transfers A and
reconstructs B = B0 + curl_h A, instead of transferring B and cleaning its
divergence?

Short answer: **No, not as the primary remap mechanism.** Transferring A gives
exactly divergence-free B and an exactly conserved mean flux for any transfer
operator, including the cheapest local ones. But it does not give better B
accuracy (it is 0 to 4 times worse per remap at the same point-location cost),
it is gauge dependent (every transfer pollutes B with the gauge part of A),
and under repeated remaps the gauge part of A grows without bound unless it is
projected out. With a Coulomb projection the growth stops, but only the
L2-projected A transfer is then stable over 400 remaps, and at p = 3 the
L2-projected-and-cleaned B transfer is 7 times more accurate. That cleaned B
transfer is mathematically the constrained least-squares A projection of the
brief, so the "hybrid constrained projection" and "remap B then clean" are
the same operator on a periodic domain.

The recommendation is in section 9.

## 1. What was built

- MFEM master (ed7e19c, v4.10.1), parallel, with hypre 2.28, METIS 5.1, GSLIB
  375eda2; the MHD-ALE fork of MFEM (842fac6, v4.8.1) and MHD-ALE (0e0564e).
  ex3p, ex4p and MHD-ALE's 2-D and 3-D Taylor-Green were reproduced before any
  change (`BUILD.md`).
- `src/vp_core`: periodic box, smooth deformations, ND_p / RT_{p-1} / L2_{p-1}
  / H1_p spaces, discrete G_h, C_h, D_h, analytic fields, diagnostics.
- `src/vp_transfer`: mesh-to-mesh transfers through GSLIB point location:
  A_pt, A_int, A_l2 (A route) and B_pt, B_int, B_l2, B_l2c (B route; "c" =
  projection onto ker D_h).
- `src/vp_gauge`: Coulomb-gauge projection (AMG-CG on G_h^T M_ND G_h) and a
  local Chebyshev-Jacobi variant.
- Drivers `experiments/exp1..exp5, exp7`, scripts, CSVs and figures in
  `results/`; MHD-ALE instrumentation patch in `mhdale/`.

## 2. Static mesh (Experiment 1)

- D_h C_h = 0 and C_h G_h = 0 hold as sparse-matrix identities (max entry
  7e-15 after scaling out the 1/detJ factor that MFEM's divergence
  interpolator carries). For every (p, N) the divergence of C_h a_h is at
  roundoff: scaled ||D_h b||/(||B||/h) <= 1.3e-14.
- B_h = C_h a_h converges at O(h^p) for p = 1..4 (observed 1.00, 2.00, 2.99,
  3.99), the RT_{p-1} rate; no order is lost by going through A. A_h itself
  converges at O(h^{p+1}) on the affine mesh. Energy and helicity converge at
  2p (p = 4, N = 16: relative energy error 5.5e-12).
- Mean flux through the six topological slice families equals B0 to 7e-15;
  the flux of the curl part is 7e-15 (discrete Stokes on a closed surface).
- The ND mass matrix is well conditioned: Jacobi-CG iterations saturate at
  about 20 to 40 independent of h.

## 3. Distorted meshes (Experiment 2)

- Over 400 valid runs (p = 1..4, geometric order q = 1..4, three deformation
  patterns, down to a normalised minimum Jacobian of 0.09) the scaled
  divergence stays <= 8.4e-15 and the flux error <= 4.7e-15. There is no trend
  with distortion. Divergence-freeness is an algebraic property of the dofs,
  not of the geometry.
- The B error grows mildly (p = 2: 2.3e-2 to 5.6e-2 from eps = 0 to 0.9) and
  keeps order p up to eps = 0.5. The A error loses one order on non-affine
  elements (p+1 to p); B is unaffected.
- A finding against the design note: a constant B0 is in RT_{p-1} only on
  affine elements. On curved elements the pointwise projection of B0 is not
  divergence free (max |D_h Pi^pt B0| up to 0.9 at p = 1), while the
  face-flux projection is (2.4e-13). All later experiments use face-flux dofs
  for B0. The "exact mean field" claim is therefore true of the flux, not of
  the pointwise constant.

## 4. One rezone (Experiment 3)

Transfer from a deformed mesh (f1, eps = 0.3) to a uniform one, N = 8, abc
field with mean field, errors relative to the floor (the error of projecting
the exact field directly onto the target mesh).

| p | A_pt | A_int | A_l2 | B_pt | B_int | B_l2 | B_l2c |
|---|---|---|---|---|---|---|---|
| 2 | 1.53 | 1.28 | 1.33 | 1.47 | 1.24 | 1.16 | 1.15 |
| 3 | 2.88 | 1.61 | 2.05 | 1.55 | 1.27 | 1.21 | 1.20 |
| 4 | 8.9 | 2.95 | 5.05 | 2.17 | 1.68 | 1.42 | 1.42 |

- A-route divergence: at roundoff for every operator (1.3e-14 scaled), mean
  flux exact to 4e-15, even for A_pt whose accuracy is the worst.
- B-route divergence without cleaning: of the same size as the B error
  itself (p = 2: B_pt 2.6e-2, B_l2 9e-3, B_int 6e-3 scaled). B_int's
  divergence is set by quadrature across source-element faces and falls only
  like 1/nq while cost grows like nq^2; no practical quadrature reaches
  roundoff. B_l2c reaches 6e-15 with a global solve. The B route also lets the
  mean flux drift (up to 1e-2 after cleaning at p = 1).
- Energy changes are at the discretisation-error level for all operators;
  none conserves energy. Helicity is perturbed at 1 to 2 orders above its
  discretisation level.
- Gauge: with a smooth O(1) gauge term added to A on the source mesh, the
  transferred B changes by 35 / 14 / 17 % (p = 1), 12.7 / 4.1 / 5.1 % (p = 2),
  2.3 / 0.8 / 1.3 % (p = 3) for A_pt / A_int / A_l2. A rough gauge gives
  order-one pollution. No transfer is gauge invariant.
- Cost: the point-location work is shared by both routes. One-shot at p = 3,
  N = 8: A_int 5 s, B_l2 4 s, B_l2c 14 s, B_int 15 s (contended machine;
  ratios only).

## 5. Repeated rezones (Experiment 4 and the gauge-fixed rerun)

Ping-pong between a uniform and a deformed mesh, N = 8, abc, 100 to 400
transfers. B error at n = 100 in units of the floor, energy drift per 100
remaps, stability:

| operator | p = 2 err | p = 3 err | energy / 100 | stable to n = 400 |
|---|---|---|---|---|
| A_pt, B_pt | diverged by n = 20-30 | diverged | - | no |
| B_int | 9.1 (diverging) | 13 | +0.3 % (p=3) | no |
| A_int (no gauge) | 6.0 | 90 | +1.3 % | no (blow-up after ~200) |
| A_l2 (no gauge) | 5.3 | 12 | +1.0 % | drifting (+58 % energy at 400) |
| A_int + Coulomb | 4.7 | 26 | +0.9 % | no (slow curl-part growth) |
| A_l2 + Coulomb | 4.8 | 19 | +0.4 % | yes (+1.9 % energy at 400) |
| B_l2 | 4.5 | 2.7 | -3.6 % | marginal (div 0.8 at 400) |
| B_l2c | 3.8 | 2.6 | -6.0 % (saturating at -11 %) | yes |

- Mechanism of the A-route drift: the gradient fraction of a_n grows
  monotonically for every A operator (A_int: 0.08 at n = 100, 0.95 at 400)
  because no transfer maps discrete gradients to discrete gradients, and
  nothing damps the gauge part; it eventually feeds back into B through the
  non-commutation. Divergence stays at roundoff throughout, i.e. "exactly
  divergence-free but wrong" is a real failure mode.
- Coulomb projection after each transfer removes the gauge part exactly
  (B unchanged to 1e-15, helicity change 4e-14, 14 to 22 AMG-CG iterations,
  bounded over 400 remaps) and costs 3 to 10 % of a B_l2c remap. A local
  Chebyshev-Jacobi smoother with 20 sweeps matched the global solve to six
  digits at n = 100, so gauge control does not need a global solve. But gauge
  control is not sufficient: A_pt still diverges and A_int still grows
  slowly; only A_l2 is stable.
- B_l2c is the only operator stable in every run with roundoff divergence,
  but it is a contraction: at p = 2 it loses 6 % energy per 100 remaps
  initially (saturating), at p = 3 only 0.03 %.
- The mean field: the A route keeps the net flux to 1e-15 over any number of
  remaps; B_l2c drifts by 6e-4 per 100 remaps at p = 2.
- Verified identity: argmin_a ||C_h a - b_target||_M equals the B_l2c output
  to 2.5e-11 in the M-norm (curl-curl solve with AMS vs Schur-complement
  clean). The constrained projection of the brief is B_l2c; the AMS route is
  5 to 10 times more expensive to solve.
- Cost per remap (p = 2 / 3, contended): A_int + Coulomb 0.6 / 2.6 s, A_l2 +
  Coulomb 0.9 / 2.6 s, B_l2c 2.3 / 4.5 s.

## 6. Smooth ALE cycle with rezones (Experiment 5)

Material motion x = X + eps s(t) f1(X), s = sin^2(pi t/T), 100 steps, frozen-in
field (ideal MHD): the dofs of A and B are exactly constant while the nodes
move, and the exact field is the Cauchy push-forward. Rezones to a fresh
uniform mesh at R equally spaced times including t = T.

- Discrete geometric conservation: with no rezone, a and b are constant,
  helicity changes by 3e-16, slice flux by 7e-16, scaled divergence stays
  <= 3e-15, and the errors at t = T equal those at t = 0 exactly. Mesh motion
  alone costs nothing in the dofs; the peak representation error at full
  deformation is +3 % (eps = 0.3) and +22 % (eps = 0.6) at p = 2. At p = 3 a
  Q2 geometry is the limiting error (3.6 to 8.7x); Q3 geometry restores
  +3 % / +14 %.
- Final B error at T (relative to the exact field; t = 0 value 2.28e-2 at
  p = 2 and 1.52e-3 at p = 3), R = 5 / R = 20 rezones:

| operator | p=2, eps=0.3 | p=2, eps=0.6 | p=3, eps=0.3 | p=3, eps=0.6 |
|---|---|---|---|---|
| A_int + Coulomb | 2.70e-2 / 2.71e-2 | 3.46e-2 / 3.75e-2 | 3.36e-3 / 2.79e-3 | 1.74e-2 / 7.39e-3 |
| A_l2 + Coulomb | 2.56e-2 / 2.38e-2 | 3.76e-2 / 3.44e-2 | 5.84e-3 / 3.09e-3 | 2.68e-2 / 1.08e-2 |
| B_l2 | 2.66e-2 / 2.64e-2 | 3.10e-2 / 3.40e-2 | 2.87e-3 / 2.66e-3 | 1.36e-2 / 6.50e-3 |
| B_l2c | 2.66e-2 / 2.63e-2 | 3.07e-2 / 3.31e-2 | 2.82e-3 / 2.64e-3 | 1.34e-2 / 5.97e-3 |

- The error right after a rezone is set by the projection floor of the new
  mesh for the deformed field; B_l2 and B_l2c land at 1.00 to 1.04x that
  floor, A_int + Coulomb at 1.02 to 1.08x, A_l2 + Coulomb at 1.3 to 1.5x.
  Information lost there is not recovered when the mesh moves back.
- The error does not grow linearly with the number of rezones: the mean jump
  per rezone falls like 1/R, and the final error is flat or decreasing in R.
  In this realistic regime the gauge-fixed A route is within 1.5 to 3 % of
  B_l2c at p = 2, eps = 0.3, and 6 to 46 % worse at p = 3 or eps = 0.6; a
  Coulomb projection on the source mesh before the transfer closes the gap to
  7 to 12 %.
- B_l2 and B_l2c have the same B error, so cleaning costs only time (2 to 3x)
  and buys divergence and (partly) flux. The A route gives div 1e-15 and exact
  slice flux for free; the B route's flux drift at T with a mean field is
  1e-6 (B_l2c), 1.6e-4 (B_l2), 8.8e-4 (B_int).
- Energy ratio stays within 1 % for every operator; helicity changes by at
  most 0.6 % in the A route.
- Structural finding: a frozen-in constant mean field deforms into a field
  whose fluctuating part is 22 % (eps 0.3) to 48 % (eps 0.6) of |B0| at full
  deformation. A plain transfer of a loses it, because a represents only the
  zero-mean part on the original mesh. A correct A-route rezone with B0 != 0
  needs a curl inversion of the deformed harmonic part (global solve) or a
  direct B transfer of that part. For cosmological boxes with zero mean field
  this does not arise; with an imposed mean field it does.

## 7. MHD-ALE integration (Experiment 6)

MHD-ALE already evolves A in ND and builds B as its discrete curl (its
"existing representation" is the A representation). Its remap is a continuous
ALE remap with GSLIB point interpolation of A (`-rma 0`) or a pseudo-time
constrained-transport evolution dA/dtau = -w x B (`-rma 2`). We instrumented
both and added a shadow direct-B RT interpolation at every remap event.
3-D Taylor-Green, orders (2,1,1,3), refinement 1, 45k ND dofs, forced remap
every 10 or 2 steps:

| run | remaps | final B error | energy jump per remap, A route / shadow B | max div, A / shadow B | A-remap share of wall |
|---|---|---|---|---|---|
| Lagrangian, t = 0.5 | 0 | 2.0e-2 | - | 4e-14 / - | - |
| every 10, `-rma 0` | 3 | 7.8e-3 | 2.9e-6 / 1.5e-4 | 5e-14 / 6.9e-2 | 0.2 % |
| every 10, `-rma 2` | 3 | 8.1e-3 | 1.0e-6 / 1.7e-4 | 5e-14 / 9.7e-2 | 72 % |
| every 2, `-rma 0` | 14 | 6.6e-3 | 3.4e-6 / 7.0e-5 | 5e-14 / 1.9e-2 | 1.2 % |

- Both A remaps keep div B at 5e-14 under forced remapping. Point
  interpolation of A is as accurate as the helicity-preserving remap at 1/200
  of the remap cost; the latter's smaller energy jumps are invisible next to
  the Lagrangian-step drift. Note that MHD-ALE's remap moves nodes by at most
  1.65 h between remaps; this is the small-displacement regime where even
  point interpolation is close to the identity, unlike the harsh ping-pong of
  Experiment 4.
- The shadow direct-B transfer has divergence at the 1e-2 level, energy jumps
  12 to 50 times larger and of one sign at low remap frequency, i.e. a
  secular energy gain if it were used.
- Taylor-Green has zero helicity and a steady exact solution, so it tests
  divergence and drift, not helicity. Remap overhead is dominated by the
  density, velocity and energy remaps, not by A.
- Total wall time with `-rma 2` is 1.7 to 2 times that with `-rma 0`.

## 8. Performance and scalability (Experiment 7)

Pending: see `results/exp7_summary.md` when available.

## 9. Assessment against the success criteria and recommendation

| criterion | A route | B route (B_l2c) |
|---|---|---|
| div_h B at roundoff | yes, for any transfer, any mesh | yes, but only after a global solve |
| O(h^p) convergence of B | yes | yes |
| robust on distorted meshes | yes (static); transfers degrade with p and eps | same |
| low error across one rezone | 1.3 to 3x floor (A_int), 9x (A_pt at p=4) | 1.15 to 1.4x floor |
| no secular loss under many remaps | gauge growth unless projected; A_l2 + Coulomb stable; A_int, A_pt not | stable; energy contraction (-6 %/100 at p=2, -0.03 % at p=3) |
| competitive with direct B | no at p >= 3; equal at p = 2 for A_l2 | - |
| cost and locality | local transfer + local gauge smoothing possible; L2 projection needs a mass solve | mass solve + Schur-complement solve (global) |
| exact mean flux | yes, 1e-15 | no, drifts |

Recommendation among the four options of the brief:

1. **Use A as the primary magnetic state: not recommended for an AREPO-style
   code.** It adds a gauge-dependent variable whose gauge part grows under
   remapping, needs gauge control, and does not improve B accuracy.
   (MHD-ALE's own choice to evolve A is a different matter: its induction
   update is built on it and its remaps move nodes only slightly.)
2. **Use A only during remapping: not recommended.** Converting B to A on the
   source mesh is itself a curl inversion (a global curl-curl or
   Schur-complement solve) and the resulting transfer is no more accurate
   than B_l2c, which is the same solve done once on the target.
3. **Retain a B representation with a compatible divergence-cleaning
   projection: recommended.** The combination "L2 projection of B into
   RT_{p-1} on the new mesh, then M-orthogonal projection onto ker D_h" is
   the most accurate and the only unconditionally stable transfer measured,
   gives roundoff divergence, and is exactly the constrained least-squares A
   projection of the brief. Its costs are one well-conditioned RT mass solve
   (bounded Jacobi-CG iterations, matrix-free) and one Poisson-like
   Schur-complement solve with AMG (bounded outer iterations). Its two
   weaknesses, energy contraction at low order and mean-flux drift, are
   addressed by using p >= 3 and by carrying the harmonic mean field B0
   separately (as the A route does), which costs nothing.
4. **Hybrid constrained projection: identical to 3** on a periodic domain;
   do not implement it through A.

Caveats that bound these conclusions: one mesh family (smooth periodic
deformations of a Cartesian hex box), smooth analytic fields, serial or
few-rank runs, N <= 32, and remap sequences that are harsher (ping-pong
between very different meshes) than a real ALE code's. The small-displacement
regime of MHD-ALE shows all operators close to the identity; there the
cheapest local A transfer was adequate for 14 remaps.
