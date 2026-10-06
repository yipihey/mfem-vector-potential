# Plan: vector-potential remapping for ALE MHD with MFEM

Goal: decide, cheaply and rigorously, whether transferring the magnetic vector
potential A (in H(curl)) across ALE rezoning events and reconstructing
B = B0 + curl_h A (in H(div)) is better than transferring B directly, for a
future AREPO + MFEM/MHD-ALE hybrid. No AREPO changes. A negative answer is a
valid outcome.

## Key fact discovered while reading MHD-ALE (changes the framing)

MHD-ALE (github.com/ruijie-xi/MHD-ALE, commit 0e0564e) already uses the
vector potential as its primary magnetic state in 3-D:

- `MHD.cpp` lines 376-405: A in `ND_FECollection(order_A)`, B in
  `RT_FECollection(order_A-1)`, with B obtained by a `ParDiscreteLinearOperator`
  + `CurlInterpolator` (`remap.cpp` 203, `tools.cpp` 131).
- Its remap options for A (`-rma`): 0 = point interpolation via GSLIB
  FindPoints (`InterpolateRemap`), 1 = DG convection (H1/L2 only),
  2 = `HPRemap` ("helicity preserving"): pseudo-time evolution
  dA/dtau = -w x B on the moving mesh, with H = M^{-1}-projection of B into ND,
  i.e. a compatible constrained-transport style remap that keeps div_h B = 0
  by construction.
- Its remap is a *continuous* ALE remap (same topology, nodes move by a mesh
  velocity), not a transfer to a new mesh.

So "MHD-ALE's existing magnetic representation" IS already A-based. The
"direct B remap" alternative therefore has to be built by us (standalone), and
the MHD-ALE integration experiment compares MHD-ALE's A-remap variants against
each other and against the standalone direct-B results. The AREPO-relevant
question becomes: for a code whose state is B (like AREPO), is it worth
introducing A (gauge-dependent) just for rezoning?

## Phases and deliverables

| Phase | Deliverable | Agent |
|---|---|---|
| 0 | Builds: GSLIB, MFEM master, MFEM fork (for MHD-ALE), MHD-ALE; reproduce ex3p/ex4p and an MHD-ALE test; `BUILD.md`, `scripts/build_all.sh` | sonnet (running) |
| 1 | `src/vp_core.{hpp,cpp}`: periodic mesh + deformation, spaces, discrete operators (G_h, C_h, D_h), analytic fields, diagnostics (energy, helicity, flux, div), transfer operators (A: ND point-interp, ND integrated, L2-projection; B: RT point-interp, RT face-flux quadrature, L2-projection, L2-projection+div cleaning). `experiments/exp1_static.cpp`, `experiments/exp2_distort.cpp` | sonnet |
| 2 | `experiments/exp3_remap.cpp` (one remap, A vs B, gauge tests), `experiments/exp4_repeat.cpp` (N remaps ping-pong), `experiments/exp5_ale_cycle.cpp` (continuous mesh motion + rezones + return) | sonnet |
| 3 | Run matrices, results CSVs in `results/`, plots via `scripts/plot_*.py` | sonnet/haiku |
| 4 | MHD-ALE instrumentation: per-remap diagnostics (div B, E_B, helicity, flux) in 3-D Taylor-Green with forced frequent remesh; compare `-rma 0` vs `-rma 2`; optional new remap mode | sonnet |
| 5 | Timing/memory/MPI scaling runs | sonnet |
| 6 | Design note + results report + recommendation (`docs/design_note.md`, `docs/results.md`) | lead |

## Discrete setting (details in docs/design_note.md)

- Periodic box [0,1]^3 (MFEM `MakePeriodic` on a Cartesian hex mesh; nodes
  live in a discontinuous L2 space, raised to geometric order q via
  `SetCurvature(q, true)`).
- ND_p (A), RT_{p-1} (B), L2_{p-1} (div B), H1_p (gauge chi). MFEM pairs
  `ND_FECollection(p)` with `RT_FECollection(p-1)` and `L2_FECollection(p-1)`.
- C_h = `DiscreteLinearOperator(ND,RT)+CurlInterpolator`,
  D_h = `DiscreteLinearOperator(RT,L2)+DivergenceInterpolator`,
  G_h = `DiscreteLinearOperator(H1,ND)+GradientInterpolator`.
  Check numerically: max|D_h C_h| and max|C_h G_h| (sparse products) and the
  action on actual vectors, on undeformed and deformed meshes.
- B_h = Pi_RT(B0) + C_h A_h; B0 constant is exactly in RT_0 subset RT_{p-1}.
  Net flux through the three topological face families (images of x_i = const
  planes) equals B0_i |cross-section| exactly for the curl part (discrete
  Stokes on a closed surface).
- Helicity H = int A_h . B_h (only reported when B0 = 0, or reported as the
  fluctuating helicity with its gauge caveat documented).

## Test fields (all 1-periodic, k = 2 pi)

1. ABC/Beltrami: A = (a sin kz + c cos ky, b sin kx + a cos kz, c sin ky + b cos kx),
   curl A = k A. Nonzero helicity, fully 3-D.
2. "Modulated" non-Beltrami: A = eps_f (sin ky cos kz, sin kz cos kx, sin kx cos ky)
   (+ its analytic curl), with a multi-mode variant (adds a k=4pi mode) to
   provide under-resolved scales.
3. Gauge perturbation: chi = g sin kx sin ky sin kz, A -> A + grad chi
   (B unchanged). Also a rough gauge: chi with k = 6 pi.
4. Mean field B0 = (0.3, -0.2, 0.5) added to selected runs.

## Mesh deformation (topology fixed, periodic displacement)

x' = x + eps f(x), f = (sin ky sin kz, sin kz sin kx, sin kx sin ky) / k and
a higher-frequency variant f2 with k2 = 4 pi. Scan eps so that min det J
goes from 1 down to ~0.1. Record min/max det J at quadrature points, and
the scaled-Jacobian / aspect metrics from the mesh-quality utilities.

## Experiments

E1 static: p = 1..4, N = 4, 8, 16 (and 32 for p <= 2); errors ||B_h - B||_L2,
||A_h - A||_L2, ||D_h B_h||_L2 and max-dof, ||D_h C_h||_max, time of C_h apply,
DOF counts, condition estimates of M_ND (for later L2 projections) via CG
iteration counts.

E2 distortion: same quantities + energy error and Jacobian stats vs eps,
for p = 1..4 at N = 8 (and 16 for p <= 2), geometric order q in {1, 2, p}.
Distinguish "div at roundoff regardless of eps" from B approximation error
growth.

E3 one remap: M0 uniform -> M1 = deformed(eps) -> M2 in {M0, deformed with
different f/eps, finer uniform}. Transfer ops:
  A-route: T_A1 ND point interpolation, T_A2 ND integrated (edge quadrature)
  interpolation, T_A3 global L2 projection into ND (CG, PA mass matrix).
  B-route: T_B1 RT point interpolation, T_B2 RT face-flux quadrature,
  T_B3 global L2 projection into RT, T_B4 = T_B3 + divergence cleaning
  (solve for phi: B <- B - grad_h phi with weak gradient; global solve).
Metrics: ||D_h B||, ||B - B_exact||, relative energy change, flux error,
helicity change, wall time, DOF counts, scaling with p, sensitivity to eps.
Gauge: repeat A-route with A + grad chi and report ||C_h T(G_h chi)|| /
||C_h T A||.

E4 repeated: M_A <-> M_B ping-pong, n = 100 (and 400 for one config), p = 1..3,
track E_B(n), ||B_n - B_0||, ||D_h B_n||, H(n), ||A_n||, with all transfer ops.

E5 smooth ALE cycle: mesh moves x(t) = x + eps s(t) f(x), s(t) = sin^2(pi t/T),
over t in [0,T]. Frozen-in (Lagrangian) field: DOFs of A (edge circulations)
and B (face fluxes) are exactly constant under the motion, so the analytic
pushed-forward field is known: B(x',t) = F B(x) / det F with F = dx'/dx.
Measure FE representation error vs t with no remap (should return exactly to
the initial state at t = T, checking geometric conservation), then insert K
rezones (to the uniform mesh and back onto the current moving mesh) during the
cycle and attribute the final error to the remaps.

E6 MHD-ALE: 3-D Taylor-Green (`testcase_TaylorGreen.cpp`, `-dim 3`), forced
remesh every few steps (`-fsr -fsri`), compare `-rma 0` vs `-rma 2`; add
per-remap diagnostics; measure E_B(t), E_K(t), ||D_h B||, helicity, spectra
if feasible, overhead. Optionally add a direct-B remap mode (RT interpolation
followed by A recovery) to compare inside MHD-ALE.

E7 cost: DOFs, memory (VmRSS), wall times of FindPoints, projection, curl,
with 1/2/4 MPI ranks and p = 1..4.

## Success criteria (from the brief)

div at roundoff; O(h^p) convergence of B; robust on distorted meshes; low
error across rezoning; no secular energy/helicity loss over many remaps;
competitive with direct B remap; acceptable cost and locality.

## Risks / things that may sink the idea (to be tested, not assumed)

- Point-interpolatory ND transfer does not commute with the gradient, so gauge
  components of A pollute B after transfer; repeated remaps may build up a
  large gauge component (||A|| drift) even if B stays bounded.
- L2 projection of A needs a global ND mass solve (well conditioned, but
  global); L2 projection of B needs an RT mass solve; div cleaning needs an
  elliptic solve. Record all of these costs explicitly.
- Direct B transfer by face-flux quadrature is "almost" div-free (only the
  quadrature error across source-element interfaces breaks it). If that error
  is tiny, the A-route's exact divergence may buy little.
- A-route may lose energy faster than B-route (curl of an interpolant is a
  lower-order approximation of B than direct interpolation of B).
