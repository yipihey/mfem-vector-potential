# Design note: a compatible vector-potential representation of B for ALE remapping

Status: design (results are reported separately in `docs/results.md`).

## 1. Setting

Periodic box Omega = [0,1]^3 (3-torus), meshed with N^3 hexahedra of
geometric order q (nodes stored in a discontinuous L2 space so that the mesh
can be periodic and curved). Mesh deformations keep the topology fixed:

    x' = Phi(x) = x + eps f(x),   f 1-periodic,   det(I + eps grad f) > 0.

## 2. Discrete de Rham complex on hexahedra

MFEM provides the exact sequence (all spaces are conforming, periodic):

    H1_p  --G_h-->  ND_p  --C_h-->  RT_{p-1}  --D_h-->  L2_{p-1}

- H1_p: continuous Q_p scalars.
- ND_p: Nedelec (first kind) edge elements of order p (`ND_FECollection(p,3)`),
  tangentially continuous, dim per element 3 p (p+1)^2.
- RT_{p-1}: Raviart-Thomas face elements (`RT_FECollection(p-1,3)`),
  normally continuous, dim per element 3 p^2 (p+1).
- L2_{p-1}: discontinuous Q_{p-1}.

G_h, C_h, D_h are the *discrete differential operators*: sparse matrices that
map coefficient vectors to coefficient vectors and are exact, i.e.
curl(sum_j a_j phi_j) = sum_i (C_h a)_i psi_i holds as functions, because
curl(ND_p) is a subspace of RT_{p-1} and the RT degrees of freedom of a
function in curl(ND_p) are computed exactly. In MFEM they are assembled by
`DiscreteLinearOperator` with `GradientInterpolator`, `CurlInterpolator`,
`DivergenceInterpolator`. Their entries are rational (for affine elements,
integers) and independent of the mesh geometry: the operators are incidence /
"topological" operators of the mesh complex, up to the fixed reference-element
basis. Consequently

    C_h G_h = 0,    D_h C_h = 0        (exactly, entries 0 after assembly),

on any valid mesh, deformed or not. This is the whole point: div_h B_h = 0 is
a property of the *algebra*, not of the geometry or of the accuracy of A_h.

Notation: coefficient vectors are lowercase (a, b), FE functions have
subscript h.

## 3. Representation of B

    B_h = Pi_RT(B0) + C_h a,      a in ND_p coefficients.

B0 is the constant mean (harmonic) field. On the torus, the first cohomology
is R^3: a constant field is divergence free but is not the curl of any
periodic A. So B0 must be carried separately. Constants are in RT_0 subset
RT_{p-1}, so Pi_RT(B0) is exact, divergence free, and its RT face fluxes are
exactly B0 . n |F|.

Mean preservation. Let S_i be a closed "slice" surface made of mesh faces that
is the image of the plane x_i = const (its face set is fixed by the topology,
so it is well defined on deformed meshes). The net RT flux of C_h a through S_i
is the sum of face fluxes = sum over faces of the circulation of a around the
face boundary (discrete Stokes); on a closed surface every edge is shared by
two faces with opposite orientation, so the sum is exactly zero. Hence the
flux of B_h through S_i is exactly that of B0, independent of a, of the mesh,
and of how a was obtained. We verify this numerically (expected: roundoff).

Energy and helicity.

    E_B = 1/2 int |B_h|^2 = 1/2 b^T M_RT b,
    H   = int A_h . B_h   = a^T M_mix b,  (M_mix: ND x RT mixed mass matrix).

Under a periodic gauge change a -> a + G_h c (c in H1_p), B_h is unchanged and
H changes by int grad(chi_h) . B_h = - int chi_h div B_h + 0 = 0 because
div B_h = 0 *as a function* (strong, since B_h is in RT with D_h b = 0, which
means div B_h = 0 in every element and normal continuity gives no face jumps).
So H is gauge invariant under periodic gauge changes, but it includes the term
int A_h . B0 which is only well defined because A is periodic; with B0 != 0
the "total" helicity (including a non-periodic potential for B0) is not
defined on the torus. We report H only for B0 = 0 runs, or label it clearly as
the fluctuating helicity.

## 4. Projection / interpolation of A onto a mesh

Given a smooth A (analytic, or a FE function on another mesh evaluated point
wise), the natural ND interpolant Pi_ND uses the ND degrees of freedom. MFEM's
default `ProjectCoefficient` for ND is *interpolatory*: dof_k = (J t_k) . A at
the dof node x_k (tangential component at points along edges / interior
nodes). `ND_HexahedronElement::ProjectIntegrated` uses integrated edge
moments instead. Only moment-based (integrated) dofs make the diagram commute
exactly, Pi_ND grad chi = G_h Pi_H1 chi; the point-interpolatory version
commutes only up to quadrature error. Both are tested, because commutation
controls gauge sensitivity of a transfer.

Error expectations: ||A - Pi_ND A||_L2 = O(h^p), and ||curl(A - Pi_ND A)||_L2
= ||B - C_h Pi_ND A||_L2 = O(h^p) (the curl of the interpolation error is
controlled by the RT interpolation error of B, by the commuting diagram:
C_h Pi_ND A = Pi_RT curl A for the moment-based interpolant). So B converges
at O(h^p) when A is in ND_p, same rate as interpolating B directly into
RT_{p-1}. No order is lost by going through A, in exact arithmetic and with
commuting interpolants.

## 5. Transfer across a rezone (M1 -> M2)

Let A1_h live on M1 (ND_p), and M2 be a new mesh of the same domain, same
topology (same N) or not. Write E1 for the "evaluation" operator: given
points y, return A1_h(y) by locating y in M1 (GSLIB FindPoints, with periodic
wrapping of y into [0,1)^3) and evaluating the FE function. Transfers:

A-route (then b2 = Pi_RT B0 + C_h a2, exactly divergence free whatever a2):
- T_A1: a2 = Pi_ND^{pt} E1 A1_h (point-interpolatory ND dofs on M2).
- T_A2: a2 = Pi_ND^{int} E1 A1_h (edge-integrated ND dofs on M2).
- T_A3: a2 = argmin ||A2_h - E1 A1_h||_L2 over ND_p(M2): global ND mass
  solve M_ND a2 = r, with r assembled by quadrature of the source. (Global but
  well conditioned, cond(M_ND) = O(1) in h for fixed p; CG with Jacobi,
  partial assembly.)

B-route (transfer B1_h = Pi_RT B0 + C_h a1 directly):
- T_B1: b2 = Pi_RT^{pt} E1 B1_h (point values of the normal component at the
  face dof nodes). No reason for D_h b2 = 0.
- T_B2: b2 = Pi_RT^{int} E1 B1_h, face fluxes by face quadrature. Because
  B1_h is exactly divergence free as a function, the exact fluxes through the
  faces of any M2 element sum to zero; the only violation is quadrature
  error on M2 faces that cut M1 element interfaces (where B1_h is only
  piecewise polynomial). D_h b2 is therefore small but not roundoff; it
  decreases with quadrature order only algebraically.
- T_B3: global L2 projection into RT_{p-1}(M2) (RT mass solve).
- T_B4: T_B3 followed by divergence cleaning: find phi in L2_{p-1}(M2) with
  (D_h M_RT^{-1} D_h^T M_L2) phi = D_h b, b <- b - M_RT^{-1} D_h^T M_L2 phi
  (the standard projection onto the discrete divergence-free subspace; a global
  saddle-point / Schur complement solve). This is the "remap B then clean"
  baseline named in the brief.

Fair comparison requires the same quadrature order and the same point
location machinery for both routes; costs are reported per route.

## 6. Gauge

A is defined modulo G_h c. For transfer operator T and gauge change G_h c on
M1,

    C_h T (a1 + G_h c) - C_h T a1 = C_h T G_h c,

which vanishes iff T maps discrete gradients to discrete gradients
(T G_h^{M1} = G_h^{M2} T_H1 for some scalar transfer). For T_A2 with exact edge
integrals of a *continuous* gradient this holds exactly (fundamental theorem of
calculus along the edge), but E1 A1_h is only tangentially continuous and
G_h c is a FE gradient of a continuous chi_h, so the edge integral of
grad chi_h along an M2 edge equals chi_h(end) - chi_h(start) exactly, and
T_A2 G_h c = G_h^{M2} (Pi_H1 chi_h) up to quadrature error of a piecewise
polynomial along an edge crossing M1 element boundaries. So gauge pollution
of B is a quadrature-level effect for T_A2, and an O(h^p) interpolation-level
effect for T_A1. We measure both, including a rough gauge.

No gauge fixing is applied unless the experiments show a need (e.g. secular
growth of the gauge component ||a_n|| under repeated remaps). Diagnostics:
||a||, the weak divergence ||G_h^T M_ND a|| (Coulomb-gauge residual), and the
split of a into gradient and gradient-orthogonal parts when cheap enough.

## 7. Optional constrained projection (only if the simple transfer fails)

    a2 = argmin_a 1/2 ||C_h a - b_target||^2_{M_RT}  + gauge condition,

normal equations C_h^T M_RT C_h a = C_h^T M_RT b_target: the curl-curl
operator, singular on gradients. Needs a gauge (e.g. Coulomb: G_h^T M_ND a = 0,
or a minimum-norm/regularised solve) and an AMS-type preconditioner; it is a
global elliptic solve, i.e. not local. We only pursue it if T_A1-T_A3 lose
too much accuracy, and then quantify its cost.

## 8. ALE motion without remap (frozen-in field)

If the mesh moves with the fluid and the field is frozen in (ideal MHD, no
resistivity), the face fluxes of B and the edge circulations of A are exactly
conserved: b and a (in a suitable gauge, a + G_h c with c = int w.A along the
motion) are *constant in time* as coefficient vectors, while the functions
change through the geometry. The pushed-forward analytic field is
B(x',t) = F B(x) / det F, A(x',t) = F^{-T} A(x). This is what makes the
"smooth ALE cycle" test clean: on return to the original geometry, a
Lagrangian code with no remap recovers the initial field exactly; every
deviation after inserting K rezones is attributable to the K transfers.

## 9. What MHD-ALE already does

MHD-ALE's 3-D state is exactly this representation (ND A, RT B = C_h A), with
B0 = 0 for its tests. Its ALE remap is continuous (nodes move, topology fixed)
with options: point interpolation of A via FindPoints (T_A1 analogue) or the
"helicity-preserving" pseudo-time evolution dA/dtau = -w x B_h (a constrained
transport remap). It has no direct B remap (its state is A). Our standalone
B-route therefore represents what an AREPO-style B-based code would do, and
the MHD-ALE experiment measures how its two A-remaps behave under frequent
forced rezoning.

## 10. Corrections and additions after the experiments

- Section 3 claimed that constants lie in RT_{p-1} so Pi_RT(B0) is exact.
  This is true only on affine elements. Under the contravariant Piola map a
  constant physical field pulls back to adj(J) B0, which is not polynomial of
  the right degree on curved elements. The pointwise projection of B0 then
  has divergence up to O(1) (Experiment 2); the face-flux projection keeps
  D_h b0 = 0 to roundoff and the exact slice flux, at the price of a
  pointwise error of order the geometric interpolation error. Always project
  B0 with face-flux dofs.
- C_h and G_h are geometry independent bitwise. D_h as assembled by MFEM's
  DivergenceInterpolator returns L2 point values of the divergence and is
  therefore a row scaling (1/detJ at the L2 nodes) of the topological
  operator; its kernel is geometry independent, its entries are not.
- A_h converges at O(h^{p+1}) on affine meshes (ND_p contains P_p), dropping
  to O(h^p) on non-affine ones; B_h converges at O(h^p) in both cases.
- Transfers of A do not commute with G_h, even the integrated ones (only up
  to the quadrature error of the kink of grad chi_h across source-element
  faces). The gauge part of A is never damped by a transfer and grows under
  repeated remaps; a Coulomb projection (AMG-CG on the H1 stiffness matrix,
  bounded iterations) or a local Chebyshev-Jacobi smoothing removes it
  without changing B (verified to 1e-15) or the helicity (4e-14).
- On the torus the constrained projection argmin ||C_h a - b||_M is the
  M-orthogonal projection of b onto ker D_h minus its harmonic part, i.e. the
  same operator as "L2-project B then clean", verified to 2.5e-11.
- A frozen-in constant mean field is not constant after a Lagrangian
  deformation; its deviation from B0 is a zero-mean curl field that a
  transfer of a alone loses. With B0 != 0 an A-route rezone therefore needs
  a curl inversion of that part or a direct B transfer of it.
