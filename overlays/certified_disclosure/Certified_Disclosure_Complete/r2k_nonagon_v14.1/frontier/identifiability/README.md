# What one boundary response cannot identify

The known-chain reconstruction has a sharp information boundary: **even perfect
one-port data at every frequency need not determine an unknown mechanical
topology**. The obstruction below uses positive unit masses, positive springs,
positive grounding, and no invisible or canceled mode. It is an explicit use of
realization symmetry, not a claim that nonidentifiability itself is new.

## An exact chain and triangle with the same transfer function

Use three unit masses and the stiffness matrices

\[
K=\begin{pmatrix}4&-1&0\\-1&3&-1\\0&-1&2\end{pmatrix},\qquad
Q=\begin{pmatrix}1&0&0\\0&12/13&5/13\\0&-5/13&12/13\end{pmatrix},
\]

\[
\widetilde K=Q^TKQ=
\begin{pmatrix}
4&-12/13&-5/13\\
-12/13&602/169&-59/169\\
-5/13&-59/169&243/169
\end{pmatrix}.
\]

For a spring network, the spring between masses i and j has stiffness
`-K_ij`; the spring from mass i to ground has stiffness equal to row i's sum.
Thus these are actual positive-spring networks:

| Property | Original network | Transformed network |
|---|---|---|
| Inter-mass edges | 1–2: 1; 2–3: 1 | 1–2: 12/13; 1–3: 5/13; 2–3: 59/169 |
| Grounding springs | 3, 1, 1 | 35/13, 387/169, 119/169 |
| Masses | 1, 1, 1 | 1, 1, 1 |
| Interconnection graph | Path | Triangle |

The graphs have different edge counts and cannot be made identical by
relabeling. Every listed spring is strictly positive. The associated energy
form is

\[
x^TKx=\sum_i g_i x_i^2+\sum_{i<j}k_{ij}(x_i-x_j)^2>0
\quad(x\ne0),
\]

so both stiffness matrices are positive definite.

Drive and observe mass 1. With `z=omega^2`, their exact compliance is

\[
H(z)=e_1^T(K-zI)^{-1}e_1
=e_1^T(\widetilde K-zI)^{-1}e_1
=\frac{z^2-5z+5}{-z^3+9z^2-24z+18}.
\]

This is an identity of rational functions, including their meromorphic
continuations; it is not a fit on selected frequencies. The denominator is
`(3-z)(z^2-6z+6)`, with positive, distinct poles `3-sqrt(3), 3, 3+sqrt(3)`.
There is no pole cancellation: the numerator and denominator have polynomial
gcd one. Independently, the Krylov matrix `[e1, K e1, K^2 e1]` has determinant
**-1**, and the transformed system has the same determinant. All three modes
are visible and controllable from port 1 in the spectral variable; the
corresponding second-order mechanical state has six first-order coordinates.

The same obstruction survives damping `C=alpha I`, for any common `alpha>=0`:
substitute `z=-s^2-alpha s` into the rational identity. Both are passive
mechanical systems. Their stored energy satisfies
`dE/dt = input_force * velocity_at_port - alpha * ||velocity||^2`.
The force-to-displacement compliance is not itself being claimed to be a
positive-real impedance.

## Why the equality is unavoidable

**Port-preserving orthogonal similarity.** Let `M=I`, let K be symmetric, and
let the columns of B specify the observed/driven ports. If Q is orthogonal and
`QB=B`, then

\[
B^T(Q^TKQ-zI)^{-1}B=B^T(K-zI)^{-1}B.
\]

Proof: `(Q^TKQ-zI)^(-1)=Q^T(K-zI)^(-1)Q`, and the factors Q drop out next to B.
This proves equality for all z outside the poles, hence as rational functions.
Orthogonal similarity also preserves eigenvalues and modal minimality. A
general hidden rotation can violate physical spring signs; positivity must be
verified separately. It is verified explicitly here.

The construction extends beyond three masses. Start with an n-node path with
uniform edge stiffness `a>0`, a spring `b>0` anchoring mass 1, unit masses, and
grounding `g>a/2` at every mass. Rotate only coordinates n-1 and n by the same
12/13–5/13 block. The last two diagonal entries become
`g+433a/169` and `g+74a/169`; their mutual spring is `59a/169`. Their previous
neighbor connects to both with springs `12a/13` and `5a/13`.

The final grounding springs become `g+218a/169` and `g-50a/169`; their previous
neighbor loses `4a/13` of its original grounding. All remain positive under the
stated bound. The path acquires a terminal triangle, changing its edge count
from n-1 to n. The rotation fixes **all first n-2 coordinates**, so their entire
multiport transfer matrix agrees, even if every one of those ports can both
drive and sense. The original path's Krylov determinant is
`(-a)^(n(n-1)/2)`, which is nonzero; similarity preserves minimality.

Consequently, taking more precise samples or more frequencies at the same
ports cannot resolve this ambiguity. Small changes of the hidden rotation
around the explicit physical example also remain physical, giving a family
of indistinguishable internal spring values.

## An added sensor separates the constructed pair

Keep driving mass 1 and also observe mass 2 in the three-mass example. At zero
frequency, the two cross-compliances are

\[
e_2^TK^{-1}e_1=\frac19,\qquad
e_2^T\widetilde K^{-1}e_1=\frac{19}{234},
\]

whose difference is exactly **7/234**. One additional static response therefore
distinguishes this pair.

For the whole path/triangle family, observe mass n-1. If
`x=K^(-1)e1`, the new response at that coordinate is
`(12/13)x_(n-1) - (5/13)x_n`, strictly smaller than `x_(n-1)` because both
coordinates are positive. Positivity follows from the connected grounded
spring network's strictly positive inverse. Every checked pair is separated
this way. This statement does **not** claim that two ports identify every
unknown network.

## A constructive remedy that does recover unknown stiffness

Assume the mass matrix M is known, coordinate labels and force calibration are
fixed, and the excitation at mass 1 is cyclic: all n mass-normalized modes are
visible. Measure the **full displacement vector**, rather than a single
compliance, at n distinct nonresonant shifts `z_j`. Set

\[
x_j=(K-z_jM)^{-1}e_1,\quad
X=[x_1\ \cdots\ x_n],\quad Z=\operatorname{diag}(z_j).
\]

Then the stiffness matrix, and therefore its spring graph, is recovered by

\[
\boxed{K=(e_1\mathbf1^T+MXZ)X^{-1}.}
\]

Proof: each observed column satisfies `K x_j=e1+z_j M x_j`; combine the columns
and right-multiply by `X^(-1)`. Under cyclicity, X is invertible: diagonalize
the mass-normalized symmetric matrix, factor X into the invertible modal
basis, a diagonal matrix of nonzero port amplitudes, and the Cauchy matrix
`[1/(lambda_i-z_j)]`. Cyclicity implies distinct eigenvalues; distinct shifts
away from the poles make this Cauchy matrix nonsingular.

This is a global known-M reconstruction under the stated assumptions, with
more sensing information than the original one-port problem. It does not
assert recovery of unknown masses, invisible modes, or unknown force scaling.
The implementation `reconstruct_known_mass` receives only the shifts, observed
full response columns, and known mass diagonal. It receives no topology,
stiffness, hidden rotation, or eigenvalues.

## Exact replay and receipts

Run with standard-library Python:

```bash
python3 network_obstruction.py
python3 verify_explicit_certificate.py
```

The main suite passes **257 physical network pairs and 6,684 exact checks**:
the explicit example plus 256 seeded fixtures, 64 each at orders 3, 4, 5 and 6.
Its extra fixtures vary rational edge, anchoring, grounding, and rotation
parameters. Each fixture checks physical spring signs, positive-definite
leading minors, exact transfer numerator/denominator equality, polynomial gcd,
minimality, observed multiport blocks, separation by the extra sensor, and
exact recovery of both stiffness matrices from full-response observations.
Adversarial controls move the observed port, use a nonphysical rotation,
repeat sample shifts, and corrupt an extra response.

The explicit JSON is checked by a second, **56-check verifier** that imports no
generator code. It uses Leibniz determinants and adjugate inverses rather than
the generator's elimination and polynomial-subset routines. Its exact
polynomial checks use n+1 evaluations with verified degree bounds, thereby
establishing coefficient identities for the stated determinants.

Files:

- `network_obstruction.py`: generator, exact challenge, and sensing remedy.
- `explicit_certificate.json`: rational matrices, spring values, shared
  transfer polynomials, minimality determinants and separating measurement.
- `verify_explicit_certificate.py`: independent explicit-certificate checker.
- `results.json`: seed, all check counts, fixture/source hashes and scope.
- `independent_verification.json`: separate verifier receipt and hashes.

The recorded run took about 2.1 seconds. These finite tests support the
implementation and explicit certificates; the all-dimension claims rely on
the algebraic derivation above. No noisy-data order-selection or physical
hardware result is claimed.
