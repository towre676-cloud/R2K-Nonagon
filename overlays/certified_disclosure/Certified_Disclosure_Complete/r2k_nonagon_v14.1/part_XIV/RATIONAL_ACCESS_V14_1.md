# Rational access certificates, Nonagon v14.1

The access question is now answered by exact rational certificates rather than
by an unsupported interpretation of a deficient finite-field nullspace.

| Experiment | Readings | Rational Jacobian rank | Nullity | First-order identifiable branch coordinates |
|---|---:|---:|---:|---:|
| IEEE 118, 20 farthest-point ports; interventions at shunt/generator buses | 13,020 | 127/179 | 52 | 74/173 |
| IEEE 118, 20-port noise-experiment design; interventions at every hidden bus | 24,780 | 178/179 | 1 | 167/173 |
| IEEE 118, six placed ports; driving-point readings after interventions at every hidden bus | 723 | 179/179 | 0 | 173/173 |
| IEEE 300, retained seven-port design; driving-point readings after interventions at every hidden bus | 2,121 | 409/409 | 0 | 408/408 |

The six IEEE 118 ports are buses **49, 87, 72, 1, 108, 55**. This replaces the
retained loop's seven-port, 847-reading completion by a six-port, 723-reading
completion: 124 fewer readings, a 14.64% reduction. It proves feasibility of that
placement, not minimum port count. Starting with the same first five ports,
adding bus 54 or 56 also gives full rank; adding bus 59 leaves rank 178/179.

## The corrected inference

Full column rank modulo a prime proves full rank over the rationals, provided
the rational entries reduce at that prime. This remains a valid local
identifiability certificate. Deficient modular rank does not prove rational rank
deficiency, and modular kernel support need not match rational kernel support.
For example, with p=2^31−1,

    J = [1, p],       ker_Q(J) = span{(-p, 1)},       J mod p = [1, 0].

Both matrices have rank one. Neither rational coordinate is individually
identifiable to first order, but the modular nullspace incorrectly suggests
that the first coordinate is. The old “undercount, never overcount” claim was
therefore false even without a rank drop. `ident14.py` now calls deficient
coordinate counts modular diagnostics. Historical `ident14_*.json` files carry
an explicit correction; their numeric fields are preserved as legacy data.

The new exact certificate confirms the original rounded IEEE 118 43% figure:
74/173 branch coordinates, and 78/179 total coordinates, are observable to first
order with the stated 20-port shunt/generator-access design. This agreement is
a result of the new rational computation, not a consequence of modular support.

## Exact certificate construction

Let J be the Jacobian of the declared rational reading map at the source model.
For each transfer-compliance reading,

    d y / d w_(a,b) = −(x_p[a]−x_p[b])(x_q[a]−x_q[b]),
    d y / d g_i     = −x_p[i] x_q[i],
    x_p = (K + Θ)^−1 e_p.

The implementation reconstructs K and G=K^−1 exactly with FLINT rationals and
checks K G=I. If G=M/d, an intervention of strength one at bus j has response
column numerator

    X_p[i] = M[i,p](d+M[j,j]) − M[i,j]M[j,p].

The omitted denominator d(d+M[j,j]) is nonzero and common to the column.
Primitive integer response columns therefore give Jacobian rows differing from
the physical rational rows only by nonzero row scalings. Rank, nullspace and
row-space coordinate tests are preserved. The reading model and intervention
strength are explicitly fixed; none of these steps uses noisy data.

A finite-field transpose RREF proposes r independent reading rows. A square
minor of those **actual reconstructed integer rows** has nonzero determinant
modulo p, proving rank_Q J≥r. Exact rational RREF of this selected subsystem
proposes U−r integer kernel generators. The programme checks Jv=0 for every
retained generator v against **every** reading, including all omitted rows.
Their free-coordinate submatrix is diagonal and nonzero, so they are independent.
Hence nullity_Q J≥U−r and rank_Q J≤r. Together these bounds certify the entire
rational kernel and rank.

A coordinate i is first-order identifiable precisely when every vector in this
complete kernel has zero coordinate i. In deficient-rank experiments, explicit
rational coefficients λ are also retained and checked to satisfy

    λ^T A = e_i^T,

where A is the selected actual integer-row subsystem. A nonzero kernel coordinate
provides the converse witness: a parameter perturbation invisible to every
linearized observation. For full rank, the verified nonzero square minor alone
certifies all coordinates, so individual inverse-row arrays are unnecessary.

Certificates retain source reading indices, selected-row and full-Jacobian
hashes, the rank minor's modular determinant, exact kernel vectors and any
coordinate row-span coefficients. Large selected integer rows are reconstructed
from the inputs during replay instead of being duplicated in JSON.

## Scope

These are exact statements about the **linearized reading map at the declared
rational model**. Full Jacobian column rank implies local identifiability there.
A deficient Jacobian, by itself, does not prove nonlinear nonidentifiability;
for example, x↦x³ is injective despite its zero derivative at zero. No deficient
coordinate count here is presented as a global uniqueness or impossibility
result. Structural harmonic-pocket theorems are separate results.

These access certificates also do not bound reconstruction error under noise.
The previous IEEE 118 20-port Gaussian-noise experiment still contains a
one-dimensional exact Jacobian kernel supported on six branches and starts its multiport estimator from a
source-derived nominal model with 10% perturbations. Its 0.366% all-branch median
error is a nominal-model-conditioned result. The new six-port full-rank design
needs its own noise/conditioning assessment before an accuracy claim is made.
Cramér–Rao values in Part XIV are analytic-Jacobian computations in binary64.

## Reproduction

From `part_XIV/`, using Python, python-flint and NumPy:

    python rational_identifiability.py case118 hw20,ports20,pocket6
    python rational_identifiability.py case300 access7
    python rational_access_replay.py
    python -m pytest -q test_rational_identifiability.py

`rational_access/replay_receipt.json` records the completed independent replay.
Replay regenerates the exact inverse and all reading rows; validates the
selected indices and their uniqueness/range; checks the matching minor; checks
every kernel dot product; checks kernel independence; and verifies the explicit
coordinate row-span witnesses. Source-derived model construction is shared,
while the retained algebraic certificates are verified without rerunning their
rational RREF search. Run without `python -O`; optimized mode is explicitly
rejected because assertion checks are part of the verifier.

The regression suite includes the modular-coordinate counterexample, an exact
intervention-vector comparison against direct rational inversion, a genuine
partial-rank example, a full-rank example, and rejection of tampered kernel and
unmatched selected-row certificates.
