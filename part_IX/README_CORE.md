# Part IX: one-port intervention tomography

Part IX reconstructs a labeled grounded spring network from scalar responses at
one fixed force/readout port. The experiment changes the system between readings:
known grounding springs are connected to individual known sites or to pairs of
sites. The exact solver uses a finite acquisition plan, not numerical derivatives.

The same mathematics describes a grounded resistor network. Replace stiffness
matrix by nodal conductance, force by injected current, displacement by voltage,
and an added grounding spring of strength t by a known shunt conductance t. The
inverse observer consumes current-normalized voltage. That electrical mapping is
a direct nodal-equation equivalence. Physical instrument calibration and error
bounds must be supplied by an acquisition adapter.

## Exact acquisition contract

There are n known, labeled sites, including port 1. The stiffness matrix K is
symmetric positive definite, has nonpositive off-diagonals and nonnegative row
sums, and its inter-node spring graph is connected. At least one grounding spring
is positive. Force and displacement units are calibrated, and every intervention
adds a known nonnegative diagonal stiffness without changing another coefficient.
Measurements are the exact unit-force static compliance at port 1.

The implemented exact artifact format restricts matrix coefficients, load
amplitudes and responses to rational numbers. The mathematical identities extend
to real coefficients. Floating-point preview exports and bounded-error physical
readings require their own observation contracts; they are not silently promoted
to exact certificates.

Put G=K^{-1}, h0=G11, aj=G1j and bj=Gjj. Connectivity and grounded spring signs
imply every entry of G is strictly positive. Choose any calibrated t>0. Acquire:

| Configuration | Number | Scalar response |
| --- | ---: | --- |
| No added load | 1 | h0 |
| Load t at each hidden site j | n-1 | hj(t) |
| Load 2t at each hidden site j | n-1 | hj(2t) |
| Load t simultaneously at each hidden pair j,k | (n-1)(n-2)/2 | hjk(t) |

There are n(n+1)/2 scalar configurations. At most two added loads are connected
simultaneously, so two programmable load channels with a site-routing matrix can
implement this schedule. Independent actuation at each labeled site is the access
requirement that replaces full-field sensing.

## Finite exact reconstruction

The rank-one inverse identity gives

    hj(t) = h0 - t aj^2 / (1+t bj).

Consequently, t/(h0-hj(t)) is affine in t. The readings at t and 2t determine its
positive intercept 1/aj^2 and positive slope bj/aj^2. Positivity fixes aj as the
positive square root. This recovers the first row and the diagonal of G.

For a hidden pair j,k write c=Gjk and

    delta = h0-hjk(t),
    A = aj^2 (1+t bk) + ak^2 (1+t bj),
    B = (1+t bj)(1+t bk).

The rank-two inverse identity is

    delta = [t A - 2 t^2 aj ak c] / [B-t^2 c^2].

The physical branch is the smaller quadratic root:

    c = [aj ak - sqrt((aj ak)^2 - delta A/t + delta^2 B/t^2)] / delta.

Branch selection follows from the model class, not an empirical sign heuristic.
After adding the two positive grounding springs, inverse positivity gives

    (1+t bk) aj - t c ak > 0,
    (1+t bj) ak - t c aj > 0.

Their product, divided by B-t^2 c^2, equals aj ak-delta c. Thus the physical root
has delta c<aj ak and is the smaller root. Every hidden off-diagonal is recovered
by one pair reading. Invert the resulting G to obtain K.

This proves uniqueness within the declared labeled network class. The scalar
configuration count equals the number of independent entries in a symmetric K.
That parameter count is an exact accounting statement. Turning it into a formal
minimality theorem for an arbitrary adaptive acquisition model additionally
requires stating its continuity and protocol assumptions.

The implementation verifies exact positive-definite LDL factors, spring signs,
grounding, connectedness, every pair root margin, the entire acquisition set and
all supplied held-out configurations. It rejects inconsistent data rather than
substituting a nearest graph. A rational square-root failure also rejects the
declared rational coefficient class, even when a real-valued model might exist.

## Observer and certificate interfaces

`intervention.reconstruct(document)` receives no hidden matrix. Its sole input is
an observation document with the known order, coordinate labels, force/readout
calibration, base load and scalar responses. `observation_document` is a synthetic
forward-acquisition helper, kept outside the observer interface.

The exact schemas are:

* `nonagon.intervention.observations.v1`
* `nonagon.intervention.certificate.v1`

Node indices in serialized load mappings are zero based. Port 1 in the mathematics
is serialized as index 0. Configuration identifiers are `baseline`,
`single_j_a`, `single_j_b`, and `pair_j_k`. A held-out configuration must differ
physically from every acquisition configuration and every other held-out load.

Each certificate contains the original scalar observations, recovered K and G,
spring and grounding values, extraction witnesses and a displacement vector for
every acquisition and held-out configuration. The latter are independently
checkable witnesses of

    (K + diag(loads)) u = e1,     u1 = recorded response.

The separate verifier is written without importing the reconstruction code. Its
own exact linear algebra checks the witness equations and physical matrix class.
Reported validation flags are metadata, not substitutes for those checks.

Run the complete synthetic campaign:

```bash
python3 campaign.py --output data
```

Reconstruct one scalar observation file:

```bash
python3 intervention.py data/twin_path.observations.json --output accepted.json
```

## Completed exact campaign

The frozen seed is 2026103008. The campaign contains 100 fresh rational connected
networks of orders 2 through 10, plus the retained path/triangle twins and one
damaged path. All 103 accepted certificates recover their hidden synthetic
stiffness matrices exactly. Their acquisition budget is 2,430 scalar readings,
followed by 309 distinct held-out readings. The core rejects 24 altered or
unsupported input controls, including calibration changes, missing acquisitions,
incorrect load strengths, malformed exact numbers and changed held-out data.
The authoritative machine-readable totals and certificate index are in
`data/campaign.json`.

The retained v7 twins have the identical unperturbed scalar transfer but different
physical graphs. Their common DC compliance is 5/18. After adding a unit grounding
spring at the same labeled node 2, the path response is 7/26 and the triangle
response is 1088/3989. The absolute response separation is 365/103714. The same physical
intervention also separates the retained nonproportionally damped twins at every
one of five tested rational real-Laplace settings. The exact similarity identities
and samples are retained in `data/damped_twins.json`.

The change experiment weakens only the path spring between nodes 1 and 2 from 1
to 1/2. Two independently acquired scalar datasets reconstruct the before/after
matrices; coefficient comparison identifies that one changed edge and no changed
grounding. The result is in `data/fault_localization.json`.

The inverse-design helper chooses the finite added grounding spring that achieves
a requested static port compliance. For a known design site j the attainable
range is

    h0-aj^2/bj < target <= h0.

For delta=h0-target, the required added stiffness is

    load = delta/(aj^2-delta bj).

The campaign verifies this design against a fresh exact forward calculation;
`data/inverse_design.json` retains the target, added load and response vector.
This is a static compliance design task, distinct from modal optimization.

## Scope and practical next step

The experiment makes accessible nodes identifiable by perturbing them. It does
not expose an arbitrary inaccessible interior. Weak port-to-site coupling can make
the exact inverse numerically sensitive; finite precision, calibration drift and
load uncertainty therefore belong in a bounded-error acquisition workflow.

A small resistor network is a direct physical implementation: one calibrated
current source, one voltage readout, a stable reference ground, and two calibrated
switchable shunts routed among known pads. The same six-configuration schedule
for the three-node fixture provides a compact first bench test, followed by
unused shunt settings for held-out prediction.
