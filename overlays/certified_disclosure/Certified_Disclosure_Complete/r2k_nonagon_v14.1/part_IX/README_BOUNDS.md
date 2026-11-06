# Bounded-error intervention tomography

`bounded_intervention.py` propagates simultaneous rational measurement intervals through the finite local-load reconstruction. It accepts one normalized response channel, known labeled sites, connected grounded nonnegative spring or conductance matrices, and independently varying calibration intervals for the added grounding loads. All calculations use `Fraction`; square roots and interval compression have exact outward dyadic bounds. No binary floating-point arithmetic supplies a verdict.

The normalized reading must be enclosed as displacement/force or voltage/current. Input-current uncertainty, sensor-gain uncertainty, repeatability and drift must already be included in that externally justified scalar interval. The optional relative load envelope `e` means that **each actual load in each acquisition**, independently, lies in `[t(1-e), t(1+e)]` around its nominal value. Pair loads need not be physically equal. Conductance uncertainty must be obtained from resistor calibration; a resistance interval `[Rlo,Rhi]` gives conductance interval `[1/Rhi,1/Rlo]`.

For unequal paired positive loads `u,v`, with `a=G[0,j]`, `b=G[0,k]`, `p=G[j,j]`, `q=G[k,k]`, `c=G[j,k]`, and measured response decrease `d`,

```
A = a^2(1+v q)/v + b^2(1+u p)/u
B = (1+u p)(1+v q)
c = (a b - sqrt((a b)^2 - d A + d^2 B/(u v))) / d.
```

This is the smaller root. Positivity of the loaded inverse fixes that branch in the stated network class. Two separated single-load calibration intervals recover each diagonal and boundary-column entry. The implementation requires positive differences and a nonnegative discriminant interval; it returns `unresolved` when these sufficient margins fail.

The compliance interval matrix has exact midpoint `Gm` and radii `E`. Let `R=Gm^-1`, `rho=||R||inf ||E||inf`. If `rho<1`, every admissible compliance is invertible and

```
||G^-1 - R||inf <= ||R||inf rho/(1-rho).
```

This produces entrywise stiffness enclosures. A spring is disclosed present only when its lower weight bound is positive. Exact absence at positive noise requires the declared gap premise: every nonzero spring is at least `gamma`. An upper bound strictly below `gamma` then discloses absence. Wide intervals remain unresolved. A positive interval entirely below `gamma` contradicts that premise.

The result is conditional on the physical class and the input enclosures; it does not prove that arbitrary uploaded data have a physically compatible realization. The exact-core certificate verifier is separate. This interval implementation currently has an exact-arithmetic containment campaign rather than a second independently written complete propagation verifier.

## Run

```bash
python bounded_campaign.py
python -O bounded_campaign.py
python bounded_intervention.py data/twin_path.observations.json --eta 1/1000000000 --relative-load-error 1/1000000000 --gap 1/4 --output data/twin_path.bounded.json
```

`infer_bounds(...)` also accepts direct interval maps keyed by canonical tuples of zero-based nominal site/load pairs. `bound_document(...)` accepts the exact-core observation format and treats those serialized responses as interval centers, with the external `eta` supplied separately.

## Retained campaign

The campaign checks 25 connected rational networks of orders 2 through 6 at four absolute response-error envelopes, and 75 additional acquisitions whose actual shunt values and reported response centers vary independently inside their declared calibration intervals. Exact zero-noise, 1e-12 and 1e-9 response envelopes enclose every random network; the conservative 1e-6 response case leaves two unresolved. At 1% load uncertainty, 21 of 25 cases remain unresolved. Every accepted matrix and edge interval contains its simulated truth; the saved record counts all predicates and preserves unresolved outcomes.

For the v7 twins, one unit load at the second displayed site separates responses `7/26` and `1088/3989`. With exact loads and a common absolute per-reading error bound `eta`, their prediction intervals are disjoint precisely when `eta < 365/207428`. This finite-pair discrimination margin is much larger than the margin needed by the conservative full-network reconstruction; the two guarantees concern different acquisition tasks.
