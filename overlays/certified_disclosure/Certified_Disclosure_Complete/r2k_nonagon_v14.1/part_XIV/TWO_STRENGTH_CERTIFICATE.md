# Two-strength response elimination with exact uncertainty

The new inversion eliminates the pair quadratic before evaluating intervals.
It replaces proximity-based root selection with a data-only rational formula
and an explicit separation gate. Symmetry and nonsingularity are required;
positive definiteness is not.

Write `u = G0i`, `v = G0j`, `A = Gii`, `B = Gjj`, and `r_l = 1/theta_l`.
The measured response difference at strength `theta_l` satisfies

\[
Q_l(x)=\delta_l x^2-2uvx+C_l=0,\qquad
C_l=u^2(r_l+B)+v^2(r_l+A)-\delta_l(r_l+A)(r_l+B).
\]

Multiplying `Q1` by `delta2` and `Q2` by `delta1` cancels the quadratic term.
Consequently every common root is

\[
x=\frac{\delta_2C_1-\delta_1C_2}
        {2uv(\delta_2-\delta_1)}.
\]

If the denominator is nonzero, at most one common root exists. The physical
response model supplies its existence. This argument also permits one of the
response differences to be zero: it never divides by `delta1` or `delta2`
individually. Nonsingularity of the loaded matrices remains part of the
physical-model assumptions.

## Computation

`two_strength_cert.py` evaluates the formula with exact rational interval
endpoints. It intersects two algebraically identical expressions, including
the cancellation-aware numerator

\[
(\delta_2-\delta_1)(u^2B+v^2A)
+(\delta_2r_1-\delta_1r_2)(u^2+v^2)
-\delta_1\delta_2(r_1-r_2)(r_1+r_2+A+B).
\]

It then contracts the enclosure using interval Newton steps for both original
quadratics. Every division requires its interval denominator to exclude zero.
Outward dyadic rounding during contraction controls rational denominator
growth; it enlarges enclosures and cannot invalidate inclusion. Residual
enclosures for both original quadratics must contain zero. This residual test
is a necessary consistency check, not a proof that an arbitrary input box has
a physical realization.

The result encloses every physically compatible input tuple in the declared
box. The inverse receives only the calibrated intervals and the measured
responses. Model truth is used afterwards to validate coverage in a synthetic
campaign.

`verify_certificate` reconstructs every exact operation from the serialized
inputs and checks the full certificate. The recorded denominator margin makes
the exclusion of zero explicit. A denominator interval crossing zero produces
a rejection, not a guessed root.

## Raw single-node calibration

The optional calibration path avoids treating inferred diagonal and port
column entries as exact. For two single-node response differences `d1,d2`
at strengths `theta,2theta`,

\[
G_{jj}=\frac{2d_1-d_2}{2\theta(d_2-d_1)},\qquad
G_{0j}^2=\frac{d_1d_2}{2\theta(d_2-d_1)}.
\]

Exact interval division and a rigorous dyadic square-root enclosure propagate
their reading errors. The sign of `G0j` must be known from calibration or from
the connected passive M-matrix assumptions. No sign is inferred from the
positive square alone.

If individual baseline and loaded readings each have deterministic error at
most `sigma*abs(G00)`, a difference has error at most `2*sigma*abs(G00)`.
Shared baseline correlations are ignored conservatively. The grid campaign
reports separately the conditional exact-calibration path and the raw
single-node uncertainty path. Its error budgets are deterministic bounds,
not Gaussian standard deviations.

## Equal-budget cancellation in the raw observations

`certify_raw_pair` improves uncertainty propagation without collecting an
additional response. Put `s_i1,s_i2` and `s_j1,s_j2` for the two single-node
differences at each endpoint, and `d1,d2` for the pair differences. Then

\[
(uv)^2=\frac{s_{i1}s_{i2}s_{j1}s_{j2}}
 {4\theta^2(s_{i2}-s_{i1})(s_{j2}-s_{j1})},
\]

\[
x=\frac{uv}{2(d_2-d_1)}
\left[
\frac{d_2(s_{i1}+s_{j1}-d_1)}{s_{i1}s_{j1}}
-\frac{d_1(s_{i2}+s_{j2}-d_2)}{s_{i2}s_{j2}}
\right].
\]

The known sign product chooses the sign of `uv`. To derive the expression,
observe that `u^2=s_i1(r+A)=s_i2(r/2+A)`, and similarly for `v`, where
`r=1/theta`. Thus `C1=(r+A)(r+B)(s_i1+s_j1-d1)` and
`C2=(r/2+A)(r/2+B)(s_i2+s_j2-d2)`. Substitution in the two-strength
elimination formula cancels the inferred diagonal entries. The interval
implementation intersects this direct raw-expression bound with the
propagated-calibration bound, so its certified interval cannot be wider.

`two_strength_raw_design.py` measures this equal-reading-budget gain on every
IEEE 118 lean pair with three intervention scales and three error budgets.
Its fixture certificates are fully serialized for exact replay.

## Four-strength acquisition

`certify_many` takes arbitrary distinct intervention strengths and intersects
every pair elimination whose separation denominator excludes zero. The
selection and intersection use only observation intervals. A new response can
therefore improve a bound even when some strength pairs are too poorly
separated to certify. A tuple inconsistent with all observations need not
exist; the guarantee always concerns every compatible physical tuple.

`two_strength_design.py` compares strengths `{1,2}` with `{0.1,1,2,10}`, all
scaled by the reciprocal median diagonal. This spends four pair responses
instead of two. It reports the calibrated-input and raw-single-calibration
paths separately; no claim of equal acquisition cost is made for this second
experiment.

## Retained numerical results

The two-grid campaign evaluates 14,020 calibrated-input boxes across 40
settings on all 207 IEEE 118 and 494 IEEE 300 lean pairs. It accepts 13,363
intervals and explicitly rejects 657 cases whose declared uncertainty is too
large for these certificate gates. Every accepted interval contains the exact
rational model value. Another 1,656 propagated raw-calibration intervals are
accepted and all contain the exact value. These are deterministic coverage
checks under the declared model and error boxes.

The equal-budget raw-cancellation experiment at normalized strength `c=1`
and `sigma=1e-6` accepts every IEEE 118 pair:

| Metric | Propagated calibration | Raw algebra cancellation |
|---|---:|---:|
| Median relative interval radius | 2.3420% | 1.0826% |
| 90th percentile relative radius | 5.0674% | 2.7474% |
| Exact model values contained | 207/207 | 207/207 |

This gain uses the same six response differences. At normalized strength
`c=10`, the number of pairs certified to a relative radius at most 1% rises
from 59 to 79 under the same error budget.

The four-response acquisition experiment at `sigma=1e-6` gives:

| Calibration supplied to the inverse | Two pair responses: median radius | Four pair responses: median radius |
|---|---:|---:|
| Exact calibrated entries | 0.028557% | 0.016381% |
| Propagated raw single-node intervals | 2.34581% | 1.76930% |

All 207 exact model values are enclosed in both rows. The separate campaigns
use different deterministic perturbation patterns, so their displayed medians
should be compared within each table.

## Reproduction

```bash
python -m pytest -q part_XIV/test_two_strength_cert.py
python part_XIV/two_strength_cert.py --campaign
python part_XIV/two_strength_raw_design.py
python part_XIV/two_strength_design.py
python part_XIV/replay_two_strength_certificates.py
```

The campaign loads the retained unmodified IEEE 118 and IEEE 300 MATPOWER
data, inverts their matrices over the rationals, and examines every lean pair
at four strength scales and five deterministic error budgets. It records
certified counts, rejected cases, truth-coverage checks, relative interval
radii, and the two earlier point-root rules on the same perturbed centres.
The point rules' numerical misses are scored at a relative tolerance of
`1e-6`; interval success is certified inclusion, so these are different
performance criteria rather than interchangeable success rates.
