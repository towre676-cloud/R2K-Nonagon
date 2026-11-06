# Independent direct-Arb certification campaign

This implementation imports **none of the Nonagon archive's Taylor/DFT
contour code**. It uses a separately written direct argument-principle engine
and Arb's adaptive integrator. Arb's function ball evaluations and integrator
remain trusted dependencies; this is numerical-method independence, not an
independent implementation of interval arithmetic.

Run `python independent_certification/run_campaign.py` from the frontier
bundle. A standard environment needs `python-flint==0.9.0`. The script also
recognizes the task's adjacent local dependency directory. Use
`--skip-lehmer` to omit the optional 60-second high-altitude challenge.

## Frozen finite campaign

The protocol freezes the entire-function choice, engine hash, 100-bit
precision, tolerance, work limits and initial **45 fixed tiles of width 2**
before numerical execution. It covers the continuous rectangle

- `10 <= Re(t) <= 100`, `-3/5 <= Im(t) <= 3/5` for Xi(t); equivalently
- `10 <= Im(s) <= 100`, `-1/10 <= Re(s) <= 11/10` for completed zeta.

This covers the whole critical strip at every height in the stated finite
range. It does not cover the range below 10 or imply a new verification-height
record. All 45 tiles succeed without adaptive splitting, with a total of **29
zeros, all simple and on the critical line**. Every empty tile is retained in
the coverage ledger, so no intervals are omitted. The root list is not input
to tile placement, counting, boundary checks, or matrix tests. A diagnostic
reference count is computed only after the completed ledger is written.

`campaign_results.json` contains the exact range, the complete attempt/leaf
ledger, certificate hashes, unresolved regions (empty in the completed run),
timings and function-evaluation counts. `tiles/` contains each full receipt.

## Independent numerical route

Each contour is handled as four complete directed sides. Before integration,
an independent dyadic subdivision of each side verifies, using function balls,
that the entire boundary is zero-free. For Xi at these positive heights, zeta
nonzero implies Xi nonzero because the completion factors are finite and
nonzero there. Accepted dyadic subintervals cover every side without gaps or
overlaps.

The engine then directly integrates the logarithmic derivative, weighted by
Chebyshev polynomials. For zeta's derivative it uses only a guarded two-term
power series, with global series cap 8. There is no high-order Taylor
remainder argument and no DFT coefficient extraction. The count enclosure must
both **contain and isolate** a nonnegative integer. The moments must enclose
real values. Positive Sylvester leading minors of the resulting Hermite
matrix prove that the enclosed zeros are real and simple. Count-one tiles use
the same symmetry consequence directly; empty tiles are vacuous successes.

All function classes currently supported are explicit: real-coefficient
polynomials and completed zeta in the t variable. The code is a reusable
contour engine for this class, not a parser that accepts arbitrary claims of
analyticity from a function name.

## Fail-closed controls and artifact checks

Five analytic controls run through the same contour engine:

| Control | Required outcome |
|---|---|
| Two distinct real roots | Real/simple certificate |
| Off-axis conjugate pair | Negative matrix direction; not all real/simple |
| Repeated root | No real/simple certificate; unresolved margin |
| Root on a contour boundary | Boundary unresolved; integration not certified |
| Two real roots separated by 1e-6 | Positive determinant and real/simple certificate |

No unresolved control is silently converted into success. Its complete reason
and any partial boundary work are retained under `controls/`.

Arb's human-readable `.str()` can deliberately coarsen wide interval balls.
Therefore each saved complex ball also includes exact dyadic midpoint/radius
fields. These are used for boundary-receipt checks; decimal strings are
diagnostics. Reconstructing a dyadic ball uses outward radius rounding and may
slightly widen it; the checker still requires the reconstructed ball to exclude
zero. `serialization_regression.py` tests this distinction explicitly.
The initial display-only run/source and its static-replay failure are retained
under `initial_display_run/` and are superseded by the current exact-dyadic
receipts.

`python independent_certification/replay.py` checks source/input hashes,
coverage, stored boundary exclusions, integer containment, positive minors and
seven deliberate corruptions. These static checks establish artifact
consistency, not the mathematical provenance of stored function evaluations.
`python independent_certification/replay.py --recompute` freshly recomputes
every tile with the independent direct-Arb engine and compares all outcomes.
The replay remains dependent on Arb and the direct engine; it does not use the
archive's Taylor kernel.

The completed fresh replay agrees on all 45 tiles. `final_results.json`
aggregates the campaign, fresh replay, five analytic controls, seven rejected
integrity corruptions, three serialization checks, and the retained Lehmer
limitation. Regenerate it with `python independent_certification/finalize_results.py`
after running the campaign, serialization checks, and full replay.

## Bounded Lehmer limitation

The initial independent direct method was also given a 60-second challenge
for the known Lehmer-pair rectangle near t=7005, without supplying root
ordinates. It did not finish within that budget. The request and timeout are
preserved. There is no completed direct-route Lehmer certificate or normalized
moment comparison from that attempt. This demonstrates a practical limit of
direct interval integration at higher height, while the corrected Taylor/DFT
route in the earlier challenge succeeds there.

All finite-height zero results here are established benchmarks reproduced by
a second numerical route, not new facts about the locations of zeta zeros or
a global proof of RH.
