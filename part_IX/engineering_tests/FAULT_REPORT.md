# Single-fault diagnosis under a calibrated electrical acquisition model

## Result and decision

A six-reading fixed-port intervention plan identified all 100 fresh single-fault specimens in the calibrated campaign at a 10 microvolt readout envelope. Every true fault's complete prediction box was disjoint from every alternative's box under that plan, so the result holds for all scalar errors and shunt variations inside the declared boxes. The same six-reading budget also allowed conventional reciprocal scans to identify all 100. At a 100 microvolt envelope, selected fixed-port interventions resolved 97/100 while a strong reciprocal scan resolved 100/100. The practical choice is governed by available physical access: the intervention protocol fixes both source and voltage readout at one pad; reciprocal scans move the electrical observation among internal pads.

The calibrated campaign uses a 1 mA source and a 1 mS conductance normalization, making a normalized response unit equal to one volt. Each shunt with nominal physical resistance R has an independent resistance envelope ±(0.0008R + 0.07 ohm), followed by an independently varying relay/contact series resistance in [0, 0.05] ohm. These are declared simulation scenarios. Current, gain and voltage error must be included in the externally supplied normalized response envelope in a bench test.

## Frozen acquisition and specimens

Two distinct frozen protocols were executed. The initial campaign contains 20 independent healthy six-node fixtures, five selected fault specimens per fixture, four response-error levels and five acquisition policies: 2,000 policy/noise trials on 100 distinct specimens. A subsequent prospective calibration extension uses a different confirmation seed and 20 new fixtures, again five faults per fixture, three response-error levels and the same five policy definitions: 1,500 trials on 100 additional specimens. The second protocol was frozen before its new confirmation fixtures were generated. It preserves the first campaign separately.

Every healthy fixture is a positive rational weighted six-cycle with zero to three extra chords and nonnegative grounding. The catalog contains the healthy fixture, each existing-edge weakening by 5% or 20%, each edge removal that preserves connectivity, an added 1/4-strength bridge at every absent edge, and a 20% increase at each ground. A confirmation specimen selects one candidate from each of the five fault categories. The healthy reference and candidate family are available to the experiment designer. The selected specimen and its matrix are hidden from the classifier.

The classifier consumes only the candidate prediction table, load/source/readout query labels, observed scalar values and error bounds. The true matrix is used by the acquisition oracle and scorer. The calibrated acquisition randomizes actual resistor and contact values independently within their envelopes; each specimen/query's actual physical signal is then shared among all policies and all noise levels that use that query. Thus a policy does not receive a more favorable physical realization.

The catalog experiment identifies a labeled finite fault alternative. Full continuous-network reconstruction is supplied by the separate Nonagon theorem and solver. A one-value exact fingerprint can distinguish a finite catalog, and that happens for the selected noiseless specimens here; it does not change the continuous reconstruction count.

## Acquisition policies and accounting

The canonical intervention plan consists of the baseline, strengths 1 and 2 at each hidden site, and simultaneous unit loads at each hidden pair: 21 scalar readings for six sites. The selected fixed intervention plan draws from 61 allowed single/pair load configurations with strengths 1/4, 1, 2 or 4. It is computed before any fault observation, repeatedly separating the largest number of as-yet-overlapping candidate prediction boxes. Its algorithm, pool and tie rule are frozen. This is fixed design conditioned on a known fixture and catalog, rather than an adaptive policy using the observed fault.

The conventional pool contains all 21 distinct reciprocal entries of the six-by-six compliance matrix. A reading of G_ij sources current at i and measures voltage at j. The natural order starts with the six responses at fixed source site 1, followed by the remaining reciprocal entries. The strong conventional comparator selects among the same 21 entries with exactly the same greedy box-separation rule used for interventions. Both have the same scalar budgets of 1, 6, 12 and 21. Per-plan records include distinct source addresses, readout addresses, shunt pads and maximum simultaneous shunts. Counting six source locations as six scalar measurements would undercount a full 21-entry scan; the retained accounting counts each scalar entry separately. No synthetic production-time or hardware-cost advantage is assigned.

## Calibrated confirmation at 10 microvolt

| Policy | Scalar readings | Identified | All-error-box guarantee | Unresolved |
|---|---:|---:|---:|---:|
| Original port | 1 | 91/100 | 91/100 | 9/100 |
| Canonical interventions | 1 | 91/100 | 91/100 | 9/100 |
| Canonical interventions | 6 | 99/100 | 99/100 | 1/100 |
| Canonical interventions | 21 | 100/100 | 100/100 | 0/100 |
| Selected fixed-port interventions | 1 | 93/100 | 93/100 | 7/100 |
| Selected fixed-port interventions | 6 | 100/100 | 100/100 | 0/100 |
| Selected fixed-port interventions | 21 | 100/100 | 100/100 | 0/100 |
| Natural reciprocal scan | 1 | 91/100 | 91/100 | 9/100 |
| Natural reciprocal scan | 6 | 100/100 | 100/100 | 0/100 |
| Natural reciprocal scan | 21 | 100/100 | 100/100 | 0/100 |
| Selected reciprocal scan | 1 | 99/100 | 99/100 | 1/100 |
| Selected reciprocal scan | 6 | 100/100 | 100/100 | 0/100 |
| Selected reciprocal scan | 21 | 100/100 | 100/100 | 0/100 |

## Calibrated confirmation at 100 microvolt

| Policy | Scalar readings | Identified | All-error-box guarantee | Unresolved |
|---|---:|---:|---:|---:|
| Original port | 1 | 62/100 | 62/100 | 38/100 |
| Canonical interventions | 1 | 62/100 | 62/100 | 38/100 |
| Canonical interventions | 6 | 90/100 | 90/100 | 10/100 |
| Canonical interventions | 21 | 94/100 | 93/100 | 6/100 |
| Selected fixed-port interventions | 1 | 69/100 | 68/100 | 31/100 |
| Selected fixed-port interventions | 6 | 97/100 | 97/100 | 3/100 |
| Selected fixed-port interventions | 21 | 97/100 | 97/100 | 3/100 |
| Natural reciprocal scan | 1 | 62/100 | 62/100 | 38/100 |
| Natural reciprocal scan | 6 | 99/100 | 99/100 | 1/100 |
| Natural reciprocal scan | 21 | 100/100 | 100/100 | 0/100 |
| Selected reciprocal scan | 1 | 68/100 | 68/100 | 32/100 |
| Selected reciprocal scan | 6 | 100/100 | 100/100 | 0/100 |
| Selected reciprocal scan | 21 | 100/100 | 100/100 | 0/100 |

At 10 microvolt, all five fault categories are resolved in all 20 examples by six selected interventions. The canonical order needs 12 readings to resolve the final case. At 100 microvolt, the selected intervention plan's three unresolved specimens remain unresolved through 21 readings. Its chosen load pool supplies insufficient box separation for those alternatives. Conventional multipoint readings provide stronger separation in that scenario.

Observed uniqueness and the all-error-box guarantee are reported separately. For each true candidate and plan prefix, the audit computes the smallest, over all alternative candidates, of the largest signed separation between their predicted scalar intervals, subtracting twice the response-error envelope. A strictly positive margin means the candidate uncertainty boxes are disjoint. At 100 microvolt the canonical protocol has 94 observed unique decisions but only 93 guaranteed decisions; a particular bounded acquisition can be identifying even when an alternative allowed acquisition remains ambiguous.

## Initial exact-load campaign

The initial campaign keeps shunts exact, allowing direct comparison to the existing one-port formulas. At normalized error 1e-5, equivalent to 1 microvolt with its separately declared 100 microampere source and 1 mS normalization, all four full protocols identify all 100 specimens by six readings. The original unloaded port resolves 86. At normalized error 1e-4, selected fixed-port interventions resolve 97/100, the canonical 21-setting plan resolves 95/100, and the selected reciprocal scan resolves 100/100 by six readings. The two campaigns use different physical current conventions and independent fixtures; their percentages should not be pooled as repeated measurements of one apparatus.

## Arithmetic checks and controls

All conductance matrices, shunt bounds, signals, prediction intervals, errors and compatibility decisions are exact rational values. Candidate intervention predictions use sequential rank-one updates after exact matrix inversion; development checks compare these with a separately executed pivoted Gaussian equilibrium solve. The calibrated campaign retains 8,200 exact checks that each actual per-query signal lies inside the corresponding candidate prediction interval. Every selected fault survives its own measurement envelope, and no trial produces a wrong unique classification.

The initial campaign additionally tests 20 two-edge changes outside the single-fault catalog. Their full noiseless canonical records are rejected as inconsistent with the catalog. These controls demonstrate that those particular excluded alternatives are not silently assigned a single-fault label. They do not turn a finite fault catalog into a general unknown-fault model. Wider component families and uncertain healthy references require a separately declared candidate or parameter model.

A separate reviewer is replaying the saved circuit predictions and acquisition histories with independent Gaussian solves. The files below retain all matrices and observations needed for that replay; any independent receipt should be read with its own stated coverage.

## Reproduce

Run `python3 fault_benchmark.py` for the initial exact-load campaign and `python3 fault_calibration.py` for the prospective calibrated campaign. Both use the Python standard library only. Existing protocol files are checked byte-for-byte and a changed protocol is refused. The outputs retain protocol/source hashes, candidate matrices, query definitions, pre-observation plans, physical signals, noisy readings, survivor histories, robust margins, and aggregate summaries.

Files: `fault_protocol.json`, `fault_results.json`, `fault_calibration_protocol.json`, `fault_calibration_results.json`. The two result archives contain 200 distinct confirmation specimens and 3,500 bounded-noise policy trials. These are electrical simulation tests; no physical instrument was connected.
