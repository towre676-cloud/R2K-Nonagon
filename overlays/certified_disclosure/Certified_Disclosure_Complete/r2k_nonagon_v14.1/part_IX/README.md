# Nonagon IX: intervention tomography

**One measurement port, two routed loads, a reconstructed network.** Part IX
expands the Nonagon from interpreting responses to designing interventions that
expose physical connections. Its first industrial fit is a grounded resistor
network: inject a calibrated current at one known pad, read voltage at that same
pad, and route calibrated shunt resistors among known internal test pads. The
result is a labeled conductance map, a diagnosis of changed branches, and a
calculated shunt that achieves a specified static response.

This is also exactly the mathematics of a grounded spring network: injected
current maps to force, voltage to displacement, nodal conductance to stiffness,
and a shunt conductance to an added grounding spring. The exact implementation
uses the spring-model names in its common mathematical observation contract. The
hardware adapter documents the electrical normalization.

## The finite reconstruction theorem

Let K be the symmetric positive-definite matrix of a connected grounded network
with nonnegative inter-node weights and nonnegative grounding. The n sites and
their labels are known. Force/readout gains and the applied grounding loads are
calibrated. A unit input and scalar readout share port 1.

For any chosen t>0, acquire the unloaded response, the responses to t and 2t at
each other site, and the response to t applied simultaneously at each pair of
other sites. These **n(n+1)/2 scalar configurations uniquely determine K** within
that class. Every configuration uses at most two added loads. The derivation
combines rank-one and rank-two inverse updates with strict positivity of K's
inverse; positivity selects the unique admissible quadratic branch.

Thus one fixed measurement port can replace full-field sensing when independently
addressable loads are available at the known sites. The theorem concerns those
sites; arbitrary inaccessible hidden nodes require a different experiment.
`README_CORE.md` gives the complete finite derivation and exact schema.

## Completed results

| Capability | Retained result |
| --- | --- |
| Exact reconstruction | 100 fresh rational connected networks, orders 2–10; all recovered exactly |
| Accepted exact artifacts | 103 network certificates; 2,430 acquisition readings and 309 held-out readings |
| Independent exact replay | 103 certificates accepted in 119,902 checks |
| Independently checked applications | Fault localization, static inverse design and retained damped twins: 1,139 checks |
| Live altered-artifact controls | 20 network-certificate alterations and nine application alterations rejected |
| Core observation controls | 24 corrupted or unsupported inputs rejected |
| Active intervention selection | 21 otherwise equivalent candidates, 15 feasible loads, 189 bounded-error trials; zero wrong identifications |
| Bounded inference | 100 response-interval evaluations plus 75 independently uncertain load/readout evaluations; 6,480 containment checks |

The retained v7 path and triangle have exactly the same original driving-point
transfer. One unit shunt at the same labeled second site changes their DC
responses to **7/26** and **1088/3989**, respectively. Intervening breaks an
ambiguity that unlimited unchanged-port frequency sampling cannot resolve.

The fault example reconstructs a spring weakened from 1 to 1/2 and isolates that
one changed edge. The inverse-design example calculates an added grounding
stiffness of 9/4 to achieve the port response 19/72, then verifies the forward
equation independently.

The active finite-library study identifies all 21 candidates with one informative
load at the tested error envelopes. Its strongest fixed load choice matches the
adaptive policy. The useful result is informative intervention design at a fixed
budget, rather than a claim that adaptation always improves the experiment.

## Bounded measurements and bench readiness

`bounded_intervention.py` accepts simultaneous intervals for normalized
voltage/current or displacement/force, with separately varying relative
calibration intervals for every applied load. Exact rational arithmetic and
outward square-root bounds propagate these inputs to conditional matrix and edge
enclosures. Wide or inconsistent margins remain unresolved. At positive error,
an absent-edge conclusion additionally needs a declared minimum positive weight.

The retained path/triangle remain distinguishable with ±1% shunt calibration
uncertainty when additional normalized response error is below approximately
0.0017108444 per reading. Full network reconstruction has stricter conditioning
requirements than distinguishing that particular pair. `README_BOUNDS.md` records
both acquisition tasks and their margins.

`hardware_plan.py` supplies instrument-neutral acquisition plans and normalized
import adapters. The three-node plan uses six configurations; the six-node plan
uses 21. Their CSV/JSON templates are in `data/`. Two calibrated programmable
load channels and a relay-routing matrix suffice for the intervention schedule.
The physical pilot should retain current, voltage, shunt-calibration and
repeatability envelopes, followed by unused load settings for prediction checks.

## Replay

Python 3.10 or later is sufficient for the exact reconstruction, independent
verifier, active finite-library study and interval campaign. Run from this folder:

```bash
python replay.py
```

The default independently checks the saved 103 certificates and three application
artifacts, reruns live corruption controls, regenerates the active finite-library
study and runs the bounded containment campaign. It does not regenerate the 100
exact acquisition fixtures. Its explicit receipt is
`data/replay_receipts.json`.

To regenerate those exact acquisitions before verification:

```bash
python replay.py --full
```

To verify behavior with Python assertions disabled:

```bash
python -O replay.py --skip-numeric-audit --output data/replay_receipts_optimized.json
```

The optional historical offset audit needs NumPy, SciPy and the earlier chapter
inputs. In the complete Nonagon package it locates the sibling `part_VI` and
`part_VII` folders automatically. For a standalone Part IX folder:

```bash
python replay.py --archive-root /path/to/nonagon-root --numeric-python /path/to/python
```

`--skip-numeric-audit` explicitly skips that optional stage. Missing numerical
dependencies or historical source inputs are recorded as **SKIPPED**, never as a
passed audit. A numerical audit that starts and fails makes replay fail.

## The BK offset repair

The historical number-variance audit finds that the archived BK and GUE reference
integrals both started at 0.0001 instead of zero. The missing contribution is
approximately **0.0002 L**: at L=16 it is 0.00319999. This accounts for the
prominent height-stable offset. The audit reproduces all six archived GUE values,
checks the sine-kernel variance against independent adaptive quadrature, and
applies the universal leading endpoint correction to all 18 archived BK rows.
After that correction, the largest residual is 1.454 archived block standard
errors. The empirical zero statistics are unchanged.

`OFFSET_ANALYSIS.md` explains the endpoint calculation, the remaining mesh error,
and the production integration repair. `offset_audit.json` retains the numerical
reproduction and binds it to the inspected source hashes. The leading endpoint
repair is distinct from a full freshly integrated BK curve.

`variance_quadrature.py` supplies the repaired zero-inclusive, endpoint-aware
integration. Its 48 regression checks compare the analytic sine kernel, adaptive
quadrature, refined integration, exact polynomial cases and invalid inputs.
`bk_multi_repaired.py` uses that integration in a replacement campaign runner;
the full BK/zero campaign is an explicit separate action. The optional numerical
replay runs the 48 lightweight regressions along with the archived offset audit.

## Files

| File | Role |
| --- | --- |
| `intervention.py` | Exact finite observer, diagnosis and static design |
| `campaign.py` | Reproducible exact synthetic acquisition campaign |
| `verify_intervention.py` | Separately written exact certificate and application verifier |
| `bounded_intervention.py` | Conditional bounded-error and load-calibration inference |
| `bounded_campaign.py` | Exact arithmetic containment campaign |
| `active_interventions.py` | Finite-library load selection and matched fixed-policy comparison |
| `hardware_plan.py` | Frozen acquisition templates and normalization adapter |
| `offset_audit.py` | Archived BK/GUE endpoint investigation |
| `variance_quadrature.py` / `test_variance_quadrature.py` | Repaired integration and 48 numerical regressions |
| `bk_multi_repaired.py` | Replacement BK campaign runner; full campaign remains explicit |
| `replay.py` | Integrated evidence replay with explicit optional-stage verdicts |
| `CHAPTER.md` | Integrated mathematical and engineering account |
| `README_CORE.md` / `README_BOUNDS.md` | Detailed exact and interval contracts |

The exact checker establishes the serialized model and its equations. The bounded
campaign establishes containment on its retained synthetic cases under the stated
input envelopes. Physical measurements are the next step for the acquisition
fixture.
