# Finite loads selected from measured pilot responses

The planner chooses **one globally shared positive base load** from a finite hardware list. Its final acquisition remains the production protocol: baseline, single-site loads \(t,2t\), and pair loads \(t,t\). No observation schema or inverse branch is changed. Two routed shunt channels still suffice.

The decision uses only calibrated pilot scalar intervals and the declared current, voltage, shunt and relay envelopes. It accepts no matrix or hidden topology. The pilot measures the baseline and \(t_0,2t_0\) at every hidden site, giving \(2n-1\) scalar readings. If those intervals cannot separate the required positive drops and intercepts, it returns the existing \(t_0\) setting explicitly as a fallback.

For \(b=G_{jj}\), \(q=G_{0j}^2\), and \(u=tb\), the relative first-order error in the recovered \(q\), under independent raw-reading errors bounded by \(\eta\), is

\[
\frac{\eta b}{q}\left(\frac4u+8+4u\right).
\]

The shared baseline is included in the three-reading Jacobian. Its finite-load optimum is **\(u=1\)**. Independent Gaussian raw readings instead give the variance factor

\[
\left(\frac{\sigma b}{q}\right)^2
\left(\frac{6.5}{u^2}+\frac{24}{u}+34+24u+8u^2\right),
\]

whose optimum is approximately \(0.9720604037\). Neither optimum is the \(1/\sqrt2\) value for the separate ideal-clamp experiment.

The deployed bounded-error diagnostic uses the pilot-enclosed upper bound \(b_+\) and lower bound \(q_-\):

\[
D_j(t)=\frac{\eta}{q_-}\left(\frac4t+8b_++4tb_+^2\right)+3r(t).
\]

Here \(r(t)\) encloses the relative calibration error of either independent shunt at every setting through \(2t\), including relay series resistance. The term \(3r\) is the local calibration contribution: independent relative load perturbations have relative-\(q\) gradient \((2,-1)\). The planner minimizes \(\max_jD_j(t)\) over the hardware list. The baseline upper bound controls the multiplicative source-current contribution to \(\eta\), since positive grounding loads decrease the response.

**This is a local sensitivity diagnostic, not a nonlinear certificate.** Optimizing single-site \(q\) sensitivity need not optimize pair-root conditioning, inversion of the full compliance matrix, or every graph decision. The unchanged rational interval inverse is the acceptance gate for all final matrix and edge claims.

## Prospective calibration benchmark

The frozen protocol contains 30 new grounded circuits, with orders 3, 4 and 6, path and complete graphs, and five conductance scales from \(1/16\) to 16. Matrix hashes are checked against the earlier frozen engineering and production fixtures. There are three declared synthetic instrument scenarios and two acquisition arms, giving **90 prospective comparisons and 180 arm cases**.

Both arms receive the same pilot. The fixed arm retains \(t_0=1\); the planned arm chooses from \(\{1/16,1/4,1,4,16\}\). Both then acquire a complete canonical plan. This isolates strength selection at equal pilot-plus-final cost. A fixed plan can be run without the pilot more cheaply; this experiment establishes no general time or cost advantage for planning.

| Order | Pilot responses | Final canonical responses | Total per arm |
|---|---:|---:|---:|
| 3 | 5 | 6 | 11 |
| 4 | 7 | 10 | 17 |
| 6 | 11 | 21 | 32 |

At the unit scale of 1 mS, nominal resistances at all candidate \(t,2t\) settings lie within the declared 20–20,000 Ω range. The current is 100 µA. Every noisy scenario includes independently varying current gain, shunt calibration and relay resistance, and bounded additive voltage error and drift. These are explicit scenarios, not combined manufacturer accuracy claims.

| Scenario, 30 circuits each | Fixed: matrix enclosed | Planned: matrix enclosed | Fixed: full graph | Planned: full graph |
|---|---:|---:|---:|---:|
| Ideal reference | 30 | 30 | 30 | 30 |
| Precision calibration | 20 | 25 | 5 | 10 |
| Readout limited | 4 | 8 | 0 | 1 |

The precision scenario uses 10 nV readout error, 1 nV drift, 1 ppm current and shunt uncertainty, and relay resistance from 0 to 5 mΩ. The readout-limited scenario changes the readout and drift bounds to 10 µV and 1 µV. There were no false edge calls or escaped true matrices. Eleven readout-limited pilots could not establish all required parameter margins and returned the fixed strength.

The planner produced five additional matrix enclosures and five additional complete graph decisions in the precision scenario, and four additional matrix enclosures and one complete graph decision in the readout-limited scenario. No case was resolved only by the fixed arm. Among cases enclosed by both arms, planned matrix intervals were narrower in 14 of 20 precision cases and 2 of 4 readout-limited cases; the remaining widths were equal. These results apply to the frozen circuits and envelopes.

Repeated baseline responses share one physical truth and can be intersected. Independently calibrated repeated loaded readings can correspond to different actual shunt values, so their scalar intervals **are not intersected**. For a repeated loaded configuration, each arm may retain the narrower individually valid interval; the load envelope covers whichever reading is used. No averaging or square-root error reduction is assumed for bounded errors. Complete graph absence decisions additionally require the declared minimum positive edge weight of \(1/32\).

No physical apparatus is connected. The records contain synthetic raw physical voltages, actual currents and shunts, calibrated scalar intervals, pilot recommendations, and final conditional interval certificates. The initial development protocol is retained as superseded provenance; protocol v2 fixes the repeat-calibration policy and uses new fixture seeds before the reported campaign.

## Run

Use Python 3.11 or later. The planner, benchmark, and independent checks need only the standard library and the sibling production interval module.

```sh
python part_IX/load_planning/benchmark.py
python part_IX/load_planning/independent_checks.py
python part_IX/load_planning/planner_input_controls.py
python part_IX/load_planning/planner.py part_IX/load_planning/results/example_pilot.json --output recommendation.json
```

For a replay without replacing the saved results:

```sh
python part_IX/load_planning/benchmark.py --output-dir replay/load_planning
```

`protocol.json` is the authoritative frozen input. `results/benchmark_receipt.json` contains aggregate counts, source and protocol hashes; `results/benchmark_cases.jsonl` retains the complete physical and inference records. A recommendation includes an executable nominal configuration sequence, the selected finite load, the measured pilot enclosures, every candidate score, and the scalar cost including its pilot. External calibration remains a prerequisite for using actual bench readings.
