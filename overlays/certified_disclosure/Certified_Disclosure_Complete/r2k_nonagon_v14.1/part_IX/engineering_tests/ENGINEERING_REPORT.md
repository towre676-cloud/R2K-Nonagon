# Nonagon Engineering Tests

## From exact reconstruction to a bench qualification

Six preselected fixed-port readings distinguish all 100 new fault specimens in the calibrated confirmation challenge, with positive separation margins across the declared measurement and shunt-error boxes. The setup uses a 1 mA source, 1 mS conductance scaling and a 10 microvolt readout-error envelope. It includes resistor-setting uncertainty and relay resistance. This is the strongest immediate result: the Nonagon now has a concrete bounded-error fault-diagnosis protocol ready for a coded resistor-fixture bench test.

Independent circuit simulation also confirms that the local-load equations translate correctly into electrical components. Across fourteen resistor fixtures, ngspice supplied 286 DC responses, and reconstruction from those scalar responses recovered the complete normalized conductance matrices to a maximum coefficient error of 2.46e-12. The circuit-simulator check and the fault challenge test different parts of the engineering chain.

The engineering campaign separates that task from unrestricted graph recovery. Distinguishing a finite set of possible faults needs less precision than enclosing every unknown conductance. The v7 path and triangle remain distinguishable across all tested calibration envelopes, while full-graph interval disclosure often needs much tighter bounds. Conventional reciprocal scans provide a strong comparison when both current injection and voltage readout can move. Fixed-port intervention tomography addresses a different access arrangement: the current and voltage measurement stay at one site, while known shunts can be routed to labeled internal sites.

All results below concern retained numerical circuit experiments and calibrated synthetic acquisitions. Physical instrument measurements in this release number zero. The test package concentrates on circuit translation, calibration, diagnosis and an executable measurement interface.

## External circuit-simulator confirmation

The SPICE fixtures comprise the two v7 twins and twelve fresh rational networks of orders three through eight. A normalized conductance unit is 1 mS and the source current is 100 microampere. Each graph is converted to actual edge resistors, grounding resistors and shunt loads. The external ngspice-39 executable performs the circuit solve; a separate floating-point inverse receives its normalized port readings and the public intervention labels.

| External check | Retained result |
| --- | ---: |
| Independent resistor fixtures | 14 |
| DC acquisition and reserved settings | 286 |
| Largest voltage difference from exact equilibrium | 3.47e-17 V |
| Largest reconstructed conductance-matrix coefficient error | 2.46e-12 normalized |
| RC transient fixtures | 4 |
| Largest transient difference from a matrix-exponential reference | 2.11e-8 V |

Every fixture receives the canonical n(n+1)/2 acquisition settings and three additional validation settings. The source matrices, physical netlists, simulator logs, response vectors and checksums are retained. Tight DC solver tolerances isolate implementation agreement; those tiny discrepancies are numerical comparisons, not an instrument-accuracy claim.

The transient test adds a 1 microfarad capacitor at each node and steps the source current. At a 1 microvolt port-settling tolerance, observed settling times range from about 6.44 to 14.76 milliseconds in the four fixtures. Conservative eigenvalue-envelope times range from 8.14 to 16.84 milliseconds. A fixed 20 millisecond delay covers these particular RC fixtures, but a real switch, fixture and instrument combination must establish its own settling criterion. The acquisition format records that criterion explicitly.

## Calibration determines the useful inference task

The frozen interval campaign contains 24 independently generated fixtures of orders three, four and six. It includes sparse and dense structures and weak remote connections of normalized strength 1/64 and 1/512. Three base shunt strengths, 1/4, 1 and 4, are tested under five declared error budgets, giving 360 acquisitions. Each physical acquisition includes source-current error, voltage error, bounded drift and independently varying shunt paths with relay series resistance. The estimator receives only the normalized response intervals, nominal load labels and calibration envelopes.

The ideal reference discloses every graph at every load scale. With a 0.08% shunt envelope, relay resistance up to 0.05 ohm, 10 ppm current uncertainty and a 0.1 microvolt readout bound plus 0.05 microvolt drift at 1 mA, the conservative reconstruction encloses only three or four of the 24 matrices at the two smaller load scales and discloses no complete graph. A much tighter precision scenario—1 ppm shunt and current errors, 10 nanovolt readout error and 1 nanovolt drift—encloses 15 to 17 matrices, with one complete graph at base strengths 1/4 and 1. These are explicit test envelopes, rather than promised specifications of an assembled commercial rack.

This behavior identifies the engineering margin that matters. The exact inverse uses small differences between a baseline response and a loaded response. Remote weak connections can make those differences tiny. Larger shunts do not automatically solve the problem: they can increase the change while reducing the margin that separates inverse branches and amplifying calibration dependence. The three tested strengths expose that tradeoff directly. Interval enclosure, complete edge disclosure and an unresolved verdict are therefore reported separately. A matrix enclosure alone does not mean every connection has been identified.

The campaign records no false edge decisions or escaped true matrices in its accepted enclosures. Its minimum-positive-edge premise is explicit: exact absence is disclosed only under a supplied gap of 1/512. The frozen v7 pair, in contrast, has a large direct prediction separation and remains guaranteed distinguishable across the retained budgets. This makes candidate diagnosis the appropriate first hardware target while the parallel mathematical work sharpens full-reconstruction bounds.

A separately frozen qualification extension adds 72 acquisitions with 120 ppm of measured voltage plus 5 microvolts absolute error, 200 ppm source-current uncertainty, 0.1% shunt uncertainty and relay resistance up to 0.05 ohm. Full-graph inference remains unresolved in all 72. The v7 candidate pair still has disjoint predicted voltage intervals with a remaining gap of 3.326 mV at 1 mA. The table reports complete-graph disclosure out of the same 24 fixtures at each base load; an unresolved entry means the current sufficient interval margin did not establish the graph.

| Declared calibration scenario | Base load 1/4 | Base load 1 | Base load 4 |
| --- | ---: | ---: | ---: |
| Ideal reference | 24/24 | 24/24 | 24/24 |
| Precision scenario | 1/24 | 1/24 | 0/24 |
| 0.08% shunt scenario, 1 mA | 0/24 | 0/24 | 0/24 |
| Combined qualification envelope | 0/24 | 0/24 | 0/24 |

## Blind fault diagnosis and a strong conventional baseline

The first fault challenge uses twenty independent six-node fixtures and five hidden specimens per fixture: 5% edge weakening, 20% weakening, an edge opening that preserves connectivity, an unexpected bridge, and 20% grounding drift. Its hundred specimens are scored against a supplied single-fault catalog. Query selection uses candidate predictions before any fault observation. It never receives the hidden specimen's identity. Every compatibility decision uses rational arithmetic and preserves all candidates whose response predictions overlap the measurement envelope.

The readout-only comparison has 2,000 trials. At normalized absolute error 1e-5, corresponding to 1 microvolt at the declared 100 microampere and 1 mS scaling, all four complete measurement protocols identify all hundred specimens by six scalar readings; the unloaded port alone identifies 86. A second prospective challenge uses twenty new fixtures and a hundred new faults. It adds calibrated resistor intervals of plus or minus (0.0008R + 0.07 ohm) and independent relay series resistance between zero and 0.05 ohm. At 1 mA and 1 mS, normalized readout error is numerically equal to voltage error. This fault comparison fixes the calibrated source current; source-gain uncertainty is exercised separately in the instrument campaign and must be included in a bench response envelope. Across its 1,500 trials, 8,200 exact checks verify the generated signals and actual load values remain in their declared envelopes.

| Protocol with shunt uncertainty | Scalar budget | Robustly distinguished at 10 microvolts | Robustly distinguished at 100 microvolts |
| --- | ---: | ---: | ---: |
| Canonical local interventions | 6 | 99/100 | 90/100 |
| Canonical local interventions | 12 | 100/100 | 93/100 |
| Canonical local interventions | 21 | 100/100 | 93/100 |
| Predesigned informative local interventions | 6 | 100/100 | 97/100 |
| Predesigned reciprocal compliance scan | 6 | 100/100 | 100/100 |

The 10 and 100 microvolt figures bound the normalized response error at the stated current and conductance scaling. Current normalization, sensor gain and drift must fit within that response envelope; shunt and relay uncertainties are included separately.

These robust counts require disjoint candidate prediction boxes in at least one selected observation for every competing catalog member. They describe all allowed errors for each retained specimen, rather than success for a favorable single noise draw. At the wider 100 microvolt envelope, canonical acquisition happens to classify 94 specimens in the realized signals, while only 93 have this stronger separation certificate. The table uses the stronger criterion throughout.

There are no wrong forced classifications; ambiguous specimens remain unresolved. Twenty compound-fault controls outside the declared catalog are rejected as inconsistent with that catalog. These controls exercise model checking, without asserting that every possible unmodeled fault must be detectable.

The conventional scan uses reciprocal compliance entries. It may change the current-injection site and voltage-readout site, while the intervention protocol keeps both at site 1 and routes at most two shunts. Both can acquire 21 independent scalar values for a six-node symmetric operator. Equal scalar counts do not imply equal switching times, contact requirements or acquisition costs. Where the richer scanning access exists, the strong scan is the better-performing baseline in this experiment. The fixed-port method's value is its supported access contract and checkable diagnosis.

## An executable path for measured CSV data

The new bench importer accepts completed acquisition sheets and converts physical readings into the normalized intervals used by the inference backend. It validates known pad addresses, the fixed source/readout port, canonical load strengths, completeness, calibration identifiers and settling declarations. Shunt conductance uncertainty is propagated through the resistor-plus-relay path; its conservative envelope is then supplied to the existing reconstruction.

Current reversal is supported. The response is derived from the difference between positive- and negative-current voltage readings, so a constant thermoelectric offset cancels. The declared error budget still includes readout uncertainty, offset change between the two polarities and drift. The independent importer tests include a constant 100 microvolt offset, current uncertainty and relay resistance. Two complete synthetic CSV examples produce conductance enclosures containing independently solved truth. Across 118 checks, all 22 malformed, unsettled, misrouted or inconsistent input controls are rejected.

A six-node canonical experiment has 21 scalar configurations. With polarity reversal it requires 42 raw voltage samples before repeat baselines and reserved validation readings are added. The CSV format counts those samples explicitly. Repeated baseline intervals are intersected and disagreement is rejected; agreement is a consistency test, not an automatic proof that the apparatus remained time-invariant.

## The first physical qualification

Begin with two coded three-node resistor fixtures implementing the v7 twins, followed by a coded six-node fault fixture. Keep the source/readout port fixed, use two calibrated programmable shunts and route them to the labeled pads. At 1 mS conductance scaling and a 1 mA source, the normalized compliance is numerically equal to the measured voltage, making the error budget easy to audit. Reverse current, record repeat baselines and reserve new load combinations for prediction.

The primary qualification should be a correct bounded diagnosis under a predeclared calibration envelope, with conventional reciprocal scans supplying an independent reference. Record every configuration, settling interval, polarity pair and source magnitude. A second source magnitude tests the assumed linear response. The complete-graph claim should be promoted only when its own edge intervals resolve the target fixture. This sequence puts the current strongest capability—reliable discrimination—on hardware immediately, while retaining a clear route to full tomography.
