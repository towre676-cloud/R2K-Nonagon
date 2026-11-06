# Instrument-scaled electrical acquisition tests

The campaign tests the existing Part IX interval reconstruction against a frozen set of 24 resistor-network fixtures. It covers orders 3, 4 and 6, with eight fixtures at each order: sparse and more strongly connected networks, followed by two remote weak-link cases. Positive conductances are at least 1/512 in normalized units. All node addresses are supplied; the connections and values are withheld from the inference function.

The physical scaling is 1 mS per normalized conductance unit. Source current is either 100 microampere or 1 milliampere. The acquisition applies two individual shunt strengths and prescribed pair loads, repeats the unmodified baseline after each five nonbaseline settings, and reserves three additional load combinations. The baseline interval supplied to inference is the intersection of the repeated baseline intervals. This includes the repeats in the measurement budget rather than treating the theorem's canonical configuration count as total bench effort.

The protocol and fixture identities were saved before generating any measurement packets. The interval algorithm was held fixed throughout the campaign. Every packet retains nominal and actual shunts, their series relay resistance, actual source current, ideal and measured voltage, additive readout and drift errors, and the scalar interval delivered to inference. Hidden circuit matrices are used only by acquisition and scoring.

## Calibration model

Each shunt's conductance varies independently inside the declared relative calibration envelope. A series relay resistance, independently between zero and the stated maximum, modifies its effective load. If the nominal normalized load is t, the conductance unit is g, the shunt relative bound is e, and maximum relay resistance is r, its actual normalized load lies between

    t(1-e) / [1+r t g(1-e)]  and  t(1+e).

The global relative load envelope passed to the existing inference code covers both single-load strengths and independently unequal pair loads. Relay resistance is therefore represented in the added branch, rather than being silently omitted or assigned to an unrelated voltage error.

Voltage error and bounded additive drift enter every readout interval. The source current can vary independently inside its declared relative envelope. Normalized response intervals enclose voltage divided by actual current, multiplied by the conductance unit. A repeated baseline narrows the shared baseline interval without assuming that source-current or readout errors average to zero.

These are synthetic, declared engineering scenarios. Their names do not assert a complete manufacturer accuracy specification or measured hardware performance. No physical apparatus was connected.

## Results

There are 360 acquisitions: 24 fixtures, five calibration scenarios, and three base-load strengths (1/4, 1 and 4). The ideal reference recovers every graph at every load scale. Under nonzero uncertainty, matrix enclosure and complete graph disclosure are separate outcomes. A matrix enclosure may still leave all or some individual edges unresolved.

The following table shows the frozen base-load-1 comparison. Each cell concerns the same 24 fixtures and 192 possible edges.

| Declared scenario | Matrix enclosed | Full graph disclosed | Edges disclosed correctly |
|---|---:|---:|---:|
| Ideal reference | 24/24 | 24/24 | 192/192 |
| 1 mA, 1 ppm source and shunt calibration, 10 nV readout, 1 nV drift, relay at most 0.005 ohm | 17/24 | 1/24 | 39/192 |
| 1 mA, 10 ppm source, 0.08% shunts, 0.1 microvolt readout, 0.05 microvolt drift, relay at most 0.05 ohm | 4/24 | 0/24 | 0/192 |
| 100 microampere, 100 ppm source, 0.08% shunts, 1 microvolt readout, 0.1 microvolt drift, relay at most 0.05 ohm | 0/24 | 0/24 | 0/192 |
| 100 microampere, 100 ppm source, 1% shunts, 10 microvolt readout, 1 microvolt drift, relay at most 0.05 ohm | 0/24 | 0/24 | 0/192 |

Across all load scales, the precision-calibration scenario encloses 15 to 17 matrices and discloses one full graph at base load 1/4 or 1. The 0.08% shunt, 1 mA scenario encloses three or four matrices at those smaller loads and none at load 4. Increasing the load is therefore not an automatic cure for this acquisition. The six weak-link fixtures remain unresolved under the nonideal base-load-1 precision scenario, while 17 of the other 18 matrices are enclosed.

All 5,893 accepted matrix and edge containment predicates pass. No matrix truth escapes its returned enclosure, and no disclosed edge verdict is false. The disconnected control returns unresolved because the inaccessible subnetwork produces no positive local response differences. An unresolved result here means that the current sufficient bounds could not disclose the requested graph; it is not an information-theoretic impossibility statement.

The order breakdown at base load 1 shows where the available margins disappear:

| Scenario | Order 3: enclosed / full graph | Order 4: enclosed / full graph | Order 6: enclosed / full graph |
|---|---:|---:|---:|
| Precision calibration | 6/8 / 1/8 | 6/8 / 0/8 | 5/8 / 0/8 |
| 0.08% shunts, 1 mA | 4/8 / 0/8 | 0/8 / 0/8 | 0/8 / 0/8 |

Among the 232 unresolved outcomes over the original 360-case campaign, 119 fail to separate the pair-root discriminant from zero, 50 fail to separate a local response difference, 44 fail the matrix-inverse perturbation margin, and 19 fail the single-load intercept margin. These diagnostics identify the actual arithmetic stages requiring a better acquisition or sharper uncertainty treatment.

## Combined qualification budget

A separate protocol was frozen before a 72-case extension using 120 ppm of the measured voltage plus 5 microvolt absolute readout error, 200 ppm source-current uncertainty, 0.1% shunt-conductance uncertainty, and relay series resistance at most 0.05 ohm. These are a declared combined engineering budget, not a manufacturer-guaranteed combined specification. The simulation varies actual source current and actual shunts independently; each retained voltage packet includes its gain-error term and its measured-voltage-dependent envelope.

At 1 mA, all 24 fixtures remain unresolved at each of the three load scales under the current full-network interval algorithm. All 2,424 normalization and effective-load checks pass. This qualification test does not support full-graph certification at that budget.

The same combined budget still separates the two frozen v7 candidate circuits. After all declared uncertainty is included, the remaining gap between their loaded-voltage intervals is **3.32562 mV at 1 mA**, or **0.323562 mV at 100 microampere**. This is the strongest immediate hardware acceptance target supported by the current acquisition: robust discrimination of specified circuits, followed by increasingly ambitious fault and graph identification tests.

## Robust candidate discrimination

The v7 path and triangle provide a practically easier task than unknown full-matrix reconstruction. For each calibration budget, the test encloses their physical loaded-voltage predictions over all declared source-current, shunt, relay, readout and drift bounds. The intervals are disjoint in all five scenarios.

| Scenario | Remaining separation between voltage intervals |
|---|---:|
| Ideal reference, 1 mA | 3.51929 mV |
| Precision calibration, 1 mA | 3.51867 mV |
| 0.08% shunts, 1 mA | 3.50529 mV |
| 0.08% shunts, 100 microampere | 0.343481 mV |
| Broad 1% shunts, 100 microampere | 0.314701 mV |

This supports an immediate experiment that distinguishes two specified concealed circuits under useful calibration envelopes. It does not transfer that error tolerance to recovery of every unknown connection.

## Acquisition burden and next bench decision

| Network order | Canonical readings | Additional baseline repeats | Reserved validation readings | Total voltage reads |
|---:|---:|---:|---:|---:|
| 3 | 6 | 1 | 3 | 10 |
| 4 | 10 | 2 | 3 | 15 |
| 6 | 21 | 4 | 3 | 28 |

The first bench experiment should therefore test robust candidate or fault discrimination before attempting a universal full-graph promise. A six-node reconstruction experiment should report both point prediction error and conditional edge disclosure, alongside a conventional scanned measurement baseline. Sharper treatment of correlated current and baseline uncertainty, additional observation ports, and stronger experimental design are concrete routes to improving graph disclosure. Any tuned method should be confirmed on fresh fixtures rather than rescored as a new protected result on this frozen set.

## Reproduction

The fixture and protocol files are already frozen. Run:

    python instrument_campaign.py

The script writes `instrument_results.json` and `instrument_case_records.jsonl`. The latter contains all 360 complete observation and interval-result packets. `instrument_protocol.json` identifies the fixture set, units, error budgets, load strengths and scoring rules. The result records SHA-256 identities for the protocol, campaign source and observation packets.

Run the separate, already frozen qualification extension with:

    python instrument_qualification.py

It writes `instrument_qualification_results.json` and 72 complete packets in `instrument_qualification_records.jsonl`. An explicit protocol amendment changes the three reserved validation settings to use at most two shunts; it leaves the canonical acquisitions, fixture set, calibration scenarios and primary inference scores unchanged. The original protocol identity is retained for traceability.
