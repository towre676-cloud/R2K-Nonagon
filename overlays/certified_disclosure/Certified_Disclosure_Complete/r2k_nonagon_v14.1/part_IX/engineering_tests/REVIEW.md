# Independent engineering review

The retained records support a direct electrical implementation of the Nonagon intervention protocol. They validate the lumped resistor model with an external circuit simulator, expose how calibration limits complete graph disclosure, and support a more immediately practical product: diagnosis against a frozen catalog of faults in a known passive assembly.

## Independent verification

`verify_engineering.py` imports no tomography, bounded-inference, SPICE-harness or fault-planning implementation. Its circuit oracle uses exact forward elimination and back substitution. It replays actual load conductances, source currents, voltage normalization, uncertainty containment, accepted inverse witnesses, candidate predictions, frozen query choices, scalar and switching costs, survivor histories, and reported summary tables. Calibrated-fault results also receive independent all-error-box margin checks: a positive margin is distinguished from success under one selected simulated readout error.

`bench_import_controls.py` exercises the importer itself using separately generated voltages and explicit current reversal. A constant 100-microvolt EMF cancels in the synthetic example. Both the exact-reference and uncertain-shunt examples produce matrix enclosures containing the independent circuit truth. Twenty-two malformed, unsettled, misrouted, incomplete or inconsistent acquisitions are rejected, including under optimized Python.

## Hardware interpretation

The count n(n+1)/2 is a count of normalized scalar responses. Current reversal doubles the raw voltage samples; repeats and reserved validation readings add further acquisitions. Conventional reciprocal scans have the same number of independent compliance entries but switch source and readout among physical addresses. No measured time or production-throughput advantage is established by counting one fixed readout channel.

The acquisition and retained validation configurations use at most two simultaneous shunts. This was explicitly checked after correcting an initial validation-only mismatch in which all nodes were loaded while two-channel hardware was being counted. The correction does not alter the inference data or its scores.

The uncertainty model includes current and voltage errors, calibrated shunt conductance, series relay resistance and declared drift. Wider envelopes legitimately leave complete reconstruction unresolved. Conditional enclosure does not prove that an arbitrary specimen is linear, reciprocal, stationary, free of hidden nodes, or accurately calibrated. Repeat consistency is a useful check, not a proof of those physical assumptions.

## Qualification conclusion

There are zero physical hardware measurements in this release. ngspice is an independent numerical circuit simulator. Its agreement establishes that the intended acquisition behaves correctly in that separate implementation of the declared electrical model; it does not certify a real test rack.

The appropriate first bench trial is a calibrated three- or six-node resistor fixture with a concealed single fault, a frozen candidate catalog, two routed shunts, current reversal, and unused validation settings. Report the smallest reliably separated fault and compare with an equally informed conventional scan. Broader unknown-graph recovery should be claimed only at calibration levels where the retained enclosures actually resolve its edges.
