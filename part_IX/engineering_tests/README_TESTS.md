# Nonagon engineering test package

Read `ENGINEERING_REPORT.md` for the results and the first physical qualification. This folder extends Part IX with external SPICE validation, instrument-scaled calibration cases, blind fault comparisons, and measured-CSV import.

The primary calibrated result is 100/100 fresh catalog faults distinguished by six preselected fixed-port readings at a 10 microvolt readout envelope, including specified resistor and relay uncertainty. A strong reciprocal compliance scan also achieves 100/100. General full-graph certification is a stricter task and frequently remains unresolved under the tested calibration envelopes.

## Folder placement and dependencies

Keep this folder inside `part_IX/engineering_tests`, with `intervention.py`, `bounded_intervention.py`, `hardware_plan.py` and the other Part IX source files in its parent. The scripts also recognize the development workspace arrangement used to produce the saved results.

Python 3.11 or later is supported by the source. Core instrument/fault/CSV tests use Python's standard library. Fresh SPICE/transient generation requires NumPy, SciPy and an independently installed ngspice console executable. Retained simulator output is included so inspection and independent replay do not require downloading the simulator.

This run used Linux and ngspice-39. No simulator binary or Debian package is redistributed. Windows users should use an official Windows ngspice console build and supply its executable path explicitly. The commands below use the Windows Python launcher; `python3` can be used on Linux.

## Inspect and independently replay retained evidence

From this folder:

```bash
py -3.11 verify_engineering.py
py -3.11 bench_import_controls.py
```

The independent engineering verifier imports no reconstruction, interval-propagation, SPICE-harness or fault-planning functions. It uses a separate rational Gaussian solver to check saved predictions, calibration envelopes and decisions. Read its `scope` fields to distinguish arithmetic replay of retained data from generating new circuit-simulator measurements. The original Part IX exact certificate replay remains available one directory above:

```bash
py -3.11 ../replay.py
```

## Regenerate engineering campaigns

The protocols and fixture hashes are retained. Do not use `--freeze` during a confirmation replay: that option creates or amends a protocol rather than merely evaluating it.

```bash
py -3.11 instrument_campaign.py
py -3.11 instrument_qualification.py
py -3.11 fault_benchmark.py
py -3.11 fault_calibration.py
py -3.11 verify_engineering.py
```

`instrument_results.json` contains the original 360 cases; `instrument_qualification_results.json` contains a separately declared 72-case combined-error extension. `fault_results.json` covers the original 100 fault specimens and 2,000 readout-only trials. `fault_calibration_results.json` covers 100 new specimens and 1,500 trials with resistor and relay uncertainty. Each campaign records its supplied assumptions and preserves unresolved outcomes.

For a fresh external circuit-simulator run:

```bash
py -3.11 -m pip install numpy scipy
py -3.11 spice_validation.py --ngspice "C:/tools/ngspice/bin/ngspice_con.exe" --output spice_rerun
```

Replace that example path with the installed console executable. The harness invokes ngspice directly; it does not replace an unavailable simulator with its own circuit solver. The output includes actual netlists, logs, transient traces, source matrices and numeric comparison results. Solver tolerances, version and binary checksum are recorded. This is numerical circuit simulation, and the package contains zero physical hardware measurements.

## Prepare a bench acquisition

Create a blank six-node sheet:

```bash
py -3.11 bench_import.py --nodes 6 --template bench_6node_template.csv
```

The template has 21 canonical scalar configurations and uses current reversal, requiring 42 raw voltage samples before repeats. Complete every reading, calibration bound, route, settling and calibration-record field. The resistance setting is the nominal shunt; its absolute conductance error and relay series-resistance interval enclose the actual path. Positive- and negative-polarity readings support cancellation of constant thermoelectric offsets, with offset-change and drift bounds retained.

Then import the measured sheet:

```bash
py -3.11 bench_import.py --nodes 6 --input measured.csv --output measured_bounds.json --gap 1/512
```

Use `--gap` only when that minimum positive conductance is an externally justified design or manufacturing premise; omit it otherwise. `enclosed` means that the program obtained conditional matrix bounds. Consult individual edge verdicts before claiming a graph is disclosed. Calibration validity and the physical reciprocal resistor model are external premises, and arbitrary imported data are not automatically proven to have that realization.

`bench_example_exact_reference.csv` and `bench_example_uncertain.csv` are complete synthetic demonstrations, not lab records. Their retained controls verify interval normalization, current reversal, routing, settling and rejection of malformed or inconsistent input.

## Measurement counts and fair comparison

A six-node fixed-port reconstruction uses 21 scalar configurations and at most two routed shunts. A complete reciprocal compliance scan also uses 21 independent scalar entries, while changing source and readout addresses. Report raw voltage samples, source positions, relay settings, settling and total elapsed time separately. Equal scalar budgets alone establish neither equal cost nor a production-speed advantage.

The fault catalogs are explicit. A unique compatible candidate is a catalog diagnosis under the recorded error model; it is not an unrestricted unknown-network reconstruction. Compound-fault controls exercise refusal behavior outside that catalog. The full mathematical observation contract remains in the parent Part IX chapter.
