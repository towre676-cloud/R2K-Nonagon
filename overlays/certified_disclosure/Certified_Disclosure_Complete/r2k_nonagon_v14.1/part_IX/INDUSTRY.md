# Existing industrial equipment for Active Nonagon tomography

Research date: 30 October 2026. The hardware mapping below uses manufacturer documentation. Detailed specifications should be confirmed for the selected module and configuration.

## Recommendation

Start with low-voltage, low-power automated testing of connected grounded resistor networks: resistor ladders, bridges, sensor-interface assemblies, and deliberately concealed passive network fixtures with known test-point addresses. This has an exact circuit model rather than a mechanical analogy. Let G be the nodal conductance matrix, with conductances as positive edge weights, nonpositive off-diagonals, and grounding conductance given by row sums. Inject calibrated current I at fixed node b and measure its voltage V_b. Then h=V_b/I=e_b^T G^{-1}e_b. A known shunt resistor r at node i adds t_i=1/r to G_ii. The Nonagon intervention theorem applies verbatim.

Two independently programmable, isolated resistor channels are enough for pairwise interventions when a relay matrix routes them to the selected node pair. All nodes need known electrical test-point access for intervention; they do not require separate simultaneous voltage readout. This should be described as one measurement channel with many addressable interventions, not arbitrary remote hidden-circuit reconstruction.

The mathematically exact model does not make hardware readings exact. The first hardware deliverable should identify a declared passive resistor fixture with calibrated shunts and bounded measurement errors, then predict reserved load settings. For general populated PCBs, active semiconductors, unknown return paths, temperature dependence, and nonlinear protection circuits require a new model. Industrial in-circuit testing already exists; the proposal is an integrated reconstruction-and-certificate workflow, not invention of PCB test.

## Actual purchasable equipment mapping

| Role | Existing equipment | Verified capability | Proposed connection |
|---|---|---|---|
| Two isolated local shunts | Pickering 40-261 PXI precision programmable resistor | Two isolated channels. Model 001: 1.5 ohm–2.9 kohm, less than 2 milliohm setting resolution. Model 002: 10 ohm–36 kohm, less than 15 milliohm setting resolution. External-DMM four-terminal calibration; Python/Windows/Linux support. Datasheet specifies typical 3 ms operation and 100 mW maximum power. | Route one channel from node i to circuit reference and the other from node j to reference. Convert calibrated resistance to conductance. |
| Automatic routing | Keysight 34980A plus 34932A | 34980A supports an internal 6.5-digit DMM, LAN/USB/GPIB. 34932A provides 128 two-wire crosspoints in dual 4-by-16, 8-by-16, or 4-by-32 arrangements. | Select test-point pairs without moving the fixed observation port. Wiring topology and relay exclusion rules must prevent unintended node shorts. |
| Fixed force-and-observation analogue | Keithley 2450 SourceMeter SMU | Official series page explicitly identifies current-source, voltage-source, DMM and measurement roles, with a specific 2450 product section. | Source a small current at b and measure V_b. Four-wire sensing at the fixture removes avoidable lead drop; calibration still includes fixture leakage and contact effects. |
| Alternative high-resolution meter | NI PXIe-4081 | 7.5-digit DMM; AC/DC voltage and current; two-/four-wire resistance. | Measure fixed-port response and calibrate shunts. A separate calibrated source is required when using it as a voltage meter. |
| Frequency-domain extension | Keysight E4990A | Official page identifies a 20 Hz–120 MHz impedance analyzer family. | Add known RC/RLC fixture dynamics after the DC reconstruction has passed; choose fixture/frequency limits compatible with relay and shunt parasitics. |

Do not confuse display resolution with absolute calibration accuracy. Pickering's 40-261 datasheet gives ±0.08% ±70 milliohm accuracy over its stated temperature range; this differs from the best-case values on the broader family page.

## Ranked equipment/application categories

1. **Automated passive resistor-network test and diagnosis — first target, exact fit.** Existing relay-matrix/SMU/PXI equipment supplies both addressable shunts and the one fixed electrical measurement. Tasks: recover resistor values, distinguish an unexpected bridge or leakage path, and identify drift while the network remains in the connected grounded class. This is the strongest place to make the full theorem operational.

2. **Resistive sensor bridges and sensor-simulation HIL fixtures — exact fit for the passive subassembly.** Pickering explicitly supplies programmable-resistance simulation for temperature, strain, pressure, load and related sensing applications. Tailor to a known-address bridge/harness or resistive network around a controller; do not claim the same theorem identifies an arbitrary active controller. The same hardware can exercise a controller with a diagnosed and calibrated sensor fixture.

3. **Passive RC/RLC filter and impedance test fixtures — strong second phase.** Existing impedance analyzers supply calibrated frequency response; programmable resistors give the interventions. DC identifies the conductance part when it remains positive definite. Dynamic capacitance/inductance inference needs its own declared pencil model and multi-frequency reconstruction. Frequency settings are limited by actual fixture parasitics, not the analyzer's headline maximum.

4. **Modal-testing rigs for spring-coupled structures — suitable controlled mechanical validation.** HBK currently lists force transducers, impedance heads, accelerometers and impact hammers for modal/structural testing. Those acquire the force/response channels. A dedicated lumped frame with calibrated adjustable springs can match K+t_i e_i e_i^T; an arbitrary flexible machine needs a measured reduced model and update localization. This is a good visual demonstration after the resistor experiment.

5. **Active floor-vibration isolation systems — real industry, later integration.** TMC's manufacturer search result states that STACIS uses piezoelectric actuators and measured floor activity for real-time vibration cancellation. The existing product therefore confirms a practical market for adaptive wave/vibration control. It does not establish that its interface exposes independently calibrated additive stiffness knobs, so the exact local-shunt theorem should not be promised for it without interface/model work.

6. **Magnetic-bearing industrial rotors and turbo blowers — high-value future adaptation.** SKF's manufacturer article documents active magnetic bearing stiffness/damping control and an actual wastewater blower application. It is dated 2014 and proves existing industrial deployment, not a fresh 2026 product catalogue. Rotating systems include feedback, gyroscopic terms and operating-point dependence, so they need a rotor-dynamic/MIMO extension rather than the grounded static theorem unchanged.

## Proposed first experiment

Use a six-node passive resistor fixture, one grounded reference, and node values in a moderate kilohm range so relay/lead resistance is a small calibrated contribution. Freeze node addresses, allowed edge-weight range, minimum positive conductance if absent-edge certification is requested, current level and error envelope. Keep the underlying resistor wiring concealed from the estimator.

Measure the original port response. Apply two calibrated shunt strengths individually at each node; apply selected pair loads with the two isolated channels. Reconstruct G with the finite-load formulas implemented in the Nonagon package. Use additional shunt combinations and different source amplitudes only for held-out prediction and linearity checks. Compare against direct multi-port voltage scans as an independent physical reference. Diagnose graph edges only when the conductance intervals support them.

The corresponding dynamic demonstration can then add grounded capacitors and frequency sweeps, so the same object acquires poles, residues, moments, Loewner models and contour certificates. This makes the resistor experiment the physical front end for the entire Nonagon rather than an isolated DC trick.

## Primary sources

- Pickering programmable resistor family, full text inspected: https://www.pickeringtest.com/products/pxi/pxi-simulation-modules/sensor-simulation-programmable-resistors
- Pickering 40-261 datasheet, full PDF text inspected: https://www.pickeringtest.com/portal/web/256/content/downloads/datasheets/40-261D.pdf
- Keysight 34980A, directly fetched HTTP 200 and inspected: https://www.keysight.com/us/en/product/34980A/multifunction-switch-measure-unit.html
- Keysight 34932A, directly fetched HTTP 200 and inspected: https://www.keysight.com/us/en/product/34932A/dual-4x16-armature-matrix.html
- Keithley 2400 Graphical Series/2450, directly fetched HTTP 200 and inspected: https://www.tek.com/en/products/keithley/source-measure-units/2400-graphical-series-sourcemeter
- NI PXIe-4081, full text inspected: https://www.ni.com/en/shop/hardware/digital-multimeters/model-pxie-4081
- Keysight E4990A, directly fetched HTTP 200 and inspected: https://www.keysight.com/us/en/product/E4990A/impedance-analyzer-20-hz-10-20-30-50-120-mhz.html
- HBK vibration transducers, directly fetched HTTP 200 and inspected: https://www.hbkworld.com/en/products/transducers/vibration
- TMC STACIS III, official search excerpt only: https://www.techmfg.com/products/stacis/stacisiii
- SKF industrial magnetic bearings, directly fetched HTTP 200 and inspected: https://evolution.skf.com/us/magnetic-approach-to-improve-the-efficiency-of-watertreatment-plants/

