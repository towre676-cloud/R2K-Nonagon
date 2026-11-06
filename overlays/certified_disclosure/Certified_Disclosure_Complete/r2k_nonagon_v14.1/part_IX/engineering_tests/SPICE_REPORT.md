# Independent circuit-simulator validation

The electrical Nonagon acquisition has now been exercised in **ngspice 39**, independently of the exact rational forward solver. Fourteen passive resistor fixtures produce 286 DC operating-point records. The suite includes the v7 path and triangle and twelve fresh connected grounded networks of orders three through eight. Every fixture is interrogated with the complete canonical single/pair shunt protocol and three additional load settings reserved for prediction.

## Results

| Quantity | Result |
|---|---:|
| Independent resistor-network fixtures | 14 |
| Canonical and held-out DC configurations | 286 |
| Largest ngspice versus exact node-voltage difference | 3.47 × 10^-17 V |
| Largest recovered normalized conductance-coefficient error | 2.46 × 10^-12 |
| Independent RC transient fixtures | 4 |
| Largest transient error versus matrix-exponential response | 2.11 × 10^-8 V |

The recovered matrices are obtained by a separate floating-point inverse that receives only ngspice's scalar port readings and public load labels. They are compared with the hidden rational matrix only afterward. This tests the complete chain from resistor netlist through a separate circuit simulator to observation-only reconstruction. Every acquisition and held-out configuration uses at most two simultaneous shunts, preserving the proposed two-channel routed bench implementation. The archived records include all node voltages, so a separate reviewer can also check Kirchhoff balance without using the reconstruction routine.

## Physical scale and circuit construction

The normalized conductance matrix K is converted to physical nodal conductance Y=(1 mS)K. Each nonzero edge conductance becomes an ordinary resistor between its two nodes, and each grounding conductance becomes a resistor to the reference node. A calibrated-current analogue of 100 microampere is applied at node 1; voltage is observed at that same fixed port. The normalized response is h=V_1(1 mS)/(100 microampere). A unit intervention is therefore a 1000-ohm shunt, and a double-strength intervention is a 500-ohm shunt.

For the v7 twins, both original responses are 27.7777778 mV. The unit node-2 shunt produces 26.9230769 mV for the path and 27.2750063 mV for the triangle. The approximately 351.929-microvolt separation is a circuit-simulator result at the declared source current. It is not a voltmeter accuracy specification.

The DC netlists specify relative tolerance 10^-12, current tolerance 10^-15 A, voltage tolerance 10^-12 V, and gmin=10^-18 S. They contain only independent sources and passive R/C elements. Optional XSPICE or OSDI models are not loaded. A supplied minimal startup file and the `-n` option make the test independent of a user's local initialization. Exact rational values are converted to 17-digit resistor values at the simulator boundary; that conversion is included in the reported numerical agreement.

## Settling and acquisition timing

The transient fixtures add one microfarad from each labeled node to reference, begin with zero capacitor voltage, and apply the DC source using ngspice's `uic` transient initialization. The independent reference solves C v'(t)+Y v(t)=I e_1 by a symmetric matrix exponential. Selected reference values are also checked with SciPy's separately implemented `expm`.

| Fixture | Port settling within 1 microvolt | Slowest RC time constant | All-node spectral-envelope time |
|---|---:|---:|---:|
| V7 path | 6.441 ms | 0.789 ms | 8.142 ms |
| V7 triangle | 6.445 ms | 0.789 ms | 8.142 ms |
| Fresh six-node fixture | 14.765 ms | 1.569 ms | 16.837 ms |
| Fresh eight-node fixture | 13.386 ms | 1.386 ms | 14.938 ms |

The envelope time uses the slowest eigenvalue and the full steady-voltage norm, so it bounds every node of the stated ideal model. The simulator settling time concerns the observed port. These times establish a useful millisecond scale for these particular R/C fixtures. Relay operation, actual fixture capacitance, meter integration time, thermal effects and calibrated uncertainty must be added to a physical acquisition budget.

## Reproduction and retained evidence

Run `python spice_validation.py --ngspice PATH_TO_NGSPICE`. The executable may be a Windows console build such as `ngspice_con.exe`. NumPy and SciPy are required by the independent reference and floating reconstruction. `--output DIRECTORY` selects the retained netlists and logs. `--library-path DIRECTORY` is available for a locally extracted Linux executable whose shared libraries are not on the system search path.

The test used the official Debian package `ngspice_39.3+ds-1_amd64.deb`, with one official `libxft2` runtime dependency. The package download URL and SHA-256 are retained as provenance; redistributing those binaries is unnecessary. The report JSON also records the actual executable hash and version banner. `netlists/` contains every circuit and raw simulator log, plus the four transient traces. `spice_results.json` binds those records to the exact fixture matrices, physical units and pass criteria.

This is independent **numerical circuit validation**. No physical apparatus was connected and no simulated voltage has been presented as a laboratory measurement.
