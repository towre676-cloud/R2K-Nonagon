# Active Spectral Tomography and Inverse Design

## The Nonagon as an experimental instrument

The next expansion of the Nonagon is a network that can be interrogated through controlled changes. Its measured response supplies poles, residues, moments and rational interpolants; its adjustable components supply experiments that expose the physical realization behind that response. The central operation is to add a known local stiffness or conductance, measure the change at one fixed port, and use the resulting resolvent identity to reconstruct the network. Once the connections are known, the same identity selects a component change that produces a desired response.

Part VII makes the need for this operation precise. A grounded path and a grounded triangle have identical driving-point response, including an example with nonproportional damping. Their common response is preserved by an orthogonal change of internal coordinates fixing the measurement port. A local intervention is tied to a physical site. Applying the same labeled intervention to both networks breaks that symmetry. The experiment therefore supplies information that an unchanged frequency sweep cannot supply.

This chapter develops that intervention into a finite acquisition theorem, a rational reconstruction program, an independent verifier, bounded-error inference, fault localization and inverse design. A programmable visual laboratory accompanies the calculation. The first industrial translation is especially direct: a passive resistor network has exactly the operator required by the theorem, and existing source-measure units, relay matrices and programmable resistor modules can perform the proposed acquisition.

## Physical model and the observation contract

Let K be the stiffness matrix of n labeled scalar displacement coordinates. An edge spring of strength w between sites i and j contributes w(e_i-e_j)(e_i-e_j)^T; a grounding spring of strength g at i contributes g e_i e_i^T. Edge and grounding strengths are nonnegative. Assume the edge graph is connected and at least one site is grounded. Then K is symmetric positive definite, its off-diagonal entries are nonpositive, and its row sums are the grounding strengths. Its inverse G=K^-1 is strictly entrywise positive.

A unit force is applied at site 1 and displacement is read at that same site. Each observation is a static equilibrium response after adding known grounding springs at selected labeled sites. The readout is a scalar, while the interventions require addressable access to every non-port site. The site set and its dimension are supplied. Unknown sites outside that set are a different inverse problem.

The dynamic extension uses D(s)=s^2 M+s C+K, with reciprocal mass and damping matrices and specified force/displacement coordinates. For a finite network, entries of D(s)^-1 are rational functions. Frequency samples, impulse moments, Loewner models and contour mode counts therefore concern the complete finite dynamic model. The production static inverse identifies K. Part VIII supplies the additional three-Laplace-point theorem for joint identification of K, C and M; the dynamic acquisition extension below makes the transient observation contract explicit. The interactive harmonic laboratory retains supplied M and C.

For an electrical implementation, replace K by the nodal conductance matrix Y. Inject a calibrated current I at site 1 and measure voltage V_1. The response is V_1/I=e_1^T Y^-1 e_1. A resistor of calibrated resistance r connected from site j to reference adds conductance 1/r to Y_jj. The theorem applies to this physical conductance matrix without a mechanical approximation.

## Exact finite acquisition

Adding a positive local load t at site j gives the rank-one resolvent identity

$$h_j(t)=h_0-\frac{t a_j^2}{1+t b_j},\qquad h_0=G_{11},\quad a_j=G_{1j}>0,\quad b_j=G_{jj}.$$

Consequently t/(h_0-h_j(t)) is affine in t, with intercept 1/a_j^2 and slope b_j/a_j^2. Two distinct positive strengths recover a_j^2 and b_j. The positive inverse fixes a_j uniquely. The implementation uses t and 2t, together with the unmodified baseline h_0.

To recover G_jk for two non-port sites, add equal positive loads t at both. Write c=G_jk, delta=h_0-h_jk(t), A=a_j^2(1+t b_k)+a_k^2(1+t b_j), and B=(1+t b_j)(1+t b_k). The two-coordinate Woodbury identity gives

$$\delta=\frac{tA-2t^2 a_j a_k c}{B-t^2c^2}.$$

This is a quadratic equation for c. Its physical root is determined, rather than guessed. The loaded network is still a connected grounded positive spring network, so both loaded cross-responses are positive. Their numerators are (1+t b_k)a_j-t c a_k and (1+t b_j)a_k-t c a_j. Multiplying these positive expressions shows delta c<a_j a_k. Thus the smaller quadratic root is the physical one:

$$c=\frac{a_j a_k-\sqrt{(a_j a_k)^2-\delta A/t+\delta^2 B/t^2}}{\delta}.$$

**Finite-load reconstruction theorem.** For the declared connected grounded network class, one baseline observation, two single-site loads at each of the n-1 non-port sites, and one simultaneous load at each pair of non-port sites determine every entry of G and hence K. The total number of scalar configurations is

$$1+2(n-1)+\binom{n-1}{2}=\frac{n(n+1)}{2}.$$

This equals the number of independent entries of a symmetric n by n matrix. Part VIII also proves minimality on the full-dimensional grounded cone for fixed and real-analytic adaptive designs; the statement here gives the implemented sufficient protocol. The protocol does not require numerical differentiation or a continuous sweep of load values. At most two shunts are active simultaneously, so a pair of programmable channels routed among known sites can implement the entire acquisition.

The theorem is over real coefficients. The portable exact implementation accepts rational observations whose reconstructed coefficients remain rational; it refuses nonrational square roots in that serialization class. That restriction belongs to the executable certificate format, rather than to the underlying uniqueness argument.

## V7 twins as a constructive starting experiment

The path has K rows (4,-1,0), (-1,3,-1), and (0,-1,2). The triangle is obtained with the Part VII internal rotation having cosine 12/13 and sine 5/13. Both have baseline compliance 5/18. A unit grounding load at site 2 gives 7/26 for the path and 1088/3989 for the triangle. The difference is 365/103714, approximately 0.00351929. If each observed loaded response has independent absolute uncertainty at most eta, those two prediction intervals are disjoint for eta<365/207428.

All six canonical scalar configurations reconstruct each complete three-site matrix. Three additional configurations are reserved for prediction, including combinations not used in extraction. The nonproportional damped twin is also retained: its unmodified quadratic-pencil responses agree, while a fixed labeled grounding change separates the responses at the tested rational Laplace settings. Those real Laplace evaluations are analytic checks, rather than a claim of sinusoidal measurements at negative or complex instrument frequencies.

## Reconstruction, certificates and uncertainty

The reconstruction function receives the declared observation document. Hidden network matrices are inputs only to the synthetic acquisition and scoring routines. The certificate contains the recovered K and G, exact extraction witnesses, the physical graph and ground strengths, and response vectors for the acquisition and validation configurations.

A separately written verifier imports no reconstruction code. It checks the acquisition contract and canonical load plan, rational types, symmetry, edge and grounding signs, connectivity, exact positive-definite LDL factors, K G=I, and every equilibrium equation. It also solves each equilibrium system independently and binds its port response to the supplied scalar observation. Extraction roots and branch inequalities are checked against the same observations. Altered calibration, missing readings and altered coefficients fail the declared exact contract.

The retained campaign reconstructs 100 fresh rational networks of orders 2 through 10, together with the two v7 twins and a damaged fixture. Across the resulting 103 certificates, 2,430 acquisition responses and 309 held-out responses agree exactly. Independent replay accepts all 103 in 119,902 explicit checks, with the same outcome under optimized Python. Separate checks verify the fault, inverse-design and damped-twin records; 29 altered actual records are rejected.

Physical measurements need uncertainty sets. The bounded route uses rational intervals and outward-enclosed square roots, propagates acquisition uncertainty through the finite-load formulas, and encloses the recovered compliance matrix. Matrix inversion uses a rational midpoint inverse and a Neumann-series bound. If the required denominator, root branch or inverse margin is not separated from its singular case, the verdict remains unresolved.

An edge is present when its conductance or stiffness interval has strictly positive lower endpoint. A zero-compatible interval does not by itself prove absence at positive error. Exact absent-edge disclosure uses an explicit manufacturing or design gap: every positive edge must have weight at least gamma. An upper endpoint below gamma then supports absence within that declared class. Current, gain, shunt calibration, drift, switching paths and repeatability must be represented in the externally justified uncertainty model before a hardware conclusion is issued.

The bounded implementation permits independently uncertain shunt strengths, including unequal pair loads. Its retained tests cover 25 random networks and 175 combinations of scalar and load uncertainty, checking 6,480 containment predicates without a false edge decision. Wider envelopes can leave the graph unresolved. Distinguishing two specified candidate circuits is easier than full reconstruction: even with independent 1% unit-shunt uncertainty, the v7 pair remains separated for additional normalized readout error below approximately 0.00171084 per reading.

The current source and readout do not become exact because their decimal values can be serialized as fractions. The hardware acquisition templates therefore distinguish recorded centers and uncertainty fields from exact synthetic proof records. The practical objective is a defensible enclosure narrow enough to decide a useful engineering question.

## Fault localization and inverse design

Reconstructing the operator before and after a change localizes component drift. For an edge spring or resistor conductance, the change has support at the two incident diagonal entries and their symmetric off-diagonal pair. Comparing recovered matrices therefore identifies the affected edge and its strength change. A disconnected region may become invisible from the chosen port; the connected-network acquisition contract detects that loss rather than reconstructing unsupported internal values.

For static inverse design, let the desired port compliance be h_target=h_0-d. The single-site formula gives the exact added load

$$t=\frac{d}{a_j^2-d b_j},\qquad 0<d<\frac{a_j^2}{b_j}.$$

The inequality is a physical feasibility condition. The largest achievable reduction at that site is a_j^2/b_j, approached as grounding tends to infinity. The program checks the feasibility margin, selects a positive load and verifies the resulting equilibrium response exactly. On a resistor fixture this changes input resistance; on a spring fixture it changes static displacement under load.

The interactive dynamic laboratory additionally scans positive grounding changes against a selected harmonic response target for a known finite M,C,K model. It displays the achieved frequency response and modal shifts. This is a numerical design experiment with declared bounds on the search, while the rational backend certifies the static target. Extending the dynamic design certificate to a physical uncertainty set is a concrete next theorem.

An acquisition companion chooses among a finite library of 21 physical rotated twins and 15 feasible load settings. Its selected intervention distinguishes the entire library in one additional reading at each retained error level. An informed fixed policy chooses the same useful load and matches that result; the gain is informative experiment selection, rather than evidence of a universal adaptive advantage. A separate Gaussian-solve verifier replays candidate predictions, policy choices and survivor histories.

## Industrial translation: equipment that already performs the operations

**Automated passive resistor-network testing is the first platform.** A grounded resistor ladder, sensor bridge or concealed passive test coupon supplies the exact nodal conductance model. A Keithley 2450 SourceMeter can supply the fixed current and voltage observation. A Pickering 40-261 provides two isolated programmable resistor channels, and a Keysight 34980A mainframe with 34932A matrix modules can route the shunts to addressed pads. Two channels suffice because the canonical acquisition uses only single and pair interventions.

Pickering documents external-DMM four-terminal calibration. Its model ranges cover the moderate resistance values suitable for an initial fixture. The plan supplied here uses a nominal conductance unit of 1 mS, a nominal source current of 100 microampere, and shunt settings of 1000 and 500 ohm. These are acquisition-design values, not a recommendation to purchase a particular configuration without matching the fixture and calibration budget. The CSV plan records switching, current, readout, shunt calibration, settling and repeatability fields.

Existing in-circuit testing and scanned multi-port measurements are the comparison. Wiring a relay matrix to every pad also permits conventional scans. The one-readout theorem does not establish an automatic production-time advantage. The proposed value is an integrated sequence of informative interventions, reconstruction, bounded component decisions, held-out prediction and independently replayable evidence. An initial bench comparison should report contact count, configuration count, settling time, total test time and the smallest reliably diagnosed defect.

**Resistive sensor and hardware-in-the-loop fixtures are a second exact application.** Programmable resistance already emulates temperature, strain, pressure and load sensors. The theorem can identify the passive bridge or harness surrounding a controller and then select a calibrated perturbation that exercises a desired controller response. Identifying an arbitrary active controller is outside the passive resistor model.

**Passive RC/RLC fixtures provide the frequency-domain extension.** Existing impedance analyzers, including the Keysight E4990A family, acquire response spectra. Adding capacitance and inductance gives a declared finite dynamic model whose modes, residues and moments are handled by the Nonagon. Fixture parasitics, shunt behavior and the chosen network model determine the usable frequency band. The analyzer's headline bandwidth does not establish the validity of a lumped model across that entire band.

**Modal-testing rigs provide the first mechanical translation.** HBK supplies force transducers, impedance heads, accelerometers and impact hammers for vibration testing. A deliberately lumped spring frame with calibrated adjustable grounding can implement the theorem. A flexible machine requires a measured reduction and a demonstrated relation between the physical adjustment and the reduced operator. This platform connects the static reconstruction to wave motion, resonance suppression and mechanical fault localization.

**Active floor isolation and magnetic-bearing machinery are later integrations.** Commercial systems already combine vibration sensing and controlled forces, and SKF documents industrial magnetic-bearing blower applications. Their existence establishes a relevant engineering market. Their control interfaces must still be mapped to the intervention model. Feedback dynamics, gyroscopic terms and operating-point changes call for an expanded reciprocal or MIMO model; they are not covered simply by relabeling K.

## The number-variance offset: isolate the common numerical operation

The three archived height windows have different log(t/(2 pi)) values but a similar large-L discrepancy. That is a sound reason to inspect operations shared across the runs. A term proportional to log(t/(2 pi))^-p changes by a factor (9.4/13.6)^p across the stated range; at p=1 the decrease is about 31%, and at p=2 about 52%. Similar observed offsets do not exclude mixtures or cancellation, but a shared numerical origin can be tested directly.

The archive starts its R2 quadrature at epsilon=10^-4. Number variance is L-2 times the integral from 0 to L of (L-r)(1-R2(r)). Since the regular pair correlation vanishes at zero, the omitted initial interval contributes approximately 2 L epsilon-epsilon^2. Omitting it raises the predicted variance. At L=16 this is approximately 0.00319999, the same size as the reported discrepancy.

Refining a grid while keeping the same positive lower endpoint tests quadrature resolution on the retained interval. It does not test the missing interval. The independent sine-kernel control exposes this directly: the archived integration gives approximately 0.63011635 at L=16, while the zero-inclusive reference is approximately 0.62693320. The difference, about 0.00318315, consists predominantly of the omitted endpoint contribution with a smaller discretization term.

An endpoint-only correction moves the three L=16 height-averaged predictions to approximately 0.32732022, 0.37153283 and 0.38486160. The corresponding observed values are approximately 0.32726574, 0.37174794 and 0.38499248. Using the archived block standard errors, residuals are about -0.10, +0.27 and +0.23. The previous large height-stable discrepancy is therefore accounted for by a specific shared numerical operation. These values are endpoint corrections to archived predictions; they are not a fresh full evaluation of the prime-product formula.

The supplied audit separates missing-interval bias, retained-grid error, exact upper-endpoint handling and unfolding conventions. The unfolding's additive 7/8 term cancels after subtracting the first unfolded coordinate and cannot produce the discrepancy. A full change of unfolding requires the underlying zero positions and is not reconstructed from a summary table. The immediate repair is zero-inclusive integration with the removable R2 endpoint handled analytically, exact inclusion of L, and a sine-kernel control at every requested L.

## The resulting program

The Nonagon now has an experimental action aligned with its spectral representations. Known local changes expose a labeled physical operator; the operator predicts unused experiments; a verified change achieves a specified response. A passive resistor fixture offers the shortest route to physical validation, while the resonator laboratory makes the same mathematics visible as wave motion and modal design. The number-variance repair follows the same principle: make the observation and integration contract explicit, then check the finite operation that carries the claimed conclusion.

## Primary equipment sources

1. Pickering, programmable resistor family and 40-261 datasheet: https://www.pickeringtest.com/products/pxi/pxi-simulation-modules/sensor-simulation-programmable-resistors ; https://www.pickeringtest.com/portal/web/256/content/downloads/datasheets/40-261D.pdf
2. Keysight 34980A and 34932A: https://www.keysight.com/us/en/product/34980A/multifunction-switch-measure-unit.html ; https://www.keysight.com/us/en/product/34932A/dual-4x16-armature-matrix.html
3. Keithley 2450 product family: https://www.tek.com/en/products/keithley/source-measure-units/2400-graphical-series-sourcemeter
4. NI PXIe-4081: https://www.ni.com/en/shop/hardware/digital-multimeters/model-pxie-4081
5. Keysight E4990A: https://www.keysight.com/us/en/product/E4990A/impedance-analyzer-20-hz-10-20-30-50-120-mhz.html
6. HBK vibration transducers: https://www.hbkworld.com/en/products/transducers/vibration
7. SKF industrial magnetic-bearing blower application: https://evolution.skf.com/us/magnetic-approach-to-improve-the-efficiency-of-watertreatment-plants/


## Mathematical and engineering integration

Part VIII is the canonical mathematical treatment: finite acquisition and its qualified minimum count, the Kron boundary for inaccessible sites, a Gaussian two-point resolution bound, and joint K, C and M recovery at three nonnegative real Laplace arguments. This chapter specifies the production observation documents and acceptance gate. The separate reference inverse in Part VIII agrees exactly on all 103 retained production matrices, using the same 2,430 acquisition readings and 309 held-out responses. The production checker remains the acceptance gate before the reference oracle runs.

The new `load_planning` extension chooses among calibrated finite common strengths using measured pilots and the shared-baseline error model. Its records include pilot overhead. The `dynamic_acquisition` extension connects finite transient observations to the real-Laplace theorem with declared error bounds. Fresh campaign results and replay instructions are provided in those folders. Archived qualification records retain their original fixtures and declared error budgets.
