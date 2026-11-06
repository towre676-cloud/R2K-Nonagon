# Independent internal review of the frontier modules

This review inspects protocols, source and finite receipts. It does not establish novelty, independently prove every analytic enclosure, or turn a locally frozen simulation into externally preregistered research. The earlier `nonagon_challenge` campaign is left unchanged.

## Early protocol review

### Adaptive acquisition

The implementation uses a fixed 1,025-frequency candidate lattice and a separate 1,024-point clean scoring lattice. The final repaired protocol uses a 33-point pilot derived as a subset of the 64-point log-uniform baseline; it counts toward the 64 total unique positive-frequency queries. Every strategy accesses a frozen candidate observation table, so querying the same candidate in different strategies returns the same noisy value. There are no repeated queries within a strategy, and inferred conjugates do not count as new observations. Earlier 17/33-pilot exploration was superseded after the generator defect described below; it must not supply confirmatory evidence.

The designer fits only observed samples and is not supplied clean responses, true poles, masses, noise fraction or scenario. Both adaptive candidates use the same final ordinary Loewner plus vector-fitting model as the fixed sampling baselines. I-optimal covariance is explicitly a design surrogate for a linearized complex model, not a calibrated confidence interval. Development selects one algorithm before fresh confirmation seeds; final freeze must bind the selected protocol, acquisition code and copied reference fitter.

The inspected I-optimal expression computes a surrogate integrated predictive-variance reduction with model-estimated relative-noise shape. Its Hermitian update is consistent with that surrogate. Physical real/conjugate parameter constraints and colored-noise covariance are not fully modeled; no statistical calibration claim should be attached to it. API isolation is not a security sandbox.

### External impedance data

The protocol fixes every cell ID before scoring, reserves IDs01–06 for development and IDs07–54 for confirmation, and fixes frequency splits, candidate model families and validation selection. Training/validation/test roles are explicit; no validation refit is performed. Positive relaxation fitting is an existing passive impedance prior, not a uniquely Nonagon algorithm or a proof that a physical cell exactly has that circuit.

Requested before confirmation: accompany any common-success comparison with a prespecified all-device failure-penalty analysis, since excluding method failures can bias the ratio. Cell bootstrap cannot account for shared batch/instrument systematic uncertainty; that limitation is already in the protocol.

### Exact one-port ambiguity

An orthogonal similarity fixing an observed port preserves its transfer function. For a mechanical-network conclusion, both stiffness matrices must additionally have nonpositive offdiagonals and nonnegative grounding row sums, with positive masses; unrestricted orthogonal similarity alone does not establish this. Extra-port measurements can distinguish a constructed pair without proving that two ports universally identify every network.

### Independent direct contour engine

The new direct integration route should guard the two-term zeta derivative series context, explicitly exclude boundary zeros, require finite moment enclosures with isolated and contained integer counts, and preserve every unresolved/failing tile in a complete coverage ledger. It shares Arb with the repaired Taylor route and must not be described as a second independent arithmetic library. A finite full-strip campaign is a local finite-range theorem, not a global RH proof.

The completed-source and receipt reviews follow below.

## Material generator discovery and independent reproduction

The first adaptive confirmation was invalidated after a noiseless sanity check exposed inaccurate physical responses from the automatic linear-solve path in the installed SciPy 1.17.0 environment. Its initially favorable numerical ratio is quarantined and must not be reported as a valid scientific gain.

`solver_reproducer.py` contains a literal 3x3 complex-symmetric, non-Hermitian matrix and right-hand side, with no generator or fitter imports. On this well-conditioned input (condition number about 14.95), this review independently obtained an absolute residual of approximately `1.7123e-3` from SciPy 1.17.0's automatic path and `4.1541e-16` from its explicit general solver. NumPy's general solve agrees with the latter. On the same literal input under the installed SciPy 1.18.1 profile, the automatic path returns a residual about `3.3307e-16`. The versioned JSON receipts record full inputs, solutions, versions and residuals. This is a particular reproducible numerical-path defect, not a blanket assertion about all automatic solvers.

The repaired generator uses explicit general solves and checks backward residuals plus agreement with NumPy's solver for every physical candidate and held-out response. Fresh development and confirmation seeds are required and are used. The repaired noiseless development sanity check is important evidence that the generated rational family is consistent before testing adaptive sampling.

The prior challenge is unchanged. Its ordinary physical generator already used an explicit general solver. This review separately audited its only default-auto branch: 3,210 solves for the ten four-mass near-collision confirmation systems, including training and test frequencies. Maximum scaled backward residual was `2.005e-16`; the observed 3x3 defect was not found there. See `prior_challenge_four_mass_audit.json` for bounded scope, source hash and details.

## Completed exact-network and passivity reviews

The path/triangle example is a valid physical construction: both matrices have unit positive masses, negative existing inter-mass offdiagonals, strictly positive grounding row sums, and positive-definite stiffness. The rational orthogonal rotation fixes the observed coordinates, so the all-frequency transfer identity follows algebraically, rather than from the three exact numerical shifts checked by the test suite. Different edge counts establish nonisomorphic graphs. Nonzero Krylov determinants rule out invisible modes. The stated all-dimension grounding bound `g>a/2` is sufficient for the displayed transformed row sums and spring signs. The added static sensor separates these constructed pairs; the stronger full-output, known-mass reconstruction formula has explicit cyclicity and nonresonance assumptions. No universal two-port uniqueness is implied.

The passivity compiler and separate verifier use exact decimal-rational signs and time constants for a restricted impedance representation. Nonnegative RC relaxation coefficients, elementary R/L/inverse-C terms, and optional nonnegative RL relaxation coefficients genuinely imply analyticity and nonnegative real part on the open right half-plane. The series composition of parallel RC or RL blocks gives the stated algebraic circuit realization. The optional RL extension was introduced after measured confirmation and has no new fitted predictive score; its theorem and monotonicity statements are sound in the stated class. The guarantee is about serialized models, not physical batteries or unique circuit interpretation. All 216 exported fits comprise 24 development and 192 confirmation models, not 216 separate devices.

## Completed external-data review

The final external protocol SHA-256 `db39d7009bc0cf56d4beca4d901d2728ee00436044b81bad708a260c18daf529` matches the preexecution protocol communicated during review and both result receipts. The scientific source hash also matches development and confirmation receipts. Confirmation comprises 48 cells, four tasks per cell, 768 selected method results and 4,608 candidate fits, with no candidate or selected failures. The all-device capped-failure analysis and common-success coverage agree because no failures occurred. Bootstrap sampling keeps all tasks for a cell together in the descriptive aggregate; per-task sampling is also by cell.

The positive relaxation model loses to the vector-fitting baseline in the aggregate (geometric error ratio about 3.48). This does not undermine its model passivity guarantee; it shows that predictive fit and certified model structure are different objectives. The external results establish no unique Nonagon predictive advantage. Intervals remain conditional on a single batch, state of charge and instrument; multiple descriptive contrasts are unadjusted.

## Final repaired adaptive confirmation and policy check

This review independently recomputed the final adaptive primary statistic from the non-quarantined raw rows. The fresh confirmation seed is `27193693`; all three frozen scientific file hashes match. There are 1,000 method/case rows, including 150 positive-noise cases in 50 independent system clusters. Adaptive disagreement divided by the fixed log-grid geometric error is **0.583096**, with the prespecified paired cluster-bootstrap 95% interval **[0.463073, 0.725213]**. Arithmetic mean NRMSE is approximately `0.035964` versus `0.056861`, 90th-percentile error `0.115995` versus `0.130044`, and the adaptive strategy wins 78% of positive-noise pairs. This is valid evidence for the specified known-order simulated family at a fixed observation budget, not hardware validation, a universal gain, or a novelty finding. Additional computation is a material cost.

The final policy uses a pilot-estimated residue-ratio rule selected only on repaired development data. Every final 64-point observation set contains all 33 pilot locations, so reading pilot values for the gate does not provide extra physical observations; reuse is from the same frozen noisy table. The rule itself uses observed model features and no hidden scenario labels. It was frozen before a separate policy confirmation seed `69377123`.

The independent policy recomputation gives the primary gated/always-log arithmetic-mean error ratio **0.975945**, with paired whole-system interval **[0.919248, 1.063477]**, over 120 positive-noise cases in 40 systems (480 total method/case rows including zero noise). This does not confirm a gate advantage. Its 90th-percentile error worsens. The rule remains experimental and the adapter defaults to the log grid. The always-adaptive control in this second fresh cohort improves its mean point estimate but has a mean-ratio interval crossing one; no universal risk-protection claim is supported.

## Final independent contour campaign

The final direct engine has SHA-256 `7c123778df29aa63bd6870d14bd41aa75343662910065ef812ce687f7cd4a228`. Its fixed 45 adjoining full-strip tiles cover exactly `10 <= Im(s) <= 100` and `-0.1 <= Re(s) <= 1.1`, with 29 zeros, all real in the t coordinate and simple, and no unresolved leaves. The finite campaign takes about 64.06 seconds in its recorded environment; a fresh direct replay agrees on all 45 tiles in about 61.49 seconds. Five analytic controls, seven stored-artifact corruptions and three serialization controls are retained. The known zero list is used only for the post-campaign diagnostic. The timed-out direct Lehmer attempt remains a limitation.

Exact dyadic midpoint/radius fields repair the loss of zero-exclusion information caused by coarse display strings. This review independently checked serialization round trips: outward reconstruction can slightly widen the ball while containing the original; it is not necessarily bit-identical. Static checks require the reconstructed boundary ball still to exclude zero. Static checks establish recorded coverage and enclosure consistency, not the provenance of saved evaluations; fresh direct replay remains the stronger check and still trusts Arb.

## Final report review

Reviewed the top-level report builder against the final valid receipt paths. Its main adaptive statistic, runtime cost, nonpromoted gate, exact physical ambiguity/remedy, finite contour scope, negative external result, and untested post-confirmation RL extension are consistent with the inspected evidence. Requested two wording refinements: adaptive acquisition uses accumulated observed samples after its pilot; test values/frequencies are not supplied to the designer rather than inaccessible through a security boundary. The impedance upturn is a fact about stored measurements with unknown measurement uncertainty; the monotonicity restriction is exact, but the diagnosis does not quantify its share of the full predictive loss.

No remaining material issue was found in the inspected scope after the generator quarantine, fresh confirmation and listed corrections. Final packaging must still check all source, input and output hashes. This review is internal independent inspection, not external peer review or a proof of every underlying theorem.
