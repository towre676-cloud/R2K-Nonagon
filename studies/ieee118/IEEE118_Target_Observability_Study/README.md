# IEEE-118 target-relative measurement study — 4 November 2026

This package extends the user's existing Nonagon known-topology intervention benchmark. It is a structural network-parameter inverse problem, not conventional AC bus-voltage state estimation and not a field-data validation.

Open `report.html` for the completed results. `results/*.json` retain exact panel indices, numerical checks, local uncertainty diagnostics, branch-decision calculations, channel-deletion checks, parameter-shift diagnostics, and nonlinear fits. `protocol.json` records the choices made before new results. `run.log` retains execution progress and `results/receipt.json` records completion.

## Reproduce

Python 3.11+ with NumPy, SciPy and matplotlib. No network is required.

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python run.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python build_report.py
```

The vendor directory retains the source benchmark modules and saved original results. The new run independently checks the bundled MATPOWER file SHA-256, parses tap-adjusted branch coefficients, compares them against the frozen fixture, and checks its vectorized forward model and derivatives against both the original implementation and centered finite differences. The original modular-rank/branch replay was executed separately; its receipt is retained in `results/original_replay.json`.

## Information and acquisition accounting

The reduced operator has dimension 117 and 179 identifiable aggregate internal-edge/reference-ground weights. Parallel physical branches are merged; separate parallel-branch parameters are not identified. Every reported measurement is one scalar port compliance after a baseline, one-node, or two-node grounding intervention. All 6,903 candidate configurations come from the same original protocol.

Design uses supplied topology, uniform nominal weights, and one baseline response to set scale and intervention strength. No actual hidden branch weights enter design scoring. The baseline is included in all budgets. The historical 180th reading was chosen adaptively from fitted candidate branches; that historical panel is labelled separately and is not treated as a fixed topology-only design.

Targets are twelve explicitly declared effective resistances, distinct from the measured port compliances. Their indices are fixed by a deterministic rule before design. These targets are illustrative scientific quantities rather than operator-validated application targets.

Every new design except random starts with a common 179-row nominal QR panel. `qr_norm` adds rows by raw log-Jacobian squared norm; `d_opt` greedily increases full-parameter log determinant; `a_opt` decreases full log-parameter covariance trace; `target_a` decreases the declared normalized target covariance trace; `target_correlated` uses the same target criterion with a common-reference/group-bias noise model. These are forward greedy heuristics from a fixed seed, not globally optimal subsets. The covariance is recomputed by square-root least-squares every twenty additions. The random comparator is one frozen nested draw; it does not estimate the distribution over random panels.

Noise has standard deviation sigma times the baseline response. The correlated model is 50% independent variance, 25% common-reference variance and 25% group variance, with group assigned by the first intervention node (a declared synthetic convention). Every scalar reading has the same marginal variance as in the iid model. Neither this grouping nor its fractions have been calibrated against instruments. Candidate-aware design accounts for correlation; GLS and a naive iid decoder are evaluated separately. Common nuisance biases have a specified distribution; unknown log-weights have no regularizing prior in the reported local uncertainty calculation.

The nominal model and chosen intervention strength use the baseline pilot without propagating its noise into the design. Results therefore condition on that pilot. A future end-to-end campaign should reacquire a noisy pilot and redesign within every replicate.

## Claims and limits

* Exact finite-field rank witnesses and original saved-branch replay are inherited mathematics/computation, explicitly replayed.
* New uncertainty and branch-error formulas are exact for their stated fixed linear/Gaussian or two-known-candidate models, subject to numerical evaluation.
* New comparisons, parameter-shift tests and nonlinear fits are exploratory synthetic experiments, not an independently confirmed superiority claim.
* A full-rank Jacobian establishes local identifiability. Separating two discovered branches does not establish global injectivity.
* A numerical full-parameter rank failure at relative singular-value threshold 1e-12 triggers the target-kernel test. When the target annihilates the numerical nullspace (gap below 1e-8), a target-only SVD decoder is evaluated without assigning finite uncertainty to unrecovered parameters. Neither test is automatically an exact-rank theorem. This qualification was added during the final audit; it corrects an overly conservative initial handling of random panels and changes no design or target choices.
* The 1,000-trial linear experiments test local decoders. Large predicted errors invalidate using the linearization as a nonlinear guarantee.
* Nonlinear fits start at the baseline-scaled uniform nominal and use 100 function evaluations at most, with each log-weight restricted to nominal +/-8. Solver termination is not recovery; errors against hidden truth are reported independently.
* Five nominal-start trials at zero noise are identical deterministic replays, not five independent starts. Distinct nonzero-noise replicates supply the stochastic comparisons.
* Parameter perturbations preserve support; no new blind-topology or outage-support recovery claim is made.
* No claims are made about PMU hardware count, communication requirements, intervention feasibility or full AC state estimation.

See `THEORY.md` and `SOURCES.md` for mathematical guarantees and connections to primary literature.
