# Closed-loop physical-frequency acquisition

The repaired, frozen confirmation found a useful acquisition gain at the same
measurement budget. The selected disagreement-driven strategy reduced geometric
mean clean held-out response error **41.7%** against a strong fixed logarithmic
grid: error ratio **0.583**, paired system-bootstrap 95% interval **0.463–0.725**.
It won 117/150 positive-noise comparisons. Arithmetic mean error and the 90th
percentile also improved in this cohort. This is synthetic, known-order
validation, with additional computation per acquisition; it is not measured
hardware or a universal identifiability result.

## Valid confirmation

Fifty independently generated physical mass–spring–damper systems, ten per
predeclared regime, were each tested at 0%, 0.1%, 1%, and 5% multiplicative complex
noise. There are 200 cases and 1,000 final fits. Every strategy gets 64 unique
positive-frequency observations from the same public 0.03–8 rad/unit-time band.
Every queried frequency has a fixed noisy value shared across strategies.

| Strategy | Geometric NRMSE | Median NRMSE | Mean NRMSE | 90th percentile | Median computation |
|---|---:|---:|---:|---:|---:|
| Strong 64-point logarithmic grid | 1.691% | 2.134% | 5.686% | 13.00% | 6.95 ms |
| Linear coverage plus common pilot | 2.519% | 3.076% | 7.403% | 16.36% | 7.61 ms |
| Random frequencies plus common pilot | 3.149% | 3.578% | 11.32% | 24.71% | 7.13 ms |
| Predicted-magnitude acquisition | 2.667% | 3.860% | 9.322% | 28.12% | 51.77 ms |
| **Selected model-disagreement acquisition** | **0.986%** | **1.024%** | **3.596%** | **11.60%** | 79.08 ms |

The comparison above uses the 150 positive-noise cases. NRMSE is the Euclidean
error divided by the clean response norm on a protected 1,024-point frequency
grid. Smaller is better. Bootstrap resamples 50 whole system clusters, retaining
their three noise levels together. The fixed log grid remains preferable when
the additional design computation is unjustified; acquisition time and hardware
cost were not measured.

The predeclared regime error ratios (adaptive/log, geometric mean) were 0.971 for
ordinary chains, 0.760 for colored noise, 0.458 for weak links, 0.210 for
near-colliding modes, and 0.950 for model mismatch. The ordinary and mismatch
intervals include one. The strongest evidence is for weakly visible and nearly
coincident modes; these subgroup comparisons are explanatory and are not a
familywise-adjusted collection of independent discoveries.

All 1,000 final fits returned finite predictions, with no recorded warnings or
exceptions. Finite output is not an accuracy guarantee: the maximum adaptive
error was 49.75% in the noisy cohort. For 40 correctly specified zero-noise
systems, log-grid median error was 3.41e-14 (maximum 5.92e-12), and adaptive
median was 2.08e-14 (maximum 2.02e-13). Mismatch cases deliberately omit an extra
oscillator from the supplied order.

## What actually adapts

The candidate lattice has 1,025 positive physical frequencies; the held-out grid
is disjoint. All strategies share a 33-point pilot that is a subset of the
strong 64-point log grid. Pilot queries count toward the total budget. The
selected design spends its remaining 31 queries in batches of four. Three
queries target disagreement between two rational models fitted to the observed
noisy data; one explores the largest remaining log-frequency gap. Models refit
after each batch. No true pole, residue, mass, noise fraction, regime label,
clean response, or held-out response is supplied to the acquisition function.

The final fitter is identical for all strategies: ordinary projected Loewner
followed by eight stable-pole vector-fitting iterations, with the same known
chain order and training-loss iterate selection. This isolates acquisition
benefit. Negative frequencies are conjugate symmetry augmentation of the same
positive observations, not extra independent measurements. Querying a frequency
again is forbidden.

Development compared a linearized pole/residue I-optimal information surrogate
with model disagreement. The disagreement design was selected on separate
repaired development systems; the information surrogate did not win. Its
covariance was a design heuristic, not a calibrated uncertainty interval. The
method uses established rational modeling components and is inspired by the
Nonagon acquisition program; these tests do not establish algorithmic novelty.

## Important invalidated precursor

An initial run was **invalidated**, and its receipts are retained under
`invalidated_auto_solver_run/`. It used SciPy 1.17.0 automatic matrix-structure
selection for a complex symmetric tridiagonal physical solve. A zero-noise
sanity check exposed large response errors for three-mass systems. Direct
residual checks showed the automatic solution could have residual around
0.0023, while the explicit general solver and NumPy agreed near machine
precision. Those initial results must not be used.

The valid run uses fresh development and confirmation seeds. Every candidate
and held-out physical response uses `scipy.linalg.solve(..., assume_a='gen')`,
checks its backward residual, and independently compares its entire solution
against `numpy.linalg.solve`. Across 102,450 physical evaluations in the valid
confirmation, maximum normalized backward residual was 1.94e-16; the independent
solutions agreed exactly in the recorded runtime. These diagnostics are retained
per generated system in every raw row. The earlier same-data identification
challenge already used the explicit general solver and is unaffected.

The valid frozen scientific files have not changed after confirmation.
`development_general_solver_only/` retains an intermediate valid development
run before the added independent NumPy check; current development receipts
include both checks. Development exploration is not external preregistration.

## Reproduce and inspect

```bash
python3 replay.py
```

`requirements.txt` records NumPy 2.3.5 and SciPy 1.17.0. Explicit solver selection
and residual checks avoid the problematic automatic path. Replay overwrites
current confirmation outputs and timings. `audit.py` verifies frozen hashes,
candidate/noise-table hashes, selected-observation hashes, unique 64-query
budgets, and the disjoint protected grid. `acquisition.py` is the main program;
`reference_core.py` is an imported fitter dependency, not the experiment entry
point. `protocol.json` and `freeze.json` record the valid design and fresh seeds.
`development_*`, `confirmation_raw.json`, `.csv`, the observation `.npz`, and
`results.json` retain development selection and all confirmation outcomes.

The next deployment question is whether an observed-pilot policy can safely
choose between adaptive and fixed acquisition. That is tested separately in
the sibling `policy/` experiment on another fresh cohort; it does not alter
these frozen results.
