# Observed-pilot acquisition policy: not confirmed for mean-risk protection

The gate did **not** establish its primary goal on a fresh cohort. Its mean
held-out NRMSE ratio against always using the fixed logarithmic grid was
**0.976**, paired system-bootstrap 95% interval **0.919–1.063**. The 90th
percentile worsened from 12.39% to 14.58%. Keep the gate experimental; the
provided JSON adapter defaults to fixed logarithmic acquisition.

The gate is interpretable: after the common 33-point noisy pilot, fit a rational
model and calculate the smallest/largest fitted residue-magnitude ratio. Use
disagreement-driven adaptive acquisition when this ratio is at most 1e-4;
otherwise complete the fixed 64-point logarithmic grid. No hidden regime label,
clean response, true pole/residue, or true physical parameter enters this rule.
All strategies use exactly 64 unique frequency measurements, including the
entire common pilot, and the same final fitter. The gate's pilot values are the
same frozen measurements included in its final query budget.

## Selection and fresh confirmation

The small fixed rule set used only the **repaired adaptive development cohort**
(ten systems, 40 noise cases). Candidates were simple thresholds on estimated
pole gap, residue ratio, damping, training residual, and model disagreement,
plus a few two-feature conjunctions. Eligibility required at least four selected
system clusters, at least 5% reduction of development mean error capped at one,
and development bootstrap upper ratio at most one. Development uncertainty was
a selection heuristic; it was not treated as independent evidence. The chosen
rule selected eight of 30 noisy development cases across five systems, with
capped-mean ratio 0.9205. All candidate scores are retained.

The rule, selection files, source code, protocol, and adaptive dependencies were
locally hash-frozen before a **new seed** generated 40 additional systems, each
at 0%, 0.1%, 1%, and 5% noise: 160 cases and 480 fits. No confirmation features
or outcomes selected or tuned the rule. The confirmation primary outcome is
uncapped mean NRMSE, with a failure penalty of 1e6; this differs deliberately
from development's capped selection objective. No final fits failed.

| Policy, noisy cases | Mean NRMSE | Geometric NRMSE | Median | 90th percentile | Median computation |
|---|---:|---:|---:|---:|---:|
| Always fixed log grid | 5.086% | 1.688% | 2.121% | 12.39% | 6.82 ms |
| Always disagreement acquisition | 4.155% | 1.121% | 1.164% | 10.53% | 77.91 ms |
| Observed-pilot gate | 4.964% | 1.296% | 1.371% | 14.58% | 15.46 ms |

The gate selected adaptive acquisition on 42/120 positive-noise cases. Always
adaptive had mean-error ratio 0.817 [0.580,1.160] against fixed log, also
inconclusive for this cohort's mean objective. Its geometric error was lower,
consistent with the sibling acquisition experiment. These observations do not
prove that pilot residue ratio safely identifies when adaptation helps or
protects against tail failures.

## Files and use

`gate.py` exposes `choose_strategy(omega, noisy_values, order, band)` and the
frozen selected rule. `policy_api.py` adds a JSON-friendly `decide(payload)`
adapter described by `policy_schema.json`; its default mode is `fixed_log`.
Explicit `experimental_gate` mode applies the tested gate. It is restricted to
the tested public band and known orders 6, 8, or 10. The frequency inputs must
be the common 33-point pilot from the sibling adaptive module.

`all_rules_development_scores.json` and `development_features.json` retain
selection evidence. `protocol.json`, `selected_rule.json` and `freeze.json`
record the frozen choice. Raw JSON/CSV and NPZ files retain all new confirmation
observations, decisions, and outcomes. `audit.py` recomputes pilot features,
checks their decision, verifies that every pilot frequency is among the 64
queries, reconstructs the exact shared noise values, and verifies frozen hashes.

Run from this directory, with the sibling `adaptive/` directory preserved:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 experiment.py confirmation
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 audit.py
```

Do not rerun `experiment.py select` to reinterpret the completed confirmation.
The source development seed is 61739311; the new policy test seed is 69377123.
This is simulation input isolation and a local freeze, not external
preregistration or measured hardware validation. The physical generator retains
the explicit general solve and independent NumPy/residual checks introduced
after invalidating the earlier automatic-solver run.
