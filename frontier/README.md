# Nonagon Frontier — acquisition, inference and independent verification

Open **report.html** for the findings and figures. `SUMMARY.json` contains the compact results. The previous challenge package is included unchanged under `previous/`.

## What advanced

- **Actual closed-loop acquisition:** on a fresh 200-case simulation, adaptive disagreement sampling reduced the prespecified geometric mean held-out error by **41.7%**, at the same 64 unique real-frequency observations. The error ratio is **0.583 [0.463, 0.725]**. Arithmetic mean and 90th percentile also improved in this cohort. Computation rose from about 7 ms to 79 ms. Model order is known; no hardware acquisition is claimed.
- **Exact inference limits and a remedy:** a grounded path and triangle can have exactly the same response, without hidden modes. The construction extends to an unchanged observed block of n−2 ports. Additional sensing separates the explicit pairs; known masses and full displacement responses at n distinct shifts reconstruct unknown stiffness/topology under cyclicity. There are **257 exact pairs / 6,684 checks**, plus 56 independent checks.
- **An independent interval route:** direct Arb integration, with no archived Taylor/DFT imports, covers the complete critical strip over heights **[10,100]**. All **45 tiles / 29 zeros** certify, and a fresh 45-tile replay agrees. This reproduces known finite-height facts; it is not a height record or global RH proof.
- **A reusable exact passive compiler:** 216 fitted model instances are certified in a supported rational class, with passive circuit exports. The expanded RC/RL verifier passes **150 fixtures / 2,400 exact checks** and **29 adversarial rejections**. Model passivity does not certify an unmeasured physical battery or its unique internal circuit.

## Results that were not promoted

The observed-data gate failed its primary risk test on 160 fresh cases: mean-error ratio **0.976 [0.919, 1.063]**, with worse upper-tail error. The research adapter defaults to fixed logarithmic acquisition.

The measured-cell experiment reserves IDs01–06 for development and IDs07–54 for confirmation. Across 48 protected cells, vector fitting wins missing-band prediction and cubic interpolation wins ordinary interpolation. The restricted positive RC model has about 3.48 times vector fitting's geometric error. All 4,608 confirmation candidate fits and 768 selected results are retained.

A post-confirmation diagnostic shows that the RC class forces monotone real impedance while the stored measurements show high-frequency upturns. The new verified RL extension can represent upturns, but it was not fitted to change the confirmation scores. Its accuracy needs a new protected test.

## Contents

| Directory | Capability and evidence |
| --- | --- |
| `adaptive/` | Frozen real-frequency acquisition comparison, raw observations, all fits, audits, invalidated precursor |
| `policy/` | Development-selected observed-data gate, fresh confirmation, adapter defaulting to fixed logarithmic acquisition and schema |
| `identifiability/` | Exact physical network construction, full derivation, added-sensor remedy, independent verifier |
| `independent_certification/` | Direct contour engine, continuous coverage ledger, exact interval serialization, replay and controls |
| `external/` | 54 measured spectra, source/license/checksum metadata, frozen design and every prediction/candidate |
| `passivity/` | Exact compiler/verifier, RC/RL theorem, fixtures, 216 certificates and circuit realizations |
| `review/` | Independent internal review and literal numerical-library failure reproducer |
| `previous/` | Previous challenge package, byte-identical |

## Evidence and numerical failures

The initial adaptive campaign is **invalidated** because SciPy 1.17's automatic solve path returned inaccurate results for a complex 3×3 physical system. `review/solver_reproducer.py` isolates the literal matrix; versioned receipts show the problem under 1.17 and its absence under 1.18.1. The repaired campaign uses explicit general solves, residual checks, independent NumPy comparisons, zero-noise sanity checks, a new local freeze and a new untouched confirmation seed. Its results are the top-level `adaptive/results.json`; invalid precursors are clearly quarantined and excluded.

The direct contour module discovered that Arb display strings can coarsen wide balls. It now stores exact dyadic midpoints and radii, preserves the old display-only run separately, and regenerates/replays the complete campaign.

All freezes are local, not external preregistration. Fitter APIs are not supplied hidden truth; this is not external-custodian blinding. The measured source represents one batch, state of charge and instrument. Subgroup intervals are unadjusted secondary comparisons. No broad novelty, physical extrapolation or universal superiority is claimed.

## Quick replay

Python's standard library is sufficient for file integrity and exact replay:

```bash
python replay.py --mode verify
python replay.py --mode quick
```

Replays operate in a fresh copy. The original saved receipts remain unchanged. A manifest checks bytes, not mathematical correctness by itself.

## Numerical replay

Use Python 3.12. Separate profiles preserve the recorded environments:

```bash
python -m venv .venv-numerics
.venv-numerics/bin/python -m pip install -r requirements-numerics.txt
python -m venv .venv-adaptive
.venv-adaptive/bin/python -m pip install -r adaptive/requirements.txt
python replay.py --mode all --numerics-python /absolute/path/.venv-numerics/bin/python --adaptive-python /absolute/path/.venv-adaptive/bin/python
```

On Windows, use each environment's `Scripts/python.exe`. Version-specific warnings, numerical last digits and wall times can vary. Recorded measured data are bundled; replay requires no data download. Installing dependencies requires internet access.

Full replay repeats frozen confirmations and direct interval verification. Module READMEs document additional generation, source download, theorem and CLI details. The expensive direct Lehmer attempt hit its 60-second budget and remains a measured limitation; it is not retried by default.

Primary sources and established mathematical context are listed in `SOURCES.json`. No messages, issues or results were posted externally.
