# External measured impedance challenge

This expansion evaluates shared Nonagon/Loewner machinery and a deliberately constrained positive-real model on an independent published measurement collection. The source contains 54 physical battery cells, rather than multiple train/test splits of one sweep. This is still one batch, one state of charge and one instrument—not a new laboratory experiment or independent calibration validation.

## Source and protection

Moertelmaier, Kasper and Clark, **EIS data of 54 21700 cells**, Zenodo version 1.0.0, DOI [10.5281/zenodo.15422339](https://doi.org/10.5281/zenodo.15422339), published May 26, 2025. The source describes separate NMC811 cylindrical 21700 cells at 30% state of charge, measured galvanostatically with a Keysight SL1007A Scienlab Battery Test System. Each CSV contains 44 complex impedance points from approximately 0.05 Hz to 10 kHz. All 54 source CSV files and their metadata are retained. Published MD5 checksums were verified; SHA-256 values are in `data/download_manifest.json`.

The record and all 54 per-cell metadata entries declare CC BY 4.0. Attribution and the separate schema license are retained in `data/ATTRIBUTION.md` and the original source metadata. No outcome-based dataset selection occurred.

IDs 01–06 are development devices. IDs 07–54 form the protected 48-cell confirmation collection. ID01 was displayed to inspect the file format before freezing. No confirmation response values were inspected or scored before the protocol was fixed. Independent review added mandatory all-device failure sensitivity before any benchmark execution. The final local protocol SHA-256 is:

`db39d7009bc0cf56d4beca4d901d2728ee00436044b81bad708a260c18daf529`

This is a recorded local freeze, not an external preregistration. Development execution was followed by confirmation with no model grid, preprocessing, metric or selection changes.

## Fixed comparison

Four predefined tasks combine 12- or 20-point training budgets with interleaved holdout or a contiguous missing band. Validation indices and test indices are disjoint from training and from one another. Every method sees identical training values in each task. Validation chooses the candidate; no model is refitted using validation or test values. Different budgets have different heldout sets, so this is not a clean monotone learning-curve experiment.

Methods:

- **Loewner:** two-sided SVD projection, fixed candidate orders 2–10, real-system conjugate symmetry. Negative-frequency conjugates are constraints, not extra independent measurements. Unstable models are recorded rather than silently discarded.
- **Vector fitting:** scikit-rf 2.2.0, stable pole relocation, fixed real/complex pole grids, 100 maximum iterations, and constant/proportional terms. Stability alone does not certify passivity.
- **Cubic interpolation:** real and imaginary components interpolated together against log frequency, with fixed not-a-knot boundary conditions.
- **Positive RC relaxation:** nonnegative least squares on fixed Debye grids with series ideal R, L and C terms; grid size and ridge strength selected on validation. This is an existing physics prior, not a unique v5 algorithm.

The constrained model is

\[
Z(s)=R_0+Ls+\frac{C^{-1}}s+\sum_j\frac{a_j}{1+s\tau_j},
\qquad R_0,L,C^{-1},a_j\ge0,\quad\tau_j>0.
\]

Every selected constrained model exports its coefficients as decimal strings with units and `s=+i*2*pi*f_hz` convention in `selected_passive_models.json`. These coefficients permit an independent exact construction check. Such a check certifies the **model's** positive-real property on the open right half-plane; it does not certify that the physical battery obeys the model outside observed frequencies.

## Fresh confirmation results

There are **48 confirmation cells × 4 tasks × 4 methods = 768 selected evaluations**, with 4,608 candidate fits. There were zero candidate exceptions and zero selected-model failures. Development adds 576 candidate fits and 96 selected evaluations, reported separately.

Median heldout complex relative RMS error across the 48 confirmation cells, expressed as percent:

| Task | Training points | Loewner | Vector fitting | Cubic log-frequency | Positive RC relaxation |
|---|---:|---:|---:|---:|---:|
| Interleaved | 12 | 0.163 | 0.105 | **0.088** | 0.714 |
| Interleaved | 20 | 0.387 | 0.096 | **0.056** | 0.743 |
| Missing band | 12 | 0.365 | **0.287** | 0.854 | 0.560 |
| Missing band | 20 | 0.421 | **0.225** | 0.424 | 0.422 |

Across all four fixed tasks, averaging log error ratios within each cell before bootstrapping cells, the positive RC model has **3.48×** the vector-fitting error (unadjusted 95% device-bootstrap interval **3.14–3.85**). Its corresponding ratio to Loewner is **1.91×** (**1.75–2.09**). These are descriptive, conditional batch comparisons. Twelve task-specific contrasts are all retained; intervals are not familywise adjusted. They do not account for shared instrument calibration errors.

The mandatory capped-failure sensitivity includes every confirmation cell; with zero failures it agrees with the common-success analysis. The score floor/cap and bootstrap seed were fixed before execution. Vector fitting emitted 37 maximum-iteration warnings and one barely-converged warning among its 1,728 confirmation candidates. All warning text is retained. The selected Loewner models have reported right-half-plane poles in 177 of 192 tasks; the selected vector-fitting models have none. No instability filter was applied after selection. Good measured-band prediction therefore does not by itself qualify a model for stable dynamical use. No passivity warning is interpreted as a physical battery diagnosis.

![Measured confirmation comparison](external_comparison.png)

## What the loss teaches

A post-confirmation structural diagnosis explains a concrete limitation of the selected passive family. On the frequency axis its real part satisfies

\[
\frac{d}{d\omega}\Re Z(i\omega)
=-\sum_j\frac{2a_j\tau_j^2\omega}{(1+\tau_j^2\omega^2)^2}\le0.
\]

All 48 confirmation curves instead have a high-frequency upturn in the stored real impedance. For example, ID07 increases by approximately 0.000663 ohm from its minimum at 3,763 Hz to 10,000 Hz. The data therefore cannot be matched exactly by this restricted family. Unknown measurement uncertainty prevents converting that statement into a claim about exact physical incompatibility or nonpassivity.

More general passive RL relaxation branches can allow a rising real part. They were **not** added after seeing the confirmation scores. A revised family needs new protected data to establish predictive improvement. The present result separates two achievements cleanly: constructing a model with a provable structural property, and choosing a model family accurate enough for real measurements.

## Replay and evidence

Use Python 3.12, NumPy 2.5.4, SciPy 1.18.1 and scikit-rf 2.2.0; Matplotlib is needed only for the plot. With retained source data, no network is needed:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 external_suite.py --phase development
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 external_suite.py --phase confirmation
python3 summarize.py
python3 diagnose_model.py
python3 plot_external.py
```

`acquire_data.py` can download the entire fixed source record again and verify its source checksums. `*_results.json` binds source, protocol and input-manifest hashes. `*_candidates.json` preserves every candidate and warning. `*_predictions.json` preserves every measured and predicted heldout value. `summary.json` and `scores.csv` provide paired comparisons and compact scores. `model_mismatch_diagnostic.json` is explicitly labeled post-confirmation. No v5 superiority, new global mathematical theorem, calibration guarantee or physical parameter identification is claimed from this experiment.
