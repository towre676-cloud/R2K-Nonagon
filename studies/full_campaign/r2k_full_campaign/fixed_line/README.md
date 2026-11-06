# Complete fixed-line sine-model covariance

The six-length covariance for S=18000 and L=.5,1,2,4,8,16 is now complete.
Start with the parent R2K_Comprehensive_Advance_Report.html. Reproduce with
build_counts_native.py (optional Linux accelerator), full_matrix.py, and
validate_full_matrix.py. All 21 entries include exact anchor-overlap weights,
spatial/separation/tail errors and rigorous sample-centering allowances.
The validator proves positivity and supplies an exact rational Loewner upper
matrix for a conservative simultaneous model bound without Gaussianity.

R2K_COUNTS_PYTHON=1 selects Python ball evaluation. An incompatible native
library also falls back to Python. Linux native execution was tested; a Windows
run is not claimed. Each matrix_cases file retains its full parameter contract.
full_matrix.py reuses completed cases and recovery prefixes; preserve then
remove that directory to request independent complete regeneration.

The original scalar proof below is retained as a separate production control.

Start with `R2K_Fixed_Line_Covariance_Report.html`.

Scope: stationary unit-intensity sine projection DPP, fixed observation span
S=18000, interval length L=4, anchors in [0,17996]. This is a conditional model
certificate, not a zeta sampling law. It does not replace the arithmetic BK mean.

The calculation certifies the fourth-order influence covariance integral and
bounds sample-centering through a projection-operator cumulant identity. It
retains exact finite-line boundary weights, overlapping count intervals, and
an analytically integrated long-distance expansion with a geometric remainder.

Reproduce in this order:

```
python exact_controls.py
python sine_covariance.py
python validate_fixed_line.py
python build_report.py
```

Dependencies: python-flint, sympy (for exact controls). No random simulation,
floating BLAS result or unverified block-independence assumption is used in the
production enclosure. The exact controls use a fixed random seed to choose
rational fixtures; every expectation and determinant is exact.

The default production run regenerates all 128 near panels. After a new run
has written its full-parameter checkpoint, `--resume` can reuse it. All published
panel balls are in `sine_fixed_line_results.json`, which is the acceptance record.
Run the independent validator before interpreting the result.

The finite-height arithmetic third/fourth-order input remains open. The scalar
report predates the now-complete six-length matrix. Do not divide BK residuals
by these sine-model standard deviations and call them calibrated significance.
