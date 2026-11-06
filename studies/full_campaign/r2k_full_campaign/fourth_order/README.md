# Fourth-order circular-model calibration

This follow-up supplies a fully specified model covariance, not a theorem about
the sampling law of zeta zeros. The primary finite-height arithmetic mean and
fixed-line anchor estimator in the parent reports remain the principal comparison.

Run `python validate_trace.py`, `python cue_covariance.py`, `python validate_covariance.py`,
`python validate_periodic.py`, and `python build_report.py` in this directory.
Python dependencies: python-flint, numpy, matplotlib. The production calculation
uses 128-bit Arb, Fourier cutoff 180000, circumference 18000, and CUE ranks
17999, 18000, 18001. No random simulation enters the covariance calculation.

The circular estimator introduces an artificial seam and a deterministic full-circle
mean. It differs from the original line-window estimator. The observed rank is used
as a model parameter; no equivalence of this conditional CUE law with zeta data
is assumed or proved. These data and model choices are exploratory.

The fourth cumulant identity is derived from the finite-rank determinantal
log determinant. All 75 ordered-partition trace spans are proved symbolically on
the cone j>=k>0; swapping j,k covers the other cone. Exact checks supplement the
derivation. Infinite covariance tails have a displayed absolute bound; Arb
propagates finite-sum rounding. A dyadic congruence and interval Gershgorin certifies positivity of each complete
6 by 6 covariance enclosure. Periodic empirical pair sums are independently
matched by 30 exact event sweeps, with propagated ordinate-coordinate error.

Primary background: Diaconis and Evans, *Linear functionals of eigenvalues of
random matrices*, Transactions AMS 353 (2001), 2615–2633:
https://statistics.berkeley.edu/sites/default/files/tech-reports/577.pdf
The finite Fourier covariance identity and bounds used here are derived in the
accompanying report; we make no novelty claim.
