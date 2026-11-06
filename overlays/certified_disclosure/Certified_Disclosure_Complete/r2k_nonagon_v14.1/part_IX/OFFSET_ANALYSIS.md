# The apparent BK offset comes from the integration endpoint

The three-window offset has a direct numerical cause. `part_VI/bk_multi.py` and `part_V/bk.py` construct their integration grid with `np.linspace(1e-4, 16, 1601)`. Their number-variance routine integrates only the retained grid points. The required lower limit is zero. The missing interval produces a positive error approximately equal to **0.0002 times L**, independently of the height window. At L = 16 its size is 0.0032, which accounts for the reported height-stable discrepancy.

This is reproducible without generating any new zeta zeros. `offset_audit.py` independently reproduces all six archived GUE values and compares them with the analytic sine-kernel variance and adaptive quadrature. It also applies the universal leading endpoint correction to all 18 archived BK comparisons. The output is `offset_audit.json`; source hashes bind it to the inspected v7 files.

## The missing term

For the normalized, distinct-point pair correlation,

$$
\Sigma^2(L)=L-2\int_0^L (L-r)\bigl(1-R_2(r)\bigr)\,dr.
$$

Starting instead at a positive value a increases the computed variance by

$$
E_a(L)=2\int_0^a (L-r)\bigl(1-R_2(r)\bigr)\,dr.
$$

Both the sine-kernel correlation and the implemented BK expression have $R_2(r)=O(r^2)$ near zero. Consequently,

$$
E_a(L)=2La-a^2+O(La^3).
$$

Here a = 0.0001. At L = 16 the leading correction is 0.00319999. The BK formula's zero limit follows from cancellation of the zeta pole, the Euler product's quadratic term, and its paired prime sum. With $z=ir$, $\ell=\log(T/(2\pi))$, $b=B(0)$ and $c=\gamma_0^2+2\gamma_1$,

$$
(\zeta'/\zeta)'(1+z)=z^{-2}-c+O(z),\quad
\zeta(1+z)\zeta(1-z)=-z^{-2}+c+O(z^2),\quad
A(z)=1-bz^2+O(z^3).
$$

After multiplication by $e^{-\ell z}$ and subtraction of $B(z)$, the real part tends to $-\ell^2/2$. The implemented normalization $1+2\Re(\cdots)/\ell^2$ therefore tends to zero. The shared finite prime cutoff preserves this cancellation because the same prime sum b occurs in both terms.

**Halving the grid spacing while keeping a fixed does not test this error.** The omitted interval survives every such mesh refinement. The previously reported small coarse-versus-fine shift measures a different error.

## Independent GUE reproduction

The analytic sine-kernel number variance is

$$
\Sigma^2_{\mathrm{sine}}(L)=
L\left(1-\frac{2}{\pi}\mathrm{Si}(2\pi L)\right)
+\frac{\log(2\pi L)+\gamma-\mathrm{Ci}(2\pi L)+1-\cos(2\pi L)}{\pi^2}.
$$

The audit evaluates this formula and independently integrates from zero with adaptive quadrature. They agree within 3e-12. It reproduces the archived legacy-grid values within 5e-13.

| L | Legacy GUE minus analytic value | Exact omitted lower integral | Remaining mesh error after endpoint correction |
|---:|---:|---:|---:|
| 0.5 | 0.0001314114 | 0.0000999900 | 0.0000314214 |
| 1 | 0.0001833218 | 0.0001999900 | -0.0000166682 |
| 2 | 0.0003833035 | 0.0003999900 | -0.0000166865 |
| 4 | 0.0007832797 | 0.0007999900 | -0.0000167103 |
| 8 | 0.0015832354 | 0.0015999900 | -0.0000167546 |
| 16 | 0.0031831475 | 0.0031999900 | -0.0000168425 |

At L = 16, the archived value is 0.6301163488, whereas the analytic value is approximately 0.6269332013. The independent GUE calibration in `part_VII/cue_calib.json` uses 0.6269332002, obtained by fine quadrature from zero. Its approximately 1e-9 difference from the analytic value is immaterial here. These two pre-existing files therefore already contain a deterministic inconsistency of approximately 0.00318 in their GUE reference values.

The final grid point below L is another detail: except at 16, the shared grid does not land exactly on the requested upper endpoint. This is especially visible at L = 0.5. At integer L, the remaining GUE mesh error is close to the expected trapezoidal contribution $-h^2/6$; it is about 0.000017 rather than 0.0032.

## Corrected BK comparison

The following table subtracts $2La-a^2$ from the archived height-averaged BK values. It does not regenerate the BK curves or claim interval-certified corrected values. Their remaining mesh error should be removed in the production rerun.

| First zero index | L | Original BK minus observed | Observed minus endpoint-corrected BK | Corrected residual / archived block SE |
|---:|---:|---:|---:|---:|
| 100,000 | 8 | 0.0019736954 | -0.0003737054 | -0.776 |
| 1,000,000 | 8 | 0.0012414379 | 0.0003585521 | 1.138 |
| 10,000,000 | 8 | 0.0014031234 | 0.0001968666 | 0.318 |
| 100,000 | 16 | 0.0032544731 | -0.0000544831 | -0.103 |
| 1,000,000 | 16 | 0.0029848769 | 0.0002151131 | 0.268 |
| 10,000,000 | 16 | 0.0030691165 | 0.0001308735 | 0.228 |

Across all 18 rows, the largest absolute endpoint-corrected residual is 1.454 archived block standard errors. This removes the previously conspicuous L = 8 and L = 16 discrepancies. These are descriptive residuals: different L values share the same zeros and counting anchors, so their standardized residuals must not be pooled as independent statistical tests.

The change is to the reference integral; none of the empirical zero statistics has been altered. The GUE and BK reference curves share the endpoint defect. Consequently the old ratio called the fraction of GUE departure explained is not protected merely because both curves used the same grid: its denominator changes under the repair.

## What the height argument does and does not establish

The inspected median heights give logarithms 9.4715086, 11.4763398 and 13.5864723. A pure fixed-coefficient $1/\log^p(T/(2\pi))$ term would decline from the first to the last window by 30.29% for p = 1, 51.40% for p = 2 and 66.12% for p = 3. Thus “at least a third for every p ≥ 1” is slightly too strong for the actual endpoints; the p = 1 decline is about 30%.

The observed stability alone would not rule out every finite-height correction. Using the archived L = 16 block errors, descriptive one-amplitude fits to the raw offsets give chi-squared values 0.099 for a constant, 1.379 for p = 1, 5.838 for p = 2 and 11.746 for p = 3, with two residual degrees of freedom if independent Gaussian errors were justified. In particular, p = 1 is not excluded by these three measurements. Mixtures of powers are still more flexible. The decisive evidence is instead the exact reproduction of the integration defect and its predicted L dependence.

## Distinguishing the other candidates

* **Lower cutoff:** produces $2aL-a^2$ to leading order. Varying a while keeping the mesh fine directly tests its predicted slope.
* **A missing constant in R2:** a change $R_2\mapsto R_2+c$ changes the variance by $cL^2$, rather than a linear term. More generally, $\Sigma^{2\prime\prime}(L)=2(R_2(L)-1)$ separates pair-correlation errors from an approximately affine integration offset.
* **An additive constant in unfolding:** the code subtracts U[0], so the +7/8 term cancels exactly. Changing that constant cannot affect the result.
* **A multiplicative unfolding error:** for a stationary process rescaled by c, the curve becomes $\Sigma^2(L/c)$. Its first-order effect follows the derivative of the variance curve, rather than necessarily the cutoff's linear law. Future empirical output should retain local mean spacing and the mean count at each L.
* **Centering of the empirical variance:** the code uses `var(count)`. Exactly, this equals `mean((count-L)^2) - (mean(count)-L)^2`. Retaining both quantities will expose any finite-window mean-count contribution.
* **Estimator bias:** the file named `cue_calib` actually simulates GUE tridiagonal matrices, with N = 40,000 and the central 20,000 points over 32 replicates. At L = 16 its estimated bias relative to the sine kernel is +0.0028175 with Monte Carlo SE 0.0014256. This is an imprecise stochastic calibration, not a deterministic explanation of the BK offset. Its block errors are generally conservative in this campaign. The unused `sigma2_gue` helper has a wrong closed-form expression; the actual saved reference uses the separate fine-quadrature function, so that dead helper did not generate the archived calibration.

## Production repair

Use a grid that contains both 0 and the exact requested L. Assign the analytic value R2(0) = 0, and evaluate near-zero cancellation through a series or adequate precision. Compare composite Simpson or adaptive integration against an independently refined grid, then vary the lower cutoff separately as a regression check. Bind the repaired BK and GUE references to the same integration convention. Retain mean counts, window lengths and cross-L covariance in the next empirical run.

No giant zero campaign is needed to establish this correction. The present audit uses only archived statistics, two independent sine-kernel evaluations, and the explicit lower-endpoint expansion. A corrected BK quadrature rerun will finish the numerical repair without changing the existing zero samples.

## Supplied corrected implementation

`variance_quadrature.py` supplies reusable `number_variance(L, x, r2)` and `integration_samples` functions. They include zero using the declared regular limit, insert the exact upper endpoint by interpolation when needed, reject extrapolation and malformed input, and support Simpson and trapezoidal integration. `sine_number_variance` supplies the corrected analytic sine-kernel formula, including a small-L series. The implementation makes no changes to v7.

`test_variance_quadrature.py` passes 48 checks, recorded in `variance_quadrature_tests.json`. These exercise independent analytic and adaptive sine values, polynomial identities, the original endpoint defect, positive-lower-grid repair, arbitrary upper endpoints, refinement, Poisson normalization, invalid inputs, and the corrected runner's counting and zero-limit interfaces.

Run the focused checks from this directory:

```bash
python test_variance_quadrature.py
```

`bk_multi_repaired.py` provides an optional production rerun. It retains the original BK Euler product, prime sum, zeta-zero acquisition, unfolding and dense-count estimator. It uses the corrected quadrature, records mean counts and block values, compares the GUE integral against the analytic reference, and distinguishes incomplete from complete output. An existing output is refused unless replacement is explicitly requested. It does not read old output as a substitute for a rerun.

```bash
python bk_multi_repaired.py --check-requirements
python bk_multi_repaired.py --indices 100000,1000000,10000000 \
  --output data/bk_multi_repaired_100000_1000000_10000000.json
```

The production campaign needs `mpmath` and `python-flint`, in addition to NumPy and SciPy. Those first two dependencies are absent in the audit environment. The runner's help, dependency report, compilation and independent components have been checked; a fresh full BK campaign has not been run. The supplied JSON containing corrected archived residuals is identified as an endpoint audit, rather than a newly acquired zero campaign.
