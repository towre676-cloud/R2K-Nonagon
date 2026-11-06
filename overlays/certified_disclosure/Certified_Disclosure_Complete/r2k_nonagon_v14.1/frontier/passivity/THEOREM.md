# A constructive positive-real certificate for a restricted rational class

This is a direct sufficient condition and an exact model certificate. It is a classical circuit fact, not a claim of a new positive-real theorem or a unique physical realization.

Let

\[
Z(s)=R_0+Ls+\frac{C_{\mathrm{inv}}}{s}
  +\sum_{j=1}^{m}\frac{a_j}{1+\tau_js}
  +\sum_{k=1}^{r}\frac{b_k\theta_ks}{1+\theta_ks},
\]

where `R0`, `L`, `Cinv`, every `a_j` and every `b_k` are nonnegative real numbers, and every `tau_j` and `theta_k` is strictly positive. A zero `Cinv` omits its term. The second relaxation sum is optional; omitting it preserves the original input schema. Then `Z` is a rational positive-real impedance: its coefficients are real, it is analytic on the open right half-plane, and `Re Z(s) >= 0` there.

For `s=x+iy`, `x>0`, every possible pole lies outside that domain: the Debye poles are `-1/tau_j<0`, and a nonzero inverse-capacitance term has only a boundary pole at `s=0`. Directly,

\[
\operatorname{Re}Z(x+iy)=R_0+Lx+
 C_{\mathrm{inv}}\frac{x}{x^2+y^2}
 +\sum_j a_j\frac{1+\tau_jx}{(1+\tau_jx)^2+(\tau_jy)^2}
 +\sum_k b_k\frac{\theta_kx(1+\theta_kx)+(\theta_ky)^2}{(1+\theta_kx)^2+(\theta_ky)^2}\ge0.
\]

Every denominator in this expression is strictly positive. Each summand is nonnegative; if any coefficient is nonzero the sum is strictly positive for `x>0`. Conjugation symmetry follows from real coefficients. The identically zero model is admitted as the nonnegative-real trivial case.

An explicit passive circuit realizes the same function. Put the following blocks in series: resistor `R0` when nonzero; inductor `L` when nonzero; capacitor `C=1/Cinv` when nonzero; and, for each `a_j>0`, a parallel RC block with resistor `R_j=a_j` and capacitor `C_j=tau_j/a_j`. That block's impedance is exactly

\[
\frac{R_j}{1+sR_jC_j}=\frac{a_j}{1+\tau_js}.
\]

Zero coefficients contribute no block. A zero model is a short circuit. All non-omitted component values are strictly positive. This construction certifies the rational model without checking a finite frequency grid, locating approximate roots, or using a sampled imaginary-axis inequality as a proxy for analyticity.

For each optional `b_k>0`, add a parallel RL block with resistor `R_k=b_k` and inductor `L_k=b_k*theta_k`. Its impedance is exactly

\[
\frac{R_ksL_k}{R_k+sL_k}=\frac{b_k\theta_ks}{1+\theta_ks}.
\]

That is another passive block in the same series construction. Its denominator has its only pole at `-1/theta_k<0`, so it does not alter analyticity in the right half-plane.

## A capacity limit and a compatible extension

Without the optional RL sum, `R0`, `Ls` and `Cinv/s` contribute no frequency-dependent real part on `s=i*omega`, `omega>0`. Therefore

\[
\frac{d}{d\omega}\operatorname{Re} Z(i\omega)
=-\sum_j\frac{2a_j\tau_j^2\omega}{(1+\tau_j^2\omega^2)^2}\le0.
\]

The RC-only relaxation class cannot represent a real-part upturn. This is a structural restriction, not an optimizer issue. The optional passive RL contribution has the opposite derivative sign:

\[
\frac{d}{d\omega}\operatorname{Re}\frac{b_k\theta_ki\omega}{1+\theta_ki\omega}
=\frac{2b_k\theta_k^2\omega}{(1+\theta_k^2\omega^2)^2}\ge0.
\]

Mixtures can therefore admit a real-part decrease followed by an upturn while preserving positive-realness. `monotonicity_results.json` gives exact rational derivative checks and an explicit mixed example. This is a mathematically checked extension of the model class. It was introduced after the measured confirmation results were inspected, and no new RL-extended fitting scores are claimed. Its predictive value needs new protected data or a prospectively frozen new experiment.

## Exact arithmetic and scope

The implementation interprets recorded finite decimal coefficient strings as exact rational numbers, including scientific notation. It checks their signs and time constants using `fractions.Fraction`, not floating point. The certificate binds the input bytes and a canonical exact-coefficient representation, and records the compiler/verifier source identities. The separate verifier recomputes those conditions and hashes. The mathematical guarantee is conditional on the stated theorem and verifier correctness; source hashes identify bytes and are not proofs by themselves.

The supported class is sufficient, not necessary, for positive-realness. A model outside this class may still be passive. Rejecting an unsupported or negative-coefficient representation is not a general proof of nonpassivity.

The guarantee concerns the **serialized fitted model** on `Re(s)>0`. It does not prove that a physical battery has the same impedance outside the observations, that its fitted circuit parameters have unique physical meaning, or that measurement calibration is correct. Circuit export is algebraically checked; no SPICE simulator execution is claimed.

## Units and frequency convention

`s` has units `1/second`; frequency-domain evaluation uses `s=+i*2*pi*f_hz`. The impedance and `R0`/`a`/`b` coefficients are in ohms, `L` in henries, `Cinv` in inverse farads, and both time constants in seconds. The JSON names both RC and RL time constants `tau`. A nonempty `rl_terms` list requires `units.b="ohm"`; old models without RL terms retain their original units mapping. Mixing Hz directly into a formula derived for `s=i*omega` would change the numerical model, so the verifier checks this convention explicitly.
