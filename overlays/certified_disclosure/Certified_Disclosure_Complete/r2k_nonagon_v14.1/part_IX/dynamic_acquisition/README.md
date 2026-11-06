# From finite step records to dynamic operator enclosures

**A finite, ordinary step response supplies the real-Laplace observations required by Part VIII.** Four fresh synthetic networks supply 32 simulated step records and 96 scalar observations. The rational interval inverse encloses all twelve matrices $P(0),P(1),P(2)$, whose linear combinations enclose all entries of $K,C,M$. No apparatus is connected.

Run with Python, NumPy and SciPy:

    python part_IX/dynamic_acquisition/replay.py
    python part_IX/dynamic_acquisition/replay.py --output-dir /absolute/output/directory
    python part_IX/dynamic_acquisition/verify_records.py /absolute/output/directory

The output consists of dynamic_receipt.json, protocol.json and acquisition_records.json. The last file preserves the exact rational test matrices, all input intervals and the reconstructed operator boxes.

## Acquisition contract

The linear, zero-initial-state model is

\[
M\ddot x+C\dot x+(K+D_\theta)x=e_0\,1_{t\geq0}.
\]

A calibrated unit step force enters port zero and displacement is measured at the same port. The unchanged static load protocol uses baseline, $t,2t$ at each hidden site, and $t,t$ at each hidden pair, with $t=1$. At most two loads are active simultaneously. For $s>0$,

\[
h_\theta(s)=s\int_0^\infty e^{-s\tau}x_0(\tau)\,d\tau
=e_0^T(K+D_\theta+sC+s^2M)^{-1}e_0.
\]

At $s=0$ a finite endpoint estimates static compliance with enclosed settling error. The same trace supplies that endpoint and both Laplace integrals. **32 step records produce 96 scalar values**; full time records contain more information than individual scalar experiments, so this does not establish a lower dynamic experiment count.

The frozen acquisition has $T=80$, $\delta=1/2048$ and 163,841 samples per trace. It assumes simultaneous absolute sample error at most $10^{-12}$ and an additional transformed numerical rounding envelope of $10^{-12}$. The calibrated step amplitude and grounding loads are exact. These are declared synthetic conditions, not manufacturer specifications. Initial-state error, drift, force uncertainty and load error need separate envelopes before applying this contract to an instrument.

## Fixed prior and analytic error bounds

The acquisition proof applies to diagonal positive $M$ and the passive grounded diagonal damper $C=2\beta M$. A **fixed external class prior**, identical for every fixture, supplies

\[
\beta=\tfrac12,\qquad
1\leq M_{ii}\leq\tfrac74,\qquad
\sum_jK_{ij}\geq2,\qquad
\|K\|_\infty\leq16.
\]

These class bounds are frozen before scalar records are generated. The envelope function receives only the prior and commanded loads; it receives no hidden matrix or true parameter. The simulator separately verifies fixture membership. A physical application must justify the prior independently. **There is no claimed tail theorem here for arbitrary grounded nonproportional damping.**

Let $m_-=1,m_+=7/4,g_-=2,k_+=16$. The spring energy identity bounds the eigenvalues of $A=M^{-1/2}(K+D_\theta)M^{-1/2}$ between

\[
\ell=g_-/m_+=8/7>\beta^2,\qquad
\Lambda=(k_++\max_i\theta_i)/m_-.
\]

Each modal step contribution is

\[
\frac{w}{\lambda}\left[1-e^{-\beta\tau}
\left(\cos\omega\tau+\frac{\beta}{\omega}\sin\omega\tau\right)\right],
\quad \omega^2=\lambda-\beta^2,\quad w\geq0,\quad\sum w=M_{00}^{-1}.
\]

For any positive rational $\underline\omega\leq\sqrt{\ell-\beta^2}$, this gives the endpoint and transfer-tail bounds

\[
E_{\rm DC}\leq
\frac{e^{-\beta T}}{m_-\ell}
\left(1+\frac{\beta}{\underline\omega}\right),
\]

\[
E_{\rm tail}\leq
\frac{e^{-sT}}{m_-\ell}
\left[1+\left(1+\frac{\beta}{\underline\omega}\right)e^{-\beta T}\right].
\]

The fourth derivative of $e^{-s\tau}x_0(\tau)$ is bounded uniformly by

\[
B_4=\frac{1}{m_-\ell}\left[
s^4+\left((s+\beta)^2+\Lambda-\beta^2\right)^2
\left(1+\frac{\beta}{\underline\omega}\right)\right].
\]

Composite Simpson integration contributes at most $sT\delta^4B_4/180$. Positive Simpson weights sum to $T$, so arbitrary bounded sample errors contribute at most $sT\varepsilon$; the DC sample contributes $\varepsilon$. No independence or averaging gain is assumed.

Envelope arithmetic is rational: an integer square root supplies the lower frequency bound, and $e^{-q}\leq2^{-\lfloor q\rfloor}$ supplies conservative rational exponential bounds. Floating modal synthesis and Simpson evaluation remain conditional on the declared transformed rounding envelope. This is not a machine-checked error proof for floating transcendental routines.

## Reconstruction and controls

Only scalar intervals, known labeled sites and exact load strengths enter the inverse observer. At $s=0,1,2$, the existing static interval inverse encloses $P(s)$. Rational interval combinations then give

\[
K=P(0),\qquad
M=\frac{P(2)-2P(1)+P(0)}2,\qquad
C=\frac{4P(1)-P(2)-3P(0)}2.
\]

The resulting boxes need not certify the exact zero pattern of $M$ or $C$. They are entrywise enclosures under the stated model; dynamic parameter-count optimality is not asserted.

An independent Fraction Gaussian forward solver checks every transform interval and every reconstructed matrix entry against exact fixture truth. Modal step samples are checked at three times per configuration against a separately constructed first-order matrix exponential. Three short-record controls require the tail term to contain the true transform. Six unsupported-model controls are rejected, and wide observation error returns unresolved. Fixture seeds are disjoint from the archived DC, SPICE, instrument and fault campaigns.

The saved campaign passes 455 acquisition and enclosure predicates. The separate record verifier imports neither the acquisition code nor the inverse code and passes 956 checks of independent Cramer responses, prior-based analytic budgets, complete plans and matrix containment. Normal and optimized Python produce byte-identical protocol, acquisition and campaign receipt files. The maximum observed transformed scalar error is $1.00018\times10^{-12}$; the largest final $K,C,M$ entry radius is $0.043993$. These widths belong to the declared synthetic prior and calibration envelope, not a general hardware precision claim.
