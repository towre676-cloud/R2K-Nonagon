from pathlib import Path
from fractions import Fraction as Q
import json,html
from flint import arb,ctx
HERE=Path(__file__).resolve().parent
ctx.prec=192
def esc(s):return html.escape(str(s))
def main():
    d=json.loads((HERE/'sine_fixed_line_results.json').read_text());v=json.loads((HERE/'validation_receipt.json').read_text())
    assert d['status']=='complete' and v['status']=='PASS'
    body=[]
    def sec(title,text):body.append('<section><h2>'+esc(title)+'</h2>'+''.join('<p>'+esc(t.strip())+'</p>' for t in text.strip().split('\n\n'))+'</section>')
    body.append('<h1>R2K Nonagon: fixed-line fourth-order covariance with sample centering</h1>')
    sec('Result and proof scope', '''The scalar fixed-line calibration proposed in the previous addendum is now complete. Under the stationary unit-intensity sine projection determinantal point process, on a fixed interval of length S=18000 with L=4 and anchors s∈[0,D], D=S−L=17996, the sampling variance of the empirical anchor variance is enclosed by 0.000020941091972389113 ± 0.000000000786. Its sampling standard deviation is enclosed by 0.00457614378843 ± 0.0000000859. These safe rounded enclosures include the effect of estimating the mean over the same anchors.

This is a finite-line result with an explicit boundary weight, not an infinite-span approximation or a circular seam. The calculation encloses the fourth-order covariance integral with 256-bit Arb, tensor Gaussian error bounds, an analytic long-distance expansion and a bound on its remainder. The remaining centering uncertainty is bounded by a projection-operator identity rather than dropped. The leading influence variance has radius below 6.145×10⁻¹⁴; the larger final enclosure comes from bounding the centering covariance.

The theorem is conditional on the sine determinantal model. It does not establish a sampling law for actual zeta zeros, and the sine-model mean is visibly different from the arithmetic finite-height mean in the present data. Only L=4 is certified here. This scalar result does not supply the full six-length fixed-line covariance matrix, a Gaussian tail approximation, or a simultaneous significance statement. The calculations below are conventional analytic arguments with validated arithmetic, not a proof-assistant formalization.''')
    body.append('<p><a href="../R2K_Nonagon_Full_Campaign_Report.html">Main campaign</a> · <a href="../fourth_order/R2K_Fourth_Order_Calibration_Report.html">Previous circular benchmark</a> · <a href="../finite_window/R2K_Finite_Window_Bridge_Report.html">Finite-window identities</a></p>')
    sec('1. Exact observation contract', '''Let N_L(s) count points in (s,s+L]. Define μ=D⁻¹∫₀ᴰN_L(s)ds, V=D⁻¹∫₀ᴰ[N_L(s)−μ]²ds, W=D⁻¹∫₀ᴰ[N_L(s)−L]²ds and Y=μ−L. The exact identity is V=W−Y². The sine kernel is K(x,y)=sin(π(x−y))/(π(x−y)), with its removable value one at x=y. It is the projection onto the Fourier band [−1/2,1/2]. Its intensity is one, so Eμ=L and EY=0.

For A=[0,L], B=[h,h+L], let X=N(A)−L and Z=N(B)−L. The influence covariance is F(h)=Cov(X²,Z²)=κ₂₂(A,B)+2 Cov(X,Z)². By stationarity and reflection, Var(W)=2D⁻²∫₀ᴰ(D−h)F(h)dh. The factor D−h is the exact finite anchor-overlap length. No large-D substitution is made.''')
    sec('2. Joint cumulant, including overlapping intervals', '''Write C=A∩B. For a cyclic word U₁⋯U_r in A,B,C, let T(U₁⋯U_r)=Tr(M_U₁ K ⋯ M_U_r K), where M_U multiplies by the interval indicator. Differentiating the determinantal log moment-generating function yields the ordered-partition cumulant formula. Combining all 75 ordered partitions for tags A,A,B,B gives

κ₂₂ = T(C)−T(AB)−2T(CC)−2T(AC)−2T(BC)+2T(ABB)+2T(AAB)+8T(ABC)−4T(AABB)−2T(ABAB).

The count covariance is T(C)−T(AB). Reversal and cyclic rotation of words preserve their traces because the sine kernel is real symmetric. exact_controls.py proves the word coefficients symbolically and checks the resulting squared-count covariance against twenty exact rational projection-DPP occupancy laws. The checks are independent controls of the operator derivation.

For h≥L, C is empty. Let A₀ be the kernel compression to one length-L interval, C_h the cross-kernel between the two disjoint intervals, and J=I−2A₀. The covariance simplifies to F(h)=−Tr(J C_h J C_h*)+2{[Tr(C_h C_h*)]²−Tr[(C_h C_h*)²]}. This avoids subtraction of many separately large traces and makes the long-distance expansion tractable.''')
    sec('3. Certified spatial and separation integration', '''The kernel traces are evaluated with positive-weight, ball-enclosed 48-node Gauss–Legendre rules on each interval piece. For 0<h<L the union is partitioned into the left-only interval, intersection, and right-only interval. For h≥L there are two disjoint intervals. Every trace is a tensor quadrature of at most four spatial variables; matrix multiplication evaluates that tensor rule efficiently.

The entire-kernel representation sinc(πz)=∫₋₁/₂¹/₂exp(2πitz)dt gives |sinc(πz)|≤exp(π|Im z|). Each spatial variable occurs in at most two kernel factors. A complex circle of radius L therefore supplies a majorant exp(2πL), while the ratio of every interval half-length to that radius is at most 1/2. Gaussian exactness through degree 2n−1 and Cauchy estimates bound a one-variable integral-plus-quadrature tail by twice its interval length times Mq^(2n)/(1−q).

Telescoping the tensor rules, and using the absolute coefficient counts (1,9,16,8) for the one-, two-, three- and four-variable terms in F, yields the uniform bound E_x=4 exp(2πL) 2^(−2n) Σ_{r=1}⁴ c_r r L^r. At n=48 and L=4 this is below 4.798×10⁻¹⁴. It is valid uniformly at every real separation. Since 2D⁻²∫₀ᴰ(D−h)dh=1, the same bound controls its contribution to Var(W).

For 0≤h≤64, separation integration uses 128 panels of width 1/2, each with a 20-node Gaussian rule. Every moving spatial node is affine in h with slope in [0,1]. Although square roots enter the real matrix representation, they cancel around every closed trace, leaving an entire expression with polynomial Jacobians. This polynomial continuation, rather than a chosen branch of those square roots, supplies the complex majorant.

On a separation circle of radius ρ=1, the total absolute quadrature weight for each interval is at most U=L+2ρ. With e=exp(πρ), a safe majorant for F is Ue+7U²e²+12U³e³+6U⁴e⁴+2(Ue+U²e²)². Multiplication by D−h costs at most D+64+ρ. Summing Gaussian Cauchy errors with ratio q=1/4 gives a variance contribution below 4.750×10⁻¹⁷. All finite evaluation and roundoff errors are enclosed by Arb.''')
    sec('4. Long-distance integral without an arbitrary truncation', '''For h≥H=64, write d=y−x, |d|≤L, and expand 1/(h+d) in eight geometric terms. The cross-kernel becomes a finite Laurent polynomial in 1/h with harmonics exp(±iπh). Its entrywise remainder has magnitude at most (L/h)^8/[π(h−L)]. The quadratic and quartic terms in F therefore have only harmonics 0, ±2π, ±4π and powers h^(−p), 2≤p≤32.

All coefficient matrices are computed in Arb/ACB. For I_p(ω)=∫_H^D exp(iωh)h^(−p)dh, the zero-frequency case uses powers and a logarithm. For ω>0, I₁=Ci(ωD)−Ci(ωH)+i[Si(ωD)−Si(ωH)], and I_p={iωI_(p−1)−[exp(iωD)/D^(p−1)−exp(iωH)/H^(p−1)]}/(p−1). Negative frequencies follow by conjugation. The weighted integral is D I_p−I_(p−1). Five independent adaptive complex integrals, including power 32, agree with these closed forms.

Here is the explicit remainder argument. Let c=L/[π(h−L)], q=(L/h)^8. The Hilbert–Schmidt norm of the cross-kernel is at most c, and that of the expansion remainder is at most cq. For the positive-weight spatial quadrature, ||J||≤1+2L. Thus the quadratic covariance error is at most (1+2L)²c²(2q+q²), and the quartic covariance error at most 4c⁴[(1+q)⁴−1]. Both are integrated explicitly after bounding h−L≥(1−L/H)h and D−h≤D. The resulting variance error is below 1.343×10⁻¹⁴.

Four direct cross-kernel controls at h=64, 64.25, 65 and 96 satisfy the displayed pointwise remainder bound. The far integral is negative, about −6.91164034×10⁻¹⁰; its sign is retained rather than replaced by an absolute truncation estimate. A failed finite enclosure is rejected. Integer powers of zero-containing real balls in the geometric coefficients are evaluated by repeated multiplication to avoid nonfinite generic-power evaluations.''')
    sec('5. Rigorous sample-centering theorem', '''A projection-DPP identity controls the centering term without evaluating a second large fourth-order integral. Let f be a bounded real test function taking values in an interval of width m, P the kernel projection and Q=I−P. Put A=PfP, B=PfQ, C=QfQ. For the centered linear statistic Y_f, β=Var(Y_f)=Tr(BB*) and

κ₄(Y_f)=||AB−BC||²_HS−2Tr[(BB*)²].

This follows by expanding the standard trace expression Tr(f⁴P)−4Tr(f³PfP)−3Tr(f²Pf²P)+12Tr(f²PfPfP)−6Tr[(fP)⁴] in the P⊕Q blocks. Subtracting the midpoint of the range of f leaves B and AB−BC unchanged and makes ||A||,||C||≤m/2. Hence ||AB−BC||_HS≤m||B||_HS, and κ₄≤m²β. Consequently Var(Y_f²)=2β²+κ₄≤2β²+m²β. These identities extend from finite projections under the finite-moment/Hilbert–Schmidt conditions satisfied by the present compactly supported trapezoidal weight. No claim of novelty is made.

For the anchor mean, f is the interval-coverage trapezoid supported on [0,S], with maximum m=L/D. The previous exact sine-bias calculation supplies β=4.961227792477709…×10⁻⁸. Thus sd(Y²)≤sqrt(2β²+m²β)<8.588×10⁻⁸. Minkowski's inequality applied to V=W−Y² gives |sd(V)−sd(W)|≤sd(Y²). This is the source of the final enclosure, and it includes the unknown covariance between W and Y² rather than assuming it vanishes.

Twenty exact finite-projection occupancy laws independently verify the block identity and the centering variance inequality. The mean is also retained: E[V]=Σ²_sine(4)−β=0.4863261172616758… . The much smaller deterministic centering bias is distinct from the standard-deviation perturbation bound.''')
    body.append('<h2>Numerical acceptance ledger</h2><table><tr><th>Quantity</th><th>Certified value or bound</th></tr>')
    for key in ['near_influence_variance','far_influence_variance','spatial_tensor_error_bound','separation_quadrature_error_bound','long_distance_polynomial_error_bound','influence_variance_enclosure','centering_bias','centering_sd_perturbation_upper','anchor_variance_sampling_variance_enclosure','anchor_variance_sampling_sd_enclosure']:
        body.append('<tr><td>'+esc(key)+'</td><td>'+esc(d[key]['ball'])+'</td></tr>')
    body.append('</table>')
    sec('6. Circular comparison and actual data', '''Comparing the full fixed-line variance enclosure with the previous rank-18000 CUE circle certificate gives |Var_line(V)−Var_CUE(V_circ)|<8.575×10⁻⁹. This bounds the total difference between these two declared model-observable pairs at this S and L. It is not an isolated estimate of a universal boundary correction: circular topology, fixed rank and mean convention all differ. The bound includes the earlier CUE certificate's Fourier-tail uncertainty.

The table retains the actual fixed-coordinate anchor variances for all five windows, together with both model means. No new zeta ordinates were acquired. In particular, the source-index-20000000 case and the 10000000 control reuse existing exploratory windows. Their selection is not retrospectively relabeled as preregistered confirmation. The BK model mean is shown for context; it is not assigned the sine-model covariance as a calibrated error bar.''')
    bk=json.loads((HERE.parent/'finite_window/fixed_continuous_height_results.json').read_text())
    lookup={r['n0']:r for r in bk['results'] if r['L']==4}
    body.append('<table><tr><th>Source index</th><th>Empirical fixed-line V</th><th>Sine-model E[V]</th><th>Arithmetic BK mean</th><th>Empirical minus BK</th></tr>')
    comparisons=[]
    for path in sorted((HERE.parent/'finite_window').glob('fixed_data_*.json')):
        data=json.loads(path.read_text());r=next(r for r in data['rows'] if r['L']==4);actual=float(Q(r['variance_exact']));pred=arb(lookup[data['n0']]['infinite_product_height_average_ball']);bkmid=float(pred.mid())
        comparisons.append({'n0':data['n0'],'L':4,'empirical_variance_exact':r['variance_exact'],'empirical':actual,'sine_mean':d['expected_empirical_anchor_variance']['mid'],'BK_mean_ball':str(pred),'residual_BK':actual-bkmid})
        body.append('<tr>'+''.join('<td>'+esc(x)+'</td>' for x in [data['n0'],f'{actual:.10f}',f"{d['expected_empirical_anchor_variance']['mid']:.10f}",f'{bkmid:.10f}',f'{actual-bkmid:+.8f}'])+'</tr>')
    body.append('</table>');(HERE/'fixed_line_comparisons.json').write_text(json.dumps(comparisons,indent=2))
    sec('7. What is now closed, and the next decisive target', '''The scalar sine-model geometry task is closed: a finite-line fourth-order variance, a rigorous sample-centering allowance, and a quantified comparison with the circular benchmark are available for L=4 and S=18000. The numerical fixed-line approximation no longer relies on block independence or a chosen block count. A complete six-length matrix remains open.

The principal scientific obstacle is now the arithmetic correlation model. The BK expression fixes a two-point mean but does not uniquely determine a compatible third/fourth-order law. A consistent four-point input is needed before assigning a justified covariance to the finite-height BK residual. The finite_window/moment_bridge.py contraction already specifies how such input must enter the exact observation geometry.

The next mathematical step is to derive and evaluate the weighted connected four-point arithmetic correction for this same scalar target, including all lower-order coincidences, and retain explicit hypotheses such as a ratios conjecture where required. An immediate check is that the proposed rho₁,…,rho₄ recover the supplied two-point BK mean and produce a nonnegative covariance. One must also establish existence/consistency of any claimed point-process model; a list of formal correlations is not automatically a stochastic law. A theorem for actual zeta zeros would require the appropriate correlation estimate or another justified transfer argument.

Further fixed-line lengths and cross-covariances can now use the same machinery, but they do not by themselves supply that arithmetic information. New external-coordinate windows with a frozen observation and decision rule remain valuable empirical confirmation. Neither this conditional sine certificate nor the small line-circle difference proves the presence or absence of an additional finite-height term.''')
    sec('Reproduction and sources', '''Run exact_controls.py, sine_covariance.py, validate_fixed_line.py and build_report.py. The production quadrature takes about two minutes on the current runtime; the validation replay is much shorter. The default calculation regenerates the near integral. The optional --resume path checks its complete parameter contract. All 128 panel enclosures and long-distance harmonic coefficients are retained, alongside independent exact controls.

Background: Alexander Soshnikov, Determinantal random point fields, Russian Mathematical Surveys 55 (2000), 923–975; and Gaussian limit for determinantal random point fields, Annals of Probability 30 (2002), 171–187. The present application uses the determinantal cumulant framework and gives its particular trace reductions and bounds explicitly. No asymptotic central-limit theorem is used as a finite-sample probability guarantee.''')
    body.append('<p><a href="https://arxiv.org/abs/math/0002099">Determinantal random point fields</a> · <a href="https://arxiv.org/abs/math/0006037">Gaussian limit</a></p>')
    for name in ['exact_controls.json','validation_receipt.json','sine_fixed_line_results.json']:
        body.append('<details><summary>'+esc(name)+'</summary><pre>'+esc((HERE/name).read_text())+'</pre></details>')
    out='<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>R2K fixed-line fourth-order certificate</title><style>body{font-family:system-ui,sans-serif;max-width:1150px;margin:36px auto;padding:0 24px;color:#17374c;line-height:1.65;background:#fbfcff}h1,h2{color:#16446d;line-height:1.25}section{margin:2em 0}table{border-collapse:collapse;width:100%;font-size:13px}td,th{border:1px solid #c9d6e4;padding:8px;text-align:left;overflow-wrap:anywhere}pre{white-space:pre-wrap;overflow-wrap:anywhere;background:#edf2f7;padding:16px;font-size:12px}summary{background:#e3edf7;padding:12px;font-weight:600;cursor:pointer}</style><main>'+''.join(body)+'</main></html>'
    assert out.count('<table>')==out.count('</table>') and out.count('<details>')==out.count('</details>')
    (HERE/'R2K_Fixed_Line_Covariance_Report.html').write_text(out)
    print('Report complete; five L=4 comparisons retained')

if __name__=='__main__':main()
