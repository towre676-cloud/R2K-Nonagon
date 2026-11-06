from pathlib import Path
from fractions import Fraction as Q
import json,html,base64
from flint import arb
HERE=Path(__file__).resolve().parent

def main():
 s=json.loads((HERE/'finite_window_summary.json').read_text());counter=json.loads((HERE/'moment_bridge_checks.json').read_text());bias=json.loads((HERE/'sine_bias_results.json').read_text());body=[]
 def sec(title,text):body.append('<section><h2>'+html.escape(title)+'</h2>'+''.join('<p>'+html.escape(x.strip())+'</p>' for x in text.strip().split('\n\n'))+'</section>')
 body.append('<h1>R2K Nonagon: the finite-window comparison and the fourth-order obstruction</h1>')
 sec('Results and scope', '''This addendum proves exact identities linking continuous-anchor variance to weighted pair gaps, verifies them on all thirty existing cases, and constructs thirty additional comparisons on fixed coordinate intervals inside the retained ordinate data. No new zeta ordinates were acquired. The original five-window campaign, its frozen primary protocol, and its results remain retained.

The new observation contract uses the same absolute smooth-unfolded rule [n₀+1000,n₀+19000] for each source index. These intervals were selected after the earlier results were visible; they are exploratory geometry checks, not newly preregistered held-out validation. Their lengths are exactly 18000 and their point counts are allowed to vary. Certified buffer checks establish that the retained consecutive ordinates cover each interval completely.

Two further mathematical results clarify what remains open. First, finite-window centering gives an exact nonnegative downward bias under a stationary unit-intensity model. This bias is evaluated in closed form for the sine-kernel DPP as a conditional benchmark. Second, even knowledge of all first-, second- and third-order correlations does not in general determine the sampling variance of the empirical number-variance statistic. An explicit stationary counterexample and a checked fourth-order moment formula identify the missing input.''')
 sec('1. Exact pair-window closure', '''Let 0≤x₁<⋯<x_n≤S and 0<L<S. Anchors range over [0,D], D=S−L. Define a_i=|[x_i−L,x_i]∩[0,D]|/D and b_ij=|[max(x_i,x_j)−L,min(x_i,x_j)]∩[0,D]|/D, with negative lengths replaced by zero. The exact mean is μ=Σa_i, and the exact variance is V=μ+Σi≠j b_ij−μ².

Writing T=Σi<j(L−x_j+x_i)₊, the clipped pair sum W equals T−E_left−E_right. With zero-based indices, E_left=Σ_{x_j<L}j(L−x_j), and E_right=Σ_{x_i>D}(n−1−i)(x_i−D). These are exact boundary terms. Pair enumeration with integer-scaled dyadic coordinates reproduces the independently event-integrated means and variances exactly on all thirty original cases and all thirty fixed-coordinate cases.

The translation-corrected pair statistic is V_tr=L−L²+2Σi<j (L−r_ij)₊/(S−r_ij). Define δ=μ−L. Then the exact identity is

V−V_tr=(1−2L)δ−δ²+2[W/D−Σi<j(L−r_ij)₊/(S−r_ij)].

This decomposes a model residual into an observed pair-functional discrepancy, a mean term and an edge term. It does not declare any one estimator to be the truth. The terms can cancel strongly, so inspecting a single component as a proposed physical correction is unsafe.

The pair denominator is certified efficiently by 1/(S−r)=Σk=0⁵r^k/S^(k+1)+remainder. All six retained gap moments are exact integers on a common dyadic grid. The positive omitted sum is at most (T/S)(L/S)⁶/(1−L/S). This supplies an enclosure, not just a floating-point pair sum. Coordinate perturbations are bounded using the Lipschitz constants of (L−r)₊/(S−r).''')
 sec('2. Fixed coordinates and density sensitivity', '''Campbell's factorial-moment identity gives E[V_tr]=Σ²(L) for a stationary unit-intensity process observed in a fixed external interval of length S, provided the required pair moments exist. It is a conditional model theorem. The previous windows instead contained exactly 20000 consecutive zeros and ended at observed points. That fixed-count selection does not automatically satisfy the fixed-window theorem.

The new coordinate intervals contain 17999–18001 points. Their means and variances are recomputed exactly, and the full local BK prediction is continuously averaged over their admissible anchor heights using the certified analytic height reduction. Each numerical target includes the ordinate/unfolding enclosure and the infinite-product model radius.

The raw pair statistic is sensitive to the observed global count. For clarity, the data also retain the density-compensated diagnostic V_comp=ρ̂L+2Σ(L−r)₊/(S−r)−ρ̂²L², where ρ̂=n/S. Under a stationary unit-intensity process in a fixed interval, E[V_comp]=Σ²(L)−L²Var[N(S)]/S². This follows by taking expectations of its three terms. Neither this identity nor the fixed-window unbiasedness claim is asserted for the old point-index-selected windows.

All estimator variants are retained. None is selected because it gives the smallest residual, and pair-statistic residuals are not divided by an anchor-variance standard error. The coordinate-window comparison changes both geometry and the included edge data, so a change in its residual alone is not an estimate of a pure selection bias.''')
 sec('3. Exact finite-window centering theorem', '''Let X be a stationary unit-intensity process with population number variance Σ²(L). Independently choose an anchor uniformly in [0,D], D=S−L. For a realization X, let μ_X be its average interval count and V_X its variance over these anchors. The law of total variance gives the exact identity E_X[V_X]=Σ²(L)−Var_X(μ_X). Thus the empirical anchor variance has a nonnegative centering bias even when the stationary model is correct.

The anchor mean is a linear statistic μ_X=Σ_{x∈X}a(x), with a the trapezoidal single-point coverage weight. Its sampling variance uses first- and second-order correlations over the full span S, rather than only correlations at distances up to L. A model known only at short separations cannot in general determine this finite-window mean correction.

For the sine-kernel DPP, the structure factor is min(|k|,1) under the exp(−2πikx) convention. Hence the bias is 2L²∫₀∞min(k,1)sinc²(Dk)sinc²(Lk)dk. The supplied closed form evaluates this integral through the cubic autocorrelation of a: Var(μ_X)=c_a(0)−2∫₀ˢc_a(h)sinc²(h)dh. Its polynomial pieces integrate exactly using Si, Ci, sine, cosine and logarithms. Independent adaptive Fourier integrals on three smaller spans enclose the same values.

At S=18000, this conditional sine-model bias is about 9.35×10⁻¹⁰ for L=0.5 and 6.84×10⁻⁷ for L=16. The small-L benchmark is far too small to explain the observed roughly 10⁻⁴–10⁻³ excursions within that null model. This is not a corresponding theorem for the finite-height BK formula or for actual zeta zeros.''')
 body.append('<table><tr><th>L</th><th>Stationary sine-model centering bias, S=18000</th></tr>')
 for r in bias['production_fixed_span_benchmarks']:body.append('<tr><td>'+str(r['L'])+'</td><td>'+f"{float(arb(r['expected_downward_bias_ball']).mid()):.12g}"+'</td></tr>')
 body.append('</table>')
 sec('4. Why pair correlation cannot calibrate the error bars', '''Consider four potential sites on a circle of circumference 3, spaced by 3/4, with independent positional jitter in [−1/16,1/16] and a common uniform rotation. The rotation makes the process stationary; the jitter gives continuous point locations, while the minimum spacing stays above L=1/2. Let b∈{0,1}⁴ be the occupancy vector and N its number of occupied sites.

Model A assigns probability 3^N/256 to each occupancy configuration. Model B assigns [3^N+(−1)^N]/256. Both distributions are positive and normalized. Summing the added parity interaction over any omitted coordinate gives zero, so all joint occupancy marginals through order three agree. Consequently the two stationary point processes have identical factorial moment measures through order three and mean intensity one.

No interval of length 1/2 contains two points. The continuous-anchor variance over the circle is therefore V=N/6−N²/36. Both models have E[V]=11/48. But Var_A(V)=11/13824, whereas Var_B(V)=1/1152=12/13824. The sampling variances differ by 1/13824 despite identical correlations through third order. Indeed, V²=(N⁴−12N³+36N²)/1296, and the parity perturbation changes E[N⁴] by 24/256 while leaving the lower moments unchanged. This gives the variance difference directly.

Periodically repeating the rotated motif gives a stationary real-line example as well. On an observation interval of length 3.5 with L=0.5, the admissible anchors cover one full period, so the same statistic and variance calculation apply. This periodic example is not a mixing or determinantal null. This proves a general obstruction, not a claim that either construction models the zeta zeros. Additional determinantal structure would determine higher correlations from a kernel, but that structure has not been established for these finite-height zeta windows. A two-point formula by itself does not provide the missing uncertainty certificate.''')
 sec('5. Exact fourth-order covariance formula', '''For fixed observation geometry, set f_L(x)=a_L(x)−a_L(x)² and g_L(x,y)=b_L(x,y)−a_L(x)a_L(y). Then V_L=Σ_i f_L(x_i)+Σ_{i≠j}g_L(x_i,x_j), using ordered pairs. Its mean is m_L=∫f_Lρ₁+∫g_Lρ₂. The formulas are written for factorial moment densities when they exist; the corresponding factorial moment measures give the same identities in general.

For two lengths L and M, the raw cross moment has four terms:

E[V_LV_M]=∫ f_L(x)f_M(x)ρ₁(x)dx
+∫[f_L(x)f_M(y)+2f_L(x)g_M(x,y)+2f_M(x)g_L(x,y)+2g_L(x,y)g_M(x,y)]ρ₂(x,y)dxdy
+∫[f_L(x)g_M(y,z)+f_M(x)g_L(y,z)+4g_L(x,y)g_M(x,z)]ρ₃(x,y,z)dxdydz
+∫g_L(x,y)g_M(z,w)ρ₄(x,y,z,w)dxdydzdw.

Subtract m_Lm_M for the covariance. The coefficients count whether the two ordered pairs share zero, one or two indices. The formula is verified against direct enumeration on twenty correlated finite occupancy laws, with exact rational arithmetic. It also reproduces the counterexample above term by term; only the fourth-order contribution changes.

The useful next target is therefore a certificate for these weighted first-through-fourth-order moment functionals, not necessarily recovery of each entire correlation function. If their error envelopes imply a mean enclosure and an upper covariance matrix C, finite-sample probability bounds can follow from Markov/Chebyshev arguments without Gaussian assumptions. For example, if the mean is m and Cov(V)≼C with C positive definite, then P[(V−m)ᵀC⁻¹(V−m)≥d/α]≤α for d reported lengths. No such higher-moment certificate is currently available for the zeta data.''')
 sec('6. Dependence diagnostics and the unresolved excursion', '''The new uncertainty diagnostic partitions the anchor coordinate, not the list of zero indices. It uses batch means of the global variance influence (N_L(s)−μ_L)². The full anchor interval is covered, including count windows that cross a batch boundary. Five, ten, twenty, forty and eighty batches are all reported.

For each partition, the six lengths are kept together. The ten-batch multiplier reference enumerates all 1024 Rademacher sign patterns. Other partitions use 32768 seeded patterns. This enumeration is computationally exact for the ten-batch reference distribution; it is not an exact significance test for zeta zeros. Independence and sign symmetry of the batches are not established.

The excursion around starting index 20000000 survives the change of coordinates. Under ten anchor batches its largest standardized residual is about 2.87, now at L=4. Across 5/10/20/40/80 batches the fixed-span maximum changes from approximately 4.03/2.87/1.46/1.23/0.97. That instability precludes a robust significance claim from the selected ten-batch result. Numerical enclosures remain much smaller than the residuals.

The proper conclusion is an unresolved empirical excursion with a clearly identified uncertainty gap. It is not an integration defect, proof of a new correction term, or an unconditional rejection of a pair-correlation conjecture. A future confirmation should fix its coordinate spans and uncertainty method before inspecting new data, and should use a justified higher-order model or explicit moment bounds.''')
 if (HERE/'finite_window_diagnostics.png').exists():body.append('<img style="width:100%" alt="Batch sensitivity and conditional sine centering bias" src="data:image/png;base64,'+base64.b64encode((HERE/'finite_window_diagnostics.png').read_bytes()).decode()+'">')
 body.append('<h2>All coordinate-window comparisons</h2><table><tr><th>Index</th><th>Points</th><th>L</th><th>Residual</th><th>Deterministic radius</th><th>10-anchor-batch nominal z</th></tr>')
 for w in s['windows']:
  if w['geometry']!='fixed_span':continue
  for r in w['rows']:
   b=arb(r['residual_ball']);vals=[w['n0'],w['point_count'],r['L'],f'{float(b.mid()):+.8g}',f"{float(arb(r['combined_deterministic_radius']).upper()):.3g}",f"{r['anchor_influence_10batch_z']:+.4f}"]
   body.append('<tr>'+''.join('<td>'+html.escape(str(v))+'</td>' for v in vals)+'</tr>')
 body.append('</table>')
 sec('Primary literature and attribution', '''The translation correction and factorial-moment partition method belong to established point-process statistics. The contribution here is the exact specialization to the Nonagon observation contract, checked finite-window implementation, conditional bias calculation and explicit counterexample.

Biscio and Lavancier, Brillinger mixing of determinantal point processes and statistical applications, arXiv:1507.06506, prove mixing and asymptotic results for stationary DPPs. Proposition 4.2 requires a positive limiting variance in addition to its process and functional assumptions; its pair estimator uses the translation overlap |W∩(W+h)|. This is relevant guidance, not a finite-sample uncertainty theorem for zeta zeros. https://arxiv.org/pdf/1507.06506

Błaszczyszyn, Yogeshwaran and Yukich, Limit theory for geometric statistics of point processes having fast decay of correlations, arXiv:1606.03988, establish expectation, variance and limit results under explicit decay and stabilization assumptions. Those assumptions cannot simply be inferred from a finite pair-correlation fit. https://arxiv.org/abs/1606.03988''')
 for name in ['fixed_preparation.json','moment_bridge_checks.json','sine_bias_results.json','finite_window_summary.json']:
  body.append('<details><summary>'+name+'</summary><pre>'+html.escape((HERE/name).read_text())+'</pre></details>')
 style='body{font-family:system-ui,sans-serif;line-height:1.62;max-width:1120px;margin:32px auto;padding:0 24px;color:#172e43;background:#fbfcff}h1,h2{color:#173f68;line-height:1.25}section{margin:2em 0}table{border-collapse:collapse;width:100%;font-size:13px}td,th{padding:8px;border:1px solid #c8d6e3;text-align:left}pre{white-space:pre-wrap;overflow-wrap:anywhere;background:#edf2f7;font-size:12px;padding:16px}summary{cursor:pointer;background:#e5eef8;padding:12px;font-weight:600}'
 (HERE/'R2K_Finite_Window_Bridge_Report.html').write_text('<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>R2K finite-window bridge</title><style>'+style+'</style><main>'+''.join(body)+'</main></html>')
 print('Finite-window report written')
if __name__=='__main__':main()
