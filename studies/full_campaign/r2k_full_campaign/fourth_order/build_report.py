from pathlib import Path
from fractions import Fraction as Q
import json,html,base64
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE=Path(__file__).resolve().parent
def esc(x):return html.escape(str(x))
def main():
    p=json.loads((HERE/'cue_covariance_results.json').read_text())
    v=json.loads((HERE/'periodic_validation.json').read_text())
    proof=json.loads((HERE/'trace_validation.json').read_text())
    cv=json.loads((HERE/'covariance_validation.json').read_text())
    assert p['status']=='complete' and v['status']==proof['status']==cv['status']=='PASS'
    m=p['models']['18000'];ls=p['lengths']
    C=np.array([[x['enclosure']['mid'] for x in r] for r in m['covariance']]);sd=np.sqrt(C.diagonal());cor=C/sd[:,None]/sd[None,:]
    gauss=np.array([m['covariance'][i][i]['Gaussian_pairing_finite']['mid'] for i in range(6)])
    fig,axes=plt.subplots(1,2,figsize=(11,4.5),layout='constrained')
    im=axes[0].imshow(cor,vmin=0,vmax=1,cmap='Blues')
    axes[0].set(xticks=range(6),yticks=range(6),xticklabels=ls,yticklabels=ls,xlabel='Interval length',ylabel='Interval length',title='CUE covariance: strong cross-length dependence')
    for i in range(6):
        for j in range(6):axes[0].text(j,i,f'{cor[i,j]:.2f}',ha='center',va='center',color='white' if cor[i,j]>.6 else '#18354a')
    axes[1].plot(ls,gauss,'o--',label='Gaussian pairing only')
    axes[1].plot(ls,sd**2,'s-',label='Full fourth-order result')
    axes[1].set(xscale='log',yscale='log',xlabel='Interval length',ylabel='Sampling variance of circular statistic',title='Fourth cumulant changes the covariance')
    axes[1].legend();axes[1].grid(alpha=.2)
    fig.suptitle('Explicit circular model, N = S = 18000; no zeta confidence claim')
    fig.savefig(HERE/'cue_covariance.png',dpi=160);fig.savefig(HERE/'cue_covariance.svg');plt.close(fig)
    body=[]
    def sec(title,text):
        body.append('<section><h2>'+esc(title)+'</h2>'+''.join('<p>'+esc(x.strip())+'</p>' for x in text.strip().split('\n\n'))+'</section>')
    body.append('<h1>R2K Nonagon: exact fourth-order calibration for a circular benchmark</h1>')
    sec('Outcome and limits', '''The missing fourth-order input can be computed rigorously for an explicitly specified circular unitary ensemble (CUE) benchmark. This addendum derives its exact trace covariance, reduces the infinite double Fourier sum to polynomial convolutions, encloses the finite arithmetic with Arb and bounds the infinite tail analytically. Three full 6 by 6 covariance matrices cover the observed ranks 17999, 18000 and 18001 on a circle of circumference 18000. No simulation or fitted block scale enters this calculation.

This is a model-conditional theorem. It does not certify the sampling law of actual zeta zeros, the finite-line sine process, or the arithmetic finite-height BK model. The circular empirical statistic introduces a seam and a deterministic full-circle mean. It is a new exploratory diagnostic, not a replacement for the frozen primary finite-window observable. Its conditional CUE mean is visibly inappropriate for these finite-height data: the known arithmetic mean correction is substantial. Adding that correction to the mean while retaining this CUE covariance would require justification.

The original source audit, the full-integral campaign and the finite-window bridge remain in the package. Read their conclusions chronologically: tasks described as pending in the historical audit were subsequently completed in the campaign and bridge. The remaining scientific gap is a justified third/fourth-order law or a separately calibrated sampling experiment, rather than more decimal places in the deterministic integral.''')
    body.append('<p><a href="../R2K_Nonagon_Full_Campaign_Report.html">Main campaign</a> · <a href="../finite_window/R2K_Finite_Window_Bridge_Report.html">Finite-window bridge</a> · <a href="../../R2K_Nonagon_Source_Audit.html">Historical source audit</a></p>')
    sec('1. Exact observable and trace identity', '''Let U be Haar distributed on U(N), with eigenangles mapped to a circle of circumference S. Let N_L(s) count eigenpoints in the circular interval (s,s+L], with 0<L<S/2 and s uniform on the full circle. Its mean over s is deterministically μ=NL/S. For T_j=Tr(U^j) and w_j(L)=[sin(πjL/S)/(πj)]², Parseval gives V_L=circular Var_s[N_L(s)]=2Σ_{j≥1}w_j(L)|T_j|². For arbitrary observed circular points the same identity follows from their Fourier series.

The second trace moment is E|T_j|²=min(j,N). If j,k>0 and a=max(j,k), b=min(j,k), the exact fourth cumulant is κ(j,−j,k,−k)=2(N−a)₊−(N−j−k)₊−(N−|j−k|)₊. Consequently Cov(|T_j|²,|T_k|²)=1_{j=k}min(j,N)²+κ(j,−j,k,−k).

The identity is derived directly from the determinantal generating function log det[I+(exp(Σt_i f_i)−1)K_N]. Differentiating four times gives the sum over ordered partitions π of {1,2,3,4}, with coefficient (−1)^(|π|−1)/|π|. Each trace is (N−range(partial block-frequency sums))₊ because K_N projects onto Fourier modes 0,…,N−1. There are 75 ordered partitions.

The supplied symbolic verifier writes j=k+d, with k,d≥0, and proves each partial-sum maximum and minimum by coefficientwise dominance on that cone. Combining the 75 terms leaves coefficients +2 at span j, −1 at span j+k, and −1 at span j−k. Symmetry covers k≥j. This is a symbolic all-frequency argument; the additional 980 integer checks are controls, not the proof of universality. Degenerate frequencies are included by continuity of the positive-part formula.

When j+k≤N, the fourth cumulant vanishes and the Gaussian trace pairing is exact at this order. For N=1, the cumulant cancels the diagonal term and every |T_j|² is deterministically one. Frequencies beyond the low-frequency Gaussian range are essential in the sharp interval statistic, so dropping this cumulant is not justified by a fixed-frequency limit theorem.''')
    sec('2. Infinite sum, convolution and explicit tail', '''The exact covariance is C_LM=4Σ_{j,k≥1}w_j(L)w_k(M)[1_{j=k}min(j,N)²+κ(j,−j,k,−k)]. All terms are absolutely summable. For a finite cutoff K, write u_j=w_j(L), v_j=w_j(M), with u_0=v_0=0. The finite sum divided by four is D+2A−H−B, where D=Σ_{j=1}^K u_jv_j min(j,N)²; A=Σ_{m=1}^{N−1}(Σ_{j≤m}u_j)(Σ_{k≤m}v_k); H=Σ_{j,k≤K}u_jv_k(N−j−k)₊; and B=Σ_{j,k≤K}u_jv_k(N−|j−k|)₊. Products of Arb polynomials supply H and the cross-correlation coefficients for B. Twenty-five independent exact rational-weight controls match direct ordered-partition double sums.

For K>N, |κ|≤N and κ is zero unless |j−k|<N. Also 0≤w_j≤1/(π²j²). In an omitted row j>K there are fewer than 2N possible k, each at least j−N. Thus the absolute omitted cumulant contribution, including both rows and columns and the outer factor four, is at most 16N²/[3π⁴K³(1−N/K)²]. The diagonal pairing tail is at most 4N²/(3π⁴K³). The sum is the certified covariance tail radius. Overlap of the omitted rows and columns is deliberately overcounted, preserving the upper bound.

The exact mean needs no infinite-tail estimate: writing ℓ=L/S and using Σ_{j≥1}w_j=ℓ(1−ℓ)/2 gives E[V_L]=Nℓ(1−ℓ)+2Σ_{j=1}^{N−1}(j−N)w_j. Every coefficient, convolution and finite sum is evaluated in 128-bit Arb. The cutoff K=180000 is fixed before inspecting residuals. No floating FFT or cutoff-agreement test supplies the proof bound.

A floating Cholesky factor proposes a preconditioning matrix only. Its entries are retained as exact dyadic numbers. Arb proves the matrix invertible, forms the complete interval covariance congruence, and verifies strictly positive Gershgorin lower margins. This certifies positive definiteness for every symmetric matrix in each enclosure. Direct interval LDL suffered excess dependency growth and is not an acceptance criterion. The float eigenvalues are retained solely as a secondary diagnostic.''')
    mx=max(c['enclosure']['radius'] for model in p['models'].values() for row in model['covariance'] for c in row)
    body.append('<p>Largest complete covariance entry radius: <strong>'+f'{mx:.10g}'+'</strong>. Rank-18000 marginal standard deviations and Gaussian-pairing sensitivity:</p><table><tr><th>L</th><th>CUE mean</th><th>Full SD</th><th>Pairing-only variance / full variance</th></tr>')
    for i,L in enumerate(ls):body.append('<tr>'+''.join('<td>'+esc(x)+'</td>' for x in [L,f"{m['mean'][i]['mid']:.10f}",f'{sd[i]:.9f}',f'{gauss[i]/C[i,i]:.6f}'])+'</tr>')
    body.append('</table><p>The pairing-only entry omits a nonzero fourth cumulant; it is a comparison, not another uncertainty certificate.</p>')
    body.append('<img style="width:100%" alt="Circular model correlation matrix and comparison of full versus Gaussian-pairing variances" src="data:image/png;base64,'+base64.b64encode((HERE/'cue_covariance.png').read_bytes()).decode()+'">')
    sec('3. Independent circular empirical checks', '''The 90000 points in the five existing fixed-coordinate intervals are reused; no new zeta ordinates were acquired. For 0<L<S/2, the exact periodic variance is NL/S−(NL/S)²+(2/S)Σ_circular pairs(L−short gap)₊. The sum includes both interior close pairs and pairs crossing the artificial seam. It is evaluated with exact rational arithmetic on the stored dyadic coordinates.

An independent event sweep integrates N_L(s) and N_L(s)² exactly around the circle and agrees with every pair calculation, thirty of thirty. Complete interval membership was already certified in the fixed-window acquisition. The coordinate perturbation bound is 4ε times the number of pairs within L+2ε, divided by S; the maximum over all thirty cases is below 6.203×10⁻¹⁰ using the conservative ε=10⁻¹¹. These are arithmetic and ordinate-location bounds, not sampling errors.

The table below deliberately uses the unmodified CUE mean, since its fourth-order law is the model actually proved. A standardized residual here is merely a model-specific scale comparison. It is not a calibrated Gaussian tail probability; no Gaussian approximation or zeta transfer is established. Observed ranks are inserted as explicit parameters, not justified as the conditional distribution of the zeta points. The geometry was selected after earlier results were known.''')
    body.append('<table><tr><th>Source index</th><th>Rank</th><th>L</th><th>Periodic empirical variance</th><th>CUE mean</th><th>Residual / CUE SD</th></tr>')
    for e in p['empirical_periodized']:
        for r in e['rows']:body.append('<tr>'+''.join('<td>'+esc(x)+'</td>' for x in [e['n0'],e['rank'],r['L'],f"{r['variance']:.10f}",f"{r['CUE_mean']:.10f}",f"{r['standardized_residual']:+.4f}"] )+'</tr>')
    body.append('</table>')
    sec('4. Reconciled project status', '''Uniform R₂ remainder: the old three-sample justification was inadequate; Part X supplied the analytic repair. The campaign extends its origin certificate to 74920≤T≤22316669 and 0<|x|≤10⁻⁴, with safe upper constant 2.987654. This is an infinite-Euler-product analytic-model bound. It is not a uniform zeta-correlation theorem, and c≤3 cannot hold at arbitrarily large height because the leading coefficient tends to π²/3>3.

Whole weighted integral: ninety quantile-level targets and thirty continuous-height comparisons have rigorous model-evaluation bounds. The largest combined deterministic radius in the primary continuous-height campaign is below 2.387×10⁻⁸. Thirty additional fixed-coordinate predictions are also certified. Finite prime cutoff, origin integration, nonzero lower cutoff, height quantile averaging and anchor discretization are now separately controlled. Agreement between finite cutoffs was not used as a tail proof.

Chebyshev–Hermite root counts: the independent full verifier replayed 46 retained contours covering 5023 segments, plus the 841-segment ten-zero fixture. Forty-nine structured moment decisions and one direct fallback are distinct from fresh contour generation. Three much higher tiles at heights 10¹⁰–10¹² remain unreplayed. A local zero count plus positive trace form certifies simple real zeros only in its declared contour; it is not a global Riemann-hypothesis proof. Conditioning and explicit reality/input acceptance checks remain engineering/mathematical hardening issues.

Constants: the old derivative estimate 2M/R is superseded by 4M/R. The printed origin circle majorant 2742.431468 is below its machine upper bound and should be outward-rounded to 2742.431538. The ideal clamp optimum 1/√2, finite (t,2t) bounded-raw-error optimum tG_jj=1, and Gaussian optimum about 0.9721 have different noise/acquisition hypotheses. These must never be interchanged.

Offset: the historical common large-L shift contains the exactly explained missing-origin term 2La−a², equal to 0.00319999 at L=16 and a=10⁻⁴. The fully repaired comparison has smaller residuals, but their statistical significance remains unresolved. Fixed-window and point-selected statistics are not interchangeable, and changing the batch count materially changes nominal standardized discrepancies. The bridge counterexample proves that even all first three correlations do not in general determine this uncertainty. This addendum supplies fourth-order calibration only for its separate circular CUE benchmark.''')
    sec('5. Most decisive next step', '''The decisive mathematical target is the covariance of the original fixed-line statistic under a finite-height arithmetic model with explicit third- and fourth-factorial correlation input and an error bound. The exact weighted contraction is already available in finite_window/moment_bridge.py. Evaluate that contraction at the frozen six lengths and five height windows; include boundary geometry and mean estimation. A theorem for the actual zeta zeros would additionally require a proved weighted correlation estimate or another justified sampling law. The two-point BK expression alone is insufficient.

The immediate computational precursor should be one completely specified fixed-line model and one scalar contraction, starting with L=4 in the source-index-20000000 coordinate window, together with a separate control window. First compute the full sine-DPP fixed-line covariance, including its sample-centering terms, and demonstrate that the circular approximation and boundary error are below a declared tolerance. Then introduce a consistent arithmetic third/fourth-order correction. Merely replacing the circular mean with the BK mean is not such a correction.

For empirical confirmation, freeze new external-coordinate windows and the unfolding, observation contract, statistic, covariance method, and joint decision rule before acquiring further data. Keep independent-window replication separate from conditional model quadrature. More neighboring blocks or more precise two-point quadrature cannot by themselves close the missing fourth-order inference. Until the required law is justified, report model residuals and deterministic enclosures without claiming a newly discovered finite-height term or its absence.''')
    sec('Reproduction and primary background', '''Run validate_trace.py, cue_covariance.py, validate_periodic.py and build_report.py. The covariance computation takes less than one minute on the present runtime; exact empirical controls take longer. The package retains all coefficients needed to reproduce the matrices, analytic tail formulas, exact dyadic moments, positivity certificates and all thirty comparisons. No randomized simulation is required.

Diaconis and Evans, Linear functionals of eigenvalues of random matrices, Transactions of the American Mathematical Society 353 (2001), 2615–2633, provides primary background on Haar-unitary trace moments. The present ordered-partition derivation is included in full and is not claimed as a new random-matrix theorem. Its use here is to make the project's conditional uncertainty calculation explicit and reproducible.''')
    body.append('<p><a href="https://statistics.berkeley.edu/sites/default/files/tech-reports/577.pdf">Primary paper: Diaconis–Evans</a></p>')
    for name in ['trace_derivation.json','trace_validation.json','periodic_validation.json','covariance_validation.json','cue_covariance_results.json']:
        body.append('<details><summary>'+esc(name)+'</summary><pre>'+esc((HERE/name).read_text())+'</pre></details>')
    out='<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>R2K fourth-order calibration</title><style>body{font-family:system-ui,sans-serif;max-width:1120px;margin:36px auto;padding:0 24px;line-height:1.65;color:#18354a;background:#fbfcff}h1,h2{line-height:1.25;color:#16446d}section{margin:2em 0}table{border-collapse:collapse;width:100%;font-size:13px}th,td{border:1px solid #c9d6e4;padding:8px;text-align:left}pre{white-space:pre-wrap;overflow-wrap:anywhere;background:#edf2f7;padding:16px;font-size:12px}summary{cursor:pointer;background:#e4eef8;padding:12px;font-weight:600}</style><main>'+''.join(body)+'</main></html>'
    assert out.count('<table>')==out.count('</table>') and out.count('<details>')==out.count('</details>')
    (HERE/'R2K_Fourth_Order_Calibration_Report.html').write_text(out)
    receipt={'status':'PASS','symbolic_trace_identity':proof['status'],'periodic_exact_matches':v['independent_exact_periodic_event_pair_matches'],
        'covariance_models':len(p['models']),'covariance_entries':3*36,'max_entry_radius':mx,
        'all_covariance_congruence_Gershgorin_margins_positive':True,'zeta_statistical_law_certified':False}
    (HERE/'validation_receipt.json').write_text(json.dumps(receipt,indent=2))
    print(json.dumps(receipt,indent=2))

if __name__=='__main__':main()
