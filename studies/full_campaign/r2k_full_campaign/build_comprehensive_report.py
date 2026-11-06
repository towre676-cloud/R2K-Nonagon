from pathlib import Path
import json,html,base64,sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from flint import arb,ctx
HERE=Path(__file__).resolve().parent;ctx.prec=192
def esc(x):return html.escape(str(x))
def read(name):return json.loads((HERE/name).read_text())

def main():
    M=read('fixed_line/full_matrix_results.json');V=read('fixed_line/full_matrix_validation.json')
    O=read('arithmetic_fourpoint/dpp_obstruction.json');J=read('arithmetic_fourpoint/ratios_jet_results.json')
    JV=read('arithmetic_fourpoint/ratios_validation.json');CB=read('arithmetic_fourpoint/connected_bridge_results.json')
    D=read('arithmetic_fourpoint/collision_diagnostics.json');H=read('root_hardening/hardening_checks.json')
    low=read('root_hardening/retained_full_replay.json');source=read('source_repairs/validation_receipt.json')
    catalogue=read('root_hardening/catalogue_validation.json');assert catalogue['status']=='PASS'
    assert M['status']==J['status']==low['status']=='complete'
    assert all(d['status']=='PASS' for d in [V,O,JV,CB,H,source])
    high=[]
    for name,fn in [('V_e10','V_e10_replay.json'),('V_e12a','V_e12a_parallel_replay.json'),('V_e12b','V_e12b_parallel_replay.json')]:
        r=read('root_hardening/'+fn);assert r['status']=='CERTIFIED' and r['segments_checked']==r['segments_recomputed']
        r['name']=name;high.append(r)
    assert len(low['results'])==46 and all(r['status']=='CERTIFIED' for r in low['results'].values())
    total_segments=sum(r['segments_recomputed'] for r in low['results'].values())+sum(r['segments_recomputed'] for r in high)
    A=np.array([[c['mid'] for c in row] for row in M['anchor_covariance']]);sd=np.sqrt(A.diagonal());cor=A/sd[:,None]/sd[None,:]
    fig,axes=plt.subplots(1,2,figsize=(11,4.4),layout='constrained')
    axes[0].imshow(cor,cmap='Blues',vmin=0,vmax=1)
    axes[0].set(xticks=range(6),yticks=range(6),xticklabels=M['lengths'],yticklabels=M['lengths'],xlabel='Interval length',ylabel='Interval length',title='Fixed-line sine-model correlation')
    for i in range(6):
        for j in range(6):axes[0].text(j,i,f'{cor[i,j]:.2f}',ha='center',va='center',color='white' if cor[i,j]>.6 else '#18354a')
    heights=[float(arb(r['T_ball']).mid()) for r in O['witnesses']];overs=[float(arb(r['R2_minus_one']).mid()) for r in O['witnesses']];errors=[float(arb(r['R2_minus_one']).rad()) for r in O['witnesses']]
    axes[1].errorbar(heights,overs,yerr=errors,fmt='o',color='#b44b27');axes[1].axhline(0,color='#46586a',lw=1)
    axes[1].set(xscale='log',xlabel='Declared model height T',ylabel='Certified BK R₂ − 1 at a retained witness',title='BK exceeds the Hermitian-DPP ceiling')
    axes[1].grid(alpha=.2);fig.suptitle('Model certificates; neither panel is a zeta sampling-law theorem')
    fig.savefig(HERE/'comprehensive_advance.png',dpi=165);fig.savefig(HERE/'comprehensive_advance.svg');plt.close(fig)
    body=[]
    def sec(title,text):body.append('<section><h2>'+esc(title)+'</h2>'+''.join('<p>'+esc(p.strip())+'</p>' for p in text.strip().split('\n\n'))+'</section>')
    body.append('<h1>R2K Nonagon: complete covariance, arithmetic four-point machinery and root-certificate closure</h1>')
    sec('Outcome and current proof status', '''This extension completes the six-length fixed-line sine-model covariance calculation, develops an independently checked finite-prime arithmetic three/four-point evaluator, proves a specific obstruction to a Hermitian determinantal completion of BK, and closes the three previously deferred high-contour replays. It also hardens the verifier and supplies reviewable source corrections. The original archive is unchanged; all earlier campaign members are retained.

The covariance certificates are conditional on the sine projection determinantal process. The arithmetic computations enclose their declared finite Euler formulas; their interpretation as infinite-prime, finite-window zeta statistics remains provisional. No calibrated significance is assigned to the surviving BK residual, and no new universal finite-height correction is claimed.''')
    body.append('<table><tr><th>Workstream</th><th>Established now or retained</th><th>Open limit</th></tr>')
    ledger=[
        ('Uniform R₂ remainder','Retained infinite-product analytic-model bound R₂(x)/x²<2.987654 on 0<|x|≤10⁻⁴, 74920≤T≤22316669','Actual zeta/model identification; no all-height c≤3 theorem'),
        ('Weighted number-variance integral','Retained ninety quantile-level targets and continuous-height certificates; primary combined deterministic radius <2.387×10⁻⁸','Sampling law and model discrepancy are separate'),
        ('Fixed-line fourth order','All 21 distinct entries at L=.5,1,2,4,8,16, span 18000; rigorous centering allowances and positive matrices','Applies to the declared sine process'),
        ('Root certificates','Fresh hardened full replay of all 49 original certificates in this extension; high tiles now included','Local contour statements; no global RH consequence'),
        ('Arithmetic completion','Five infinite-product witnesses exclude exact scalar Hermitian-DPP completion of BK','Does not exclude general processes or weak/asymptotic approximations'),
        ('Arithmetic three/four point','Finite-prime ratios compiler, cutoffs 31–10007, independent RMT and pair controls','Infinite-prime derivative tails, collision-uniform quadrature and ratios/short-window errors'),
        ('General covariance bridge','Exact connected-measure identity, all fifteen collision patterns, 21 explicit error-transfer targets','Requires justified arithmetic cumulant input')]
    for row in ledger:body.append('<tr>'+''.join('<td>'+esc(x)+'</td>' for x in row)+'</tr>')
    body.append('</table><img style="width:100%" alt="Six-length fixed-line model correlation and certified BK values above one" src="data:image/png;base64,'+base64.b64encode((HERE/'comprehensive_advance.png').read_bytes()).decode()+'">')
    sec('1. Complete fixed-line covariance, without a seam', '''For length L let D_L=S−L and W_L=D_L⁻¹∫₀ᴰᴸ[N_L(s)−L]²ds. Let Y_L be the anchor mean minus L. The actual empirical anchor variance is V_L=W_L−Y_L². All six lengths use the same fixed observation span S=18000; their admissible anchor ranges differ.

Center the count intervals. For A=[−L/2,L/2], B=[r−M/2,r+M/2], define F_LM(r)=Cov((N(A)−L)²,(N(B)−M)²). Reflection makes F even. Write δ=|L−M|/2 and R=S−(L+M)/2. Then Cov(W_L,W_M)=2/(D_LD_M)∫₀ᴿ w(r)F_LM(r)dr, where w(r)=min(D_L,D_M) on [0,δ] and w(r)=R−r on [δ,R]. This trapezoidal overlap weight is exact.

The joint cumulant reduction in the scalar report applies without requiring L=M. Spatial quadrature uses separate Gaussian grids on A, B and their intersection; cross-kernel block products evaluate the traces. The analytic circle radius for a domain of length ℓ is ρ=n/π. The per-variable tensor error factor is 2 exp(2n)[πℓ/(2n)]^(2n)/[1−πℓ/(2n)]. Orders are selected to put this factor below 10⁻²⁰. The actual positive-weight ball nodes and all matrix operations are evaluated with FLINT/Arb.

Separation panels include every containment, touching and overlap-weight transition. The 20-node separation rule has half-width at most 1/4 and Cauchy radius one. Beyond H=max(64,ceil(4(L+M))), the same eight-term cross-kernel geometric expansion gives harmonics 0, ±2π, ±4π, integrated analytically to R. The explicit remainder uses |y−x|≤(L+M)/2 and the factor 2/max(D_L,D_M). Every entry retains its panels, harmonic coefficients and three distinct deterministic error bounds.

The centering theorem proved in the scalar report gives q_L²=2β_L²+(L/D_L)²β_L as an upper bound for Var(Y_L²), where β_L is the certified sine-model anchor-mean variance. Cauchy–Schwarz therefore gives |Cov(V_L,V_M)−Cov(W_L,W_M)|≤sd(W_L)q_M+sd(W_M)q_L+q_Lq_M. It includes the unknown cross-covariances rather than setting them to zero. All 21 entries now carry this allowance.

Thirty independent Python/native count controls agree, including containment, partial overlap, disjoint intervals and the largest length. The earlier separately computed L=4 production certificate overlaps the new diagonal entry. Exact-dyadic preconditioning followed by interval Gershgorin proves positivity of both complete covariance enclosures. A numerical Cholesky factor proposes the basis only; it does not decide acceptance. The optional Linux native accelerator can be replaced by the slower Python ball implementation.''')
    body.append('<table><tr><th>L</th><th>Mean empirical anchor variance, sine model</th><th>Sampling SD midpoint</th><th>Variance enclosure radius</th></tr>')
    for i,l in enumerate(M['lengths']):body.append('<tr>'+''.join('<td>'+esc(x)+'</td>' for x in [l,f"{M['expected_anchor_variances'][i]['mid']:.10f}",f'{sd[i]:.9f}',M['anchor_covariance'][i][i]['radius']])+'</tr>')
    body.append('</table><h3>Complete covariance midpoints (multiply entries by 10⁻⁶)</h3><table><tr><th>L / M</th>'+''.join('<th>'+esc(l)+'</th>' for l in M['lengths'])+'</tr>')
    for i,l in enumerate(M['lengths']):body.append('<tr><th>'+esc(l)+'</th>'+''.join('<td>'+f'{A[i,j]*1e6:.7f}'+'</td>' for j in range(6))+'</tr>')
    body.append('</table><p>Midpoints are a display aid. The complete entry balls and centering bounds are in <a href="fixed_line/full_matrix_results.json">full_matrix_results.json</a>; acceptance is in <a href="fixed_line/full_matrix_validation.json">full_matrix_validation.json</a>.</p>')
    sec('2. A simultaneous model bound without a Gaussian approximation', '''The validator constructs an exact rational matrix U that dominates every covariance in the certified enclosure in Loewner order. It takes a dyadic midpoint matrix and adds the row sums of the complete entry-error bounds to its diagonal. Symmetric diagonal dominance proves the padding dominates every allowed perturbation; interval arithmetic independently proves U positive definite. The exact rational entries are retained, so this is a specified matrix rather than an ambiguous rounded display.

For a future fixed observation obeying this sine model, Markov's inequality gives P{(V−EV)ᵀU⁻¹(V−EV)≥6/α}≤α. At α=.05 the quadratic threshold is 120. This finite-model statement needs neither Gaussianity nor independent lengths. It is deliberately conservative. It supplies no post-selection coverage for the exploratory zeta windows and no confidence region for a BK mean paired with a sine covariance.''')
    sec('3. Root-count replay and acceptance hardening', '''The independent full verifier recomputes the zero-free contour segments and moment balls, then tests the Chebyshev–Hermite form. This extension freshly replays all forty-nine original Part VI certificates with the hardened checker, instead of relying only on the earlier forty-six-contour receipt. It also freshly replays the ten-zero stress fixture (841 segments, inertia (10,0)) and the five-zero union fixture (147 segments, inertia (5,0)). These two controls are separate from the forty-nine-certificate catalogue.

The three previously deferred tiles now have full analytic replay. The high computations parallelize independent segment enclosures with atomic checkpoints; they do not substitute stored balls for recomputation. Every final moment uses the verifier's own enclosures. Fresh is the default replay mode; --resume is explicitly a recovery path and must not be described as independent regeneration of its stored prefix.

The checker now rejects nonreal full-contour moment enclosures, negative zero counts, nonfinite or negative-radius recorded balls, empty/boolean moment indices, malformed flags, invalid evaluator selections, invalid orders and insufficient DFT length. Sixteen adversarial inputs are rejected, forty-nine retained light decisions remain accepted, and positive/empty count controls pass. Light acceptance remains stored-evidence checking; full acceptance requires analytic recomputation. The original verifier, all forty-nine sealed source certificates and a reviewable patch are included.

These are local statements within the exact binary64 rectangles encoded by the certificates. A successful real-root criterion establishes simple critical-line zeros inside those rectangles, not at other heights or outside the displayed real-part range. The Davenport–Heilbronn negative control remains a valid off-line example; success of the catalogue does not mean every supported function has only real zeros.''')
    body.append('<table><tr><th>Newly closed tile</th><th>Recomputed segments</th><th>Count</th><th>Inertia</th><th>Decision</th></tr>')
    for r in high:body.append('<tr>'+''.join('<td>'+esc(x)+'</td>' for x in [r['name'],r['segments_recomputed'],r['count'],r['inertia'],r['status']])+'</tr>')
    body.append('</table><p>Fresh full catalogue total: <strong>'+str(total_segments)+'</strong> segments; including the two separate fixtures, 51 full replays and 6434 segments. The two adjacent 10¹² tiles contain 3+1=4 simple line zeros. Their common boundary is zero-free by the component certificates.</p>')
    sec('4. A proved obstruction to an exact Hermitian-DPP BK kernel', '''For a scalar Hermitian determinantal process with unit intensity, the pair density is ρ₂(x,y)=1−|K(x,y)|²≤1. Five retained BK-model witnesses, one in each production height window, instead have R₂>1. Their balls include rigorous infinite-prime tails at cutoff 2000000. Thus no exact scalar Hermitian determinantal completion can have this BK pair function at those points. Stationarity is not needed for the inequality when the diagonal intensity is one.

For the Euler factor A, every omitted prime contributes |u_p|≤4/(p−1)². With the retained certified sums S₀=Σ_{p>P}(p−1)⁻² and S₂=Σ_{p>P}(log p)²/(p−1)², the log-A tail is at most 4S₀/(1−4/P²), and the B tail at most S₂. These propagate through the exact regularized two-point evaluator. The overshoots below are orders of magnitude larger than the resulting numerical radii.

This also gives a quantitative pointwise approximation obstruction: any unit-intensity Hermitian-DPP pair function has sup error at least the displayed overshoot on a domain containing that witness. It does not rule out a general process, non-Hermitian representation, suitable mixture with varying intensity, leading-order approximation, or agreement of selected weighted integrals. In particular, it does not refute the BK asymptotic conjecture for zeta zeros.''')
    body.append('<table><tr><th>Source index</th><th>Unfolded separation</th><th>Certified R₂ − 1</th></tr>')
    for r in O['witnesses']:body.append('<tr><td>'+esc(r['n0'])+'</td><td>'+esc(r['x'])+'</td><td>'+esc(r['R2_minus_one'])+'</td></tr>')
    body.append('</table>')
    sec('5. Arithmetic three/four-point evaluator: what is actually computed', '''Conrey and Snaith (2008), Correlations of eigenvalues and Riemann zeros, gives higher correlations conditional on the ratios conjecture and RH. We use its local small-shift version as a computational starting point. It does not provide an unconditional short-window theorem or an explicit numerical remainder for this campaign. The paper also distinguishes the exact functional-equation factor from its exponential approximation. Our local evaluator uses exp(−ℓΣ shifts), consistent with the present BK two-point model; the omitted functional-equation correction is not included in its finite-arithmetic balls.

The implemented compiler uses square-zero mixed jets with up to four independent variables: ε_i²=0. Sixteen coefficients encode every distinct mixed derivative through order four. It differentiates only numerator shifts; denominator shifts are held fixed. Zeta poles in reciprocal factors are handled through H(z)=zζ(1+z), so 1/ζ(1+z)=z/H(z) remains analytic at zero. Balanced swaps supply the oscillatory exponential terms.

At each prime, the local arithmetic factor is the constant Laurent coefficient of ∏_c(1−u_c z)∏_d(1−v_d/z)/[∏_a(1−u_a z)∏_b(1−v_b/z)], divided by its universal first-order Euler normalization. With equal numerator/denominator sizes, this coefficient is the sum of the residue at zero and the simple residues at the B roots. Independent power-series constant-term sums through degree 256, with a Cauchy bound for all omitted coefficients, validate all sixteen mixed coefficients at three primes. Nine separate interval comparisons match the CUE determinant for orders two, three and four, including a rank-two, order-four zero determinant. Five independent elementary BK pair evaluations agree at the tested prime cutoffs.

For order n, the correlation assembly sums disjoint assignments K,L,M of the n point labels, with weight ℓ^{|M|} times the mixed logarithmic-derivative term for K and L, divided by ℓⁿ. At n=4 this includes twelve pair terms, twenty-four triple terms and fourteen genuine four-derivative assignments, besides the unit-density term. This implementation retains the arithmetic factors rather than imposing a sine determinant.

The table uses unfolded points (0,.37,1.13,2.41), their initial subsets for lower orders, and the fixed local height near 9503636.190163. All displayed balls rigorously enclose the stated finite-prime formula only. Prime-cutoff changes are evidence about convergence, not bounds for the infinite product or for the ratios conjecture. The four-point cutoff changes are nonmonotone, so reporting only a favorable pair of cutoffs would be misleading. Symmetry and sampled positivity do not establish existence of a joint point process; consistency remains an additional requirement for a probabilistic interpretation.''')
    body.append('<p><a href="https://arxiv.org/abs/0803.2795">Primary source and hypotheses: Conrey–Snaith</a></p><table><tr><th>Prime cutoff</th><th>Finite-input ρ₂</th><th>Finite-input ρ₃</th><th>Finite-input ρ₄</th></tr>')
    for P in [31,97,379,997,10007]:
        vals=[next(r for r in J['arithmetic_results'] if r['prime_cutoff']==P and r['order']==n) for n in [2,3,4]]
        body.append('<tr><td>'+str(P)+'</td>'+''.join('<td>'+f"{float(arb(r['finite_Euler_correlation_ball']).mid()):.12f}"+'</td>' for r in vals)+'</tr>')
    body.append('</table>')
    sec('6. Collision evidence and exact general covariance contraction', '''At cutoff 97, translation and point-permutation controls agree. With points (0,ε,1.13,2.41), ε=.1,.01,.001, the four-point value decreases approximately quadratically; its ratio to ε² approaches a stable scale. These are pointwise evidence for the expected collision cancellation. They are not a uniform removable-singularity certificate over the integration cells, and they do not bound the infinite-prime derivative tail.

A new exact identity supplies the proper covariance target without assuming a determinantal model. Let a_L(x) be the coverage fraction of anchors containing x, b_L(x,y) the joint coverage fraction, and g_L(x,y)=b_L(x,y)−a_L(x)a_L(y). Then ∫g_L(x,y)dy=0, because ∫b_L(x,y)dy=L a_L(x) and ∫a_L(y)dy=L. With η=N−dx and unit intensity, V_L=∫g_L(x₁,x₂)η(dx₁)η(dx₂), including its diagonal atoms.

It follows exactly that Cov(V_L,V_M)=∫g_L12 g_M34 [κ₄^ordinary+κ₂^ordinary(13)κ₂^ordinary(24)+κ₂^ordinary(14)κ₂^ordinary(23)]. Ordinary fourth cumulants include every coincidence pattern. They convert to factorial cumulants by summing over the fifteen set partitions of four labels and identifying all positions in a block. The block counts by factorial order are 1,7,6,1. Thus the compact identity still needs third-order factorial information; that information is contained in the collapsed fourth ordinary cumulant.

The conversion from factorial densities is explicit: k₂=ρ₂−1; k₃=ρ₃−Σ_pairsρ₂+2; k₄=ρ₄−Σ_(3,1)ρ₃−Σ_(2,2)ρ₂ρ₂+2Σ_(2,1,1)ρ₂−6. Twenty exact correlated finite occupancy laws verify the complete contraction, all fifteen collision patterns and the pairing term. The connected formulation cancels disconnected contributions algebraically before numerical quadrature. It incorporates sample centering exactly and therefore supplies an arithmetic target without borrowing the projection-DPP centering theorem.

For bounded factorial-cumulant density errors ε₂,ε₃,ε₄ at fixed unit intensity, the package gives an explicit uniform error-transfer inequality. Put m_L=L/D_L, m_M=M/D_M and d=max(D_L,D_M). Define C₂=LM+4m_LM²+4m_ML²+4min(m_LM²,m_ML²), C₃=2LM²+2ML²+16L²M²/d and C₄=4L²M². If the reference |k₂|≤B₂, then |ΔCov|≤C₂ε₂+C₃ε₃+C₄ε₄+16L²M²(B₂+1/d)ε₂+8L²M²ε₂².

The proof uses ||g_L||₁≤2L², ||g_L||∞≤m_L, ∫|g_L(x,y)|dy≤2L a_L(x), and the fifteen collision contractions. The last three terms bound the change of the two pairing products, including their delta measures. All twenty-one length-pair constants are retained exactly as rational numbers. At L=M=4, C₄=1024: a crude uniform fourth-cumulant error alone would need to be below about 9.77×10⁻¹² to put its covariance contribution below 10⁻⁸. Weighted error bounds can be sharper, but a few stable point samples cannot meet this acceptance target by themselves.''')
    sec('7. Source constants, repaired status labels and surviving gaps', '''The source overlay corrects the printed circle majorant from 2742.431468 to the safe outward value 2742.431538 in both the main Part X paper and its source. The machine upper value 2742.431537628173828125 is unchanged. The existing degree-eighteen tail and origin theorem remain valid. A Part VIII README note redirects its historical open status to Part X and the later campaign; its old three-sample evidence is not relabeled as a uniform proof. Original and replacement hashes and reviewable patches are included.

The other constants keep distinct hypotheses: the corrected inner-disk derivative estimate is 4M/R; the ideal-clamp optimum is 1/√2; the finite (t,2t) protocol with bounded raw-reading errors has tG_jj=1; the Gaussian optimum is about .9721. Their differing values are not interchangeable contradictions. The finite-height origin bound c≤3 cannot be extended indefinitely, since its leading coefficient tends to π²/3>3.

The historical common offset includes the proved missing-origin contribution 2La−a², equal to .00319999 at L=16 and a=10⁻⁴. The complete repaired integral and empirical anchor calculation already control the deterministic numerical effects. A smaller observed residual remains. The present sine covariance neither proves that residual significant under the arithmetic model nor proves it absent.

The decisive next mathematical/computational target is now a collision-regularized, weighted connected arithmetic covariance contraction, first at L=4 in the same fixed observation geometry. It needs a certified infinite-prime derivative tail, uniform cancellation near colliding arguments, and a justified treatment of the exact functional-equation and height-window factors. Any ratios-conjecture/RH or short-window assumptions must remain visible. The finite-input compiler and the new error-transfer constants make this target concrete, but do not finish it.

New frozen external-coordinate windows remain an empirical confirmation route. They should preserve the unfolding, observation contract, covariance method and joint rule before acquiring data. The existing exploratory windows are not newly preregistered tests. More block partitions, more pointwise digits or an exact Hermitian BK kernel cannot substitute for the missing arithmetic correlation control.''')
    sec('Reproduction and trust boundary', '''For the full matrix run fixed_line/build_counts_native.py optionally, then fixed_line/full_matrix.py and fixed_line/validate_full_matrix.py. The native route was exercised on the present Linux runtime; the Python ball fallback is supplied, not claimed as a tested Windows execution. All 21 retained case contracts, panel balls and harmonic coefficients are included.

For arithmetic run arithmetic_fourpoint/dpp_obstruction.py, ratios_jet.py, validate_ratios.py, collision_diagnostics.py and connected_bridge.py. The ratios script resumes completed finite-input rows; remove or preserve its output before demanding a genuinely fresh entire run. For roots run root_hardening/check_hardening.py, replay_retained_contours.py and full_multiprocess_replay.py; replay V_e10 with r2kverify_hardened.py certs/V_e10.json --full. Default full replay recomputes its declared target; --resume explicitly reuses recovery data.

All mathematical status claims depend on the displayed arguments and correct implementation/runtime arithmetic. Exact controls, interval certificates, historical receipts, finite-prime approximation evidence and conjectural zeta transfer are kept separate. These are not end-to-end Lean proofs. The final archive retains all earlier package members and verifies every manifest hash and ZIP CRC.''')
    body.append('<details><summary>Current validation summary</summary><pre>'+esc(json.dumps({'fresh_root_contours':49,'fresh_root_segments':total_segments,'high_tiles':high,'fixed_line_entries':21,
           'matrix_positive_definite':True,'max_covariance_radius':V['maximum_anchor_entry_radius'],
           'DPP_obstruction_witnesses':len(O['witnesses']),'finite_input_arithmetic_values':len(J['arithmetic_results']),
           'RMT_controls':len(JV['RMT_determinant_controls']),'local_prime_mixed_controls':48,
           'BK_pair_controls':len(JV['finite_prime_BK_pair_controls']),'connected_exact_law_checks':20,
           'zeta_sampling_law_certified':False},indent=2))+'</pre></details>')
    out='<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>R2K comprehensive mathematical advance</title><style>body{font-family:system-ui,sans-serif;max-width:1180px;margin:36px auto;padding:0 24px;color:#17364c;background:#fbfcff;line-height:1.65}h1,h2,h3{color:#17486f;line-height:1.25}section{margin:2.2em 0}table{border-collapse:collapse;width:100%;font-size:12px}td,th{border:1px solid #c6d4e2;padding:8px;text-align:left;overflow-wrap:anywhere}pre{white-space:pre-wrap;overflow-wrap:anywhere;padding:16px;background:#eaf1f8;font-size:12px}summary{font-weight:600;cursor:pointer;padding:12px;background:#e3edf7}</style><main>'+''.join(body)+'</main></html>'
    assert out.count('<table>')==out.count('</table>') and out.count('<details>')==out.count('</details>')
    (HERE/'R2K_Comprehensive_Advance_Report.html').write_text(out)
    receipt={'status':'PASS','fresh_original_catalogue_full_replays':49,'fresh_segments':total_segments,
             'separate_full_fixture_replays':2,'all_fresh_segments_including_fixtures':6434,
             'high_counts':{r['name']:r['count'] for r in high},'distinct_fixed_line_covariance_entries':21,
             'certified_Hermitian_DPP_obstruction_witnesses':5,'finite_prime_arithmetic_values':len(J['arithmetic_results']),
             'open':['infinite-prime 3/4-point derivative tails','uniform coincidence quadrature','functional-equation correction in local model','ratios/RH and short-window transfer'],
             'zeta_sampling_law_certified':False}
    (HERE/'comprehensive_validation.json').write_text(json.dumps(receipt,indent=2));print(json.dumps(receipt,indent=2))

if __name__=='__main__':main()
