#!/usr/bin/env python3
"""Render a bounded research report from the final, non-quarantined receipts."""
import html
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parent
def read(name):return json.loads((ROOT/name).read_text())
def table(headers,rows):
    return '<table><thead><tr>'+''.join('<th>'+html.escape(str(x))+'</th>' for x in headers)+'</tr></thead><tbody>'+''.join('<tr>'+''.join('<td>'+html.escape(str(x))+'</td>' for x in row)+'</tr>' for row in rows)+'</tbody></table>'

def build():
    adaptive=read('adaptive/results.json');policy=read('policy/results.json')
    external=read('external/summary.json');ident=read('identifiability/results.json')
    exact=read('identifiability/explicit_certificate.json');passive=read('passivity/results.json')
    direct=read('independent_certification/results.json')
    primary=next(x for x in adaptive['comparisons'] if x['scenario']=='ALL' and x['strategy']=='adaptive_disagreement' and x['comparator']=='log_uniform')
    base=next(x for x in adaptive['summary'] if x['scenario']=='ALL' and x['noise']=='POSITIVE' and x['strategy']=='log_uniform')
    advanced=next(x for x in adaptive['summary'] if x['scenario']=='ALL' and x['noise']=='POSITIVE' and x['strategy']=='adaptive_disagreement')
    gate=next(x for x in policy['comparisons'] if x['scenario']=='ALL' and x['strategy']=='gated')
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    scenarios=['ordinary','colored_noise','weak_link','near_collision','model_mismatch']
    rows=[next(x for x in adaptive['comparisons'] if x['scenario']==s and x['strategy']=='adaptive_disagreement' and x['comparator']=='log_uniform') for s in scenarios]
    fig,axes=plt.subplots(1,2,figsize=(13.3,4.8),constrained_layout=True)
    vals=np.array([x['geometric_error_ratio'] for x in rows]);low=np.array([x['ratio_ci95'][0] for x in rows]);high=np.array([x['ratio_ci95'][1] for x in rows])
    axes[0].errorbar(vals,np.arange(5),xerr=np.array([vals-low,high-vals]),fmt='o',color='#14556c',capsize=4)
    axes[0].axvline(1,color='#a96e3c',ls='--',lw=1)
    axes[0].set_yticks(np.arange(5),['Ordinary','Colored noise','Weak link','Close modes','Model mismatch']);axes[0].invert_yaxis()
    axes[0].set(xlabel='Adaptive / fixed-log geometric error ratio',title='Fresh adaptive confirmation · 64 observations')
    axes[0].grid(axis='x',alpha=.18)
    tasks=[('interleaved',12),('interleaved',20),('missing_band',12),('missing_band',20)]
    methods=['loewner','vector_fitting','cubic_log_frequency','positive_real_relaxation']
    for j,(m,label,color) in enumerate(zip(methods,['Loewner','Vector fitting','Cubic','Passive RC'],['#64899c','#153f52','#b59b78','#bf703b'])):
        y=[100*next(x['median_capped_error'] for x in external['per_task_method'] if x['scenario']==s and x['budget']==b and x['method']==m) for s,b in tasks]
        axes[1].bar(np.arange(4)+(j-1.5)*.2,y,width=.18,label=label,color=color)
    axes[1].set_xticks(np.arange(4),['Interleave\n12','Interleave\n20','Missing band\n12','Missing band\n20'])
    axes[1].set(ylabel='Median held-out relative RMS error (%)',title='Protected confirmation · 48 measured cells')
    axes[1].legend(frameon=False,fontsize=8);axes[1].grid(axis='y',alpha=.15)
    fig.savefig(ROOT/'frontier_comparison.png',dpi=180,facecolor='white');plt.close(fig)

    fig,axes=plt.subplots(1,2,figsize=(11,3.7),constrained_layout=True)
    coords=np.array([[0,0],[1,0],[2,0]],float)
    for ax,title,edges,grounds in [(axes[0],'Grounded path',exact['original_inter_mass_springs'],exact['original_grounding_springs']),(axes[1],'Grounded triangle',exact['transformed_inter_mass_springs'],exact['transformed_grounding_springs'])]:
        xy=coords.copy()
        if 'triangle' in title:xy[1,1]=.7
        for u,v,k in edges:
            p,q=xy[u-1],xy[v-1];ax.plot([p[0],q[0]],[p[1],q[1]],color='#466b7c',lw=3)
            midpoint=(p+q)/2;ax.text(midpoint[0],midpoint[1]+.08,k,ha='center',fontsize=10,bbox={'facecolor':'white','edgecolor':'none','alpha':.9,'pad':2})
        ax.scatter(xy[:,0],xy[:,1],s=600,c=['#bf703b','#14556c','#14556c'],zorder=3)
        for j,p in enumerate(xy):
            ax.text(*p,str(j+1),ha='center',va='center',color='white',weight='bold',zorder=4)
            ax.text(p[0],-.45,'Ground spring '+grounds[j],ha='center',fontsize=8)
        ax.text(0,-.2,'Drive and observe',ha='center',fontsize=9)
        ax.set_title(title);ax.set(xlim=(-.4,2.4),ylim=(-.65,1.05));ax.axis('off')
    fig.savefig(ROOT/'same_response_different_networks.png',dpi=180,facecolor='white');plt.close(fig)

    metricrows=[]
    for label,key in [('Geometric mean error','geometric_nrmse'),('Arithmetic mean error','mean_capped_nrmse'),('90th percentile error','p90_nrmse')]:metricrows.append([label,f'{100*base[key]:.3f}%',f'{100*advanced[key]:.3f}%'])
    metricrows.append(['Median computation',f'{1000*base["median_seconds"]:.2f} ms',f'{1000*advanced["median_seconds"]:.2f} ms'])
    measuredrows=[]
    for s,b in tasks:
        measuredrows.append([s.replace('_',' '),b,*[f'{100*next(x["median_capped_error"] for x in external["per_task_method"] if x["scenario"]==s and x["budget"]==b and x["method"]==m):.3f}%' for m in methods]])
    reduction=100*(1-primary['geometric_error_ratio']);ci=primary['ratio_ci95']
    text=f'''<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Nonagon Frontier</title><style>body{{margin:0;background:#f1eee6;color:#17313e;font:17px/1.65 system-ui,sans-serif}}main{{max-width:1080px;margin:auto;padding:48px 24px 90px}}h1{{font-size:clamp(2.3rem,5vw,4rem);line-height:1.06;margin:20px 0}}h2{{font-size:1.55rem;margin-top:42px}}.eyebrow{{font-size:12px;letter-spacing:.15em;text-transform:uppercase;color:#4d737c}}.lead{{font-size:1.22rem;max-width:950px}}.callout{{padding:18px 24px;background:#fffbef;border-left:5px solid #b8753e;margin:24px 0}}table{{display:block;overflow:auto;border-collapse:collapse;width:100%;background:white;font-size:14px}}th,td{{padding:10px 13px;border-bottom:1px solid #d8e0df;text-align:left}}th{{background:#173f51;color:white}}img{{width:100%;height:auto;margin:18px 0}}code{{background:#e0e5e3;padding:2px 5px;font-size:.9em}}.formula{{font:19px/1.7 Georgia,serif;padding:17px 24px;background:white;overflow:auto}}.small{{font-size:14px;color:#4b6470}}a{{color:#1c607c}}</style><main>
<div class="eyebrow">Rank-Two Nonagon · Frontier campaign · October 29, 2026</div>
<h1>Choose better observations.<br>Know what they determine.<br>Verify the result.</h1>
<p class="lead">The next level is a measured increase in acquisition capability, an exact boundary on physical inference, an independent interval route, and a reusable passive-model compiler. The experiment also preserves failed methods and a numerical-library failure that invalidated an early run.</p>
<h2>1. Adaptive sampling produces a real improvement in the frozen simulation</h2>
<p>The prespecified primary result is <strong>{reduction:.1f}% lower geometric mean held-out response error</strong> at the same 64 unique real-frequency observations. The adaptive/fixed-log error ratio is {primary['geometric_error_ratio']:.4f}, with paired system-cluster 95% interval [{ci[0]:.4f}, {ci[1]:.4f}]. The same final Loewner-plus-vector-fitting algorithm is used for every acquisition strategy, so this comparison isolates where observations were collected.</p>
{table(['Metric, positive-noise cases','Fixed log grid','Adaptive disagreement'],metricrows)}
<p>Fresh confirmation: 200 cases and 1,000 fits, including 150 positive-noise cases clustered by 50 systems. Scenarios cover ordinary chains, correlated noise, weak links, close modes and model mismatch. The pilot uses 33 of the 64 observations. The selector uses accrued observations, known order and the public band. The candidate-frequency noise table is common across strategies; clean test frequencies are disjoint and are not supplied to the acquisition designer.</p>
<p>This is a simulation result with known model order, not an executed hardware experiment. Computation rises by approximately {advanced['median_seconds']/base['median_seconds']:.1f}×. Subgroup intervals are descriptive, unadjusted secondary contrasts. A second fresh policy cohort does not establish a universally safer default.</p>
<img src="frontier_comparison.png" alt="Adaptive acquisition ratios in five synthetic regimes and all four measured battery tasks">
<div class="callout"><strong>A gate was tested and not promoted.</strong> A simple rule based only on pilot-estimated residue visibility was trained on development cases, frozen, and tested on 160 new cases. Its primary arithmetic-mean error ratio was {gate['mean_error_ratio']:.4f} [{gate['ratio_ci95'][0]:.4f}, {gate['ratio_ci95'][1]:.4f}], and the 90th percentile worsened. It remains experimental; the deployable adapter defaults to the fixed log design.</div>
<h2>2. Exact data cannot always identify the hidden machine</h2>
<p>A three-mass path and a differently connected three-mass triangle have identical driving-point responses at every frequency. Both have unit masses, strictly positive grounding and inter-mass springs, positive-definite stiffness, and every mode visible at the observed port.</p>
<img src="same_response_different_networks.png" alt="A grounded path and a grounded triangle with exact positive spring constants and identical one-port response">
<div class="formula">H(z) = (z² − 5z + 5) / (18 − 24z + 9z² − z³), &nbsp; H(z)=e₁ᵀ(K−zI)⁻¹e₁.</div>
<p>One additional static displacement measurement distinguishes the displayed pair: 1/9 versus 19/234, an exact difference of <strong>7/234</strong>. The construction extends to arbitrary order: even a complete transfer block on the first n−2 ports can remain identical after changing the hidden connections. No increase in numerical precision or frequency density at those unchanged ports resolves this exact ambiguity.</p>
<p>The constructive remedy uses known masses and full displacement coordinates at n distinct, nonresonant shifts:</p>
<div class="formula">K = (e₁1ᵀ + MXZ)X⁻¹.</div>
<p>Here the columns of X are measured displacement responses and Z contains the shifts. Under cyclicity, X is nonsingular and the unknown stiffness and graph are recovered. The single added sensor separates the constructed pairs; it is not a universal two-port reconstruction theorem. Validation: <strong>{ident['fixture_count']} exact physical network pairs and {ident['exact_checks']:,} checks</strong>, plus 56 checks in a separate determinant/adjugate verifier. The broader nonuniqueness principle is established realization theory; this is a checked constructive capability for the toolkit.</p>
<h2>3. Independent continuous full-strip verification</h2>
<p>A direct Arb contour program, importing no archived Taylor or DFT code, covers the entire rectangle <strong>10 ≤ Im(s) ≤ 100, −0.1 ≤ Re(s) ≤ 1.1</strong>. All {direct['campaign_leaves']} adjoining tiles certify: {direct['campaign_zero_count']} zeros, all simple and on the critical line, with no unresolved regions. Generation took {direct['campaign_seconds']:.2f} seconds; a fresh 45-tile replay independently reran this direct implementation and agreed throughout.</p>
<p>Each tile retains complete zero-free boundary partitions, moment enclosures and positive Hermite minors. Boundary interval midpoints and radii are stored as exact dyadics because a human-readable Arb string can deliberately coarsen a wide ball and lose its zero-exclusion information. Five analytic controls and seven ledger/trace corruptions exercise failure behavior; three serialization regressions check the repaired representation.</p>
<p class="small">This independently implemented numerical route still trusts Arb's function evaluations and adaptive integrator. Static hashes and ledger checks alone do not prove the interval evaluations. The range is a known finite RH fact, not a new verification-height record or global proof; heights below 10 are not covered by this campaign. A separate attempted direct Lehmer calculation hit a 60-second budget and is retained as a limitation.</p>
<h2>4. Measured batteries expose a model-class limitation</h2>
<p>The external experiment uses 54 physically measured NMC811 cells at 30% state of charge, recorded with a Keysight SL1007A. Six devices support development; the other 48 are protected confirmation. There are 4,608 confirmation candidate fits and 768 selected evaluations, with zero fitting exceptions. All data, published checksum checks and CC BY attribution are retained.</p>
{table(['Task','Training samples','Loewner','Vector fitting','Cubic','Passive RC'],measuredrows)}
<p>The positive RC-plus-ideal-RLC model loses to vector fitting by about 3.48× in geometric error. Post-confirmation diagnosis identifies a structural limitation: this class forces the real impedance to be nonincreasing with frequency, while all 48 stored confirmation spectra show a high-frequency upturn. Measurement uncertainty is not available here, and the upturn does not quantify the entire prediction loss. The model family is restrictive; this observation does not diagnose nonpassive batteries. No model was refitted after confirmation to change these scores.</p>
<h2>5. A reusable exact passive-model compiler</h2>
<div class="formula">Z(s)=R₀+Ls+C⁻¹/s+Σ aⱼ/(1+τⱼs)+Σ bⱼτⱼs/(1+τⱼs).</div>
<p>For nonnegative coefficients and positive time constants, every term is analytic with nonnegative real part on Re(s)&gt;0. The independent verifier interprets serialized decimal coefficients as exact rationals, checks the supported class and SI convention, and certifies this model-level property without a frequency grid. The compiler exports passive circuit realizations.</p>
<p><strong>216 fitted models verified</strong> (24 development plus 192 confirmation task models, not 216 devices). The expanded RC/RL verifier passes 150 fixtures, 2,400 exact checks and 29 adversarial rejections, including optimized-Python attacks. RL relaxations permit real-part upturns while retaining the positive-real guarantee. All 216 fitted models contain zero RL terms: the new representation is mathematically checked, and its fitting accuracy remains untested on new protected data. Circuit exports are algebraically checked, not simulator-tested.</p>
<p class="small">The restricted certificate is sufficient, not necessary for all positive-real functions. It certifies the serialized model, not battery behavior outside the observed range, a unique internal circuit, or hardware calibration. Ordinary Loewner fits frequently have right-half-plane poles despite low finite-band prediction error; prediction and structural validity are separate questions.</p>
<h2>Numerical failure and evidence discipline</h2>
<p>A zero-noise sanity check exposed a SciPy 1.17 automatic complex-tridiagonal solve failure on a well-conditioned 3×3 physical system. A separate literal-matrix reproducer shows a residual near 0.0017 versus machine precision with an explicit general solve; SciPy 1.18.1 does not reproduce that path defect. The affected initial adaptive experiment is invalidated and quarantined. The final experiment uses explicit general solves, residual checks and a separate NumPy solve comparison, followed by a fresh protocol freeze and untouched confirmation seed. The previous challenge was audited for this issue and left unchanged.</p>
<p>Protocol freezes are local, not external preregistration. The fitter is not supplied simulation truth; this is not external-custodian blinding. The battery study represents one batch, state of charge and instrument, not independent laboratories. Losses, timeouts, invalidated receipts, warnings and scope restrictions are included.</p>
<h2>Use the package</h2>
<p>Start with <code>README.md</code> and <code>SUMMARY.json</code>. Run <code>python replay.py --mode verify</code> for file integrity or <code>--mode quick</code> for exact tests in a fresh copy. Module READMEs document numerical replay. <code>previous/</code> includes the original challenge package unchanged. See <code>SOURCES.json</code> for primary data and established mathematical foundations.</p>
</main></html>'''
    (ROOT/'report.html').write_text(text)
    summary={'adaptive':{'confirmation_cases':200,'fits':1000,'same_unique_query_budget':64,'primary_geometric_error_ratio':primary['geometric_error_ratio'],'ratio_ci95':ci,'geometric_error_reduction_percent':reduction,'base':base,'adaptive':advanced},'policy':{'confirmation_cases':160,'fits':480,'primary_gated_mean_error_ratio':gate['mean_error_ratio'],'ratio_ci95':gate['ratio_ci95'],'promoted':False},'identifiability':{'exact_pairs':ident['fixture_count'],'exact_checks':ident['exact_checks'],'independent_checks':56,'shared_H_numerator_ascending':exact['shared_transfer_numerator_ascending'],'shared_H_denominator_ascending':exact['shared_transfer_denominator_ascending'],'added_sensor_difference':exact['exact_extra_sensor_response_difference']},'independent_contour':read('independent_certification/final_results.json'),'external':{'confirmation_cells':48,'candidate_fits':4608,'selected_evaluations':768,'passive_RC_superiority_established':False},'passivity':{'certified_fitted_models':216,'exact_fixtures':150,'exact_checks':2400,'adversarial_rejected':29,'RL_extension_empirically_tested':False},'invalidated_experiments_excluded':True,'previous_package_unchanged':True}
    (ROOT/'SUMMARY.json').write_text(json.dumps(summary,indent=2)+'\n')
    print('Report and figures built from final receipt paths.')

if __name__=='__main__':build()
