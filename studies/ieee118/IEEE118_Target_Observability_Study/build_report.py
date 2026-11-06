#!/usr/bin/env python3
import json,sys,base64,html,hashlib,csv
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parent;O=ROOT/'results'
def load(n):return json.loads((O/n).read_text())
def fmt(v):
    if v is None:return 'unresolved'
    if isinstance(v,str):return html.escape(v)
    return f'{v:.4g}'
def table(headers,rows):
    return '<table><thead><tr>'+''.join('<th>'+html.escape(x)+'</th>' for x in headers)+'</tr></thead><tbody>'+''.join('<tr>'+''.join('<td>'+fmt(v)+'</td>' for v in row)+'</tr>' for row in rows)+'</tbody></table>'
def fig(name):return '<figure><img alt="'+name+'" src="data:image/png;base64,'+base64.b64encode((ROOT/name).read_bytes()).decode()+'"></figure>'
def main():
    receipt=load('receipt.json');linear=load('linear.json');non=load('nonlinear.json');branch=load('branch.json');drop=load('dropout.json');shift=load('shift.json');checks=load('checks.json');follow=load('followup.json')
    methods=['qr_norm','d_opt','a_opt','target_a','target_correlated'];labels={'qr_norm':'QR + row norm','d_opt':'D-optimal greedy','a_opt':'Parameter A-optimal','target_a':'Target A-optimal','target_correlated':'Target + correlation','historical':'Historical adaptive','full':'Full 6,903'}
    colors=['#8b91a3','#8c6ad9','#ca8a2c','#167f8c','#d04c67']
    def row(m,b,n):return next(x for x in linear if x['method']==m and x['budget']==b and x['noise']==n)
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'figure.facecolor':'white','axes.titleweight':'bold'})
    f,axs=plt.subplots(1,2,figsize=(12,4.2))
    for ax,n in zip(axs,['iid','correlated']):
      for method,col in zip(methods,colors):
        rr=sorted([x for x in linear if x['method']==method and x['noise']==n],key=lambda x:x['budget'])
        ax.plot([x['budget'] for x in rr],[x['target_rms_per_sigma'] for x in rr],'-o',ms=4,label=labels[method],color=col)
      full=row('full',6903,n)['target_rms_per_sigma'];ax.axhline(full,ls='--',color='#333',label='Full-panel reference (6,903)')
      ax.set(xlabel='Scalar readings',ylabel='Relative target RMS / noise fraction',title='Independent noise' if n=='iid' else 'Common-reference + group noise',yscale='log');ax.grid(alpha=.2)
    axs[0].legend(fontsize=8);f.tight_layout();f.savefig(ROOT/'budget_frontier.png',dpi=170);plt.close(f)
    f,axs=plt.subplots(1,2,figsize=(12,4.2))
    for ax,metric,title in zip(axs,['target_rms_per_sigma','parameter_rms_per_sigma'],['Declared target precision','Full-parameter precision']):
      vals=[row(m,379,'iid')[metric] for m in methods]
      ax.bar(np.arange(len(methods)),vals,color=colors);ax.set_xticks(np.arange(len(methods)),['QR','D','A','Target','Target+R']);ax.set(ylabel='RMS / noise fraction',title=title);ax.grid(axis='y',alpha=.2)
    f.tight_layout();f.savefig(ROOT/'tradeoff.png',dpi=170);plt.close(f)
    f,axs=plt.subplots(1,2,figsize=(12,4.2))
    noise=[x['sigma'] for x in branch['rows']]
    for key,label in [('max_gap_error','Largest disagreement'),('random_expected_error','Random unused reading'),('qr_extra_error','QR extra reading')]:
      axs[0].semilogx(noise,[x[key] for x in branch['rows']],'-o',label=label)
    axs[0].set(xlabel='Noise standard deviation / baseline',ylabel='Equal-prior decision error',title='Two fixed known candidate branches',ylim=(-.02,.53));axs[0].legend(fontsize=8);axs[0].grid(alpha=.2)
    for method,col in zip(['qr_norm','target_a','target_correlated'],[colors[0],colors[3],colors[4]]):
      med=[]
      for b in [180,250,379]:
        zz=[x['target_relative_rms'] for x in non if x['method']==method and x['budget']==b and x['sigma']==1e-5];med.append(np.median(zz)*100)
        axs[1].scatter([b]*len(zz),np.array(zz)*100,s=10,alpha=.3,color=col)
      axs[1].plot([180,250,379],med,'-o',color=col,label=labels[method])
    axs[1].set(xlabel='Scalar readings',ylabel='Target relative RMS (%)',title='Actual nonlinear fits; noise fraction 1e-5',yscale='log');axs[1].legend(fontsize=8);axs[1].grid(alpha=.2)
    f.tight_layout();f.savefig(ROOT/'branch_and_nonlinear.png',dpi=170);plt.close(f)
    baseline=row('qr_norm',379,'iid')['target_rms_per_sigma'];tar=row('target_a',379,'iid')['target_rms_per_sigma'];aa=row('a_opt',379,'iid')['target_rms_per_sigma']
    red=100*(1-tar/baseline);ared=100*(1-tar/aa)
    # QA: exact declared rows and numerical Monte Carlo/analytic consistency.
    ratios=[x['target_mc_rms_per_sigma']/x['target_rms_per_sigma'] for x in linear if x['target_rms_per_sigma'] is not None]
    assert min(ratios)>.8 and max(ratios)<1.2
    assert receipt['nonlinear_fits']==135 and len(follow['rows'])==33
    qa={'status':'PASS','linear_monte_carlo_to_analytic_ratio_range':[min(ratios),max(ratios)],'nonlinear_fit_count':len(non),'followup_fit_count':len(follow['rows']),'target_improvement_vs_qr_percent':red,'target_improvement_vs_parameter_a_percent':ared}
    (O/'report_qa.json').write_text(json.dumps(qa,indent=2)+'\n')
    budgetrows=[[labels[x['method']],x['budget'],x['target_rms_per_sigma'],x['parameter_rms_per_sigma']] for x in linear if x['noise']=='iid' and x['budget'] in (180,379,6903) and x['method']!='random']
    nr=[]
    for method in ['qr_norm','target_a','target_correlated']:
     for b in [180,250,379]:
      for sigma in [0,1e-7,1e-5]:
       zz=[x for x in non if x['method']==method and x['budget']==b and x['sigma']==sigma]
       nr.append([labels[method],b,sigma,float(np.median([x['target_relative_rms'] for x in zz]))*100,float(np.median([x['weight_relative_l2'] for x in zz]))*100,f"{sum(x['solver_success'] for x in zz)}/{len(zz)}"])
    fr=[]
    for method in ['qr_norm','target_a','target_correlated']:
     for sigma in [0,1e-7,1e-5]:
      zz=[x for x in follow['rows'] if x['method']==method and x['sigma']==sigma]
      fr.append([labels[method],zz[0]['budget'],sigma,float(np.median([x['target_relative_rms'] for x in zz]))*100,float(np.median([x['weight_relative_l2'] for x in zz]))*100])
    dr=[]
    for method in ['historical','target_a']:
      dd=[x for x in drop if x['method']==method];dr.append([labels[method],len(dd),sum(x['target_recoverable_at_tolerance'] for x in dd),sum(not x['target_recoverable_at_tolerance'] for x in dd)])
    sr=[]
    for method in ['qr_norm','target_a','target_correlated']:
     for sd in [.03,.1]:
      zz=[x['target_rms_per_sigma'] for x in shift if x['method']==method and x['budget']==379 and x['logweight_sd']==sd]
      sr.append([labels[method],sd,min(zz),float(np.median(zz)),max(zz)])
    corrrows=[]
    for method in methods:
      x=row(method,379,'correlated');corrrows.append([labels[method],x['target_rms_per_sigma'],x['target_mc_rms_per_sigma'],x['naive_iid_decoder_mc_rms_per_sigma']])
    secs=sum(load('panels.json')['selection_seconds'].values())
    body=f'''<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1"><title>IEEE-118 · target-relative measurement study</title><style>
    :root{{--ink:#172438;--muted:#526278;--accent:#167f8c}}*{{box-sizing:border-box}}body{{margin:0;background:#f4f6f8;color:var(--ink);font:16px/1.65 system-ui,sans-serif}}main{{max-width:1120px;margin:auto;padding:45px 30px 80px;background:white}}h1{{font-size:clamp(30px,5vw,48px);line-height:1.12;letter-spacing:-1.5px}}h2{{margin-top:48px;font-size:25px}}p{{max-width:95ch}}.eyebrow{{letter-spacing:2px;text-transform:uppercase;font-size:12px;color:var(--accent);font-weight:750}}.lead{{font-size:21px}}.callout{{border-left:4px solid var(--accent);padding:16px 22px;background:#eef7f7}}.warning{{border-color:#bc6b25;background:#fff6eb}}table{{border-collapse:collapse;width:100%;font-size:13px;margin:20px 0;display:block;overflow:auto}}th,td{{padding:9px 12px;text-align:left;border-bottom:1px solid #dce3eb;white-space:nowrap}}th{{background:#edf2f6}}figure{{margin:24px -8px}}img{{max-width:100%;height:auto}}a{{color:#096d80}}.small{{font-size:13px;color:var(--muted)}}code{{background:#eef1f4;padding:2px 4px}}footer{{margin-top:40px;border-top:1px solid #dce3eb;padding-top:20px;color:var(--muted);font-size:13px}}</style><main>
    <div class="eyebrow">Research experiment · 4 November 2026 · IEEE-118</div>
    <h1>Fewer readings, chosen for the target—and checked for ambiguity</h1>
    <p class="lead">At 379 scalar readings, target-focused selection reduced local predicted target RMS by <strong>{red:.1f}%</strong> versus QR plus row-norm selection and <strong>{ared:.1f}%</strong> versus full-parameter A-optimal selection. The benefit trades away precision in nuisance parameters.</p>
    <div class="callout warning"><strong>The decisive distinction:</strong> full local rank and solver success do not guarantee the correct nonlinear branch. Nonadaptive 180-row designs can return the saved wrong branch even with zero noise. The inherited adaptive 180th reading replays correctly; a separate 181-row follow-up combines precision selection with that known separator.</div>
    <p class="small">Completed: 1,000-trial local decoder checks for each feasible panel/noise case; 135 primary nonlinear fits; 33 exploratory follow-up fits; 360 single-deletion checks; 90 parameter-shift comparisons. These are simulations, not field measurements or a global injectivity theorem.</p>
    <h2>1 · What is being recovered</h2>
    <p>The unknowns are 179 aggregate internal-edge and reference-ground weights in a 117-dimensional grounded DC operator. Measurements are scalar one-port compliances after calibrated grounding interventions. The full candidate protocol has 6,903 readings. All designs use the supplied topology and a uniform nominal weighting scaled by one baseline reading; the actual MATPOWER weights enter evaluation only.</p>
    <p>The twelve targets are effective resistances at deterministic bus pairs listed in <code>results/panels.json</code>. Selection optimizes their relative first-order uncertainty while retaining a full-parameter rank-completing seed. This is not a proof of the minimum number of readings required for these twelve targets.</p>
    <h2>2 · Independent replay and implementation checks</h2>
    <p>The original 179/180-reading result passed replay: both saved exact-fit branches, the branch separator, final saved solutions, and both finite-field rank determinants. The unchanged MATPOWER case file matches its frozen SHA-256. Independently parsed tap-adjusted coefficients agree with the frozen fixture.</p>
    {table(['Check','Result'],[['Response difference / baseline',checks['vectorized_response_relative_max']],['Log-Jacobian relative difference',checks['vectorized_logjac_relative_error']],['Worst centered-difference discrepancy',max(checks['finite_difference_relative_errors'])],['MC / analytic target RMS minimum',min(ratios)],['MC / analytic target RMS maximum',max(ratios)]])}
    <p>Numerical precision uses QR-based square-root least squares. Inverting normal equations would square the condition number. A full-parameter rank failure triggers a separate target-kernel test and target-only SVD decoder. Only panels failing the target test are left unresolved. These are numerical diagnoses, not exact-rank impossibility theorems.</p>
    <h2>3 · Local measurement-budget frontier</h2>
    {fig('budget_frontier.png')}
    <p>The vertical quantity multiplies the noise fraction sigma. For example, a coefficient 205.4 at sigma=1e-5 predicts relative target RMS 0.2054% in the local model. The full 6,903-row iid coefficient is {row('full',6903,'iid')['target_rms_per_sigma']:.2f}. Target A-optimal selection at 379 rows uses {6903/379:.2f} times fewer readings, with {tar/row('full',6903,'iid')['target_rms_per_sigma']:.2f} times the full-panel local target error.</p>
    {table(['Method','Readings','Target RMS / sigma','Parameter RMS / sigma'],budgetrows)}
    <p><strong>Recovering the target can require less information than recovering every weight.</strong> The frozen random panels at 200, 250 and 379 readings have numerical parameter ranks 174, 174 and 175, yet all twelve targets annihilate the corresponding numerical nullspaces. Their target-only decoders remain valid at the declared tolerances, although noise amplification is worse. The 179/180-row random panels fail the target test. This qualification was caught during the final audit; targets and panels were unchanged.</p>
    {table(['Random readings','Parameter rank','Target kernel gap','IID target RMS / sigma'],[[x['budget'],x['rank_rtol1e12'],x.get('target_kernel_gap'),x['target_rms_per_sigma']] for x in linear if x['method']=='random' and x['noise']=='iid'])}
    {fig('tradeoff.png')}
    <p>Full-parameter A-optimality is the strongest parameter-error comparator here. Target selection improves the declared target but worsens unimportant coordinates. This is the intended tradeoff, and rules out describing the target design as uniformly better.</p>
    <h2>4 · Correlated noise and decoder choice</h2>
    <p>Correlated noise has 50% independent variance, 25% common-reference variance, and 25% first-intervention-node group variance. Every reading retains unit marginal variance. The grouping and fractions are synthetic, not instrument-calibrated. A covariance-aware selector is not a minimax-robust selector.</p>
    {table(['379-row design','GLS analytic target RMS / sigma','GLS MC target RMS / sigma','Naive iid decoder MC RMS / sigma'],corrrows)}
    <p>Correlation can reduce error in differences when shared offsets cancel; it need not make every target worse. These results establish behavior under this specified covariance only. All uncertainty calculations condition on the baseline pilot used to set nominal scale; its noise is not propagated through redesign.</p>
    <h2>5 · Local rank versus the correct nonlinear branch</h2>
    {fig('branch_and_nonlinear.png')}
    <p>The two retained branches differ by {100*branch['branch_relative_weight_distance']:.4f}% in relative weight L2. Their strongest unused response separation is {branch['normalized_max_gap']:.6g} times the baseline. Both absolute-gap and historical relative-gap selection choose index {branch['max_absolute_gap_index']}.</p>
    {table(['Noise fraction','Best-gap decision error','Random unused decision error','QR extra decision error'],[[x['sigma'],x['max_gap_error'],x['random_expected_error'],x['qr_extra_error']] for x in branch['rows']])}
    <p>This is an exact equal-prior Gaussian decision formula for two fixed known candidates using one new scalar reading. Tiny probabilities can underflow to displayed zero. It excludes uncertainty in fitting the candidates and any unknown additional branches. It is not an end-to-end adaptive error bound.</p>
    <h2>6 · Actual nonlinear recovery from the nominal start</h2>
    <p>Each fit starts from baseline-scaled uniform weights, uses the declared correlated-noise model and GLS loss, and receives at most 100 evaluations. Log-weights are bounded to nominal +/-8. Target and weight errors are measured against hidden truth. Solver success means numerical termination, not correct recovery. Zero-noise rows contain five identical deterministic replays, not five independent starts.</p>
    <div class="callout warning"><strong>Target accuracy does not mean weight recovery:</strong> at 379 readings and noise fraction 1e-5, target A-optimal selection gives median target error 0.150%, versus 0.306% for QR+norm, but its median full-weight error is 29.5%, versus 2.14%. Four of its five solvers terminate successfully. Weak nuisance directions and nonlinear fitting remain material limitations.</div>
    {table(['Design','Readings','Noise fraction','Median target RMS %','Median weight L2 %','Solver successes'],nr)}
    <h2>7 · Concrete repair: add the known separator</h2>
    <p>After inspecting the primary results, a separately labelled follow-up adds the inherited separator to each 180-reading design. These panels contain 181 rows. This uses an already discovered ambiguity; an operational system would need to discover competing candidates from its own observations. Reproduce separately with <code>python followup.py</code>.</p>
    {table(['Design + separator','Readings','Noise fraction','Median target RMS %','Median weight L2 %'],fr)}
    <h2>8 · Single-reading loss and parameter shifts</h2>
    {table(['180-row panel','Deletions checked','Target kernel test passes','Target kernel test fails'],dr)}
    <p>Deletion counts use the same singular-value tolerance and a target-nullspace gap threshold 1e-8. Passing nominal/full rank does not imply all channels are dispensable. No exact adversarial-corruption theorem is claimed by these numerical checks.</p>
    {table(['379-row design','Logweight perturbation SD','Minimum target RMS / sigma','Median','Maximum'],sr)}
    <p>These ten parameter-shift draws preserve support. They test local sensitivity away from frozen weights, not nonlinear convergence under topology changes. Runtime: {receipt['seconds']:.1f} s for the primary run; total greedy selection work {secs:.1f} s on this host. Repeated zero-noise fits are included in runtime. Timings are implementation- and machine-specific.</p>
    <h2>9 · What the evidence supports</h2>
    <ul><li>Target covariance is a useful design objective when the scientific target is fixed in advance.</li><li>Rank completion, noise conditioning and nonlinear ambiguity require separate checks.</li><li>The 181-row follow-up is a concrete combination of noise-oriented selection with a known branch separator; it remains exploratory.</li><li>Full raw indices, source code and results are retained. No globally minimum target panel, universal robust-noise guarantee or full AC state-estimation result is claimed.</li></ul>
    <p>Next protected experiment: freeze the combined policy, reacquire a noisy baseline within every replicate, discover branches without saved candidates, and test untouched weight/topology cases with application-chosen targets. Additional literature methods—quadratic target design, worst-case covariance design, MISDP, interval estimation, GNNs and PowerBench transfer—have not been reproduced in this package.</p>
    <h2>Primary sources and exact derivations</h2>
    <p>The accompanying <code>SOURCES.md</code> distinguishes literature contributions from algorithms implemented here. <code>THEORY.md</code> proves the kernel criterion and states noise, branch-decision and corruption guarantees with their assumptions.</p>
    <p><a href="https://arxiv.org/abs/2307.08923">Functional observability and sensor placement</a> · <a href="https://doi.org/10.1137/24M1667543">Robust A-optimal design</a> · <a href="https://arxiv.org/abs/2411.07532">Quadratic goal-oriented design</a> · <a href="https://arxiv.org/abs/2502.03338">DAE/Kalman PMU placement</a> · <a href="https://arxiv.org/abs/2510.16693">Bounded-uncertainty grid estimation</a> · <a href="https://zenodo.org/records/15400564">PowerBench</a></p>
    <footer>Protocol SHA-256: {receipt['protocol_sha256']}<br>Exploratory computational study · external MATPOWER fixture · no hardware experiment</footer></main></html>'''
    (ROOT/'report.html').write_text(body)
    # Plain-text summary stays usable without a renderer.
    (ROOT/'RESULTS.md').write_text(f'# Completed IEEE-118 study\n\nAt 379 readings, target A-optimal selection reduces local iid target RMS by {red:.2f}% versus QR+norm and {ared:.2f}% versus full-parameter A-optimal selection.\n\nThe full protocol has 6,903 readings: the 379-row target design is {6903/379:.2f}x smaller but has {tar/row("full",6903,"iid")["target_rms_per_sigma"]:.2f}x the local target RMS. Parameter precision is worse than with parameter A-optimal selection.\n\nPrimary nonlinear fits: 135; post-inspection follow-up fits: 33. See report.html and results/nonlinear.json for truth errors independent of solver success.\n\nOriginal 179/180 branch and modular witnesses replayed. Nonadaptive 180-row panels can converge to the wrong branch. The known separator is evaluated separately in 181-row follow-ups. No global injectivity claim.\n')
    for name,data in [('linear',linear),('nonlinear',non),('followup',follow['rows'])]:
      with (O/(name+'.csv')).open('w',newline='') as out:
       writer=csv.DictWriter(out,fieldnames=list(dict.fromkeys(k for x in data for k in x)));writer.writeheader();writer.writerows(data)
    print(json.dumps(qa,indent=2))
if __name__=='__main__':main()
