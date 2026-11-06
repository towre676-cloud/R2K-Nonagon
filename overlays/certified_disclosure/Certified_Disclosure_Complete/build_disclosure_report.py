"""Build the offline mathematical report from the verified extension ledgers."""
from __future__ import annotations
import base64
import html
import io
import json
import re
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
from matplotlib.mathtext import math_to_image
from matplotlib.font_manager import FontProperties

HERE=Path(__file__).resolve().parent
P=HERE/'possq-0.4.0'
R=HERE/'r2k_nonagon_v14.1'
OUT=(HERE.parent/'deliverables' if HERE.name=='work' else HERE)/'Certified_Disclosure.html'
load=lambda path:json.loads(path.read_text())
spectral=load(P/'results/resolvent_campaign.json')
cycle=load(P/'results/frontier_cycle.json')
star=load(P/'results/frontier_star.json')
hcycle=load(P/'results/frontier_horizons_cycle.json')
hstar=load(P/'results/frontier_horizons_star.json')
access=[load(p) for p in sorted((R/'part_XIV/rational_access').glob('rational_*.json'))]
captions=load(HERE/'figures/captions.json')
equation_number=0
figure_number=0


def eq(tex):
    global equation_number
    equation_number+=1
    tex=tex.replace(r'\mathsf T',r'\mathsf{T}').replace(r'\hbox',r'\mathrm').replace(r'\frac1{',r'\frac{1}{')
    for letter in ['Q','Z','E']:
        tex=tex.replace('\\mathbb '+letter,'\\mathbb{'+letter+'}')
    buf=io.BytesIO()
    math_to_image('$'+tex+'$',buf,prop=FontProperties(size=18),format='svg',dpi=150,color='#142435')
    svg=buf.getvalue().decode()
    svg=svg[svg.find('<svg'):]
    return '<div class="equation" role="img" aria-label="'+html.escape(tex,quote=True)+'">'+svg+f'<span>({equation_number})</span></div>'


def fig(name,caption):
    global figure_number
    path=HERE/'figures'/f'{name}.svg'
    if not path.exists():return ''
    figure_number+=1
    caption=re.sub(r'^Figure \d+\.\s*','',caption)
    caption=f'Figure {figure_number}. '+caption
    svg=path.read_text();svg=svg[svg.find('<svg'):]
    return '<figure class="scientific">'+svg+'<figcaption>'+caption+'</figcaption></figure>'


def table(headers,rows):
    return '<div class="tablewrap"><table><thead><tr>'+''.join('<th>'+h+'</th>' for h in headers)+'</tr></thead><tbody>'+''.join('<tr>'+''.join('<td>'+str(v)+'</td>' for v in row)+'</tr>' for row in rows)+'</tbody></table></div>'


def probability(d):
    return f"{d['lower_float']:.12e}"


best=spectral['campaigns'][-1]['adaptive'][-1]
certrows=[]
for c in spectral['campaigns']:
    a=c['adaptive'][-1]
    certrows.append([c['noise_relative_h0'],a['readings'],f"{a['lower_float']:.9f}",f"{a['upper_float']:.9f}",
                     f"{a['width_float']:.9f}"])
accessrows=[]
for a in access:
    accessrows.append([a['case'].replace('case','IEEE '),a['scheme'],len(a['ports']),f"{a['readings']:,}",
                       f"{a['rational_rank']}/{a['unknowns']}",a['rational_nullity'],
                       f"{a['first_order_identifiable_branches']}/{a['branches']}"])
copyrows=[]
for target,c in cycle['threshold_crossings'].items():
    s=star['threshold_crossings'][target]
    copyrows.append([target,c['minimal_copies'],probability(c['probability']),probability(s['probability'])])

body='''
<header>
 <div class="edition">NONAGON 14.1 / POSSQ 0.4 · 31 OCTOBER 2026</div>
 <h1>Certified<br><em>Disclosure</em></h1>
 <p class="deck">Exact access, spectral risk, and counting laws beyond 10<sup>−100</sup></p>
 <p class="lead">A hidden system need not be reconstructed completely to answer a consequential question about it. This work develops that principle into executable mathematics: rational access certificates establish which parameter directions the readings determine; positive spectral measures yield target bounds directly from uncertain responses; and sparse validated propagation resolves rare stochastic events while retaining every omitted-path error.</p>
 <div class="hero-rule"></div>
</header>
<nav aria-label="Contents"><a href="#spectral">I · Spectral disclosure</a><a href="#access">II · Exact access</a><a href="#roots">III · Intervention inversion</a><a href="#counting">IV · Counting frontier</a><a href="#verification">V · Reproduction</a></nav>
<main>
<section class="abstract"><h2>Abstract</h2><p>The extension joins two complementary forms of inference. Nonagon identifies network parameters from interventions; possq certifies functionals over all positive measures compatible with uncertain observations. A new rational-resolvent basis connects them without numerical differentiation. On an exact 32-node grounded network, twelve response samples certify a small low-decay participation band even though a nontrivial relabelling preserves every port response and prevents unique recovery of labelled branch weights. The grid-access analysis is repaired by complete rational nullspace and row-span witnesses, and a six-port IEEE 118 design reduces the reading count from 847 to 723. A sparse Arb-polynomial uniformization engine extends Cyclecraft from ten-copy risks near 10<sup>−16</sup> to seventy-copy risks below 10<sup>−100</sup>, with relative enclosure widths below 10<sup>−22</sup>. Chronological kernel composition separately resolves an uninterrupted 700-period experiment, retaining its state correlations.</p></section>
<section id="spectral"><div class="section-number">I</div><h2>Spectral decisions from the readings themselves</h2>
<p>A collocated response of a symmetric positive-definite first-order system is a Stieltjes transform of a positive spectral measure. For a calibrated positive storage matrix C, whitening replaces K by C<sup>−1/2</sup>KC<sup>−1/2</sup> and b by C<sup>−1/2</sup>b. The identity below is written in those coordinates.</p>
'''
body+=eq(r'h_b(s)=b^{\mathsf T}(K+sI)^{-1}b=\int_I\frac{d\mu_b(\lambda)}{\lambda+s},\qquad \mu_b(I)=\|b\|^2')
body+='''<p>The instrument supplies finitely many values ĥ<sub>j</sub> with declared deterministic errors ε<sub>j</sub>. The spectral support enclosure I and the mass normalization are supplied independently. The new method consumes those values in their native response basis. It neither estimates derivatives nor converts them into a potentially ill-conditioned sequence of high-order moments.</p>'''
body+=eq(r'q(\lambda)=c_0+\sum_{j=1}^m\frac{c_j}{\lambda+s_j},\qquad D(\lambda)=\prod_{j=1}^m(\lambda+s_j)>0\quad(\lambda\in I)')
body+=eq(r'U(c)=c_0M+\sum_{j=1}^m c_j\widehat h_j+|c_0|\varepsilon_M+\sum_{j=1}^m|c_j|\varepsilon_j')
body+='''<div class="theorem"><div class="theorem-label">Theorem 1 · Raw-response spectral certificate</div><p>If q ≥ 1 on a target band S and q ≥ 0 on its complement in I, every compatible positive spectral measure satisfies μ<sub>b</sub>(S) ≤ U(c). Reversing the sandwiches gives a lower bound with the uncertainty term subtracted. Rational input data and rational coefficients admit an exact verifier.</p><p><strong>Proof.</strong> Positivity gives μ<sub>b</sub>(S) ≤ ∫q dμ<sub>b</sub>. Expanding q and applying the independent response-error intervals gives U(c). Writing q = N/D, positivity of D reduces the continuous side conditions to polynomial signs of N and N − D. The verifier decides strict positivity by exact rational Sturm sequences, accepting the identically zero polynomial separately. A constant shift repairs an approximate optimizer proposal, and the repaired objective is recomputed exactly. ∎</p></div>
<p>The optimizer is therefore a proposal mechanism. The final certificate consists of raw-basis rational coefficients, the exact repaired objective, the target intervals and a replayable polynomial-sign decision. A separate rational atomic witness is accepted only when all response and mass intervals hold exactly. Its target mass demonstrates what the observations still permit.</p>
<div class="theorem"><div class="theorem-label">Theorem 2 · Refinement preserves certificates</div><p>Suppose new reading intervals refine the old ones and additional shifts are acquired. Embed the old coefficient vector in the larger raw basis, placing zero coefficients on new readings. Its spectral side conditions are unchanged; its upper objective cannot increase and its lower objective cannot decrease. Exact interval-containment checks and a fresh Sturm replay establish both claims. Retaining the strongest verified old or new certificate makes the implemented bounds monotone even when a numerical proposal deteriorates. ∎</p></div>
<h3>An exactly soluble network with permanently ambiguous branches</h3>
<p>The benchmark is a five-dimensional weighted cube: 32 vertices, 80 labelled branches, uniform ground conductance γ = 1/64, and axis conductances (1/32, 1/16, 1/8, 1/4, 1/2). The port is vertex zero. Its Walsh eigenvectors give a fully rational reference measure.</p>'''
body+=eq(r'K=\gamma I+\sum_{j=1}^{5}w_j(I-X_j),\qquad \lambda_A=\gamma+2\sum_{j\in A}w_j,\qquad \mu_0=\frac1{32}\sum_{A\subseteq\{1,\ldots,5\}}\delta_{\lambda_A}')
body+='''<p>Exchanging two coordinate axes fixes the port and conjugates K by a permutation. It changes the labelled branch weights while preserving h<sub>0</sub>(s) at every admissible s. The package verifies the matrix conjugacy exactly and checks four direct rational matrix inverses against the closed response law. Complete labelled reconstruction is therefore impossible from this port, even with infinitely many exact frequencies. The low-decay participation mass remains a well-defined observable target.</p>
<div class="instrument"><div class="instrument-top"><span>Exact cube / coordinate conjugacy</span><button id="swap">Exchange hidden axes</button><button id="motion">Pause motion</button></div><canvas id="network" width="1200" height="560" aria-label="Rotating projection of the actual 32-node, 80-edge weighted cube"></canvas><p>One collocated port is highlighted. Axis exchange changes the labelled network; the complete port response remains identical by exact permutation conjugacy.</p></div>
<p>The target is S = [1/64, 1/8], whose reference mass is exactly 1/16 = 6.25%. The table reports exact-verified response certificates at twelve selected shifts. Reading errors are absolute bounds ε = ηh(0). Acquisition selection is scored with the nominal benchmark model; each resulting certificate uses only its selected observations, declared errors and support enclosure.</p>
'''
body+=table(['Relative reading budget η','Readings','Certified lower','Certified upper','Interval width'],certrows)
body+=fig('spectral_disclosure','Figure 1. Direct-response spectral bounds as the acquisition budget increases. Solid and dashed traces distinguish nominally selected and wide-band designs; the horizontal reference is the exactly known mass 1/16. Certificate reuse preserves every earlier selected-design bound. These are deterministic uncertainty intervals, not confidence intervals.')
body+=fig('spectral_witnesses',captions['spectral_witnesses'])
body+='''<p>At η = 10<sup>−10</sup>, the twelve-reading certificate is <strong>6.1226001% ≤ μ(S) ≤ 8.2168456%</strong>. A distinct, exactly compatible atomic measure places <strong>8.2150349%</strong> in the same band. The worst-case upper target is therefore bracketed to a relative gap of <strong>0.0221%</strong>. This witness establishes an information limit: the declared reading boxes alone cannot support an upper certificate below its mass, even though the reference network has mass 6.25%.</p>'''
body+='''<div class="live-band"><label for="readings">Inspect a retained certificate at η = 10<sup>−8</sup></label><input id="readings" type="range" min="0" max="5" step="1" value="3"><output id="bandText"></output><canvas id="band" width="1200" height="240"></canvas><p>The slider selects discrete, independently verified certificates. It does not interpolate between them.</p></div>
<p>At the smallest budget in the table, the lower certificate lies close to the true mass and the upper certificate excludes most competing slow-mode mass. This is the intended engineering use: establish a target spectral statement while leaving irrelevant parameter ambiguity intact. The spectral measure records participation at the port. It is not a count of every hidden eigenmode; additional ports are needed for modes orthogonal to the original excitation.</p></section>
<section id="access"><div class="section-number">II</div><h2>Exact rational access, with witnesses on both sides</h2>
<p>The access question is local to a declared rational network model. Let J be the Jacobian of the acquired response vector with respect to the unknown branch and ground weights. Full column rank proves local injectivity. In deficient cases, a coordinate is determined to first order precisely when every rational kernel vector vanishes in that coordinate.</p>'''
body+=eq(r'\delta y=J\delta w,\qquad JN=0,\qquad \mathrm{rank}_{\mathbb Q}J+\dim_{\mathbb Q}\ker J=U')
body+=eq(r'\delta w_i=0\ \hbox{for every }\delta w\in\ker J\quad\Longleftrightarrow\quad e_i^{\mathsf T}=\alpha_i^{\mathsf T}J')
body+='''<p>A nonzero minor modulo a prime is a valid rational rank lower bound. A modular kernel’s zero-coordinate pattern, by itself, is not a rational coordinate certificate: the one-row matrix [1, p] becomes [1, 0] modulo p although its rational null vector (−p, 1) moves both coordinates. The new verifier closes that logical gap with actual rational witnesses.</p>
<div class="theorem"><div class="theorem-label">Theorem 3 · Complete first-order access certificate</div><p>Retain r genuine Jacobian rows with a nonzero r × r minor modulo p. Construct U − r independent rational kernel vectors and verify JN = 0 against every acquired reading row. Then the rational rank is exactly r and those vectors span the complete kernel. A coordinate with a nonzero kernel component is not determined to first order; an explicit row-span witness proves each claimed determined coordinate. Full-rank certificates need only the verified minor. ∎</p></div>
<p>The implementation clears inverse-column denominators before constructing Jacobian rows. It reconstructs these integer rows directly from the MATPOWER rational inputs, verifies their hashes, checks the rank minor and streams every kernel product. Certificate row indices must be distinct, in range and matched to the actual experiment. Assertion-disabled Python execution is rejected.</p>
'''
body+=table(['Model','Access scheme','Ports','Readings','Rational rank','Nullity','Branch coordinates'],accessrows)
body+=fig('access_design','Figure 2. Exact access and placement. All full-rank and deficient-rank conclusions shown here are backed by rational witnesses. The new six-port IEEE 118 design uses 723 readings rather than the inherited seven-port design’s 847: a 14.64% reading reduction and one fewer port. No claim of globally minimum port count is needed for this improvement.')
body+='''<p>The restricted-switching IEEE 118 result retains the earlier numerical percentage—74 of 173 branch coordinates—but now establishes it over the rationals. Its 52-dimensional kernel is checked against all 13,020 readings. The improved full-access design places its extra port inside the identified harmonic pocket. This is an explicit connection between a structural obstruction and a cheaper acquisition scheme. Deficient-rank statements concern first-order access; they do not assert a global nonlinear parameterization of the entire observation fiber.</p></section>
<section id="roots"><div class="section-number">III</div><h2>Intervention inversion without guessing a quadratic root</h2>
<p>For a pair of nodes, write u = G<sub>0i</sub>, v = G<sub>0j</sub>, A = G<sub>ii</sub>, B = G<sub>jj</sub> and x = G<sub>ij</sub>. A pair response at strength θ yields a quadratic constraint. Two strengths eliminate the squared term exactly.</p>'''
body+=eq(r'Q_\ell(x)=\delta_\ell x^2-2uvx+C_\ell,\qquad C_\ell=u^2(\theta_\ell^{-1}+B)+v^2(\theta_\ell^{-1}+A)-\delta_\ell(\theta_\ell^{-1}+A)(\theta_\ell^{-1}+B)')
body+=eq(r'x=\frac{\delta_2C_1-\delta_1C_2}{2uv(\delta_2-\delta_1)}')
body+='''<p>The formula depends only on the observations and calibrated or enclosed single-node quantities. Exact rational interval arithmetic propagates all supplied errors. The denominator must exclude zero; otherwise the measurement set is rejected as unresolved. Algebraic cancellation precedes interval evaluation, and interval-Newton contractions sharpen the enclosure against the original quadratics.</p>
<p>With more than two pair strengths, every admissible two-strength enclosure contains the true value. Intersecting these enclosures preserves soundness and can only narrow the result. The accompanying benchmark compares two and four pair strengths, both with independently calibrated inputs and with uncertainty propagated from raw single-node readings. The truth is used only to score coverage and radius after inversion. The raw cancellation path also requires known port-column signs; these are positive for the passive IEEE 118 M-matrix and otherwise form part of the calibration input.</p>
'''
body+=fig('two_strength_intervals',captions['two_strength_intervals'])
body+='''<h3>Sharper bounds from the same raw observations</h3><p>Additional readings are not the only route to a tighter enclosure. The raw single-node differences already determine the squared product (uv)<sup>2</sup>. Substituting this identity into the two-strength expression cancels the inferred diagonal variables before interval evaluation. This removes repeated appearances of the same uncertain quantities and reduces dependency inflation without adding any observations.</p>'''
body+=eq(r'(uv)^2=\frac{s_{i1}s_{i2}s_{j1}s_{j2}}{4\theta^2(s_{i2}-s_{i1})(s_{j2}-s_{j1})}')
body+='''<p>On all 207 IEEE 118 lean pairs at normalized strength one and deterministic raw-reading budget 10<sup>−6</sup>h(0), the median relative radius falls from <strong>2.3420% to 1.0826%</strong>. Every accepted interval contains the exact rational reference. At normalized strength ten, the number of pairs with radius below 1% rises from 59 to 79 under the same budget.</p>'''
body+=fig('raw_cancellation',captions['raw_cancellation'])
body+=fig('two_strength_strengths',captions['two_strength_strengths'])
body+='''</section>
<section id="counting"><div class="section-number">IV</div><h2>The counting frontier: preserve the law, account for every omission</h2>
<p>The Cyclecraft extension concerns the declared five-state, three-hold dyadic Markov protocols. Each signed counted transition shifts a polynomial coefficient index. A count window kills paths that leave it; a finite positive Poisson series omits additional paths. Every retained coefficient is therefore a submeasure contribution, and the missing total mass bounds all omitted events together.</p>'''
body+=eq(r'e^{Q\tau}=e^{-\Lambda\tau}\sum_{k=0}^{\infty}\frac{(\Lambda\tau)^k}{k!}P^k,\qquad P=I+Q/\Lambda,\qquad \Lambda\geq\max_j(-Q_{jj})')
body+='''<p>The engine applies this uniformization through sparse state transitions and Arb polynomial shifts, rather than exponentiating a large dense lifted matrix. All generator entries and hold times are exact rationals; Arb encloses exponentials, products and series weights. The periodic start is rebuilt at 512 bits with a Doeblin residual certificate. Finite count windows, truncated Poisson sums and initial-law error all enter a single explicit lost-mass budget.</p>'''
body+=eq(r'\nu(E)\leq\Pr(E)\leq\nu(E)+L,\qquad L\leq 1-\nu(\mathbb Z)+\delta_{\mathrm{init}}')
body+='''<p>For independent ten-period experiments, scalar convolution propagates the law of their sum. For one uninterrupted trajectory, the state at the period boundary must remain in the computation. Joint state/count kernels compose by matrix multiplication and count convolution.</p>'''
body+=eq(r'A^{(n+m)}_r=\sum_k A^{(m)}_{r-k}A^{(n)}_k')
body+='''<p>These are different experiments, and the extension retains both. The independent-copy certificates determine the smallest copy count meeting each stated risk target. Every smaller positive count is excluded by a certified lower bound, rather than by an assumed monotonicity rule.</p>'''
body+=table(['Risk target','Minimum copies, both networks','Cycle probability at crossing','Star probability at crossing'],copyrows)
body+=fig('counting_frontier','Figure 4. Certified wrong-direction probabilities for independent ten-period experiments and uninterrupted trajectories. The separate horizontal axes preserve the different information and dependence structures. At seventy independent copies, the cycle risk is 6.0345877770 × 10⁻¹⁰¹; at 700 consecutive periods it is 2.3569749568 × 10⁻¹⁰².')
body+=table(['Experiment at total process-time 4200','Cycle','Star'],[
    ['70 independent ten-period experiments',probability(cycle['copies']['70']),probability(star['copies']['70'])],
    ['One uninterrupted 700-period trajectory',probability(hcycle['horizons']['700']['probability']),probability(hstar['horizons']['700']['probability'])]])
body+='''<p>At equal total process-time, uninterrupted operation has approximately 25.60 times lower wrong-direction probability on the cycle and 25.85 times lower on the star. This comparison follows from the retained state correlations and periodic restarts in these specific protocols. It is a model result, not a universal ordering of continuous and restarted systems.</p>
<p>The cycle’s seventy-copy probability has relative enclosure width below 9.474 × 10<sup>−24</sup>; its 700-period probability has width below 2.426 × 10<sup>−22</sup>. The omitted-mass budgets remain far below the probabilities being reported. The calculation therefore resolves relative risk at these scales, rather than merely printing small floating-point values.</p>
<h3>A bound for every horizon</h3><p>A rational tilt selected by numerical optimization is certified separately through a positive-vector Collatz–Wielandt inequality. A componentwise bound M(z)v ≤ ρv gives the following all-horizon statement.</p>'''
body+=eq(r'\Pr(N_n\leq0)\leq\mathbb E[z^{N_n}]\leq C\rho^n\qquad(0<z<1,\ n\geq0)')
body+='''<p>The new cycle factor is ρ ≤ 0.7192219152797379, improving the earlier grid choice 0.7198899414. The corresponding all-horizon certificate is below 10<sup>−100</sup> from 701 periods onward; the star certificate reaches that level from 694 onward. These exponential bounds complement the narrower finite-horizon probability enclosures. They do not replace them.</p></section>
<section id="verification"><div class="section-number">V</div><h2>Executable mathematics and independent replay</h2>
<p>The release retains exact rational response coefficients, atomic witness laws, rank minors, complete kernel generators, coordinate row-span witnesses and dyadic probability endpoints. Each reported spectral bound replays without running an optimizer. The rational-access replay rebuilds the experiment from the original grid inputs. The counting laws are reproducible from the declared rational generators and hold protocols.</p>
<p>A separate control programme checks rational atomic measures against the spectral inequalities, compares sparse uniformization with independently constructed dense Arb exponentials, verifies scalar convolution against exact rational laws, checks joint-kernel composition coefficient by coefficient, and encloses independently generated compatible intervention roots. Its 3,145 controls pass. The access certificates receive an additional independent replay. The package test suite includes tampered-certificate rejection, denominator and pole rejection, interval-refinement checks, and the modular-coordinate counterexample.</p>
<div class="reproduction"><code>python -m pip install -r requirements.txt<br>python replay_extension.py<br>python -m possq.resolvent_campaign<br>python -m possq.cyclecraft_frontier results cycle star<br>python -m possq.cyclecraft_frontier results cycle star --horizons</code></div>
<p>The proofs and bounds are specialized to their declared classes: positive measures for spectral disclosure, symmetric first-order models for the Stieltjes bridge, rational reading maps for access certificates, and the specified Markov protocols for counting. No positive-measure assumption is imposed on Nonagon’s separate indefinite contour-zero machinery. The exact access and spectral tools can now answer a common operational question: which further observation changes the certified target answer enough to justify acquiring it?</p>
<h3>Mathematical lineage</h3><p>The generalized moment optimization viewpoint connects to Bertsimas and Popescu, <em>Optimal Inequalities in Probability Theory: A Convex Optimization Approach</em>, SIAM Journal on Optimization 15 (2005), 780–804. The positive spectral-measure representation belongs to the matrix-moment and quadrature framework developed by Golub and collaborators. The present extension supplies a raw-resolvent certificate implementation with exact sign verification, complete rational access witnesses, observation-level intervention inversion, and a sparse validated counting engine.</p>
</section>
</main><footer>NONAGON 14.1 × POSSQ 0.4 · Exact coefficients and dyadic endpoints accompany the source packages. This document is self-contained and works offline.</footer>
'''

css='''
*{box-sizing:border-box}html{scroll-behavior:smooth}body{margin:0;color:#142435;background:#f5f3ed;font-family:Georgia,'Times New Roman',serif;font-size:18px;line-height:1.72}header{padding:80px max(7vw,24px) 50px;background:#102638;color:#f5f3ed;overflow:hidden}.edition{font:12px/1.5 ui-monospace,monospace;letter-spacing:.17em;color:#8dd2c7}h1{font-weight:400;font-size:clamp(70px,12vw,156px);line-height:.91;letter-spacing:-.065em;margin:40px 0}h1 em{color:#b9d9ce;font-weight:400}.deck{font-size:clamp(23px,3vw,38px);line-height:1.28;max-width:840px}.lead{max-width:880px;color:#d5dedf;font-size:21px;line-height:1.7}.hero-rule{height:2px;background:linear-gradient(90deg,#b47d4d,#8dd2c7,transparent);margin-top:48px}nav{display:flex;flex-wrap:wrap;gap:14px 28px;padding:18px max(7vw,24px);border-bottom:1px solid #ccd2ce;font:12px/1.5 ui-monospace,monospace}a{color:#176b65;text-decoration:none}a:hover{text-decoration:underline}main{max-width:1120px;padding:0 46px;margin:0 auto}section{position:relative;padding:64px 0;border-bottom:1px solid #ccd2ce}section.abstract{padding:38px 0;font-size:19px}.section-number{font:15px ui-monospace,monospace;color:#8b5d37;letter-spacing:.2em;margin-bottom:12px}h2{font-size:clamp(31px,4vw,46px);font-weight:400;line-height:1.2;letter-spacing:-.035em;margin:0 0 28px;max-width:900px}h3{font-size:26px;line-height:1.35;font-weight:400;margin:40px 0 15px}p{max-width:900px;margin:0 0 21px}strong{font-weight:700}.equation{display:flex;align-items:center;gap:20px;overflow:auto;padding:24px 10px;max-width:100%;margin:22px 0;background:#eeeee7;border-left:2px solid #77a79c}.equation svg{max-width:calc(100% - 35px);height:auto;flex-shrink:0}.equation span{margin-left:auto;font:12px ui-monospace,monospace;color:#6d777a}.theorem{border-left:3px solid #b47d4d;padding:24px 26px;margin:30px 0;background:#fffdf7}.theorem p:last-child{margin-bottom:0}.theorem-label{font:12px/1.6 ui-monospace,monospace;text-transform:uppercase;letter-spacing:.08em;color:#8a552c;margin-bottom:16px}.tablewrap{overflow:auto;margin:28px 0}table{border-collapse:collapse;min-width:650px;width:100%;font:14px/1.65 ui-sans-serif,system-ui,sans-serif;font-variant-numeric:tabular-nums}th{text-align:left;border-top:2px solid #142435;border-bottom:1px solid #142435;padding:11px 12px;color:#1c5e59;font-weight:600}td{padding:11px 12px;border-bottom:1px solid #d7dbd5}tr:nth-child(even){background:#ebece5}.scientific{margin:40px -12px}.scientific>svg{width:100%;height:auto;display:block}.scientific figcaption{font:13px/1.7 ui-sans-serif,system-ui,sans-serif;color:#506069;padding:14px 14px 0;max-width:930px}.instrument{background:#102638;color:#d5e7e1;margin:30px 0;padding:14px 18px}.instrument-top{display:flex;gap:12px;align-items:center;flex-wrap:wrap;font:12px ui-monospace,monospace}.instrument-top span{margin-right:auto}button{font:12px ui-monospace,monospace;border:1px solid #689d99;background:transparent;color:#dbefea;padding:10px 12px;cursor:pointer}button:hover{background:#234338}.instrument canvas{width:100%;height:auto;display:block}.instrument p,.live-band p{font:12px/1.7 ui-sans-serif,system-ui,sans-serif;margin:10px 0}.live-band{margin:28px 0;padding:24px;background:#fffdf7;border-top:2px solid #176b65}.live-band label{font:15px ui-sans-serif,system-ui,sans-serif}input[type=range]{width:100%;accent-color:#176b65;display:block;margin:20px 0}output{display:block;font:15px/1.6 ui-monospace,monospace;color:#176b65}.live-band canvas{width:100%;height:auto}.reproduction{font-size:13px;background:#e8ece5;padding:22px 24px;overflow:auto;line-height:2;color:#194b46;margin:28px 0}footer{padding:42px max(7vw,24px);font:12px/1.7 ui-monospace,monospace;color:#53666a}@media(max-width:650px){body{font-size:17px;line-height:1.7}header{padding:48px 22px 34px}.lead{font-size:18px}main{padding:0 22px}section{padding:44px 0}.equation{margin:20px -8px}.equation svg{max-width:none}.scientific{margin:30px -14px}.theorem{padding:19px}.edition{font-size:10px;letter-spacing:.11em}nav{padding:16px 22px;font-size:11px;gap:10px 16px}.instrument{padding:12px;margin:26px -10px}h1{margin:35px 0}.deck{font-size:25px}}@media(prefers-reduced-motion:reduce){html{scroll-behavior:auto}}@media print{body{background:white;font-size:11pt}header{background:white;color:#142435;padding:20px 0}.lead{color:#142435}.edition{color:#176b65}h1{font-size:65pt}main{padding:0}nav,.instrument,.live-band{display:none}section{page-break-inside:auto}figure,.theorem{page-break-inside:avoid}.scientific{margin:20px 0}}
'''

series=next(c for c in spectral['campaigns'] if c['noise_relative_h0']=='1/100000000')['adaptive']
interactive={'series':[{'readings':a['readings'],'lower':a['lower_float'],'upper':a['upper_float']} for a in series],
             'eigenvalues':[float(__import__('fractions').Fraction(x)) for x in spectral['spectrum']]}
js=r'''
'use strict';
const DATA=__DATA__;
const canvas=document.getElementById('network'),ctx=canvas.getContext('2d');
let swapped=false,moving=!matchMedia('(prefers-reduced-motion: reduce)').matches,angle=.7;
document.getElementById('swap').onclick=()=>{swapped=!swapped;drawNetwork()};
document.getElementById('motion').onclick=()=>{moving=!moving;document.getElementById('motion').textContent=moving?'Pause motion':'Resume motion'};
function drawNetwork(){
 const W=canvas.width,H=canvas.height;ctx.clearRect(0,0,W,H);
 const axes=[[1,.15,.3],[.1,1,-.2],[-.3,.1,1],[.6,-.6,.4],[-.5,-.4,.6]];
 const pos=[];const ca=Math.cos(angle),sa=Math.sin(angle);
 for(let i=0;i<32;i++){let x=0,y=0,z=0;for(let a=0;a<5;a++){const v=(i>>a&1)?1:-1;x+=v*axes[a][0];y+=v*axes[a][1];z+=v*axes[a][2]};const rx=ca*x+sa*z,rz=-sa*x+ca*z;pos.push([W/2+rx*145,H/2+y*92+rz*18,rz])}
 const weights=swapped?[2,1,4,8,16]:[1,2,4,8,16];
 for(let i=0;i<32;i++)for(let a=0;a<5;a++){let j=i^(1<<a);if(i<j){ctx.strokeStyle=a<2?(swapped?'#d7a56c':'#78cbbd'):'#477282';ctx.lineWidth=.8+weights[a]/9;ctx.globalAlpha=.55;ctx.beginPath();ctx.moveTo(...pos[i].slice(0,2));ctx.lineTo(...pos[j].slice(0,2));ctx.stroke()}}
 ctx.globalAlpha=1;for(let i=0;i<32;i++){ctx.fillStyle=i===0?'#f0bd7d':'#bfdbd5';ctx.beginPath();ctx.arc(pos[i][0],pos[i][1],i===0?9:3.6,0,Math.PI*2);ctx.fill()}
 ctx.fillStyle='#d4e9e3';ctx.font='18px monospace';ctx.fillText(swapped?'K′ = P K Pᵀ':'K : original axis assignment',30,H-30);ctx.fillText('h₀(s) unchanged for every s',W-400,H-30);
}
let last=0;function frame(t){if(moving&&t-last>40){angle+=.004;drawNetwork();last=t};requestAnimationFrame(frame)}drawNetwork();requestAnimationFrame(frame);
const slider=document.getElementById('readings'),band=document.getElementById('band'),bc=band.getContext('2d');slider.max=DATA.series.length-1;
function drawBand(){let row=DATA.series[Number(slider.value)];document.getElementById('bandText').textContent=`${row.readings} readings · certified participation ${(100*row.lower).toFixed(5)}% to ${(100*row.upper).toFixed(5)}% · reference 6.25%`;
 const W=band.width,H=band.height;bc.clearRect(0,0,W,H);const x=v=>60+v/.45*(W-110);const y=112;
 bc.strokeStyle='#b4bfba';bc.lineWidth=1;bc.beginPath();bc.moveTo(x(0),y);bc.lineTo(x(.45),y);bc.stroke();
 for(let n=0;n<=.45001;n+=.05){bc.fillStyle='#64746c';bc.font='17px monospace';bc.fillText((n*100).toFixed(0)+'%',x(n)-14,y+42);bc.beginPath();bc.moveTo(x(n),y-5);bc.lineTo(x(n),y+5);bc.stroke()}
 bc.fillStyle='#8cc9b9';bc.fillRect(x(row.lower),y-24,x(row.upper)-x(row.lower),48);bc.strokeStyle='#176b65';bc.lineWidth=3;bc.strokeRect(x(row.lower),y-24,x(row.upper)-x(row.lower),48);
 bc.strokeStyle='#b47746';bc.lineWidth=3;bc.beginPath();bc.moveTo(x(.0625),y-54);bc.lineTo(x(.0625),y+54);bc.stroke();bc.fillStyle='#8c562d';bc.fillText('exact reference',x(.0625)+10,y-54);
}slider.oninput=drawBand;drawBand();
'''.replace('__DATA__',json.dumps(interactive))
doc='<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><meta name="description" content="Certified spectral disclosure, rational network access and stochastic counting laws beyond 1e-100."><title>Certified Disclosure · Nonagon 14.1 × possq 0.4</title><style>'+css+'</style></head><body>'+body+'<script>'+js+'</script></body></html>'
OUT.parent.mkdir(exist_ok=True,parents=True)
OUT.write_text(doc)
(HERE/'report_visualization.js').write_text(js)
print(OUT,OUT.stat().st_size,'bytes',equation_number,'equations')
