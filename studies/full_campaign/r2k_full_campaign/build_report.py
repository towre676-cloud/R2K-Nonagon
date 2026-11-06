from pathlib import Path
import json,html,base64,hashlib
from flint import arb
HERE=Path(__file__).resolve().parent

def esc(s):return html.escape(str(s))
def main():
 d=json.loads((HERE/'summary_results.json').read_text());rows=d['cases'];complete=d['status']=='complete';N=len({r['n0'] for r in rows})*20000
 body=[]
 def section(title,text):
  body.append('<section><h2>'+esc(title)+'</h2>'+''.join('<p>'+esc(p.strip())+'</p>' for p in text.strip().split('\n\n'))+'</section>')
 body.append('<h1>R2K Nonagon: full variance certification and independent-window checks</h1>')
 body.append('<p class="lead">Fresh ordinates: '+str(N)+'. Completed comparisons: '+str(len(rows))+'. Status: '+esc(d['status'])+'. The certificates concern the declared infinite-product BK model; the identification of that model with finite-height zeta correlations remains conjectural.</p>')
 if (HERE/'comprehensive_validation.json').exists():
  section('Current comprehensive advance', 'The complete six-length fixed-line covariance is now certified under the sine model. The latest work also provides arithmetic three/four-point finite-prime computations, an exact connected covariance contraction, a proved Hermitian-DPP obstruction for BK, full hardened replay of the complete original contour catalogue, and reviewable manuscript corrections. Infinite-prime arithmetic covariance and zeta sampling-law transfer remain open.')
  body.append('<p><a href="R2K_Comprehensive_Advance_Report.html">Start with the comprehensive report, current proof ledger and all completed workstreams.</a></p>')
 if (HERE/'fixed_line/validation_receipt.json').exists():
  section('Latest advance: fixed-line fourth-order certificate', 'For L=4 on span 18000, the sine-model fourth-order covariance is now certified in the original line-window geometry, including a rigorous sample-centering allowance and a quantified comparison with the circular benchmark. This closes the proposed scalar geometry calculation. Arithmetic finite-height third/fourth-order correlations and the full six-length fixed-line matrix remain open.')
  body.append('<p><a href="fixed_line/R2K_Fixed_Line_Covariance_Report.html">Open the fixed-line theorem, analytic error bounds and current next-step target.</a></p>')
 if (HERE/'fourth_order/validation_receipt.json').exists():
  section('Fourth-order calibration and current proof ledger', 'The latest addendum derives and certifies complete covariance matrices for an explicitly circular CUE benchmark, including the fourth cumulant and all six-length cross covariances. It verifies thirty new circular empirical identities and reconciles the project proof status and constants. This conditional benchmark does not supply a finite-height zeta sampling law or replace the original observation contract.')
  body.append('<p><a href="fourth_order/R2K_Fourth_Order_Calibration_Report.html">Open the fourth-order calibration, proofs, matrices and current status ledger.</a> · <a href="../R2K_Nonagon_Source_Audit.html">Historical source audit (pending tasks superseded by these follow-ups)</a></p>')
 if (HERE/'finite_window/validation_receipt.json').exists():
  section('Finite-window follow-up', 'The follow-up verifies exact pair/anchor identities, adds thirty exploratory fixed-coordinate comparisons, derives an exact stationary centering-bias theorem and supplies a fourth-order obstruction to calibrating uncertainty from pair correlation alone. The coordinate-window residual persists, while the batch-based uncertainty varies substantially. These additions do not change the frozen primary campaign below.')
  body.append('<p><a href="finite_window/R2K_Finite_Window_Bridge_Report.html">Open the complete finite-window bridge report, proofs and all comparison tables.</a></p>')
 section('What changed', '''The origin-only repair has been extended to a certified calculation of the entire weighted number-variance integral. The calculation uses regularized origin coefficients, ball-enclosed Gauss–Legendre nodes and weights, and analytic Cauchy remainder bounds on every panel. It then adds rigorous infinite-prime tail bounds, rather than treating agreement between two finite prime cutoffs as a proof.

The original three 20000-zero windows have been reacquired. Continuous-anchor moments are computed exactly for the stored dyadic unfolded ordinates. A perturbation theorem propagates the enclosing ordinate balls and numerical unfolding error to the empirical variance. The final calculation evaluates the smooth unfolding in Arb and subtracts the first value before rounding displacements to binary64, avoiding loss from subtracting large floating-point values.

The new windows at starting indices 20000000 and 50000000 use the frozen formula, lengths, equal three-quantile prediction and ten-block primary diagnostic. No additive constant or other term is fitted. The complete table reports every length, including unfavorable comparisons. Alternate block partitions and overlapping-batch estimates are diagnostic sensitivity analyses; none is promoted to a statistical confidence certificate.''')
 body.append('<h2>Numerical error budget</h2><table><tr><th>Campaign</th><th>Comparisons</th><th>Maximum model radius</th><th>Maximum empirical ordinate/unfolding bound</th><th>Maximum combined radius</th><th>Largest |nominal z|</th></tr>')
 for kind in ['replication','held_out','all']:
  if kind not in d['summaries']:continue
  s=d['summaries'][kind];body.append('<tr>'+''.join('<td>'+esc(v)+'</td>' for v in [kind,s['cases'],f"{s['max_model_radius']:.4g}",f"{s['max_empirical_ordinate_error']:.4g}",f"{s['max_combined_deterministic_radius']:.4g}",f"{s['max_absolute_primary_nominal_z']:.5f}"])+'</tr>')
 body.append('</table><p>Each combined radius bounds the computational uncertainty in the exact smooth-unfolded finite-window variance minus the declared three-quantile infinite-product prediction. It excludes model discrepancy, nonstationarity and statistical extrapolation. All reported targets are checked against 10⁻⁶.</p>')
 section('The full finite-product integral', '''Set ℓ=log(T/(2π)), α=2π/ℓ and r=αx. The finite-product model is evaluated using the same arithmetic formula as the supplied Nonagon sources, with primes through P=2000000. Expensive Euler factors are evaluated once at a shared r-node and reused across all quantile heights. Every frozen target r=αL is an explicit panel boundary.

The origin contribution uses exact vanishing constant and linear coefficients, even coefficients through degree forty, and an analytic tail from a larger circle. Away from zero, each Gauss panel has a zero-free analytic neighborhood certified by ball enclosures of the zeta function over a covering square. A failed lower bound triggers refinement, rather than acceptance.

If a function is analytic on |r−c|≤ρ with bound M, n-node Gaussian quadrature is exact through degree 2n−1. With panel half-width h and q=h/ρ<1, Cauchy bounds the unintegrated Taylor tail by Mq^(2n)/(1−q). The integral and positive-weight quadrature of that tail each contribute at most 2h times that bound, giving total error at most 4hMq^(2n)/(1−q). The first moment rR(r) uses the corresponding majorant (|c|+ρ)M. These errors are included in the returned balls.

The native accelerator only performs scalar Arb/ACB arithmetic through exported FLINT functions; it does not replace ball arithmetic with floating-point estimates. Separate Python loops agree with it at six small-cutoff controls and three production-cutoff controls. The Linux accelerator is optional: a slower pure-Python Arb path is included. These are rigorous computational certificates conditional on the displayed analytic arguments and correct runtime arithmetic, rather than proof-assistant formalizations.''')
 section('New integrated B-tail theorem', '''The tail of B(z)=Σp(log p)²/(p^(1+z)−1)² admits an absolutely convergent expansion on the real r path:

B(ir)=Σp Σm≥2 (m−1)(log p)² p^(−m) exp(−imr log p).

Integrating the cosine terms against the triangular variance weight cancels their logarithmic frequency factors. The omitted-prime contribution is exactly

δΣ²_B(L,T)=−π⁻² Σp>P Σm≥2 [(m−1)/m²]p^(−m)[1−cos(mαL log p)].

Since (m−1)/m²≤1/4 and 0≤1−cos≤2, its absolute value is at most (2π²)⁻¹Σp>P1/[p(p−1)]. Summing over all integers instead gives the universal bound 1/(2π²P), about 2.54×10⁻⁸ at the production cutoff.

The executed bound is sharper: exact ball sums over primes from P to 100000000, plus an all-integer residue-wheel envelope beyond that endpoint, bound S₀=Σp>P(p−1)⁻². The B-tail bound S₀/(2π²) is then used. The wheel 2310 has 480 allowable residues. For decreasing f, each allowable progression above Q contributes at most f(Q)+2310⁻¹∫Q∞f; this gives an elementary certified tail without a prime-density conjecture.''')
 section('New contour-shift A-tail theorem', '''Let A_P(z) be the finite Euler product and A_>(z) its omitted tail. Put U_P(z)=ζ(1+z)ζ(1−z)A_P(z)[A_>(z)−1] and Fℓ(r)=exp(−iℓr)U_P(ir). The apparent pole at zero is removable because A_>(z)−1=O(z²). Thus Fℓ is analytic throughout the rectangle needed to move the r integral downward by iσ; the executed calculation uses σ=1/2.

The A-tail error is 4/(ℓ²α) times the real part of ∫₀ᴿ(L−r/α)Fℓ(r)dr, where R=αL. Deform the path to 0→−iσ→R−iσ→R. The bottom path gains exp(−ℓσ). On the left edge Fℓ(−it) is real, so the term proportional to L contributes only an imaginary integral. Its real contribution has weight t/α. On the right edge the triangular weight vanishes at R and also leaves weight it/α. Both vertical bounds therefore gain a factor t and avoid the original large L weight.

Explicitly, the absolute real integral is bounded by α⁻¹∫₀σt|Fℓ(−it)|dt + ∫₀ᴿ[L−u/α+σ/α]|Fℓ(u−iσ)|du + α⁻¹∫₀σt|Fℓ(R−it)|dt.

For Re z=t≥0, |p^(−z)−1|≤min(|z|log p,1+p^(−t)). With S₂=Σp>P(log p)²/(p−1)² and K(t)=(1+P^(−t))²S₀, the product difference is bounded by exp(4S₀)min(|z|²S₂,K(t)). Writing H(z)=zζ(1+z) removes the pole in the bounding expression as well. All three path integrals are bounded by outward ball enclosures over complete rectangles. They are upper sums over intervals, not sampled-point maxima.

This is the step that reduces the conservative real-path A-tail envelope to a useful full-integral bound. It is a theorem about the specified analytic Euler-product model. It does not prove a finite-height pair-correlation theorem for actual zeta ordinates.''')
 section('Exact counting and ordinate perturbations', '''For each length L, N_L(s) changes only at s=u_j and s=u_j−L. Sorting those events integrates N_L and N_L² over all admissible anchors. For dyadic stored coordinates, the integrals and variance are exact rational numbers. Both Var(N_L) and E[(N_L−L)²] are retained; their difference is (E[N_L]−L)².

To propagate ordinate errors, normalize anchors to y∈[0,1]. Each point contributes two threshold events u_j/D and (u_j−L)/D, where D is the admissible anchor span. The sum S of all certified event movements bounds ∫₀¹|N−N′|dy. If both count functions are at most M, then |Var(N)−Var(N′)|≤4MS. The code bounds M by an exact sliding-window count with the coordinate-error envelope included. Independent rational perturbation examples verify the estimate against exact recomputed variances.

The original binary64 unfolding evaluated values of order the zero index and then subtracted two large numbers. At high index that loses precision unnecessarily. The final pipeline evaluates and centers in Arb first; the remaining binary64 displacement error is then explicitly propagated. Original-arithmetic replication records are retained separately.''')
 section('Statistical interpretation and remaining gaps', '''The common 0.00319999 discrepancy at L=16 has the known omitted-origin explanation. The fresh results do not restore it as a new correction law. Computational error bars are now much smaller than the observed finite-window residuals. Such residuals need not vanish: the numerical model is a conjectural statistical prediction, not an exact equality for each deterministic window.

The primary z column divides by the predeclared ten-block nominal standard error. Five- and twenty-block versions can differ considerably. For example, the original highest replication window at L=2 appears roughly −1.47 ten-block standard errors away, while a five-block estimate makes it about −6.46. Overlapping-batch diagnostics at lengths 250, 500, 1000, 2000 and 4000 show why choosing one favorable partition would be misleading. No column here is a calibrated simultaneous p-value.

The primary contract rigorously encloses the equal average of three specified height-quantile predictions. The secondary calculation reported below removes this three-quantile approximation for the continuous average of the local BK variance over all admissible anchor heights. That still does not prove that finite-window zeta statistics equal an average of local stationary predictions: a local nonstationary correlation theorem with an explicit comparison remainder remains necessary. The within-window prediction spread is reported separately, and can be much larger than the numerical radius.

The exact smooth unfolding is also a declared convention; numerical enclosure of that convention does not prove stationarity or remove all finite-height density corrections. The detailed BK/ratios formula remains conditional in its number-theoretic interpretation. The most decisive remaining mathematical problem is therefore a justified finite-window/height-averaging comparison with an explicit remainder, together with a defensible statistical uncertainty model. Further refinement of the already tiny origin-only constant is not the bottleneck.''')
 body.append('<h2>All comparisons</h2><table><tr><th>Set</th><th>Starting index</th><th>L</th><th>Empirical variance</th><th>Residual</th><th>Combined deterministic radius</th><th>Primary nominal z</th><th>5-block z</th><th>20-block z</th><th>Quantile model spread</th></tr>')
 for r in rows:
  vals=[r['kind'],r['n0'],r['L'],f"{r['empirical_variance_float']:.10f}",f"{r['residual_midpoint']:+.8g}",f"{float(arb(r['combined_deterministic_radius']).upper()):.3g}",f"{r['primary_nominal_z']:+.4f}",f"{r['nominal_5block_z']:+.3f}",f"{r['nominal_20block_z']:+.3f}",f"{r['height_quantile_model_spread']:.3g}"]
  body.append('<tr>'+''.join('<td>'+esc(v)+'</td>' for v in vals)+'</tr>')
 body.append('</table>')
 if (HERE/'residuals.png').exists():body.append('<h2>Residuals and error scales</h2><img alt="Residuals and numerical versus sampling error" style="width:100%" src="data:image/png;base64,'+base64.b64encode((HERE/'residuals.png').read_bytes()).decode()+'">')
 if (HERE/'extended_origin_results.json').exists():
  o=json.loads((HERE/'extended_origin_results.json').read_text())['continuous_height_certificate']
  section('Extended continuous-height origin certificate', 'The original Part X certificate was rerun over the interval from '+str(o['T_interval'][0])+' to '+str(o['T_interval'][1])+'. The infinite-product ratio R₂(x)/x², for 0<|x|≤10⁻⁴, is enclosed between '+o['uniform_R2_over_x2_lower']+' and '+o['uniform_R2_over_x2_upper']+'. The c≤3 gate returns '+str(o['certified_c_le_3'])+'. This is an extension of the stated finite height interval, not an all-heights assertion.')
 if (HERE/'continuous_height_results.json').exists():
  avg=json.loads((HERE/'continuous_height_results.json').read_text())
  if avg.get('status')=='complete':
   section('Certified continuous-height averaging', '''This secondary comparison uses the exact smooth-unfolded anchor range rather than three representative heights. Write N̄(T)=exp(ℓ)(ℓ−1)+7/8, with ℓ=log(T/(2π)). For a window from its first to its last ordinate, the admissible upper height solves N̄(T_b)=N̄(T_last)−L. The normalization is D=N̄(T_last)−N̄(T_first)−L, and dN̄=ℓexp(ℓ)dℓ.

Let Dζ(r)=(log ζ)″(1+ir), C(r)=ζ(1+ir)ζ(1−ir)A_P(ir), R_a=2πL/ℓ_a and b(r)=min(ℓ_b,2πL/r). Define K_q(r) as the difference, at b(r) and ℓ_a, of exp(qℓ)[L/q−r(ℓ/q−1/q²)/(2π)]. Fubini's theorem gives the continuous finite-product height average as L + 2/(πD) Re ∫₀ᴿᵃ{[Dζ(r)−B_P(ir)]K₁(r)+C(r)K_(1−ir)(r)}dr.

The grouped expression is regular at the origin; its singular-looking components must not be integrated separately there. The implementation integrates its even origin polynomial analytically in height, using elementary primitives of exp(ℓ) times powers of ℓ. The two analytic branches are split at r=2πL/ℓ_b. Ball-enclosed Gaussian quadrature and Cauchy tails cover the remaining panels. A uniform-in-height contour-shift bound covers the infinite-prime tail. This certifies the specified continuously averaged local model, not a number-theoretic local-stationarity assertion.''')
   body.append('<table><tr><th>Starting index</th><th>L</th><th>Continuous-height model</th><th>Change from three quantiles</th><th>Model radius</th><th>Nominal residual z</th></tr>')
   for ar in avg['results']:
    rr=next(r for r in rows if r['n0']==ar['n0'] and r['L']==ar['L']);model=arb(ar['infinite_product_height_average_ball']);primary=arb(rr['model_infinite_ball']);delta=model-primary;res=rr['empirical_variance_float']-float(model.mid())
    vals=[ar['n0'],ar['L'],f'{float(model.mid()):.10f}',f'{float(delta.mid()):+.5g}',f'{float(model.rad().upper()):.3g}',f'{res/rr["nominal_block_se"]:+.4f}']
    body.append('<tr>'+''.join('<td>'+esc(v)+'</td>' for v in vals)+'</tr>')
   body.append('</table>')
 section('What was not claimed', '''No unconditional zeta pair-correlation theorem, global Riemann-hypothesis theorem, formal proof-assistant certificate, calibrated simultaneous confidence region, or proof that a new finite-height term is absent follows from this campaign. The three very high contour tiles deferred in the earlier audit were not required to acquire these lower-height production windows and are not relabeled as freshly replayed.

The original v15 source archive was preserved. The earlier printed circle-majorant mismatch should still be corrected to an outward-rounded value such as 2742.431538. Existing intervention optima retain their separate acquisition/noise hypotheses. The new sources and receipts are an additive research package, not an unreviewed replacement of all older project files.''')
 section('Reproduction', '''Install the pinned dependencies, build the optional native accelerator on Linux with build_native.py, and run the programs from the extracted package root. The reference origin and quadrature modules are supplied under audit_work. The numerical campaign reuses the retained 20000-ordinate files and model curves. For a genuinely fresh finite-product analytic replay, pass --fresh to certified_campaign.py; the default can resume a checked local prefix.

Run proof_checks.py, native_euler.py, numerical_campaign.py, certified_campaign.py --fresh, shifted_tail.py, ordinate_error.py, overlap_uncertainty.py and summarize.py. Held-out runs set R2K_INDICES, R2K_CERT_TAG and R2K_HEIGHTS_FILE as documented by run_heldout.py. Raw ordinate balls, dyadic coordinates, exact rational moments, all model curves, all Gaussian panel bounds, the frozen held-out protocol and every comparison are retained. Programs must run without Python -O.''')
 for name in ['summary_results.json','proof_checks.json','native_production_checks.json','tail_constants.json','heldout_protocol.json','extended_origin_results.json','height_identity_checks.json','validation_receipt.json']:
  p=HERE/name
  if p.exists():body.append('<details><summary>'+esc(name)+'</summary><pre>'+esc(p.read_text())+'</pre></details>')
 style='body{font-family:system-ui,sans-serif;line-height:1.62;color:#172b3a;background:#fbfcff;max-width:1180px;margin:36px auto;padding:0 24px}h1,h2{line-height:1.25;color:#173f68}section{margin:2.2em 0}.lead{font-size:1.15em;padding:18px;background:#e9f2fb}table{border-collapse:collapse;width:100%;font-size:12px}th,td{border:1px solid #c9d6e4;padding:7px;text-align:left}pre{white-space:pre-wrap;overflow-wrap:anywhere;font-size:12px;background:#edf2f7;padding:16px}summary{cursor:pointer;font-weight:600;background:#e5eef8;padding:12px}'
 report='<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>R2K full variance campaign</title><style>'+style+'</style><main>'+''.join(body)+'</main></html>'
 (HERE/'R2K_Nonagon_Full_Campaign_Report.html').write_text(report);print('Report written:',len(report),'characters; status',d['status'])
if __name__=='__main__':main()
