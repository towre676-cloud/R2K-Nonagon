from pathlib import Path
import time,json,subprocess,sys,hashlib,zipfile
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
p=HERE/'continuous_height_results.json'
while not p.exists() or json.loads(p.read_text()).get('status')!='complete':time.sleep(2)
for name in ['overlap_uncertainty.py','summarize.py','validate_deliverable.py','plot_results.py','build_report.py','finite_window/validate_bridge.py','finite_window/plot_bridge.py','finite_window/build_report.py','fourth_order/validate_trace.py','fourth_order/validate_covariance.py','fourth_order/validate_periodic.py','fourth_order/build_report.py']:
 r=subprocess.run([sys.executable,str(HERE/name)],capture_output=True,text=True)
 (HERE/(name+'.receipt.log')).write_text(r.stdout+r.stderr)
 if r.returncode:raise RuntimeError(name+' failed: '+r.stderr)
 print(name,'PASS',flush=True)
report=(HERE/'R2K_Nonagon_Full_Campaign_Report.html').read_text()
assert report.count('<table>')==report.count('</table>') and report.count('<details>')==report.count('</details>')
assert 'Certified continuous-height averaging' in report
readme=HERE/'README.md';s=readme.read_text().replace('The secondary continuous-height averaging calculation is checkpointed while running.','The secondary continuous-height averaging calculation is complete for all thirty cases, with rigorous infinite-prime tail bounds.')
addition='''\n## Final acceptance\n\nAll 100000 acquired ordinate balls and all thirty exact continuous-anchor moment\nidentities pass validation. Ninety quantile-level infinite-product model targets\nand thirty continuous-height model averages have certified deterministic bounds.\nvalidation_receipt.json records the actual maximum errors and descriptive residuals.\nNo fitted correction term is introduced. Nominal block/OBM uncertainties are not\nrigorous statistical confidence intervals.\n\ncontinuous_height.py certifies the continuously averaged local stationary model\nusing the analytic height-integral reduction in the report. This does not prove\nthat the local stationary model equals the true finite-window zeta correlations.\nheight_identity_checks.py supplies independent exact-identity controls.\n\nThe original-arithmetic numerical_results_original_unfolding.json and empirical\nfiles are retained; the final numerical_results.json and heldout_numerical_results.json\nuse the more accurate centered Arb unfolding.\n\nFor a fresh quantile-level analytic replay, use certified_campaign.py --fresh.\nTo recompute continuous-height quadrature, preserve then remove its existing\ncheckpoint JSON before running continuous_height.py. Retained tails are independently\nrecomputable from the exact prime sum and the displayed residue-wheel estimate.\n'''
if '## Final acceptance' not in s:s+=addition
if '## Finite-window follow-up' not in s:s+='\n## Finite-window follow-up\n\nThe finite_window directory contains exact pair closure, thirty externally specified\ncoordinate-window comparisons, conditional sine-model centering bias, the fourth-order\ncounterexample, and descriptive joint multiplier diagnostics. Start with its HTML report.\nThese coordinate-window checks are exploratory, not newly preregistered held-out tests.\nNo unconditional significance level is claimed. finite_window/validation_receipt.json\nrecords the completed analytic and arithmetic checks.\n'
if '## Fourth-order circular benchmark' not in s:s+='\n## Fourth-order circular benchmark\n\nThe fourth_order directory supplies an all-frequency CUE trace covariance identity,\nthree complete Arb covariance matrices with analytic Fourier tails, thirty exact\nperiodic data controls, a current proof-status ledger and a detailed next-step target.\nThe circular model changes geometry and mean estimation. It does not certify\nthe actual zeta sampling law or justify combining its error bars with a BK mean.\nThe historical source audit at the package root is retained as provenance; its\npending tasks are superseded where the later reports provide completed receipts.\n'
readme.write_text(s)
extra=['R2K_Nonagon_Source_Audit.html','continuous_anchor_stats.py','audit_work/part_X/origin_cauchy.py','audit_work/part_IX/variance_quadrature.py','audit_work/part_VI/bk_multi_100000_1000000_10000000.json']
files=[p for p in HERE.rglob('*') if p.is_file() and '__pycache__' not in p.parts and p.suffix!='.tmp' and p.name not in ['MANIFEST.json','finalize.log','bench.py']]
files += [ROOT/f for f in extra]
manifest={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
(HERE/'MANIFEST.json').write_text(json.dumps({'schema':'r2k-full-campaign-manifest/1','files':manifest},indent=2))
files.append(HERE/'MANIFEST.json')
out=ROOT/'R2K_Nonagon_Full_Campaign.zip'
with zipfile.ZipFile(out,'w',zipfile.ZIP_DEFLATED) as z:
 for p in files:z.write(p,p.relative_to(ROOT).as_posix())
with zipfile.ZipFile(out) as z:
 assert z.testzip() is None
 for name,digest in manifest.items():assert hashlib.sha256(z.read(name)).hexdigest()==digest
print('FINAL_ARTIFACT',out.stat().st_size,'bytes; hashes and CRC verified',flush=True)
