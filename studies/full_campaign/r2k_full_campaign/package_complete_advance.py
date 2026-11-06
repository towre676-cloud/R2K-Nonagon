"""Validate accepted workstreams and package an additive, hash-checked update."""
from pathlib import Path
import subprocess,sys,json,hashlib,zipfile,platform
import flint,numpy
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent

def main():
    prior=json.loads((HERE/'MANIFEST.json').read_text());assert len(prior['files'])>=175
    # These are receipt/content validators. Expensive numerical and contour
    # generators already ran; packaging does not silently replace their output.
    for name in ['root_hardening/validate_catalogue_replays.py','build_comprehensive_report.py','build_report.py']:
        r=subprocess.run([sys.executable,str(HERE/name)],capture_output=True,text=True)
        if r.returncode:raise RuntimeError(name+'\n'+r.stdout+r.stderr)
        print(name,'PASS',flush=True)
    for name in ['fixed_line/full_matrix_validation.json','arithmetic_fourpoint/ratios_validation.json',
                 'arithmetic_fourpoint/connected_bridge_results.json','arithmetic_fourpoint/dpp_obstruction.json',
                 'source_repairs/validation_receipt.json','source_repairs/patch_validation.json','comprehensive_validation.json']:
        assert json.loads((HERE/name).read_text())['status']=='PASS',name
    readme=HERE/'README.md';s=readme.read_text()
    if '## Comprehensive completion and arithmetic advance' not in s:
        s='''## Comprehensive completion and arithmetic advance

Start with R2K_Comprehensive_Advance_Report.html. It supersedes pending-task
statements in earlier stage reports where completed receipts are now supplied.
The six-length fixed-line sine covariance is complete. All 49 original contours
and the two separate fixtures have fresh hardened full replay (6434 segments).
Arithmetic three/four-point finite-input evaluation and a general connected
covariance bridge are supplied. A certified BK R2>1 obstruction excludes exact
scalar Hermitian-DPP completion of that pair model. Infinite-prime arithmetic
covariance and actual zeta sampling-law transfer remain open.

For packaging use package_complete_advance.py. Source review overlays and two
cleanly applicable patches are in source_repairs and root_hardening. The original
archive is unchanged. All previous package members are retained.

'''+s
        readme.write_text(s)
    provenance={'python':sys.version,'platform':platform.system(),'machine':platform.machine(),
         'python_flint':flint.__version__,'numpy':numpy.__version__,
         'parent_manifest_sha256':hashlib.sha256((HERE/'MANIFEST.json').read_bytes()).hexdigest(),
         'source_archive_name':'r2k_nonagon_v15.0.zip',
         'source_archive_sha256':hashlib.sha256((ROOT.parent/'original_sources/r2k_nonagon_v15.0.zip').read_bytes()).hexdigest(),
         'scope':'Linux native execution tested; portable Python ball fallback included; no formal proof-assistant claim'}
    (HERE/'complete_advance_provenance.json').write_text(json.dumps(provenance,indent=2))
    names=set(prior['files']);excluded={'MANIFEST.json','run.log','retry.log','near_checkpoint.json','finalize.log'}
    for p in HERE.rglob('*'):
        if p.is_file() and '__pycache__' not in p.parts and p.suffix!='.tmp' and p.name not in excluded:names.add(p.relative_to(ROOT).as_posix())
    missing=[n for n in prior['files'] if not (ROOT/n).is_file()];assert not missing,missing
    hashes={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in sorted(names)}
    retention={'status':'PASS','previous_members':len(prior['files']),'previous_members_retained':len(set(prior['files'])&names),
         'previous_members_updated':[n for n,h in prior['files'].items() if hashes[n]!=h],
         'new_members':sorted(names-set(prior['files']))}
    receipt=HERE/'complete_advance_retention.json';receipt.write_text(json.dumps(retention,indent=2))
    n=receipt.relative_to(ROOT).as_posix();names.add(n);hashes[n]=hashlib.sha256(receipt.read_bytes()).hexdigest()
    (HERE/'MANIFEST.json').write_text(json.dumps({'schema':'r2k-full-campaign-manifest/1','files':hashes},indent=2))
    dest=ROOT.parent/'delivery_complete';dest.mkdir(exist_ok=True);out=dest/'R2K_Nonagon_Full_Campaign.zip';tmp=out.with_suffix('.tmp')
    with zipfile.ZipFile(tmp,'w',zipfile.ZIP_DEFLATED) as z:
        for n in sorted(names):z.write(ROOT/n,n)
        z.write(HERE/'MANIFEST.json','r2k_full_campaign/MANIFEST.json')
    with zipfile.ZipFile(tmp) as z:
        assert z.testzip() is None
        for n,h in hashes.items():assert hashlib.sha256(z.read(n)).hexdigest()==h,n
        main=z.read('r2k_full_campaign/R2K_Nonagon_Full_Campaign_Report.html').decode()
        assert 'R2K_Comprehensive_Advance_Report.html' in main
    tmp.replace(out)
    print('PACKAGE',out,out.stat().st_size,'bytes;',len(hashes),'verified hashes;',len(prior['files']),'earlier members retained',flush=True)

if __name__=='__main__':main()
