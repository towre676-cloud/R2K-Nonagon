"""Add the accepted fixed-line result while retaining every earlier package member."""
from pathlib import Path
import json,hashlib,zipfile,subprocess,sys
HERE=Path(__file__).resolve().parent;CAMPAIGN=HERE.parent;ROOT=CAMPAIGN.parent

def main():
    prior=json.loads((CAMPAIGN/'MANIFEST.json').read_text())
    for script in [HERE/'exact_controls.py',HERE/'validate_fixed_line.py',HERE/'build_report.py',CAMPAIGN/'build_report.py']:
        r=subprocess.run([sys.executable,str(script)],capture_output=True,text=True)
        if r.returncode:raise RuntimeError(str(script)+'\n'+r.stdout+r.stderr)
        print(script.name,'PASS',flush=True)
    readme=CAMPAIGN/'README.md';s=readme.read_text()
    if '## Fixed-line scalar fourth-order certificate' not in s:
        s+='\n## Fixed-line scalar fourth-order certificate\n\nThe fixed_line directory completes the sine-model covariance calculation for\nL=4, span 18000, in the original fixed-line anchor geometry. The report proves\nthe trace reduction, all spatial/separation/tail bounds, and a projection-DPP\nsample-centering inequality. Start with its HTML report. This does not certify\nthe arithmetic finite-height zeta sampling law. The earlier suggested scalar\ngeometry computation is now complete; the arithmetic correlation input and\ncomplete six-length fixed-line matrix remain open.\n'
        readme.write_text(s)
    # All previously manifested members remain, even if a later generic
    # scratch exclusion would otherwise have hidden one.
    names=set(prior['files'])
    for p in HERE.rglob('*'):
        if p.is_file() and '__pycache__' not in p.parts and p.name not in ['run.log','retry.log','near_checkpoint.json']:
            names.add(p.relative_to(ROOT).as_posix())
    manifest={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in sorted(names)}
    retention={'status':'PASS','previous_manifested_members':len(prior['files']),
          'previous_members_retained':len(set(prior['files'])&set(names)),
          'new_members':sorted(set(names)-set(prior['files'])),
          'updated_previous_members':[n for n,h in prior['files'].items() if manifest[n]!=h]}
    receipt=HERE/'package_retention_receipt.json';receipt.write_text(json.dumps(retention,indent=2))
    names.add(receipt.relative_to(ROOT).as_posix());manifest[receipt.relative_to(ROOT).as_posix()]=hashlib.sha256(receipt.read_bytes()).hexdigest()
    (CAMPAIGN/'MANIFEST.json').write_text(json.dumps({'schema':'r2k-full-campaign-manifest/1','files':manifest},indent=2))
    dest=ROOT.parent/'delivery_fixed_line';dest.mkdir(exist_ok=True);out=dest/'R2K_Nonagon_Full_Campaign.zip';tmp=out.with_suffix('.tmp')
    with zipfile.ZipFile(tmp,'w',zipfile.ZIP_DEFLATED) as z:
        for n in sorted(names):z.write(ROOT/n,n)
        z.write(CAMPAIGN/'MANIFEST.json','r2k_full_campaign/MANIFEST.json')
    with zipfile.ZipFile(tmp) as z:
        assert z.testzip() is None
        for n,h in manifest.items():assert hashlib.sha256(z.read(n)).hexdigest()==h
        assert 'fixed_line/R2K_Fixed_Line_Covariance_Report.html' in z.read('r2k_full_campaign/R2K_Nonagon_Full_Campaign_Report.html').decode()
    tmp.replace(out)
    print('PACKAGE',str(out),out.stat().st_size,'bytes;',len(manifest),'hashes verified;',len(prior['files']),'previous members retained',flush=True)

if __name__=='__main__':main()
