#!/usr/bin/env python3
"""Replay canonical theory controls without mpmath; HP samples remain explicitly archival.

NumPy and SciPy are needed only for the descriptive stability campaign. Run from
any working directory. Original historical sources live in research/claude_v8.
"""
from pathlib import Path
import argparse
import hashlib
import json
import subprocess
import sys
import tempfile

HERE=Path(__file__).resolve().parent

def require(condition,message):
    if not condition:raise RuntimeError(message)

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=HERE/'theory_receipt.json')
    args=parser.parse_args()
    tasks=[]
    with tempfile.TemporaryDirectory(prefix='nonagon-theory-') as unrelated_cwd:
        for filename in ['intervene.py','origin.py','stability.py','finite_design_theory.py']:
            result=subprocess.run([sys.executable,str(HERE/filename)],cwd=unrelated_cwd,capture_output=True,text=True,timeout=120)
            require(result.returncode==0,f'{filename} failed: {result.stderr[-2000:]}')
            tasks.append({'script':filename,'exit_code':0,'unrelated_working_directory':True})
        verifier_results=[]
        for optimize in [False,True]:
            command=[sys.executable]+(['-O'] if optimize else [])+[str(HERE/'verify_intervention.py'),'80','1']
            result=subprocess.run(command,cwd=unrelated_cwd,capture_output=True,text=True,timeout=120)
            require(result.returncode==0,'independent algebraic campaign failed')
            parsed=json.loads(result.stdout)
            require(parsed['all_recovered_exactly_with_n(n+1)/2_measurements'],'failed exact recovery')
            verifier_results.append({'optimized_python':optimize,**parsed})
        require(verifier_results[0]['trials']==verifier_results[1]['trials']==80,'wrong trial count')
        (HERE/'verify_intervention.out').write_text(json.dumps({k:v for k,v in verifier_results[0].items() if k!='optimized_python'})+'\n')
        assembled=subprocess.run([sys.executable,str(HERE/'assemble8.py')],cwd=unrelated_cwd,capture_output=True,text=True,timeout=30)
        require(assembled.returncode==0,'chapter assembly failed: '+assembled.stderr)
        tasks.append({'script':'assemble8.py','exit_code':0,'unrelated_working_directory':True})
    exact=json.loads((HERE/'intervene.json').read_text())
    require(all(x['exact_recovery'] for x in exact['random_networks']),'primary exact campaign failure')
    require(exact['root_selection']['plus_root_ever_correct']==0,'wrong physicalroot')
    require(exact['kron_ambiguity']['h_equal_on_samples'],'Kron fixture failed')
    require(all(all(v.values()) for v in exact['damped_extension'].values() if isinstance(v,dict)),'damped extension failed')
    origin=json.loads((HERE/'origin.json').read_text())
    require(len(origin['rows'])==18 and origin['max_abs_residual_in_archived_se']<1.454,'archive correction regression')
    require(origin['uniform_remainder_certified'] is False,'overclaimed remainder')
    finite=json.loads((HERE/'finite_design_theory.json').read_text())
    require(finite['status']=='PASSED','finite theory controls failed')
    receipt={'status':'PASSED','tasks':tasks,'independent_algebraic_campaigns':verifier_results,
        'primary_exact_networks':exact['random_summary']['networks'],
        'physical_pairs':exact['root_selection']['pairs_checked'],
        'finite_design_controls':finite,
        'origin_correction_rows':18,'max_corrected_archived_se':origin['max_abs_residual_in_archived_se'],
        'hp_samples_recomputed':False,'uniform_origin_remainder_certified':False,
        'reference_is_external_certificate_validator':False,
        'source_sha256':{x.name:hashlib.sha256(x.read_bytes()).hexdigest() for x in sorted(HERE.glob('*.py'))}}
    args.output.write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({'status':receipt['status'],'primary_exact_networks':receipt['primary_exact_networks'],
        'reference_trials_per_mode':80,'physical_pairs':receipt['physical_pairs'],
        'finite_design_predicates':finite['predicates'],'origin_rows':18,'output':str(args.output)}))

if __name__=='__main__':main()
