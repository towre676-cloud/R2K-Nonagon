#!/usr/bin/env python3
"""Verify all files, then replay selected experiments in a fresh copy."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

ROOT=Path(__file__).resolve().parent

def verify():
    rows=json.loads((ROOT/'MANIFEST.json').read_text())['files']
    for row in rows:
        path=ROOT/row['path']
        if not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest()!=row['sha256']:
            raise ValueError('File integrity mismatch: '+row['path'])
    return len(rows)

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--mode',choices=['verify','quick','all'],default='verify')
    p.add_argument('--numerics-python',default=sys.executable)
    p.add_argument('--adaptive-python',default=sys.executable)
    p.add_argument('--destination',type=Path)
    args=p.parse_args()
    print(f'Verified {verify()} file identities.',flush=True)
    if args.mode=='verify':return
    target=args.destination or Path(tempfile.mkdtemp(prefix='nonagon-frontier-'))/'campaign'
    if target.exists():raise ValueError('Replay destination must not exist')
    shutil.copytree(ROOT,target,ignore=shutil.ignore_patterns('__pycache__','*.pyc'))
    env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
    env.pop('PYTHONOPTIMIZE',None)
    commands=[('identifiability/network_obstruction.py',sys.executable,[]),
              ('identifiability/verify_explicit_certificate.py',sys.executable,[]),
              ('passivity/test_passivity.py',sys.executable,[]),
              ('passivity/certify_external.py',sys.executable,[])]
    if args.mode=='all':
        commands += [('adaptive/replay.py',args.adaptive_python,[]),
                     ('adaptive/audit.py',args.adaptive_python,[]),
                     ('policy/experiment.py',args.adaptive_python,['confirmation']),
                     ('policy/audit.py',args.adaptive_python,[]),
                     ('external/external_suite.py',args.numerics_python,['--phase','confirmation']),
                     ('external/summarize.py',args.numerics_python,[]),
                     ('independent_certification/replay.py',args.numerics_python,['--recompute'])]
    for script,python,extra in commands:
        path=target/script
        print('Replaying '+script,flush=True)
        subprocess.run([python,str(path),*extra],cwd=path.parent,env=env,check=True)
    print('Fresh replay outputs: '+str(target.resolve()))

if __name__=='__main__':main()
