"""Replay the immutable selected confirmation design with normal Python."""
import hashlib,json,os,subprocess,sys
from pathlib import Path
root=Path(__file__).resolve().parent
freeze=json.loads((root/'freeze.json').read_text())
for name,wanted in freeze['sha256'].items():
    if hashlib.sha256((root/name).read_bytes()).hexdigest()!=wanted:
        raise ValueError('frozen input changed: '+name)
env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
subprocess.run([sys.executable,str(root/'acquisition.py'),'confirmation'],cwd=root,env=env,check=True)
subprocess.run([sys.executable,str(root/'audit.py')],cwd=root,env=env,check=True)
