from pathlib import Path
import time,json,subprocess,os,sys
HERE=Path(__file__).resolve().parent
indices=[20000000,50000000]
while not all((HERE/'data'/f'empirical_{n}.json').exists() for n in indices):time.sleep(2)
heights={'heights':[]}
for n in indices:
 e=json.loads((HERE/'data'/f'empirical_{n}.json').read_text());heights['heights'].append({'n0':n,'T_quantiles':e['T_quantiles']})
p=HERE/'heldout_heights.json';p.write_text(json.dumps(heights,indent=2))
env={**os.environ,'R2K_INDICES':'20000000,50000000','R2K_CERT_TAG':'heldout_','R2K_HEIGHTS_FILE':str(p),'OPENBLAS_NUM_THREADS':'1'}
procs=[]
for file in ['numerical_campaign.py','certified_campaign.py']:
 log=open(HERE/('heldout_'+file+'.log'),'w')
 procs.append(subprocess.Popen([sys.executable,str(HERE/file)],env=env,stdout=log,stderr=subprocess.STDOUT))
for proc in procs:
 rc=proc.wait()
 if rc:raise RuntimeError('Heldout subprocess failed '+str(rc))
with open(HERE/'heldout_shifted_tail.log','w') as log:
 subprocess.run([sys.executable,str(HERE/'shifted_tail.py')],env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
subprocess.run([sys.executable,str(HERE/'ordinate_error.py')],check=True)
print('Both held-out campaigns and bounds completed.',flush=True)
