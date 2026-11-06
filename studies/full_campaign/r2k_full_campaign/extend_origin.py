from pathlib import Path
import sys,json,time,math
from flint import ctx,arb
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'audit_work/part_X'))
import origin_cauchy as oc
while not (HERE/'heldout_heights.json').exists():time.sleep(2)
d=json.loads((HERE/'heldout_heights.json').read_text());top=max(T for h in d['heights'] for T in h['T_quantiles'])
datasets=[json.loads(p.read_text()) for p in (HERE/'data').glob('zeros_*.json')]
low=min(math.floor(v['ordinates'][0]) for v in datasets);top=max(math.ceil(v['ordinates'][-1]) for v in datasets)
assert all(arb(v['ordinate_balls'][0])>low and arb(v['ordinate_balls'][-1])<top for v in datasets)
ctx.prec=192;ctx.cap=64;oc.N=20
A,B,dt=oc.compute_products(2000000);AI,BI,tail=oc.complete_euler_tail(A,B,2000000)
out={'status':'complete','precision_bits':ctx.prec,'prime_cutoff':2000000,'infinite_euler_tail':tail,
 'continuous_height_certificate':oc.continuous_range(AI,BI,low,top)}
(HERE/'extended_origin_results.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
