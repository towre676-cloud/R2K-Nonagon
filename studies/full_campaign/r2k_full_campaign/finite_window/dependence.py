"""Anchor-coordinate influence batches and joint multiplier diagnostics.
Reference sign enumeration is computationally exact for 10 batches, but is
not an exact test for zeta zeros. No symmetry/independence assumption is proved.
"""
from pathlib import Path
import json,math
import numpy as np
from fractions import Fraction as Q
HERE=Path(__file__).resolve().parent;BASE=HERE.parent
LENGTHS=[.5,1,2,4,8,16];BATCH_COUNTS=[5,10,20,40,80]

def field(u,S,L,mu):
 D=S-L;events={0.:0,D:0}
 for x in u:
  if 0<x<D:events[x]=events.get(x,0)-1
  if 0<x-L<D:events[x-L]=events.get(x-L,0)+1
 nodes=np.array(sorted(events));counts=sum(0<x<=L for x in u)+np.cumsum([events[t] for t in nodes]);width=np.diff(nodes)
 G=np.r_[0.,np.cumsum(width*(counts[:-1]-mu)**2)]
 return nodes,G,D

def main():
 out=[];references={}
 for kind in ['point_selected','fixed_span']:
  for n0 in [100000,1000000,10000000,20000000,50000000]:
   if kind=='point_selected':
    e=json.loads((BASE/'data'/f'empirical_stable_{n0}.json').read_text());u=np.array([float.fromhex(x) for x in e['unfolded_float_hex']]);S=u[-1]
    means=[r['continuous']['mean_count']['float'] for r in e['rows']];variances=[r['continuous']['variance']['float'] for r in e['rows']]
   else:
    e=json.loads((HERE/f'fixed_data_{n0}.json').read_text());u=np.array([float.fromhex(x) for x in e['points_hex']]);S=18000.
    means=[float(Q(r['mean_exact'])) for r in e['rows']];variances=[float(Q(r['variance_exact'])) for r in e['rows']]
   fields=[field(u,S,L,mu) for L,mu in zip(LENGTHS,means)];tiers=[]
   for b in BATCH_COUNTS:
    columns=[]
    for (nodes,G,D),v in zip(fields,variances):
     edges=np.arange(b+1)*D/b;mb=np.diff(np.interp(edges,nodes,G))*b/D
     assert abs(float(np.mean(mb))-v)<1e-9
     columns.append(mb)
    matrix=np.column_stack(columns);centered=matrix-matrix.mean(axis=0);cov=centered.T@centered/(b*(b-1));se=np.sqrt(np.diag(cov));assert np.all(se>0)
    if b==10:
     signs=2*((np.arange(1<<b)[:,None]>>np.arange(b))&1)-1;mode='all 1024 sign patterns'
    else:
     signs=np.random.default_rng(70000+n0+b+(0 if kind=='point_selected' else 100)).choice([-1.,1.],size=(32768,b));mode='32768 seeded Rademacher patterns'
    maxT=np.max(np.abs((signs@centered/b)/se),axis=1)
    key=f'{kind}_{n0}_{b}';references[key]=maxT
    tiers.append({'batches':b,'influence_batch_means':matrix.tolist(),'mean_covariance_estimate':cov.tolist(),'nominal_se':se.tolist(),
     'reference_mode':mode,'reference_90':float(np.quantile(maxT,.9)),'reference_95':float(np.quantile(maxT,.95)),
     'reference_99':float(np.quantile(maxT,.99)),'reference_npz_key':key,'reference_draw_count':len(maxT)})
   out.append({'geometry':kind,'n0':n0,'span':float(S),'lengths':LENGTHS,'empirical_variances':variances,'tiers':tiers});print('dependence',kind,n0,flush=True)
 np.savez_compressed(HERE/'multiplier_reference.npz',**references)
 (HERE/'dependence_results.json').write_text(json.dumps({'status':'complete','windows':out,
 'scope':'descriptive joint block-multiplier reference; does not certify statistical coverage or transfer a DPP central limit theorem to zeta zeros'},indent=2))
if __name__=='__main__':main()
