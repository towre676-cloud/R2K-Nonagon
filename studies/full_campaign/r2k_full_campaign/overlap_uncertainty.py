"""Overlapping-batch means diagnostic for continuous-anchor influence counts.
This estimates sampling uncertainty under stationarity/mixing assumptions;
it is not a confidence certificate. Batch locations are integrated as piecewise
quadratics, avoiding an arbitrary grid of batch origins.
"""
from pathlib import Path
import json,math
import numpy as np
HERE=Path(__file__).resolve().parent
BLOCKS=[250.,500.,1000.,2000.,4000.]

def profile(u,L,mu):
 D=u[-1]-L;events={0.:0,D:0}
 for x in u:
  if 0<x-L<D:events[x-L]=events.get(x-L,0)+1
  if 0<x<D:events[x]=events.get(x,0)-1
 nodes=np.array(sorted(events));delta=np.array([events[x] for x in nodes]);n0=sum(0<x<=L for x in u)
 counts=n0+np.cumsum(delta);g=(counts[:-1]-mu)**2
 prefix=np.r_[0.,np.cumsum(np.diff(nodes)*g)];mean=prefix[-1]/D;out=[]
 for B in BLOCKS:
  stop=D-B
  br=np.unique(np.r_[0.,stop,nodes[(nodes>0)&(nodes<stop)],(nodes-B)[((nodes-B)>0)&((nodes-B)<stop)]])
  vals=(np.interp(br+B,nodes,prefix)-np.interp(br,nodes,prefix))/B-mean
  integral=np.sum(np.diff(br)*(vals[:-1]**2+vals[:-1]*vals[1:]+vals[1:]**2)/3)
  se=math.sqrt(B/(D-B)*integral/(D-B))
  out.append({'batch_length':B,'nominal_obm_se':se,'batch_origin_breakpoints':len(br)})
 return mean,out

def main():
 hs=[]
 for n0 in [100000,1000000,10000000,20000000,50000000]:
  p=HERE/'data'/f'empirical_{n0}.json'
  stable=HERE/'data'/f'empirical_stable_{n0}.json'
  if stable.exists():p=stable
  if not p.exists():continue
  e=json.loads(p.read_text());u=np.array([float.fromhex(x) for x in e['unfolded_float_hex']]);rows=[]
  for r in e['rows']:
   mean,obm=profile(u,r['L'],r['continuous']['mean_count']['float'])
   assert abs(mean-r['continuous']['variance']['float'])<1e-9
   rows.append({'L':r['L'],'influence_mean_roundoff_difference':mean-r['continuous']['variance']['float'],'profiles':obm})
  hs.append({'n0':n0,'rows':rows})
 out={'status':'complete','scope':'descriptive overlapping-batch means under sampling assumptions, no rigorous statistical interval','batch_lengths':BLOCKS,'heights':hs}
 (HERE/'overlap_uncertainty_results.json').write_text(json.dumps(out,indent=2));print('completed OBM diagnostics',len(hs),'windows')
if __name__=='__main__':main()
