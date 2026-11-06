"""Independent exact event sweep of every new circular empirical observable."""
from pathlib import Path
from fractions import Fraction as Q
from collections import defaultdict
from flint import arb,ctx
import json
HERE=Path(__file__).resolve().parent
ctx.prec=128

def sweep(points,L,S):
    events=defaultdict(int)
    for x in points:
        events[(x-L)%S]+=1;events[x]-=1
    n=sum(0<x<=L for x in points);start=Q(0);m1=m2=Q(0)
    for x,d in sorted(events.items()):
        m1+=(x-start)*n;m2+=(x-start)*n*n;n+=d;start=x
    m1+=(S-start)*n;m2+=(S-start)*n*n
    return m1/S,m2/S-(m1/S)**2

def near_pairs(xx,L,S):
    pairs=0;n=len(xx)
    for i,x in enumerate(xx):
        j=i+1
        while j<n and xx[j]-x<L:pairs+=1;j+=1
    for i,x in enumerate(xx):
        j=0
        while j<i and S-x+xx[j]<L:pairs+=1;j+=1
    return pairs

def main():
    p=json.loads((HERE/'cue_covariance_results.json').read_text());checks=0;maxerr=0
    controls=[]
    for e in p['empirical_periodized']:
        d=json.loads((HERE.parent/'finite_window'/('fixed_data_%s.json'%e['n0'])).read_text())
        xx=[Q(float.fromhex(x)) for x in d['points_hex']];eps=Q(1,10**11)
        assert arb(d['max_coordinate_error'])<arb(1)/10**11
        for row in e['rows']:
            L=Q(row['L']);mu,v=sweep(xx,L,p['span'])
            assert mu==Q(row['mean_exact']) and v==Q(row['variance_exact'])
            bound=4*eps*near_pairs(xx,L+2*eps,p['span'])/p['span']
            maxerr=max(maxerr,float(bound));checks+=1
            controls.append({'n0':e['n0'],'L':float(L),'coordinate_variance_error_bound_exact':str(bound)})
    out={'status':'PASS','independent_exact_periodic_event_pair_matches':checks,
         'maximum_coordinate_variance_error_bound':maxerr,'controls':controls}
    (HERE/'periodic_validation.json').write_text(json.dumps(out,indent=2))
    print('PASS',checks,'exact identities; max coordinate error',maxerr)

if __name__=='__main__':main()
