"""Exact continuous-anchor statistics for stored rational/dyadic ordinates.
Input JSON on stdin: {"points": [...], "length": L}. JSON numbers are read
as their decimal rationals; Python callers passing floats retain their exact
dyadic values. Zero/unfolding uncertainty is outside this finite computation.
"""
from fractions import Fraction as Q
from collections import defaultdict
import json,sys

def statistics(points,length):
    u=[Q(x) for x in points];L=Q(length)
    if len(u)<2 or any(a>=b for a,b in zip(u,u[1:])) or L<=0:
        raise ValueError('strictly increasing ordinates and positive length required')
    A=u[0];B=u[-1]-L
    if B<=A:raise ValueError('no positive-length anchor domain')
    events=defaultdict(int)
    for x in u:
        if A<x-L<B:events[x-L]+=1
        if A<x<B:events[x]-=1
    count=sum(A<x<=A+L for x in u);last=A;m1=Q(0);m2=Q(0)
    for x,delta in sorted(events.items()):
        gap=x-last;m1+=gap*count;m2+=gap*count*count
        count+=delta;last=x
    gap=B-last;m1+=gap*count;m2+=gap*count*count
    mean=m1/(B-A);second=m2/(B-A);variance=second-mean*mean
    return {'anchor_lower':A,'anchor_upper':B,'length':L,'mean_count':mean,
            'second_moment_count':second,'variance':variance,
            'mean_squared_deviation_from_L':second-2*L*mean+L*L,
            'mean_count_bias_squared':(mean-L)**2}

def selftest():
    import random
    r=statistics([0,1,3,4],1)
    assert r['mean_count']==Q(2,3) and r['variance']==Q(2,9)
    assert r['mean_squared_deviation_from_L']==Q(1,3)
    rng=random.Random(379)
    for _ in range(100):
        u=sorted(Q(x,7) for x in rng.sample(range(100),10));L=Q(rng.randint(1,10),7)
        rr=statistics(u,L);A=u[0];B=u[-1]-L
        breaks=sorted({A,B}|{x for z in u for x in [z,z-L] if A<x<B})
        m1=Q(0);m2=Q(0)
        for left,right in zip(breaks,breaks[1:]):
            s=(left+right)/2;n=sum(s<=x<s+L for x in u)
            m1+=(right-left)*n;m2+=(right-left)*n*n
        assert rr['mean_count']==m1/(B-A)
        assert rr['second_moment_count']==m2/(B-A)
        assert rr['mean_squared_deviation_from_L']==rr['variance']+rr['mean_count_bias_squared']
    return {'exact_fixture':'2/3 mean, 2/9 variance, 1/3 mean squared deviation',
            'independent_midpoint_enumeration_cases':100,'status':'PASS',
            'production_zero_windows_evaluated':False}

if __name__=='__main__':
    if '--selftest' in sys.argv:print(json.dumps(selftest(),indent=2))
    else:
        data=json.load(sys.stdin,parse_float=Q)
        print(json.dumps({k:str(v) for k,v in statistics(data['points'],data['length']).items()},indent=2))
