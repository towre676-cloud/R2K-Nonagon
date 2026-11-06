"""Exact centered-measure covariance contraction and explicit error-transfer constants.
The identity applies to general unit-intensity processes; no determinantal law
is needed. Exact finite-law fixtures test every collision partition.
"""
from pathlib import Path
from itertools import product
from functools import lru_cache
from fractions import Fraction as Q
import math,json,random
HERE=Path(__file__).resolve().parent

def partitions(n):
    if n==0:return [[]]
    out=[]
    for p in partitions(n-1):
        out.append(p+[[n-1]])
        for j in range(len(p)):
            q=[b[:] for b in p];q[j].append(n-1);out.append(q)
    return out

def main():
    rng=random.Random(20261105);checks=[];n=4;bits=list(product([0,1],repeat=n))
    for trial in range(20):
        weights=[rng.randrange(1,30) for _ in bits];den=sum(weights);law={b:Q(w,den) for b,w in zip(bits,weights)}
        E=lambda fn:sum((p*fn(b) for b,p in law.items()),Q(0))
        intens=[E(lambda b,i=i:b[i]) for i in range(n)]
        @lru_cache(None)
        def rho(ix):
            if len(set(ix))<len(ix):return Q(0)
            return E(lambda b:math.prod(b[i] for i in ix))
        @lru_cache(None)
        def kfac(ix):
            return sum((Q((-1)**(len(p)-1)*math.factorial(len(p)-1))*math.prod(rho(tuple(ix[j] for j in block)) for block in p) for p in partitions(len(ix))),Q(0))
        def matrix():
            mat=[[Q(0) for _ in range(n)] for _ in range(n)]
            for _ in range(2):
                a=[Q(rng.randrange(-4,5)) for i in range(n-1)];a.append(-sum(a[i]*intens[i] for i in range(n-1))/intens[-1])
                for i in range(n):
                    for j in range(n):mat[i][j]+=a[i]*a[j]
            assert all(sum(mat[i][j]*intens[j] for j in range(n))==0 for i in range(n))
            return mat
        G=matrix();H=matrix();fn=lambda b,M:sum(M[i][j]*b[i]*b[j] for i in range(n) for j in range(n))
        direct=E(lambda b:fn(b,G)*fn(b,H))-E(lambda b:fn(b,G))*E(lambda b:fn(b,H))
        terms=[Q(0)]*4
        for p in partitions(4):
            r=len(p);assign={i:j for j,block in enumerate(p) for i in block}
            for xx in product(range(n),repeat=r):
                ix=[xx[assign[i]] for i in range(4)]
                terms[r-1]+=G[ix[0]][ix[1]]*H[ix[2]][ix[3]]*kfac(xx)
        C=[[E(lambda b,i=i,j=j:b[i]*b[j])-intens[i]*intens[j] for j in range(n)] for i in range(n)]
        pairing=2*sum((G[i][j]*H[k][l]*C[i][k]*C[j][l] for i,j,k,l in product(range(n),repeat=4)),Q(0))
        assert sum(terms)+pairing==direct
        checks.append({'trial':trial,'direct_covariance':str(direct),'collision_orders':[str(x) for x in terms],'pairing':str(pairing)})
    LS=[Q(1,2),Q(1),Q(2),Q(4),Q(8),Q(16)];S=Q(18000);bounds=[]
    for i,L in enumerate(LS):
        for M in LS[i:]:
            dl=S-L;dm=S-M;ml=L/dl;mm=M/dm;dd=max(dl,dm)
            c2=L*M+4*ml*M*M+4*mm*L*L+4*min(ml*M*M,mm*L*L)
            c3=2*L*M*M+2*M*L*L+16*L*L*M*M/dd
            c4=4*L*L*M*M
            bounds.append({'L':str(L),'M':str(M),'factorial_cumulant_2_coefficient':str(c2),
                'factorial_cumulant_3_coefficient':str(c3),'factorial_cumulant_4_coefficient':str(c4),
                'pairing_linear_epsilon2_delta_coefficient':str(16*L*L*M*M/dd),
                'pairing_linear_epsilon2_B2_coefficient':str(16*L*L*M*M),
                'pairing_quadratic_epsilon2_coefficient':str(8*L*L*M*M)})
    out={'status':'PASS','exact_correlated_law_covariance_checks':checks,'collision_partitions':len(partitions(4)),
         'uniform_error_transfer_constants':bounds,
         'identity':'Cov(VL,VM)=integral gL12 gM34 [ordinary cumulant_4 + cumulant_2(13)cumulant_2(24)+cumulant_2(14)cumulant_2(23)]',
         'scope':'exact when intensity is the declared reference measure and each g has zero reference row integral; useful for general arithmetic correlations, without assuming a DPP'}
    (HERE/'connected_bridge_results.json').write_text(json.dumps(out,indent=2));print('PASS 20 exact correlated laws, all 15 collision partitions, 21 explicit error-transfer targets')

if __name__=='__main__':main()
