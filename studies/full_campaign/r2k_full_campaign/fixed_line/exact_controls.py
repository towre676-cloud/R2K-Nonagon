"""Exact finite-projection checks of joint-count and centering cumulants."""
from pathlib import Path
from itertools import combinations
from fractions import Fraction as F
from collections import defaultdict
import sys,json,random
import sympy as s
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'fourth_order'))
from trace_derivation import PARTS

def canonical(word):
    rev=word[::-1];return min(tuple(w[i:]+w[:i]) for w in [word,rev] for i in range(len(w)))

def symbolic_joint():
    tags=['A','A','B','B'];out=defaultdict(F)
    for p in PARTS:
        word=[]
        for block in p:
            t={tags[i] for i in block};word.append(next(iter(t)) if len(t)==1 else 'C')
        out[canonical(word)]+=F((-1)**(len(p)-1),len(p))
    expected={('C',):F(1),('A','B'):F(-1),('C','C'):F(-2),('A','C'):F(-2),('B','C'):F(-2),
              ('A','B','B'):F(2),('A','A','B'):F(2),('A','B','C'):F(8),('A','A','B','B'):F(-4),('A','B','A','B'):F(-2)}
    assert {k:v for k,v in out.items() if v}==expected
    return {''.join(k):str(v) for k,v in out.items() if v}

def main():
    sym=symbolic_joint();rng=random.Random(20261105);checked=0;bounds=0
    n=5;rank=2;Id=s.eye(n)
    for trial in range(20):
        vec=s.Matrix([rng.randrange(1,5) for _ in range(n)]);O=Id-2*vec*vec.T/(vec.T*vec)[0]
        P=O*s.diag(1,1,0,0,0)*O.T;Q=Id-P
        law={sub:P.extract(sub,sub).det() for sub in combinations(range(n),rank)}
        assert sum(law.values())==1 and all(p>=0 for p in law.values())
        aa=[rng.randrange(2) for _ in range(n)];bb=[rng.randrange(2) for _ in range(n)]
        A=s.diag(*aa);B=s.diag(*bb);C=A*B;X=A*P;Y=B*P;Z=C*P
        kap=s.trace(Z)-s.trace(X*Y)-2*s.trace(Z*Z)-2*s.trace(Z*X)-2*s.trace(Z*Y)+2*s.trace(X*Y*Y)+2*s.trace(Y*X*X)+8*s.trace(Z*X*Y)-4*s.trace(X*X*Y*Y)-2*s.trace(X*Y*X*Y)
        ma=s.trace(X);mb=s.trace(Y)
        ex=lambda fn:sum(p*fn(sub) for sub,p in law.items())
        u=lambda sub:sum(aa[i] for i in sub)-ma
        v=lambda sub:sum(bb[i] for i in sub)-mb
        cov=ex(lambda sub:u(sub)*v(sub))
        direct=ex(lambda sub:u(sub)**2*v(sub)**2)-ex(lambda sub:u(sub)**2)*ex(lambda sub:v(sub)**2)
        assert direct==kap+2*cov**2;checked+=1
        ff=[s.Rational(rng.randrange(6),7) for _ in range(n)];f=s.diag(*ff);mean=s.trace(f*P)
        yy=lambda sub:sum(ff[i] for i in sub)-mean
        var=ex(lambda sub:yy(sub)**2);k4=ex(lambda sub:yy(sub)**4)-3*var**2
        Ap=P*f*P;Bp=P*f*Q;Cp=Q*f*Q;T=Ap*Bp-Bp*Cp
        rhs=s.trace(T*T.T)-2*s.trace((Bp*Bp.T)**2)
        assert k4==rhs and var==s.trace(Bp*Bp.T)
        m=max(ff)-min(ff)
        assert k4<=m*m*var
        assert ex(lambda sub:(yy(sub)**2-var)**2)<=2*var**2+m*m*var;bounds+=1
    out={'status':'PASS','symbolic_75_partition_joint_cumulant':sym,'exact_joint_count_covariance_checks':checked,
         'exact_projection_block_identity_and_centering_bounds':bounds,'scope':'controls of the displayed all-degree operator identities, not a substitute for their analytic proof'}
    (HERE/'exact_controls.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))

if __name__=='__main__':main()
