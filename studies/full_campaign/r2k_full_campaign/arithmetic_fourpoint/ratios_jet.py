"""Mixed-derivative ratios compiler through order four.

Arithmetic mode encloses ONLY its stated finite Euler product. No infinite-prime
tail, ratios-conjecture remainder or short-window transfer is certified here.
RMT mode provides independent exact determinantal controls of the compiler.
"""
from itertools import combinations,product
from pathlib import Path
import math,time,json,sys
from flint import arb,acb,acb_series,acb_mat,ctx
HERE=Path(__file__).resolve().parent

class Jet:
    def __init__(self,n,c=0):self.n=n;self.a=[acb(0) for _ in range(1<<n)];self.a[0]=acb(c)
    def coerce(self,x):return x if isinstance(x,Jet) else Jet(self.n,x)
    def __add__(self,x):
        x=self.coerce(x);z=Jet(self.n);z.a=[a+b for a,b in zip(self.a,x.a)];return z
    __radd__=__add__
    def __neg__(self):z=Jet(self.n);z.a=[-a for a in self.a];return z
    def __sub__(self,x):return self+-self.coerce(x)
    def __rsub__(self,x):return self.coerce(x)+-self
    def __mul__(self,x):
        x=self.coerce(x);z=Jet(self.n)
        for mask in range(1<<self.n):
            sub=mask
            while True:
                if self.a[sub]!=0 and x.a[mask^sub]!=0:z.a[mask]+=self.a[sub]*x.a[mask^sub]
                if sub==0:break
                sub=(sub-1)&mask
        return z
    __rmul__=__mul__
    def inv(self):
        assert self.a[0].abs_lower()>0,'noninvertible jet constant'
        z=Jet(self.n,1/self.a[0])
        for mask in range(1,1<<self.n):
            sub=mask;v=acb(0)
            while sub:
                v+=self.a[sub]*z.a[mask^sub];sub=(sub-1)&mask
            z.a[mask]=-v/self.a[0]
        return z
    def __truediv__(self,x):return self*self.coerce(x).inv()
    def __rtruediv__(self,x):return self.coerce(x)*self.inv()
    def exp(self):
        base=self.a[0];v=self-Jet(self.n,base);out=Jet(self.n,1);term=Jet(self.n,1)
        for k in range(1,self.n+1):term=term*v/k;out+=term
        return out*base.exp()

def analytic(z,coeff):
    delta=z-Jet(z.n,z.a[0]);v=Jet(z.n,1);out=Jet(z.n)
    for a in coeff:out+=v*a;v*=delta
    return out

def zetaH(z):
    # H(z)=z*zeta(1+z), analytic and equal to one at z=0.
    c=acb_series([1+z.a[0],1],prec=z.n+1).zeta(deflate=True)
    return 1+z*analytic(z,[c[k] for k in range(z.n+1)])

def zfunc(z,mode,inverse=False):
    if mode=='rmt':
        val=1+ -(-z).exp()
        return val if inverse else val.inv()
    h=zetaH(z)
    return z/h if inverse else h/z

def prod(values,one):
    result=one
    for x in values:result=result*x
    return result

def Y(A,B,C,D,mode):
    n=A[0].n if A else B[0].n;out=Jet(n,1)
    for aa,bb in [(A,B),(C,D)]:
        for a in aa:
            for b in bb:out*=zfunc(a+b,mode)
    for aa,bb in [(A,D),(C,B)]:
        for a in aa:
            for b in bb:out*=zfunc(a+b,mode,True)
    return out

def local_arithmetic(A,B,C,D,p):
    n=A[0].n if A else B[0].n;one=Jet(n,1);lp=arb(p).log()
    root=lambda a:(-(a+arb(1)/2)*lp).exp()
    u=[root(a) for a in A];v=[root(b) for b in B];c=[root(a) for a in C];d=[root(b) for b in D]
    # Exact constant term by residues at zero and the B roots. Sizes B,D agree.
    ct=prod(d,one)/prod(v,one)
    for j,b in enumerate(v):
        numerator=prod((1-a*b for a in c),one)*prod((b-a for a in d),one)
        denominator=b*prod((1-a*b for a in u),one)*prod((b-a for k,a in enumerate(v) if k!=j),one)
        ct+=numerator/denominator
    norm=prod((1-a*b for a in u for b in v),one)*prod((1-a*b for a in c for b in d),one)
    norm/=prod((1-a*b for a in u for b in d),one)*prod((1-a*b for a in c for b in v),one)
    return ct*norm

def primes(P):
    a=bytearray(b'\1')*(P+1);a[:2]=b'\0\0'
    for p in range(2,math.isqrt(P)+1):
        if a[p]:a[p*p::p]=b'\0'*((P-p*p)//p+1)
    return [i for i in range(2,P+1) if a[i]]

def jstar(aa,bb,ell,mode,ps):
    if not aa or not bb:return acb(1 if not aa and not bb else 0)
    n=len(aa)+len(bb);A=[];B=[]
    for i,a in enumerate(aa+bb):
        z=Jet(n,a);z.a[1<<i]=1
        (A if i<len(aa) else B).append(z)
    C=[Jet(n,z.a[0]) for z in A];D=[Jet(n,z.a[0]) for z in B];out=Jet(n)
    for k in range(min(len(A),len(B))+1):
        for si in combinations(range(len(A)),k):
            for ti in combinations(range(len(B)),k):
                Ap=[z for i,z in enumerate(A) if i not in si]+[-B[i] for i in ti]
                Bp=[z for i,z in enumerate(B) if i not in ti]+[-A[i] for i in si]
                term=(-ell*(sum((A[i] for i in si),Jet(n))+sum((B[i] for i in ti),Jet(n)))).exp()*Y(Ap,Bp,C,D,mode)
                if mode=='arithmetic':
                    for p in ps:term*=local_arithmetic(Ap,Bp,C,D,p)
                out+=term
    val=out.a[-1];assert val.is_finite();return val

def correlation(xs,ell,mode='arithmetic',P=97):
    alpha=2*arb.pi()/ell;n=len(xs);ps=primes(P) if mode=='arithmetic' else [];memo={};out=acb(0);counts={}
    for assignment in product(range(3),repeat=n):
        A=tuple(i for i,a in enumerate(assignment) if a==0);B=tuple(i for i,a in enumerate(assignment) if a==1);m=sum(a==2 for a in assignment)
        if (not A or not B) and (A or B):continue
        key=(A,B);conj=False
        if B<A:key=(B,A);conj=True
        if key not in memo:
            aa=[acb(0,alpha*xs[i]) for i in key[0]];bb=[acb(0,-alpha*xs[i]) for i in key[1]]
            memo[key]=jstar(aa,bb,ell,mode,ps)
        v=memo[key].conjugate() if conj else memo[key]
        out+=ell**m*v;counts[len(A)+len(B)]=counts.get(len(A)+len(B),0)+1
    out/=ell**n;assert out.imag.contains(0)
    return out.real,counts,len(memo)

def rmt_controls():
    ctx.prec=192;ctx.cap=12;checks=[]
    for N in [2,5,9]:
        for n in [2,3,4]:
            xs=[arb(x) for x in ['0','0.37','1.13','2.41'][:n]]
            rho,counts,calls=correlation(xs,arb(N),'rmt')
            K=acb_mat([[acb(1) if i==j else acb(((xs[i]-xs[j])*arb.pi()).sin()/(N*((xs[i]-xs[j])*arb.pi()/N).sin())) for j in range(n)] for i in range(n)])
            det=K.det().real;assert rho.overlaps(det),(N,n,rho,det)
            checks.append({'N':N,'order':n,'ratios_compiler':str(rho),'independent_CUE_determinant':str(det),'assignment_counts_by_differentiated_order':counts,'distinct_J_calls':calls})
    return checks

def main():
    ctx.prec=192;ctx.cap=12;t=time.time();controls=rmt_controls();print('RMT controls PASS',len(controls),flush=True)
    xs=[arb(x) for x in ['0','0.37','1.13','2.41']];ell=(arb('9503636.19016305477494')/(2*arb.pi())).log();rows=[]
    result_path=HERE/'ratios_jet_results.json'
    if result_path.exists():
        rows=json.loads(result_path.read_text()).get('arithmetic_results',[])
        for row in rows:
            assert row['order'] in [2,3,4] and arb(row['ell']).overlaps(ell),'cached height contract mismatch'
            assert len(row['points'])==row['order'] and all(arb(a).overlaps(b) for a,b in zip(row['points'],xs[:row['order']])),'cached point contract mismatch'
    for P in [31,97,379,997,10007]:
        for n in [2,3,4]:
            if any(r['prime_cutoff']==P and r['order']==n for r in rows):continue
            value,counts,calls=correlation(xs[:n],ell,'arithmetic',P)
            rows.append({'prime_cutoff':P,'order':n,'points':[str(x) for x in xs[:n]],'ell':str(ell),
                         'finite_Euler_correlation_ball':value.str(45,radius=True,more=True),'assignments':counts,'distinct_J_calls':calls})
            print('ARITHMETIC',P,n,str(value),'seconds',round(time.time()-t,1),flush=True)
            result_path.write_text(json.dumps({'status':'running','RMT_controls':controls,'arithmetic_results':rows},indent=2))
    out={'status':'complete','scope':'mixed-derivative assembly verified against exact CUE; arithmetic balls enclose finite Euler inputs only; prime tails and ratios/short-window errors are not certified',
         'RMT_controls':controls,'arithmetic_results':rows,'elapsed_seconds':time.time()-t,'zeta_sampling_law_certified':False}
    (HERE/'ratios_jet_results.json').write_text(json.dumps(out,indent=2));print('COMPLETE',flush=True)

if __name__=='__main__':main()
