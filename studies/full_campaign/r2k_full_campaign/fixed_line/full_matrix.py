"""All 21 fixed-line sine-DPP influence covariances, with centering enclosures.
Native evaluation is optional and uses FLINT ball arithmetic, not floating BLAS.
"""
from pathlib import Path
from collections import defaultdict
import sys,json,math,time,ctypes,os
from flint import arb,acb,arb_mat,acb_mat,ctx
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'finite_window'))
from sine_bias import bias
from sine_covariance import integer_power,txt,enc
LS=[.5,1.,2.,4.,8.,16.];S=18000;J=8;NH=20;BITS=128;NATIVE_BITS=96

def order(l):
    n=8
    while True:
        rho=arb(n)/arb.pi();q=arb(l)/(2*rho)
        if q<1 and 2*(2*n*arb(1)).exp()*q**(2*n)/(1-q)<arb('1e-20'):return n
        n+=1

def qrule(l):
    n=order(l);q=[arb.legendre_p_root(n,k,weight=True) for k in range(n)]
    assert all(w>0 and abs(x)<1 for x,w in q)
    return q

def make_points(l,q,center=arb(0)):
    return [center+x*arb(l)/2 for x,w in q],[(w*arb(l)/2).sqrt() for x,w in q]

def kernel(x,w,y,v):return arb_mat([[(b-a).sinc_pi()*ww*vv for b,vv in zip(y,v)] for a,ww in zip(x,w)])
def inner(A,B):return sum((a*b for a,b in zip(A.entries(),B.entries())),arb(0))

class Counts:
    def __init__(self,l,m,native=True):
        self.l=arb(l);self.m=arb(m);self.qa=qrule(l);self.qb=qrule(m);self.qc=qrule(min(l,m))
        self.x,self.w=make_points(l,self.qa);self.y,self.v=make_points(m,self.qb)
        self.A=kernel(self.x,self.w,self.x,self.w);self.B=kernel(self.y,self.v,self.y,self.v)
        self.lib=None
        if native and (HERE/'native_counts.so').exists() and os.environ.get('R2K_COUNTS_PYTHON')!='1':
            try:lib=ctypes.CDLL(str(HERE/'native_counts.so'))
            except OSError:return
            P=ctypes.POINTER(ctypes.c_char_p)
            lib.counts_config.argtypes=[ctypes.c_long]*4+[ctypes.c_char_p]*2+[P]*6;lib.counts_config.restype=ctypes.c_int
            lib.counts_F.argtypes=[ctypes.c_char_p,ctypes.c_int];lib.counts_F.restype=ctypes.c_void_p
            lib.counts_release.argtypes=[ctypes.c_void_p]
            arrays=[]
            for q in [self.qa,self.qb,self.qc]:
                for k in [0,1]:arrays.append((ctypes.c_char_p*len(q))(*[txt(pair[k]).encode() for pair in q]))
            assert lib.counts_config(len(self.qa),len(self.qb),len(self.qc),NATIVE_BITS,txt(self.l).encode(),txt(self.m).encode(),*arrays)==0
            self.lib=lib

    def __call__(self,h):
        if self.lib:
            p=self.lib.counts_F(txt(h).encode(),int(h<(self.l+self.m)/2))
            if not p:raise ArithmeticError('native count enclosure failed')
            result=arb(ctypes.string_at(p).decode());self.lib.counts_release(p);assert result.is_finite();return result
        return self.python(h)

    def python(self,h):
        y=[p+h for p in self.y];C=kernel(self.x,self.w,y,self.v)
        P=self.A*C;Q=C*self.B;T=C*C.transpose();t2=inner(C,C)
        f=-t2+2*inner(P,C)+2*inner(Q,C)-4*inner(P,Q)-2*inner(T,T)
        lo=(-self.l/2).max(h-self.m/2);hi=(self.l/2).min(h+self.m/2)
        c=hi-lo
        if c>0:
            z,sw=make_points(c,self.qc,(lo+hi)/2)
            AC=kernel(self.x,self.w,z,sw);BC=kernel(y,self.v,z,sw);CC=kernel(z,sw,z,sw)
            f+=c-2*inner(CC,CC)-2*inner(AC,AC)-2*inner(BC,BC)+8*inner(C*BC,AC)
            v=c-t2
        else:v=t2
        return f+2*v*v

def spatial_bound(l,m):
    lens={'A':arb(l),'B':arb(m),'C':arb(min(l,m))};eps={}
    for key,v in lens.items():
        n=order(float(v));rho=arb(n)/arb.pi();q=v/(2*rho)
        eps[key]=2*(arb(2)*n).exp()*q**(2*n)/(1-q)
    out=arb(0)
    for word,c in [('C',1),('AB',1),('CC',4),('AC',2),('BC',2),('ABB',2),('AAB',2),('ABC',12),('AABB',4),('ABAB',4)]:
        vol=arb(1)
        for k in word:vol*=lens[k]
        out+=c*vol*sum((eps[k] for k in word),arb(0))
    return out

def harmonic(H,R,freq,maxp):
    if freq==0:return {p:acb((arb(R)/H).log() if p==1 else (arb(H)**(1-p)-arb(R)**(1-p))/(p-1)) for p in range(1,maxp+1)}
    if freq<0:return {p:z.conjugate() for p,z in harmonic(H,R,-freq,maxp).items()}
    o=arb(freq)*arb.pi();a=o*H;b=o*R
    out={1:acb(b.ci()-a.ci(),b.si()-a.si())};eR=acb(0,b).exp();eH=acb(0,a).exp()
    for p in range(2,maxp+1):out[p]=(acb(0,o)*out[p-1]-(eR/arb(R)**(p-1)-eH/arb(H)**(p-1)))/(p-1)
    return out

def far_integral(l,m,H):
    with ctx.workprec(256):
        C=Counts(l,m,native=False);x,w=C.x,C.w;y,v=C.y,C.v;na=len(x);nb=len(y)
        IA=arb_mat(na,na);IB=arb_mat(nb,nb)
        for i in range(na):IA[i,i]=1
        for i in range(nb):IB[i,i]=1
        JA=acb_mat(IA-2*C.A);JB=acb_mat(IB-2*C.B);E={}
        for r in range(1,J+1):
            mat=acb_mat([[integer_power(-(b-a),r-1)*ww*vv/(2*arb.pi())*acb((b-a).sin_pi(),-(b-a).cos_pi()) for b,vv in zip(y,v)] for a,ww in zip(x,w)])
            E[r,1]=mat;E[r,-1]=mat.conjugate()
        keys=list(E);Z={k:JA*E[k]*JB for k in keys}
        left=acb_mat([[Z[k][i,j] for i in range(na) for j in range(nb)] for k in keys]);right=acb_mat([[E[k][i,j] for k in keys] for i in range(na) for j in range(nb)])
        gram=left*right;coef=defaultdict(acb);ds={}
        for i,(r,a) in enumerate(keys):
            for j,(s,b) in enumerate(keys):
                coef[r+s,a+b]-=gram[i,j];key=(r+s,a+b);mat=E[r,a]*E[s,b].transpose()
                ds[key]=ds[key]+mat if key in ds else mat
        dk=list(ds);tr=[ds[k].trace() for k in dk]
        left=acb_mat([[ds[k][i,j] for i in range(na) for j in range(na)] for k in dk]);right=acb_mat([[ds[k][j,i] for k in dk] for i in range(na) for j in range(na)])
        gram=left*right
        for i,(p,a) in enumerate(dk):
            for j,(q,b) in enumerate(dk):coef[p+q,a+b]+=2*(tr[i]*tr[j]-gram[i,j])
        R=S-(l+m)/2;d1=S-l;d2=S-m
        ints={o:harmonic(H,R,o,4*J) for o in [-4,-2,0,2,4]}
        total=sum((c*(R*ints[o][p]-ints[o][p-1]) for (p,o),c in coef.items()),acb(0))*2/(d1*d2)
        assert total.imag.contains(0) and total.is_finite()
        a=arb(l+m)/2;den=1-a/H;lm=arb(l*m);HH=arb(H)
        e2=(1+2*l)*(1+2*m)*lm/(arb.pi()**2*den**2)*(2*a**J/((J+1)*HH**(J+1))+a**(2*J)/((2*J+1)*HH**(2*J+1)))
        e4=4*lm**2/(arb.pi()**4*den**4)*sum((math.comb(4,r)*a**(J*r)/((3+J*r)*HH**(3+J*r)) for r in range(1,5)),arb(0))
        err=2*(e2+e4)/max(d1,d2)
        coeff=[{'power':p,'harmonic_pi':o,'real':txt(c.real),'imag':txt(c.imag)} for (p,o),c in sorted(coef.items())]
        return total.real,err,coeff

def run_pair(l,m):
    ctx.prec=BITS;t=time.time();H=max(64,math.ceil(4*(l+m)));R=S-(l+m)/2;delta=abs(l-m)/2;d1=S-l;d2=S-m
    C=Counts(l,m);points={float(k)/2 for k in range(2*H+1)}|{delta,(l+m)/2}
    cuts=sorted(x for x in points if 0<=x<=H);qh=[arb.legendre_p_root(NH,k,weight=True) for k in range(NH)]
    path=HERE/'matrix_cases'/('L%s_M%s.json'%(l,m));path.parent.mkdir(exist_ok=True)
    records=[];prefix=arb(0)
    if path.exists():
        old=json.loads(path.read_text())
        assert old['contract']==[l,m,S,H,NH,BITS,NATIVE_BITS,order(l),order(m)]
        if old.get('status')=='complete':return old
        records=old['panels'];prefix=arb(old['prefix'])
    contract=[l,m,S,H,NH,BITS,NATIVE_BITS,order(l),order(m)]
    for idx in range(len(records),len(cuts)-1):
        lo=arb(cuts[idx]);hi=arb(cuts[idx+1]);c=(lo+hi)/2;half=(hi-lo)/2;value=arb(0)
        for x,w in qh:
            h=c+half*x;weight=arb(min(d1,d2)) if h<delta else R-h
            value+=half*w*weight*C(h)
        prefix+=value;records.append({'lo':txt(lo),'hi':txt(hi),'weighted_integral':txt(value)})
        if idx%16==15:
            path.write_text(json.dumps({'status':'running','contract':contract,'prefix':txt(prefix),'panels':records},indent=2))
            print('PAIR',l,m,'through',float(hi),'seconds',round(time.time()-t,1),flush=True)
    near=2*prefix/(d1*d2);far,tail,coeff=far_integral(l,m,H);space=spatial_bound(l,m)
    U=arb(max(l,m)+2);e=arb.pi().exp();MF=U*e+7*U**2*e**2+12*U**3*e**3+6*U**4*e**4+2*(U*e+U**2*e**2)**2
    # Account for shorter panels at containment/touching transitions by using
    # their actual half-width; all are at most 1/4.
    he=arb(0)
    for lo,hi in zip(cuts,cuts[1:]):
        half=arb(hi-lo)/2;q=half
        he+=2/(d1*d2)*4*half*(S+H+1)*MF*q**(2*NH)/(1-q)
    whole=near+far+arb(0,(space+he+tail).upper())
    out={'status':'complete','contract':contract,'L':l,'M':m,'span':S,'H':H,'rho_spatial':'n/pi','near':enc(near),'far':enc(far),
         'spatial_error':enc(space),'separation_error':enc(he),'tail_error':enc(tail),'influence_covariance':enc(whole),
         'panels':records,'harmonic_coefficients':coeff,'seconds':time.time()-t,'native_used':C.lib is not None}
    path.write_text(json.dumps(out,indent=2));print('COMPLETE PAIR',l,m,enc(whole),'seconds',round(time.time()-t,1),flush=True);return out

def assemble():
    rows=[]
    for i,l in enumerate(LS):
        row=[]
        for j,m in enumerate(LS):
            path=HERE/'matrix_cases'/('L%s_M%s.json'%(min(l,m),max(l,m)))
            if not path.exists():return
            d=json.loads(path.read_text())
            if d.get('status')!='complete':return
            row.append(arb(d['influence_covariance']['ball']))
        rows.append(row)
    q=[];means=[]
    for l in LS:
        B,_=bias(S,l);q.append((2*B*B+(arb(l)/(S-l))**2*B).sqrt());x=2*arb.pi()*l
        pop=l-2*l/arb.pi()*x.si()+(1-x.cos()+arb.const_euler()+x.log()-x.ci())/arb.pi()**2
        means.append(enc(pop-B))
    final=[];errors=[]
    for i in range(6):
        a=[];er=[]
        for j in range(6):
            bound=rows[i][i].sqrt()*q[j]+rows[j][j].sqrt()*q[i]+q[i]*q[j]
            a.append(enc(rows[i][j]+arb(0,bound.upper())));er.append(enc(bound))
        final.append(a);errors.append(er)
    out={'status':'complete','scope':'stationary sine projection DPP, fixed-line anchor variances; not a zeta law',
         'span':S,'lengths':LS,'influence_covariance':[[enc(c) for c in row] for row in rows],
         'centering_sd_bounds':[enc(x) for x in q],'centering_covariance_bounds':errors,
         'anchor_covariance':final,'expected_anchor_variances':means,'zeta_sampling_law_certified':False}
    (HERE/'full_matrix_results.json').write_text(json.dumps(out,indent=2));print('ASSEMBLED SIX-LENGTH MATRIX',flush=True)

if __name__=='__main__':
    ctx.prec=BITS
    for i,l in enumerate(LS):
        for m in LS[i:]:run_pair(l,m)
    assemble()
