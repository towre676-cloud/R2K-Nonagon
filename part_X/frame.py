"""Executable algebra of the ported elimination calculus.
Exact rational checks, random binary64 replay and 80-digit control.
Run: python frame.py. Dependencies: numpy, scipy, sympy, mpmath.
"""
from pathlib import Path
import json, itertools, time
import numpy as np
import sympy as s
import mpmath as mp

def network(n,rng):
    W=np.zeros((n,n))
    for j in range(1,n):
        i=int(rng.integers(j)); W[i,j]=W[j,i]=rng.uniform(.25,2)
    for i in range(n):
        for j in range(i+1,n):
            if rng.random()<.3:W[i,j]=W[j,i]=rng.uniform(.25,2)
    return np.diag(W.sum(1)+rng.uniform(.2,1.4,n))-W

def rational_frame():
    K=s.Matrix([[5,-2,-1,0],[-2,6,-1,-1],[-1,-1,5,-1],[0,-1,-1,4]])
    G=K.inv(); b=s.eye(4)[:,0]; U=s.eye(4)[:,1:]
    a=U.T*G*b; S=U.T*G*U; h=(b.T*G*b)[0]
    z=s.symbols('z1:4'); D=s.diag(0,*z)
    tau=s.Poly((K+D).det(),*z)
    rho=s.Poly((K+D)[1:,1:].det(),*z)
    assert all(all(v<=1 for v in m) for m,c in tau.terms()+rho.terms())
    response=rho.as_expr()/tau.as_expr()
    zero=dict.fromkeys(z,0)
    grad=[s.diff(response,x).subs(zero) for x in z]
    Hess=s.Matrix([[s.diff(response,x,y).subs(zero) for y in z] for x in z])
    assert grad==[-a[i]**2 for i in range(3)]
    assert Hess==2*s.diag(*a)*S*s.diag(*a)
    # The observed response and its two-jet recover all inverse entries.
    Grec=s.zeros(4);Grec[0,0]=h
    for i in range(3):Grec[0,i+1]=Grec[i+1,0]=s.sqrt(-grad[i])
    for i in range(3):
        for j in range(3):Grec[i+1,j+1]=Hess[i,j]/(2*Grec[0,i+1]*Grec[0,j+1])
    assert Grec==G
    strengths=[s.Rational(1,2),s.Rational(3,2),s.Rational(2)]
    Theta=s.diag(*strengths)
    hf=h-(a.T*(Theta.inv()+S).inv()*a)[0]
    assert hf==(b.T*(K+U*Theta*U.T).inv()*b)[0]
    # Eliminate site 3, retaining port and two controllable sites.
    A=K[:3,:3];B=K[:3,3:];C=K[3:,3:]
    Q=A-B*C.inv()*B.T
    zz=s.diag(0,z[0],z[1])
    assert s.cancel((Q+zz).inv()[0,0]-(K+s.diag(0,z[0],z[1],0)).inv()[0,0])==0
    # A nontrivial congruence transforms the directions along with the matrix.
    T=s.Matrix([[1,0,0,0],[0,1,1,0],[0,0,1,0],[0,0,0,2]])
    assert (T.T*b).T*(T.T*K*T).inv()*(T.T*b)==s.Matrix([[h]])
    # Three Schur steps exactly give a scalar continued fraction for a path.
    chain=s.Matrix([[4,-1,0,0],[-1,5,-2,0],[0,-2,6,-1],[0,0,-1,3]])
    cf=1/(s.Rational(4)-1/(s.Rational(5)-4/(s.Rational(6)-1/s.Rational(3))))
    assert cf==chain.inv()[0,0]
    # Stieltjes-Jacobi fractions at complex arguments; polynomial identity.
    x=s.symbols('x');diag=[4,5,6,3];off=[1,2,1]
    den=x+diag[-1]
    for i in range(2,-1,-1):den=x+diag[i]-s.Rational(off[i]**2)/den
    assert s.cancel(1/den-(chain+x*s.eye(4)).inv()[0,0])==0
    vals=[float(v) for v in np.linalg.eigvalsh(np.array(chain,dtype=float))]
    # Jacobi eigenvalues interlace after deleting the first coordinate.
    poles=np.linalg.eigvalsh(np.array(chain,dtype=float));zeros=np.linalg.eigvalsh(np.array(chain[1:,1:],dtype=float))
    assert all(poles[i]<zeros[i]<poles[i+1] for i in range(3))
    return dict(K=str(K),G=str(G),h=str(h),a=[str(x) for x in a],tau=str(tau.as_expr()),rho=str(rho.as_expr()),
                exact_two_jet=True,exact_elimination=True,exact_congruence=True,
                chain_response=str(cf),chain_poles=poles.tolist(),chain_zeros=zeros.tolist())

def random_replay():
    rng=np.random.default_rng(20261030)
    mx=dict(frame=0.,jet_recovery=0.,schur=0.,principal_minor=0.,positive_real=0.)
    cases=0; bounds=[]
    for n in range(2,13):
        for rep in range(24):
            K=network(n,rng);G=np.linalg.inv(K);g=G[:,0];p=n-1
            t=rng.uniform(.05,3,p);S=G[1:,1:];a=g[1:]
            h=G[0,0]-a@np.linalg.solve(np.diag(1/t)+S,a)
            Hd=np.linalg.inv(K+np.diag(np.r_[0,t]))[0,0]
            mx['frame']=max(mx['frame'],abs(h-Hd)/abs(Hd))
            gradient=-a*a; Hess=2*a[:,None]*S*a[None,:]
            recovered=Hess/(2*np.sqrt(-gradient[:,None])*np.sqrt(-gradient[None,:]))
            mx['jet_recovery']=max(mx['jet_recovery'],np.max(np.abs(recovered-S))/np.max(abs(S)))
            keep=list(range(max(1,n//2)));elim=list(range(max(1,n//2),n))
            A=K[np.ix_(keep,keep)]; B=K[np.ix_(keep,elim)];C=K[np.ix_(elim,elim)]
            Q=A-B@np.linalg.solve(C,B.T)
            mx['schur']=max(mx['schur'],np.max(abs(np.linalg.inv(Q)-G[np.ix_(keep,keep)]))/np.max(abs(G)))
            # Principal minor generating polynomial, all subsets at n<=7.
            if n<=7:
                ratio=0.
                for r in range(p+1):
                    for I in itertools.combinations(range(p),r):
                        ratio+=(np.linalg.det(S[np.ix_(I,I)]) if I else 1)*np.prod(t[list(I)])
                direct=np.linalg.det(K+np.diag(np.r_[0,t]))/np.linalg.det(K)
                mx['principal_minor']=max(mx['principal_minor'],abs(ratio-direct)/abs(direct))
            # Passive mobility: energy identity for full second-order pencil.
            C=network(n,rng);M=np.diag(rng.uniform(.5,2,n));z=.3+1.7j
            A=K+z*C+z*z*M;x=np.linalg.solve(A,np.eye(n)[:,0]);Y=z*x[0]
            # Re Y = Re(s)||x||_K² + |s|²||x||_C² + |s|² Re(s)||x||_M².
            rhs=z.real*np.vdot(x,K@x).real+abs(z)**2*np.vdot(x,C@x).real+abs(z)**2*z.real*np.vdot(x,M@x).real
            mx['positive_real']=max(mx['positive_real'],abs(Y.real-rhs)/abs(Y.real))
            assert Y.real>0
            cases+=1
    assert max(mx.values())<1e-11
    return dict(cases=cases,orders=[2,12],max_relative_errors=mx,seed=20261030)

def hp_control():
    mp.mp.dps=80
    K=mp.matrix([[5,-2,-1,0],[-2,6,-1,-1],[-1,-1,5,-1],[0,-1,-1,4]])
    G=K**-1;T=mp.diag([mp.mpf('.5'),mp.mpf('1.5'),mp.mpf('2')]);a=mp.matrix([G[i,0] for i in range(1,4)]);S=mp.matrix([[G[i,j] for j in range(1,4)] for i in range(1,4)])
    h=G[0,0]-(a.T*(T**-1+S)**-1*a)[0]
    D=mp.diag([0,mp.mpf('.5'),mp.mpf('1.5'),mp.mpf('2')]);direct=(K+D)**-1
    err=abs(h-direct[0,0]);assert err<mp.mpf('1e-75')
    return dict(decimal_digits=80,frame_absolute_error=mp.nstr(err,12),response=mp.nstr(h,80))

def run():
    out=dict(schema='nonagon.part_x.frame.v1',exact=rational_frame(),binary64=random_replay(),high_precision=hp_control(),verdict='PASS')
    Path(__file__).with_name('frame_results.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))
if __name__=='__main__':run()
