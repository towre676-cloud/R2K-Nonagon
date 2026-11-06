import sys,json
from pathlib import Path
import numpy as np
P=Path('/mnt/data/nonagon_work/r2k_nonagon_v10.2/external_validation/matpower'); sys.path.insert(0,str(P))
import target_relative_panel as T
from fixtures import CASES
PRIMES=(1000000007,1000000009)

def inv_mod_matrix(A,p):
    A=np.asarray(A,dtype=np.int64)%p; n=A.shape[0]
    aug=np.concatenate([A,np.eye(n,dtype=np.int64)],axis=1)
    for c in range(n):
        piv=None
        for r in range(c,n):
            if int(aug[r,c])%p:
                piv=r; break
        if piv is None: raise ValueError('singular mod p')
        if piv!=c: aug[[c,piv]]=aug[[piv,c]]
        inv=pow(int(aug[c,c]),p-2,p); aug[c]=(aug[c]*inv)%p
        for r in range(n):
            if r==c: continue
            f=int(aug[r,c])
            if f: aug[r]=(aug[r]-f*aug[c])%p
    return aug[:,n:]

def det_mod(A,p):
    A=np.asarray(A,dtype=np.int64)%p; n=A.shape[0]; det=1
    for c in range(n):
        piv=None
        for r in range(c,n):
            if int(A[r,c])%p:
                piv=r; break
        if piv is None: return 0
        if piv!=c:
            A[[c,piv]]=A[[piv,c]]; det=(-det)%p
        pv=int(A[c,c]); det=(det*pv)%p; inv=pow(pv,p-2,p)
        # Normalize pivot row to keep all products < p^2.
        A[c,c:]=(A[c,c:]*inv)%p
        for r in range(c+1,n):
            f=int(A[r,c])
            if f: A[r,c:]=(A[r,c:]-f*A[c,c:])%p
    return int(det%p)

def small_inv(M,p):
    if M.shape==(1,1):
        a=int(M[0,0])%p
        return np.array([[pow(a,p-2,p)]],dtype=np.int64)
    a,b,c,d=map(lambda x:int(x)%p,[M[0,0],M[0,1],M[1,0],M[1,1]])
    det=(a*d-b*c)%p; invd=pow(det,p-2,p)
    return np.array([[d*invd%p,(-b)*invd%p],[(-c)*invd%p,a*invd%p]],dtype=np.int64)

def certificate_for_prime(Bi,Kunit,G,pats,saved,p):
    g0=G[:,0].copy(); inv2=pow(2,p-2,p); rows=[]
    for ii in saved:
        typ,S,loads=pats[ii]
        if not S:
            z=g0.copy()
        else:
            S=np.asarray(S,dtype=int)
            M=G[np.ix_(S,S)].copy()%p
            for a,load in enumerate(loads):
                dinv=1 if float(load)==1.0 else inv2
                M[a,a]=(int(M[a,a])+dinv)%p
            Mi=small_inv(M,p)
            coeff=(Mi @ (g0[S]%p))%p
            z=(g0 - (G[:,S]@coeff)%p)%p
        q=(Bi.T@z)%p
        rows.append((-q*q)%p)
    return det_mod(np.asarray(rows,dtype=np.int64),p)


def main():
    buses,K,B,w,keys=T.edge_parameterization(CASES['case118'])
    saved=json.load(open(P/'target_relative_results.json'))['ieee118_adaptive']['local_kernel_panel']['selected_indices']
    Bi=np.rint(B).astype(np.int64); Kunit=Bi@Bi.T; pats=T.patterns(len(buses),1.0)
    certs=[]
    for p in PRIMES:
        G=inv_mod_matrix(Kunit,p)
        d=certificate_for_prime(Bi,Kunit,G,pats,saved,p)
        certs.append({'prime':p,'determinant_mod_prime':d,'nonzero':bool(d)})
    out={
      'schema':'nonagon.target-relative-exact-rank.v1','case':'case118',
      'selected_count':len(saved),'parameter_count':len(keys),'certificates':certs,
      'witness_point':'all 179 aggregate weights = 1; single shunts 1 and 2; pair shunts 1,1',
      'meaning':'a nonzero determinant modulo either prime proves the selected 179 response rows have exact full Jacobian rank over Q at this rational witness point; two primes are retained as an implementation cross-check',
      'does_not_prove':'rank at every positive weight/load choice or global injectivity of the nonlinear response map'
    }
    (P/'target_relative_exact_rank.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))

if __name__=='__main__': main()
