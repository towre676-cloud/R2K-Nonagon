"""Exact first-order identifiability certificates for Part XIV.

Finite-field rank is a one-sided rational rank lower certificate.  A deficient
finite-field nullspace does not certify rational nullspace support.  This module
computes exact rational kernel and row-span witnesses and checks every reading.
It makes no assertion of nonlinear/global identifiability from deficient rank.
"""
from __future__ import annotations
import argparse, hashlib, json, math, sys, time
from collections import deque
from pathlib import Path
import flint
import numpy as np

if not __debug__:
    raise RuntimeError("Exact certificate verification requires assertions; run without python -O")

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT / "part_XI"))
import gridlib
P = 2147483647


def load_exact(case):
    rec = json.loads((ROOT / "part_XIII" / ("record_" + case + ".json")).read_text())
    nodes0, edges0, ground0, meta = gridlib.network(ROOT / "part_XI" / "data" / (case + ".m"))
    adj0 = {i: set() for i in range(len(nodes0))}
    for a, b in edges0:
        adj0[a].add(b); adj0[b].add(a)
    start = nodes0.index(rec["port_bus"]); keep = {start}; todo = [start]
    for a in todo:
        for b in adj0[a]:
            if b not in keep: keep.add(b); todo.append(b)
    keep = sorted(keep); ix = {v: k for k, v in enumerate(keep)}
    nodes = [nodes0[i] for i in keep]; n = len(nodes)
    edges = {(ix[a], ix[b]): w for (a, b), w in edges0.items() if a in ix and b in ix}
    ground = [ground0[i] for i in keep]; el = sorted(edges); gs = [i for i,g in enumerate(ground) if g]
    K = flint.fmpq_mat(n,n)
    for (a,b),w in edges.items():
        q=flint.fmpq(w.numerator,w.denominator)
        K[a,a]+=q;K[b,b]+=q;K[a,b]-=q;K[b,a]-=q
    for i,g in enumerate(ground):
        if g: K[i,i]+=flint.fmpq(g.numerator,g.denominator)
    G=K.inv(); M,d=G.numer_denom()
    assert K*G==flint.fmpq_mat([[int(i==j) for j in range(n)] for i in range(n)])
    adj={i:set() for i in range(n)}
    for a,b in el: adj[a].add(b);adj[b].add(a)
    buses,_=gridlib.parse(ROOT / "part_XI" / "data" / (case+".m"));idx={b:i for i,b in enumerate(nodes)}
    hw=sorted(idx[int(float(b[0]))] for b in buses if int(float(b[0])) in idx and (float(b[5])!=0 or int(float(b[1]))==2))
    return dict(case=case,n=n,nodes=nodes,el=el,gs=gs,port=nodes.index(rec["port_bus"]),adj=adj,hw=hw,
                M=[[int(x) for x in row] for row in M.tolist()],d=int(d),U=len(el)+len(gs))


def farthest_ports(D,k):
    ports=[D["port"]]; dmin=[10**9]*D["n"]
    while len(ports)<k:
        dist=[10**9]*D["n"];dist[ports[-1]]=0;q=deque([ports[-1]])
        while q:
            a=q.popleft()
            for b in D["adj"][a]:
                if dist[b]>dist[a]+1:dist[b]=dist[a]+1;q.append(b)
        dmin=[min(a,b) for a,b in zip(dmin,dist)]
        ports.append(max(range(D["n"]),key=lambda i:dmin[i]))
    return ports


def primitive(values):
    g=0
    for a in values:g=math.gcd(g,abs(a))
    return [a//g for a in values] if g>1 else values


def exact_columns(D,ports,j):
    M,d=D["M"],D["d"]
    if j is None:return [primitive([M[i][p] for i in range(D["n"])]) for p in ports]
    h=d+M[j][j]
    if not h:raise ValueError("singular intervened model")
    return [primitive([M[i][p]*h-M[i][j]*M[j][p] for i in range(D["n"])]) for p in ports]


def exact_rows(D,ports,where,driving_only=False):
    for j in [None]+list(where):
        X=exact_columns(D,ports,j)
        pairs=[(a,b) for a in range(len(ports)) for b in range(a,len(ports))] if j is None or not driving_only else [(a,a) for a in range(len(ports))]
        for a,b in pairs:
            x,y=X[a],X[b]
            row=[(x[i]-x[k])*(y[i]-y[k]) for i,k in D["el"]]+[x[i]*y[i] for i in D["gs"]]
            yield primitive(row)


def modular_rows(D,ports,where,driving_only=False):
    n=D["n"];di=pow(D["d"]%P,-1,P)
    G=np.array([[v%P*di%P for v in row] for row in D["M"]],dtype=np.int64)
    ea=np.array([a for a,b in D["el"]],dtype=np.int64);eb=np.array([b for a,b in D["el"]],dtype=np.int64)
    rows=[]
    for j in [None]+list(where):
        X=G[:,ports]
        if j is not None:
            f=pow(int(1+G[j,j])%P,-1,P)
            X=(X-(G[:,j,None]*((G[j,ports]*f)%P)[None,:])%P)%P
        pairs=[(a,b) for a in range(len(ports)) for b in range(a,len(ports))] if j is None or not driving_only else [(a,a) for a in range(len(ports))]
        pa=np.array([a for a,b in pairs]);pb=np.array([b for a,b in pairs])
        dp=(X[ea][:,pa]-X[eb][:,pa])%P;dq=(X[ea][:,pb]-X[eb][:,pb])%P
        rows.extend(np.vstack([(dp*dq)%P,(X[D["gs"]][:,pa]*X[D["gs"]][:,pb])%P]).T.tolist())
    return rows


def hash_matrix(rows):
    h=hashlib.sha256()
    for row in rows:h.update((",".join(map(str,row))+"\n").encode())
    return h.hexdigest()


def rational_kernel(A):
    R,r=A.rref();U=A.ncols();piv=[]
    for i in range(r):piv.append(next(j for j in range(U) if R[i,j]))
    free=[j for j in range(U) if j not in piv];cols=[]
    for f in free:
        v=[flint.fmpq(0) for _ in range(U)];v[f]=flint.fmpq(1)
        for i,c in enumerate(piv):v[c]=-R[i,f]
        den=1
        for x in v:den=math.lcm(den,int(x.denominator))
        cols.append(primitive([int(x*den) for x in v]))
    return R,r,piv,free,cols


def certify(D,ports,where,driving_only=False,include_witnesses=True):
    t=time.time(); U=D["U"]
    Jmod=flint.nmod_mat(modular_rows(D,ports,where,driving_only),P)
    Rt,r=Jmod.transpose().rref()
    print("modular rank",r,"of",U,flush=True)
    rowids=[next(j for j in range(Rt.ncols()) if Rt[i,j]) for i in range(r)]
    chosen=set(rowids);rows=[None]*r;pos={v:i for i,v in enumerate(rowids)};m=0
    for i,row in enumerate(exact_rows(D,ports,where,driving_only)):
        if i in chosen:rows[pos[i]]=row
        m+=1
    A=flint.fmpq_mat(rows)
    print("selected integer rows ready",len(rows),"max bits",max(abs(x).bit_length() for row in rows for x in row),flush=True)
    if r==U:
        R=None;rr=r;piv=list(range(U));free=[];Ncols=[]
    else:
        R,rr,piv,free,Ncols=rational_kernel(A)
    print("rational kernel ready",len(Ncols),flush=True)
    if rr!=r:raise ValueError("modular rank below rational subset rank; change prime")
    # The selected square minor is nonzero modulo p after exact row clearing.
    minor=flint.nmod_mat([[rows[i][j]%P for j in piv] for i in range(r)],P)
    minor_det=int(minor.det());assert minor_det
    sparseN=[[(i,a) for i,a in enumerate(v) if a] for v in Ncols]
    checks=0;Jhash=hashlib.sha256()
    for row in exact_rows(D,ports,where,driving_only):
        Jhash.update((",".join(map(str,row))+"\n").encode())
        for v in sparseN:
            if sum(row[i]*a for i,a in v)!=0:raise ValueError("rational nullspace failed on omitted reading")
            checks+=1
    print("all exact kernel products verified",checks,flush=True)
    # N free-coordinate block is diagonal and nonzero, proving U-r independent columns.
    assert all(Ncols[k][f]!=0 and all(Ncols[l][f]==0 for l in range(len(free)) if l!=k) for k,f in enumerate(free))
    ident=[i for i in range(U) if all(v[i]==0 for v in Ncols)]
    witnesses=[]
    if include_witnesses and ident and r<U:
        B=flint.fmpq_mat([[rows[i][j] for j in piv] for i in range(r)])
        E=flint.fmpq_mat(r,len(ident))
        for col,i in enumerate(ident):E[piv.index(i),col]=1
        L=B.transpose().solve(E).transpose();T=L*A
        for k,i in enumerate(ident):
            lam=[L[k,j] for j in range(r)]
            assert all(T[k,j]==int(i==j) for j in range(U))
            witnesses.append(dict(coordinate=i,coefficients=[str(x) for x in lam]))
    result=dict(case=D["case"],ports=[D["nodes"][p] for p in ports],intervention_buses=[D["nodes"][j] for j in where],
                driving_point_only_after_intervention=driving_only,readings=m,unknowns=U,
                rational_rank=r,rational_nullity=U-r,first_order_identifiable_unknowns=len(ident),
                first_order_identifiable_branches=sum(i<len(D["el"]) for i in ident),branches=len(D["el"]),
                identifiable_coordinate_indices=ident,rank_witness_prime=P,rank_witness_rows=rowids,
                rank_witness_columns=piv,rank_witness_minor_determinant_mod_p=minor_det,
                kernel_free_coordinates=free,kernel_generators=Ncols,exact_kernel_dot_products_checked=checks,
                all_integer_jacobian_rows_sha256=Jhash.hexdigest(),selected_integer_rows_sha256=hash_matrix(rows),
                selected_rows_reconstructed_from_exact_input=True,coordinate_row_span_witnesses=witnesses,
                coordinate_classification_via_complete_kernel=not include_witnesses,
                full_rank_proves_all_coordinates_without_individual_row_span_witnesses=r==U,
                certificate_scope="exact linearized reading-map identifiability at the declared rational model; no nonlinear/global claim for deficient rank",
                seconds=time.time()-t)
    return result


def modular_counterexample():
    # Same rank over Q and F_p does not rescue individual coordinate support.
    J=flint.fmpq_mat([[1,P]]);v=flint.fmpq_mat([[-P],[1]])
    assert J*v==flint.fmpq_mat([[0]])
    assert flint.nmod_mat([[1,P]],P).rank()==J.rank()==1
    return {"prime":P,"J":[[1,P]],"rational_kernel":[-P,1],"rational_first_order_identifiable_coordinates":[],"mod_p_apparently_identifiable_coordinates":[0]}


def main():
    ap=argparse.ArgumentParser();ap.add_argument("case");ap.add_argument("schemes",help="e.g. kron20,hw20,ports20,access7")
    ap.add_argument("--output-dir",type=Path,default=HERE/"rational_access")
    ap.add_argument("--no-row-span-witnesses",action="store_true",help="use verified complete-kernel duality for coordinate classification")
    args=ap.parse_args();args.output_dir.mkdir(exist_ok=True);t=time.time();D=load_exact(args.case)
    print("exact inverse ready",D["n"],"nodes",time.time()-t,flush=True)
    for scheme in args.schemes.split(","):
        prefix="".join(c for c in scheme if c.isalpha());k=int("".join(c for c in scheme if c.isdigit()))
        if prefix=="pocket":
            if args.case!="case118" or k!=6:raise ValueError("pocket6 is the declared IEEE118 design")
            ports=farthest_ports(D,5)+[D["nodes"].index(55)]
        elif prefix=="access":
            old=json.loads((HERE/("portloop14_"+args.case+".json")).read_text());ports=[D["nodes"].index(b) for b in old["port_buses"][:k]]
        else:ports=farthest_ports(D,k)
        where=[] if prefix=="kron" else (D["hw"] if prefix=="hw" else list(range(D["n"])))
        out=certify(D,ports,where,driving_only=prefix in ("access","pocket"),include_witnesses=not args.no_row_span_witnesses)
        out["scheme"]=scheme;out["modular_counterexample"]=modular_counterexample()
        fn=args.output_dir/("rational_"+args.case+"_"+scheme+".json")
        fn.write_text(json.dumps(out,separators=(",",":")))
        print(json.dumps({k:v for k,v in out.items() if k in ["case","scheme","readings","rational_rank","rational_nullity","first_order_identifiable_branches","branches","seconds"]}),flush=True)

if __name__=="__main__":main()
