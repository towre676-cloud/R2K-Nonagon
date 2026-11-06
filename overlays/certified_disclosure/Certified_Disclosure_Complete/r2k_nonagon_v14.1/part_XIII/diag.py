import sys, json, subprocess
sys.path[:0]=['../p11']
import gridlib, numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spl
from fractions import Fraction as Fr
case=sys.argv[1]; p=int(sys.argv[2]) if len(sys.argv)>2 else 2147483647
nodes,edges,ground,meta=gridlib.network('../p11/data/%s.m'%case); n=len(nodes); idx={b:k for k,b in enumerate(nodes)}
rec=json.load(open('record_%s.json'%case)); port=idx[rec['port_bus']]; extra=[idx[b] for b in rec['extra_source_buses']]
rows,cols,vals=[],[],[]
for (a,c),w in edges.items():
    w=float(w); rows+=[a,c,a,c]; cols+=[a,c,c,a]; vals+=[w,w,-w,-w]
for i,g in enumerate(ground):
    if g: rows.append(i); cols.append(i); vals.append(float(g))
A=sp.csc_matrix((vals,(rows,cols)),shape=(n,n))
perm=spl.splu(A,permc_spec='MMD_AT_PLUS_A').perm_c.tolist(); iperm={v:k for k,v in enumerate(perm)}
E=sorted(edges,key=lambda e:min(iperm[e[0]],iperm[e[1]])); Gs=[i for i in range(n) if ground[i]!=0]
def res(x): x=Fr(x); return x.numerator%p*pow(x.denominator%p,-1,p)%p
L=['%d %d %d %d %d'%(p,n,port,len(E),len(Gs)),' '.join(map(str,perm))]
L+=['%d %d %d'%(a,c,res(edges[(a,c)])) for a,c in E]+['%d %d'%(i,res(ground[i])) for i in Gs]+[str(len(extra)),' '.join(map(str,extra)),str(len(extra)),' '.join(str(rec.get('widened_sources',{}).get(str(nodes[x]),1)) for x in extra)]
r=subprocess.run(['./tomo_k'],input='\n'.join(L)+'\n',capture_output=True,text=True)
free=[int(x.split()[1]) for x in r.stdout.split('\n') if x.startswith('free')]
deg=np.zeros(n,int)
for a,c in edges: deg[a]+=1; deg[c]+=1
import collections
for line in r.stdout.split('\n'):
    if line.startswith('kern'):
        head,cnt=line.split('|'); ids=[int(t) for t in head.split()[2:]]
        print('kernel vector support size',cnt.strip())
        for f in ids[:40]:
            if f<len(E): a,c=E[f]; print('  edge',nodes[a],nodes[c],round(float(edges[E[f]]),4))
            else: i=Gs[f-len(E)]; print('  ground',nodes[i],round(float(ground[i]),4))
for f in []:
    if f<len(E):
        a,c=E[f]; print('edge',nodes[a],nodes[c],'w',float(edges[E[f]]),'deg',deg[a],deg[c],'grounds',float(ground[a]),float(ground[c]))
    else:
        i=Gs[f-len(E)]; print('ground',nodes[i],float(ground[i]),'deg',deg[i])
