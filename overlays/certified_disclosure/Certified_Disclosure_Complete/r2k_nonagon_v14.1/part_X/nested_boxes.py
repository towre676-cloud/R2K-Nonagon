#!/usr/bin/env python3
"""Glue adjacent true contour boxes and localize their moment form.
Requires experiments/union_7005_contour.json for an independently integrated union.
"""
import json,math,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'.optional-local-deps'))
sys.path.insert(0,str(Path(__file__).resolve().parent))
from flint import arb,arb_mat,ctx
from certification import moments,power_y,nested,sub,sign,r2kverify
ctx.prec=192
V9=ROOT if (ROOT/'part_VI').is_dir() else ROOT/'source/v9/r2k_nonagon_v9.0'

def affine(m,a,b,K):
 return [sum((math.comb(k,j)*a**j*b**(k-j)*m[j] for j in range(k+1)),arb(0)) for k in range(K+1)]
def root_poly(m,n):
 e=[arb(1)]
 for k in range(1,n+1):e.append(sum(((-1)**(j-1)*e[k-j]*m[j] for j in range(1,k+1)),arb(0))/k)
 return [(-1)**(n-k)*e[n-k] for k in range(n+1)]
def main():
 children=[json.load(open(V9/'part_VI/certs'/('III_scan7005_%02d.json'%i))) for i in range(2)]
 child_verdicts=[r2kverify.verify(c,full=True,evaluator='series') for c in children]
 assert all(v['status']=='CERTIFIED' and v['count']==n for v,n in zip(child_verdicts,[3,2]))
 union=json.load(open(Path(__file__).resolve().parent/'union_7005_contour.json'))
 full_verdict=r2kverify.verify(union,full=True,evaluator='series')
 assert full_verdict['status']=='CERTIFIED' and full_verdict['count']==5
 assert full_verdict['segments_recomputed']==147
 ns=[3,2];N=sum(ns);K=2*N-2
 lo=arb(float.fromhex(union['rectangle']['lo']));hi=arb(float.fromhex(union['rectangle']['hi']))
 cm=[]
 for c,n in zip(children,ns):
  mm=moments(c);yy=power_y(mm,K);yy[0]=arb(n)
  cl=arb(float.fromhex(c['rectangle']['lo']));ch=arb(float.fromhex(c['rectangle']['hi']))
  cm.append(affine(yy,(ch-cl)/(hi-lo),(cl-lo)/(hi-lo),K))
 mu=[a+b for a,b in zip(*cm)];flat=power_y(moments(union),K)
 assert all(a.overlaps(b) for a,b in zip(mu,flat))
 H=arb_mat([[mu[i+j] for j in range(N)] for i in range(N)])
 A,B=[root_poly(m,n) for m,n in zip(cm,ns)]
 cols=[]
 for p,n in [(B,ns[0]),(A,ns[1])]:
  for k in range(n):cols.append([arb(0)]*k+p+[arb(0)]*(N-k-len(p)))
 V=arb_mat([[cols[j][i] for j in range(N)] for i in range(N)])
 assert sign(V.det())
 C=V.transpose()*H*V
 cross=[C[i,j] for i in range(ns[0]) for j in range(ns[0],N)]
 assert all(v.contains(0) for v in cross)
 ai,_=nested(sub(C,range(ns[0]),range(ns[0])))
 bi,_=nested(sub(C,range(ns[0],N),range(ns[0],N)))
 fi,_=nested(H)
 assert fi==[ai[0]+bi[0],ai[1]+bi[1]]
 r={'independent_full_child_verifiers':child_verdicts,'independent_full_union_verifier':full_verdict,'child_boxes':[[6990,6992],[6992,6994]],'child_counts':ns,'parent_count':N,'moment_order':K,
 'direct_parent_moments_overlap_child_sum':True,'basis_determinant_ball':V.det().str(30,radius=True),
 'transformed_cross_blocks_contain_zero':True,'maximum_cross_block_abs_upper':max(float(v.abs_upper()) for v in cross),
 'child_weighted_inertias':[ai,bi],'nested_parent_inertia':fi,
 'all_parent_moment_balls':[v.str(25,radius=True) for v in mu],
 'all_direct_union_moment_balls':[v.str(25,radius=True) for v in flat],
 'proof':'Exact child polynomials from Newton identities vanish on their child roots. The basis [p_B y^k; p_A y^l] has exact zero cross blocks and nonzero determinant; the interval computation encloses its exact coefficients.'}
 (Path(__file__).resolve().parent/'nested_boxes_results.json').write_text(json.dumps(r,indent=2));print(json.dumps(r,indent=2))
if __name__=='__main__':main()
