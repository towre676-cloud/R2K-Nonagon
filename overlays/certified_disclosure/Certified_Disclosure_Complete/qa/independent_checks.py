"""Independent exact/interval controls for the disclosure extension.

This script does not use the proposal optimizers, generated ledgers, or replay
receipts as truth. It checks externally specified rational atomic measures,
independently built dense killed generators, and interval root boxes.
"""
from pathlib import Path
from fractions import Fraction as F
import json, random, sys, time

ROOT = Path(__file__).resolve().parents[1]
WORK = ROOT/'work' if (ROOT/'work/possq-0.4.0').is_dir() else ROOT
sys.path.insert(0,str(WORK/'possq-0.4.0'))
sys.path.insert(0,str(WORK/'r2k_nonagon_v14.1/part_XIV'))

from flint import arb,arb_mat,arb_poly,ctx,fmpq
from possq.resolvent import ResolventData,certify_resolvent_candidate,verify_resolvent_certificate
from possq.cyclecraft_frontier import uniformized_stage,A,convolve,wrong_direction,compose_kernel
from possq.cyclecraft_exact import lifted_generator
from possq.cyclecraft import TAU,N
from two_strength_cert import Interval,pair_delta,certify_pair


def spectral_controls():
 rng=random.Random(91822);checks=0
 for case in range(40):
  xs=[F(rng.randint(11,198),100) for _ in range(7)]
  nums=[rng.randint(1,20) for _ in xs];total=sum(nums);ws=[F(x,total) for x in nums]
  shifts=[F(1,10),F(1,3),F(1),F(3)]
  truth=[sum((w/(x+s) for x,w in zip(xs,ws)),F(0)) for s in shifts]
  eps=[F(1,10**5)]*4
  centers=[v+F(rng.randint(-3,3),4)*e for v,e in zip(truth,eps)]
  d=ResolventData((F(1,10),F(2)),shifts,centers,eps)
  S=[(F(1,5),F(3,5)),(F(11,10),F(7,5))]
  mass=sum((w for x,w in zip(xs,ws) if any(a<=x<=b for a,b in S)),F(0))
  for side in ['upper','lower']:
   proposal=[F(rng.randint(-20,20),7) for _ in range(5)]
   cert=certify_resolvent_candidate(d,S,proposal,side)
   assert verify_resolvent_certificate(d,cert)
   assert (mass<=cert.value) if side=='upper' else (cert.value<=mass)
   checks+=1
 return {'random_atomic_measures':40,'upper_and_lower_certificates':checks,'all_exact_mass_inequalities':True}


def uniformization_controls():
 ctx.prec=256; checks=0; max_tail=arb(0)
 for graph in ['cycle','star']:
  for stage in range(3):
   lo,hi=-3,4;width=hi-lo+1
   K=lifted_generator(graph,stage,0,lo,hi)
   E=(arb_mat([[A(v) for v in row] for row in K])*A(TAU)).exp()
   for start in [0,2,4]:
    vector=[arb_poly([0]*(-lo)+[int(i==start)]) for i in range(N)]
    out,j,tail=uniformized_stage(vector,graph,stage,width,arb('1e-45'))
    for level in range(width):
     for state in range(N):
      dense=E[level*N+state,(-lo)*N+start]
      sparse=out[state][level]
      # Full held-stage kernel is sparse retained sum plus the omitted
      # positive Poisson tail, entrywise at most the scalar tail.
      assert sparse.lower()<=dense.upper()
      assert dense.lower()<=(sparse+tail).upper()
      checks+=1
    max_tail=max(max_tail,tail)
 return {'graphs':2,'stages':3,'start_states_per_stage':3,'dense_vs_sparse_entry_checks':checks,'max_poisson_tail':max_tail.str(20)}


def independent_convolution_controls():
 ctx.prec=256
 rng=random.Random(6105); checks=0
 for _ in range(25):
  xs=list(range(-2,4));ws=[rng.randint(1,9) for x in xs];tot=sum(ws);w=[F(v,tot) for v in ws]
  p={x:(A(y),A(y)) for x,y in zip(xs,w)}
  acc=p
  exact=dict(zip(xs,w))
  for copies in range(2,7):
   acc=convolve(acc,p)
   exact2={}
   for x,y in exact.items():
    for z,u in zip(xs,w):exact2[x+z]=exact2.get(x+z,F(0))+y*u
   exact=exact2
   exact_event=sum((q for x,q in exact.items() if x<=0),F(0))
   lo,hi=wrong_direction(acc,arb(0))
   assert lo<=A(exact_event)<=hi
   checks+=1
 return {'rational_convolution_tail_checks':checks,'all_exact_probabilities_enclosed':True}


def kernel_composition_controls():
 ctx.prec=256
 rng=random.Random(891);checks=0
 def random_kernel(shifts):
  data={k:[[F(0) for j in range(N)] for i in range(N)] for k in shifts}
  for j in range(N):
   weights=[rng.randint(1,9) for _ in range(len(shifts)*N)];total=sum(weights)
   for x,k in enumerate(shifts):
    for i in range(N):data[k][i][j]=F(weights[x*N+i],total)
  return data
 for _ in range(12):
  left,right=random_kernel([-1,1]),random_kernel([0,1])
  first=[[arb_poly([A(left.get(k,[[F(0)]*N]*N)[i][j]) for k in range(-1,2)]) for j in range(N)] for i in range(N)]
  second=[[arb_poly([A(right.get(k,[[F(0)]*N]*N)[i][j]) for k in range(-1,2)]) for j in range(N)] for i in range(N)]
  product=compose_kernel(first,second,-1,-2,2)
  for k in range(-2,3):
   for i in range(N):
    for j in range(N):
     exact=sum((right[b][i][a]*left[c][a][j] for b in right for c in left if b+c==k for a in range(N)),F(0))
     ball=product[i][j][k+2]
     ml,el=ball.lower().man_exp();mh,eh=ball.upper().man_exp()
     assert F(int(ml))*F(2)**int(el)<=exact<=F(int(mh))*F(2)**int(eh)
     checks+=1
 return {'correlated_hidden_state_entry_checks':checks,'all_exact_joint_kernels_enclosed':True}


def root_box_controls():
 rng=random.Random(711);checks=0;rejected=0
 for case in range(80):
  # Generate a family of genuinely compatible rational input points inside one
  # declared box; verify its roots all survive interval elimination/Newton.
  base=[F(4,5),F(3,4),F(5,2),F(7,3),F(1,6)]
  radius=F(1,10**7); pts=[]
  for j in range(9):
   pts.append([v+radius*F(rng.randint(-5,5),5) for v in base])
  ds=[]
  for u,v,A0,B,x in pts:
   ds.append((pair_delta(F(1,3),u,v,A0,B,x),pair_delta(F(5,3),u,v,A0,B,x)))
  inp=[Interval(min(p[i] for p in pts),max(p[i] for p in pts)) for i in range(4)]
  d1=Interval(min(d[0] for d in ds),max(d[0] for d in ds))
  d2=Interval(min(d[1] for d in ds),max(d[1] for d in ds))
  try: cert=certify_pair(F(1,3),F(5,3),*inp,d1,d2)
  except (ZeroDivisionError,ValueError):rejected+=1;continue
  xlo,xhi=map(F,cert['enclosure'])
  assert all(xlo<=p[4]<=xhi for p in pts)
  checks+=len(pts)
 return {'independent_consistent_points_enclosed':checks,'boxes_rejected':rejected}


def main():
 t=time.perf_counter()
 result={'schema':'nonagon_possq.independent_controls.v1','spectral':spectral_controls(),
         'uniformization':uniformization_controls(),'convolution':independent_convolution_controls(),
         'joint_kernel':kernel_composition_controls(),'root_boxes':root_box_controls(),'verdict':'PASS','seconds':time.perf_counter()-t}
 (ROOT/'qa/independent_receipt.json').write_text(json.dumps(result,indent=2))
 print(json.dumps(result,indent=2))

if __name__=='__main__':main()
