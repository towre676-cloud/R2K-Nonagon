from itertools import product
from fractions import Fraction as Q
from collections import defaultdict
import json
from pathlib import Path

def ordered_partitions(n):
 for m in range(1,n+1):
  for assignment in product(range(m),repeat=n):
   if set(assignment)==set(range(m)):
    yield [[i for i,a in enumerate(assignment) if a==j] for j in range(m)]
PARTS=list(ordered_partitions(4))
def span(seq):
 x=0;seen=[0]
 for k in seq:x+=k;seen.append(x)
 assert x==0
 return max(seen)-min(seen)
def cumulant(N,j,k):
 freqs=[j,-j,k,-k];out=Q(0)
 for p in PARTS:
  w=span([sum(freqs[i] for i in b) for b in p]);m=len(p);out+=Q((-1)**(m-1),m)*max(N-w,0)
 return out

def main():
 coeff=defaultdict(Q)
 for p in PARTS:
  f=[7,-7,3,-3];w=span([sum(f[i] for i in b) for b in p]);m=len(p);coeff[w]+=Q((-1)**(m-1),m)
 print('partitions',len(PARTS),'coefficients for j=7,k=3',dict(coeff))
 # Recover linear positive-part expression for j>=k.
 labels={0:'0',3:'k',7:'j',10:'j+k',4:'j-k'}
 print({labels[k]:str(v) for k,v in coeff.items() if v})
 for N in [1,2,3,5,10]:
  for j in range(1,15):
   for k in range(1,15):
    a=max(j,k);b=min(j,k)
    form=sum((v*max(N-w,0) for w,v in [(0,coeff[0]),(b,coeff[3]),(a,coeff[7]),(a+b,coeff[10]),(a-b,coeff[4])]),Q(0))
    assert cumulant(N,j,k)==form
    cov=form+(min(N,j)**2 if j==k else 0)
    if N==1:assert cov==0
 out={'status':'PASS','ordered_partitions':len(PARTS),'positive_part_coefficients':{labels[k]:str(v) for k,v in coeff.items() if v},'exact_integer_cases':5*14*14}
 Path(__file__).with_name('trace_derivation.json').write_text(json.dumps(out,indent=2))
if __name__=='__main__':main()
