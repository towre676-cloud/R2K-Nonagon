from pathlib import Path
import ctypes,time,json,sys,os
import numpy as np
from flint import arb,acb,ctx
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'audit_work/part_X'))
from origin_cauchy import sieve
class Euler:
 def __init__(self,Pmax,bits=128):
  self.primes=sieve(Pmax).astype(np.uint64)
  if os.environ.get('R2K_PURE_PYTHON')=='1' or not Path(__file__).with_name('euler_native.so').exists():
   self.lib=None;return
  self.lib=ctypes.CDLL(str(Path(__file__).with_name('euler_native.so')))
  self.lib.init_euler.argtypes=[ctypes.POINTER(ctypes.c_uint64),ctypes.c_size_t,ctypes.c_long];self.lib.init_euler.restype=ctypes.c_int
  self.lib.get_ab.argtypes=[ctypes.c_char_p,ctypes.POINTER(ctypes.c_void_p)];self.lib.get_ab.restype=ctypes.c_int
  self.lib.release_string.argtypes=[ctypes.c_void_p]
  self.primes=sieve(Pmax).astype(np.uint64)
  rc=self.lib.init_euler(self.primes.ctypes.data_as(ctypes.POINTER(ctypes.c_uint64)),len(self.primes),bits)
  assert rc==0,rc
 def __call__(self,r):
  if self.lib is None:
   z=acb(0,arb(r));A=acb(1);B=acb(0)
   for p in self.primes:
    p=int(p);lp=arb(p).log();q=(-z*lp).exp();A*=1-(q-1)**2/arb(p-1)**2;B+=lp**2/(p/q-1)**2
   return A,B
  out=(ctypes.c_void_p*4)();rc=self.lib.get_ab(arb(r).str(60,radius=True).encode(),out);assert rc==0,rc
  vals=[]
  for p in out:
   vals.append(arb(ctypes.string_at(p).decode()));self.lib.release_string(p)
  return acb(vals[0],vals[1]),acb(vals[2],vals[3])
def test():
 ctx.prec=128;start=time.time();E=Euler(20000)
 checks=[]
 for rs in ['0','0.01','0.5','1','3','10']:
  r=arb(rs);a,b=E(r);A=acb(1);B=acb(0);z=acb(0,r)
  for p in E.primes:
   p=int(p);lp=arb(p).log();q=(-z*lp).exp();A*=1-(q-1)**2/arb(p-1)**2;B+=lp**2/(p/q-1)**2
  assert a.overlaps(A) and b.overlaps(B)
  checks.append({'r':rs,'A':a.str(25,radius=True),'B':b.str(25,radius=True),'independent_python_overlap':True})
 print(json.dumps({'status':'PASS','checks':checks,'seconds':time.time()-start},indent=2))
if __name__=='__main__':test()
