"""Build the optional Linux-64 FLINT arithmetic accelerator."""
from pathlib import Path
import glob,flint,subprocess,sys
base=Path(flint.__file__).resolve().parent.parent
libs=list((base/'python_flint.libs').glob('libflint*.so*'))
if not libs:raise SystemExit('Linux FLINT shared library not found. Use R2K_PURE_PYTHON=1 for the slower portable evaluator.')
lib=libs[0];here=Path(__file__).resolve().parent
subprocess.run(['cc','-O3','-shared','-fPIC',str(here/'euler_native.c'),str(lib),'-o',str(here/'euler_native.so'),'-Wl,-rpath,'+str(lib.parent)],check=True)
print(here/'euler_native.so')
