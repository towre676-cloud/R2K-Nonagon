from pathlib import Path
import flint,subprocess
HERE=Path(__file__).resolve().parent
libs=Path(flint.__file__).resolve().parent.parent/'python_flint.libs'
lib=next(libs.glob('libflint-*.so.*'))
subprocess.run(['gcc','-O3','-shared','-fPIC',str(HERE/'native_counts.c'),str(lib),'-Wl,-rpath,'+str(libs),'-o',str(HERE/'native_counts.so')],check=True)
print('Built optional FLINT scalar/matrix accelerator')
