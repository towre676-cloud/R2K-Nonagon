from pathlib import Path
import subprocess,sys
if sys.flags.optimize:raise RuntimeError("Run without -O")
root=Path(__file__).resolve().parent
for sub in ["possq-0.4.0","r2k_nonagon_v14.1"]:
 subprocess.run([sys.executable,"replay_extension.py"],cwd=root/sub,check=True)
subprocess.run([sys.executable,str(root/"qa/independent_checks.py")],cwd=root,check=True)
print("All extension checks passed.")
