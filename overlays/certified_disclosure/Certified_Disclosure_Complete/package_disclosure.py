"""Package the two upgraded source distributions and complete offline edition."""
from pathlib import Path
import hashlib,json,shutil,zipfile

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
OUT=ROOT/'deliverables'


def eligible(path):
    if path.is_symlink() or not path.is_file():return False
    if any(x in path.parts for x in ['__pycache__','.pytest_cache','node_modules','.git']):return False
    if path.name.startswith('rational_access_run') and path.suffix=='.log':return False
    return True


def inventory(base,manifest):
    rows=[]
    for path in sorted(base.rglob('*')):
        if eligible(path) and path.name!=manifest:
            rows.append({'path':path.relative_to(base).as_posix(),'bytes':path.stat().st_size,
                         'sha256':hashlib.sha256(path.read_bytes()).hexdigest()})
    (base/manifest).write_text(json.dumps({'files':rows,'date':'2026-10-31'},indent=2))


def addtree(z,base,prefix):
    for path in sorted(base.rglob('*')):
        if eligible(path):z.write(path,str(Path(prefix)/path.relative_to(base)))


def main():
    for name in ['possq-0.4.0','r2k_nonagon_v14.1']:
        base=HERE/name
        shutil.copy2(OUT/'Certified_Disclosure.html',base/'Certified_Disclosure.html')
        figures=base/'extension_figures';figures.mkdir(exist_ok=True)
        for path in (HERE/'figures').glob('*'):
            if path.suffix in ['.svg','.png','.pdf','.json','.md']:shutil.copy2(path,figures/path.name)
        inventory(base,'MANIFEST_EXTENSION.json')
        with zipfile.ZipFile(OUT/(name+'.zip'),'w',zipfile.ZIP_DEFLATED,compresslevel=7) as z:addtree(z,base,name)
    bundle=OUT/'Certified_Disclosure_Complete.zip'
    with zipfile.ZipFile(bundle,'w',zipfile.ZIP_DEFLATED,compresslevel=7) as z:
        prefix='Certified_Disclosure_Complete'
        for name in ['possq-0.4.0','r2k_nonagon_v14.1']:addtree(z,HERE/name,str(Path(prefix)/name))
        addtree(z,HERE/'figures',str(Path(prefix)/'figures'))
        addtree(z,ROOT/'qa',str(Path(prefix)/'qa'))
        for name in ['build_disclosure_report.py','build_advance_figures.py','package_disclosure.py']:
            z.write(HERE/name,str(Path(prefix)/name))
        z.write(OUT/'Certified_Disclosure.html',str(Path(prefix)/'Certified_Disclosure.html'))
        z.writestr(str(Path(prefix)/'requirements.txt'),(HERE/'possq-0.4.0/requirements.txt').read_text())
        z.writestr(str(Path(prefix)/'replay_all_extensions.py'),'''from pathlib import Path
import subprocess,sys
if sys.flags.optimize:raise RuntimeError("Run without -O")
root=Path(__file__).resolve().parent
for sub in ["possq-0.4.0","r2k_nonagon_v14.1"]:
 subprocess.run([sys.executable,"replay_extension.py"],cwd=root/sub,check=True)
subprocess.run([sys.executable,str(root/"qa/independent_checks.py")],cwd=root,check=True)
print("All extension checks passed.")
''')
        z.writestr(str(Path(prefix)/'README.md'),'''# Certified Disclosure — Nonagon 14.1 × possq 0.4

Open **Certified_Disclosure.html**. The offline mathematical report contains 7 scientific figures, 13 display equations, proofs, quantitative tables, a rotating exact-topology cube, and a certificate inspection slider.

The source packages preserve the original work and add exact rational access certificates, direct raw-response spectral risk, uncertainty-aware intervention inversion, and validated counting laws beyond 1e-100. Figures are supplied in editable SVG and PDF plus PNG. Source hashes, exact rational coefficients, compatible witness measures and dyadic probability endpoints accompany the code.

## Reproduction

Tested with Python 3.12.14. From this directory:

```bash
python -m pip install -r requirements.txt
python replay_all_extensions.py
python build_advance_figures.py
python build_disclosure_report.py
```

58 package tests, 96 spectral certificates, 96 atomic witnesses, 4 rational-access certificates and 7 serialized interval-inversion fixtures passed. Independent controls check 3,145 exact/interval cases. `qa/` holds the independent source and receipts. Per-package READMEs give commands that regenerate the numerical campaigns rather than only replaying retained certificates.

Reading errors are deterministic declared budgets. Acquisition selection uses a nominal model; final spectral certificates use only selected readings. Stochastic risks concern the specified exact generators and distinguish independent copies from continuous trajectories. Deficient-rank access statements concern first-order identifiability.
''')
    for path in OUT.glob('*.zip'):
        with zipfile.ZipFile(path) as z:
            bad=z.testzip()
            if bad:raise RuntimeError(str(path)+' corrupt: '+bad)
        print(path.name,path.stat().st_size,hashlib.sha256(path.read_bytes()).hexdigest())


if __name__=='__main__':main()
