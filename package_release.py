#!/usr/bin/env python3
"""Make a self-contained ZIP and SHA-256 file manifest of this release."""
from pathlib import Path
import argparse, hashlib, json, zipfile
ROOT=Path(__file__).resolve().parent
SKIP={'__pycache__','.git','.DS_Store'}
def members():
 return sorted(p for p in ROOT.rglob('*') if p.is_file() and not any(x in SKIP for x in p.relative_to(ROOT).parts) and p.suffix not in {'.pyc','.pyo'} and p.name!='BUNDLE_MANIFEST.json')
def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();target=a.output.resolve()
 if target==ROOT or ROOT in target.parents:raise ValueError('Write archive outside its own source tree')
 paths=members(); manifest={'schema':'nonagon.release.sha256-manifest.v1','release':'9.0','root':'r2k_nonagon_v9.0','physical_hardware_measurements':0,'scope':'All distributed files except this manifest; cached bytecode and Git metadata excluded','files':[{'path':p.relative_to(ROOT).as_posix(),'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}for p in paths]}
 (ROOT/'BUNDLE_MANIFEST.json').write_text(json.dumps(manifest,indent=2)+'\n')
 target.parent.mkdir(parents=True,exist_ok=True)
 with zipfile.ZipFile(target,'w',zipfile.ZIP_DEFLATED,compresslevel=6)as z:
  for p in paths+[ROOT/'BUNDLE_MANIFEST.json']:z.write(p,'r2k_nonagon_v9.0/'+p.relative_to(ROOT).as_posix())
 with zipfile.ZipFile(target)as z:
  bad=z.testzip()
  if bad:raise RuntimeError('Corrupt ZIP member '+bad)
 print(json.dumps({'status':'PASSED','archive':str(target),'files':len(paths)+1,'bytes':target.stat().st_size,'sha256':hashlib.sha256(target.read_bytes()).hexdigest()}))
if __name__=='__main__':main()
