"""Replay v14.1 exact access and intervention inversion certificates."""
from pathlib import Path
import hashlib,json,platform,subprocess,sys,time

if sys.flags.optimize:raise RuntimeError('Certificate replay requires Python without -O')
HERE=Path(__file__).resolve().parent


def main():
    start=time.perf_counter()
    commands=[['-m','pytest','part_XIV/test_rational_identifiability.py','part_XIV/test_two_strength_cert.py','-q'],
              ['part_XIV/rational_access_replay.py'],['part_XIV/replay_two_strength_certificates.py']]
    for cmd in commands:subprocess.run([sys.executable,*cmd],cwd=HERE,check=True)
    names=['part_XIV/rational_identifiability.py','part_XIV/rational_access_replay.py','part_XIV/two_strength_cert.py']
    receipt={'verdict':'PASS','python':platform.python_version(),'seconds':time.perf_counter()-start,
             'source_sha256':{n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in names}}
    (HERE/'part_XIV/extension_replay_receipt.json').write_text(json.dumps(receipt,indent=2))
    print(json.dumps(receipt,indent=2))


if __name__=='__main__':main()
