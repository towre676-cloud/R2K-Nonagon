"""Replay every fully serialized retained pair certificate (no model truth read)."""
from pathlib import Path
import hashlib
import json

from two_strength_cert import (verify_certificate, verify_many_certificate,
                               verify_raw_certificate)


def main():
    here=Path(__file__).resolve().parent
    items=[]
    d=json.loads((here/'two_strength_design_results.json').read_text())
    for k,f in enumerate(d['fixtures']):
        items.append((f'design_fixture_{k}_two',verify_certificate(f['two'])))
        items.append((f'design_fixture_{k}_four',verify_many_certificate(f['four'])))
    d=json.loads((here/'two_strength_raw_design_results.json').read_text())
    for k,f in enumerate(d['fixtures']):
        items.append((f'raw_fixture_{k}',verify_raw_certificate(f['certificate'])))
    assert items and all(ok for _,ok in items)
    names=['two_strength_cert.py','test_two_strength_cert.py','two_strength_design.py',
           'two_strength_raw_design.py','two_strength_cert_results.json',
           'two_strength_design_results.json','two_strength_raw_design_results.json']
    receipt=dict(schema='nonagon.two_strength_replay.v1',status='PASS',
                 fixture_certificates=len(items),fixture_results=dict(items),
                 sha256={n:hashlib.sha256((here/n).read_bytes()).hexdigest() for n in names},
                 scope='serialized certificates replayed without consulting true_Gij fields; full campaigns separately reproducible')
    (here/'two_strength_replay_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({k:v for k,v in receipt.items() if k!='sha256'},indent=2))


if __name__=='__main__':main()
