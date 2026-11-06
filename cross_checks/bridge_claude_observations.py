"""Lossless bridge from engineering certificates to Claude's recovery oracle.

No source files are overwritten. The engineering verifier remains the strict
certificate gate. Claude's independent Bareiss/cofactor implementation then
recovers the operator from exactly the same scalar observations and evaluates
all held-out responses using Cramer's rule.
"""
from __future__ import annotations
import argparse
from fractions import Fraction as F
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
import time


def load(path,name):
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec);sys.modules[name]=module
    spec.loader.exec_module(module)
    return module


def rational(value):
    if isinstance(value,(bool,float)):raise ValueError('Only exact rational serialization is admissible')
    return F(value)


def adapt(cert):
    doc=cert['observations'];n=doc['contract']['network_order']
    t=rational(doc['plan']['base_load']);measurements={};serialized=[]
    for row in doc['observations']:
        theta=tuple(sorted((int(i),rational(v)) for i,v in row['loads'].items()))
        if theta in measurements:raise ValueError('Duplicate measurement configuration')
        value=rational(row['response']);measurements[theta]=value
        serialized.append({'id':row['id'],'theta':[{'site':i,'weight':str(v)} for i,v in theta],
                           'response':str(value)})
    if len(measurements)!=n*(n+1)//2:raise ValueError('Unexpected scalar acquisition count')
    return n,t,measurements,serialized


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    here=Path(__file__).resolve().parent;root=here.parent
    production=root/'part_IX'
    reference=root/'part_VIII'
    ap.add_argument('--engineering',type=Path,default=production)
    ap.add_argument('--claude',type=Path,default=reference)
    ap.add_argument('--output',type=Path,default=here/'cross_schema_receipt.json')
    ap.add_argument('--records',type=Path,default=here/'bridge_records.jsonl')
    args=ap.parse_args();started=time.monotonic()
    strict_path=args.engineering/'verify_intervention.py';oracle_path=args.claude/'verify_intervention.py'
    strict=load(strict_path,'nonagon_engineering_certificate_gate')
    oracle=load(oracle_path,'nonagon_claude_bareiss_oracle')
    files=sorted((args.engineering/'data').glob('*.certificate.json'))+sorted((args.engineering/'data/cases').glob('*.certificate.json'))
    result={'schema':'nonagon.integration.cross-schema.v1',
            'mapping':{'port':'Both use zero-based site0','our_loads':'single(t,2t), pair(t,t)',
                       'claude_parameters':'t1=t,t2=2t,tp=t','acquisition_values':'unchanged exact rational scalar readings',
                       'strict_gate':'engineering verifier.verify(cert)',
                       'independent_recovery':'Claude recover(n,measurement_dictionary,t1,t2,tp)',
                       'heldout':'Claude Cramer/Bareiss response(K,theta)'},
            'source_sha256':{'engineering_verifier':hashlib.sha256(strict_path.read_bytes()).hexdigest(),
                             'claude_verifier':hashlib.sha256(oracle_path.read_bytes()).hexdigest()},
            'cases':[]}
    checks=0;heldouts=0;scalar_count=0;certificate_count=0;strict_checks=0
    args.records.parent.mkdir(parents=True,exist_ok=True)
    with args.records.open('w') as records:
        for index,path in enumerate(files):
            cert=json.loads(path.read_text());gate=strict.verify(cert)
            if not gate.get('valid'):raise ValueError('Strict engineering certificate gate rejected '+str(path))
            strict_checks+=gate['checks'];n,t,meas,serialized=adapt(cert)
            recovered=oracle.recover(n,meas,t,2*t,t)
            expected=[[rational(x) for x in row] for row in cert['recovered']['K']]
            if recovered!=expected:raise ValueError('Cross-implementation operator mismatch '+str(path))
            checks+=n*n;certificate_count+=1;scalar_count+=len(meas)
            for row in cert['observations']['observations']+cert['observations'].get('heldout',[]):
                loads={int(i):rational(v) for i,v in row['loads'].items()}
                predicted=oracle.response(recovered,loads)
                if predicted!=rational(row['response']):raise ValueError('Independent Cramer response mismatch '+str(path))
                checks+=1
            heldouts+=len(cert['observations'].get('heldout',[]))
            row={'certificate':str(path.relative_to(args.engineering)),
                 'certificate_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
                 'n':n,'scalar_readings':len(meas),'heldout_readings':len(cert['observations'].get('heldout',[])),
                 'operator_recovered_exactly':True,'all_cramer_predictions_exact':True}
            result['cases'].append(row)
            records.write(json.dumps({'schema':'nonagon.integration.scalar-measurements.v1',
                                      'source':row['certificate'],'source_sha256':row['certificate_sha256'],
                                      'n':n,'port_zero_based':0,'t1':str(t),'t2':str(2*t),'tp':str(t),
                                      'measurements':serialized},separators=(',',':'))+'\n')
            if (index+1)%20==0:print(json.dumps({'cross_schema_progress':index+1}),flush=True)
    if (certificate_count,scalar_count,heldouts,strict_checks,checks)!=(103,2430,309,119902,6994):
        raise ValueError('Retained cross-schema campaign totals differ from the integration contract')
    result.update(status='passed',certificates=certificate_count,scalar_observations=scalar_count,
                  heldout_observations=heldouts,strict_certificate_checks=strict_checks,
                  independent_operator_and_response_checks=checks,
                  reverse_mapping='Claude builtin(1,3,pair2) observations cannot be relabeled as our fixed(t,2t,pairt) protocol; generalize schema or acquire new readings explicitly.',
                  records_sha256=hashlib.sha256(args.records.read_bytes()).hexdigest(),
                  elapsed_seconds=round(time.monotonic()-started,3))
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('cases','mapping','reverse_mapping')},indent=2))


if __name__=='__main__':main()
