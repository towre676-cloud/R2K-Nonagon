#!/usr/bin/env python3
"""Download the entire fixed public Zenodo record; preserve source bytes/checksums."""
import concurrent.futures, hashlib, json, time, urllib.request
from pathlib import Path
ROOT=Path(__file__).resolve().parent
DATA=ROOT/'data'

def get(url):
    with urllib.request.urlopen(url,timeout=60) as response:return response.read()

def run():
    DATA.mkdir(exist_ok=True)
    record_bytes=get('https://zenodo.org/api/records/15422339')
    (DATA/'zenodo_record.json').write_bytes(record_bytes)
    record=json.loads(record_bytes)
    selected=[x for x in record['files'] if x['key'].startswith('ID') and x['key'].endswith('.csv') or x['key'] in ['metadata.jsonld','keysight_eis_csvw_schema_0.1.0.json']]
    def download(entry):
        path=DATA/entry['key']
        if path.exists():payload=path.read_bytes()
        else:
            payload=get(entry['links']['self'])
            path.write_bytes(payload)
        md5=hashlib.md5(payload).hexdigest()
        if entry['checksum']!='md5:'+md5:raise ValueError('source checksum mismatch: '+entry['key'])
        return {'name':entry['key'],'source_url':entry['links']['self'],'source_checksum':entry['checksum'],'sha256':hashlib.sha256(payload).hexdigest(),'size_bytes':len(payload)}
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
        receipts=list(pool.map(download,selected))
    if len([x for x in receipts if x['name'].startswith('ID')])!=54:raise ValueError('missingcellfiles')
    manifest={'source':'https://zenodo.org/records/15422339','doi':'10.5281/zenodo.15422339','record_sha256':hashlib.sha256(record_bytes).hexdigest(),'metadata_license':record['metadata'].get('license'),'files':receipts}
    (DATA/'download_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print('Downloaded/verified',len(receipts),'files',flush=True)
if __name__=='__main__':run()
