"""Check preservation of the two input branches and their distinct interfaces."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent


def run(output):
    manifest = json.loads((HERE / 'source_preservation.json').read_text())
    checked = 0
    for name, entry in manifest['files'].items():
        path = (ROOT / name).resolve()
        if not path.is_relative_to(ROOT.resolve()):
            raise ValueError('Preservation manifest path escapes release')
        if not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest() != entry['sha256']:
            raise ValueError('Retained input changed: ' + name)
        checked += 1
    production = ROOT / 'part_IX/verify_intervention.py'
    reference = ROOT / 'part_VIII/verify_intervention.py'
    if production.resolve() == reference.resolve():
        raise ValueError('Production gate and research oracle resolve to one file')
    if 'def verify(' not in production.read_text() or 'def recover(' not in reference.read_text():
        raise ValueError('Distinct production/research entry points are missing')
    result = {'schema': 'nonagon.integration.provenance-replay.v1', 'status': 'passed',
              'preserved_files': checked, 'separate_verifier_namespaces': True,
              'scope': manifest['scope'], 'archive_sha256': manifest['archive_sha256'],
              'source_sha256': {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
                                for p in (production, reference)}}
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({'status': 'passed', 'preserved_files': checked,
                      'separate_verifier_namespaces': True}))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=HERE / 'provenance_receipt.json')
    run(parser.parse_args().output)
