"""Bounded real-data replay in an isolated tree, from an unrelated cwd."""
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]


class RelocatedPartXIII(unittest.TestCase):
    def test_case118_generation_and_two_independent_ranks(self):
        with tempfile.TemporaryDirectory() as directory:
            base = Path(directory)
            checkout = base / 'relocated'
            here = checkout / 'part_XIII'
            (here / 'data').mkdir(parents=True)
            (checkout / 'part_XI').mkdir()
            cwd = base / 'unrelated'
            cwd.mkdir()
            for name in ('record.py', 'verify13.py', 'verify13s.py', 'reconstruction.py', 'tomo.c'):
                shutil.copy2(ROOT / 'part_XIII' / name, here / name)
            shutil.copy2(ROOT / 'part_XI/gridlib.py', checkout / 'part_XI/gridlib.py')
            shutil.copy2(ROOT / 'part_XIII/data/case118.m', here / 'data/case118.m')
            expected = json.loads((ROOT / 'part_XIII/record_case118.json').read_text())
            subprocess.run(['gcc', '-O2', str(here / 'tomo.c'), '-o', str(here / 'tomo')],
                           check=True, capture_output=True, text=True, timeout=60)
            env = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1', MKL_NUM_THREADS='1')
            # No inherited Python path or replay tuning; preserve PYTHONOPTIMIZE.
            for key in ('PYTHONPATH', 'TAG', 'ORDER', 'PRIME_TOP'):
                env.pop(key, None)
            for script in ('record.py', 'verify13.py', 'verify13s.py'):
                result = subprocess.run([sys.executable, str(here / script), 'case118'],
                    cwd=cwd, env=env, capture_output=True, text=True, timeout=60)
                self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            fresh = json.loads((here / 'record_case118.json').read_text())
            self.assertTrue(fresh['exact_rational_recovery'])
            self.assertEqual(fresh['rational_mismatches'], 0)
            self.assertEqual(fresh['solution_sha256'], expected['solution_sha256'])
            for prefix in ('verify_', 'verify2_'):
                checked = json.loads((here / (prefix + 'case118.json')).read_text())
                self.assertTrue(checked['agrees_with_record'])
                self.assertTrue(checked['identifiable'])
                self.assertEqual(checked['rank'], fresh['unknowns'])
            # The verifier scripts emit a disagreement with exit 0; consumers must inspect JSON.
            fresh['identifiable'] = False
            (here / 'record_case118.json').write_text(json.dumps(fresh))
            for script, prefix in (('verify13.py', 'verify_'), ('verify13s.py', 'verify2_')):
                subprocess.run([sys.executable, str(here / script), 'case118'], cwd=cwd,
                    env=env, check=True, capture_output=True, text=True, timeout=60)
                checked = json.loads((here / (prefix + 'case118.json')).read_text())
                self.assertFalse(checked['agrees_with_record'])
            self.assertEqual(list(cwd.iterdir()), [])


if __name__ == '__main__':
    unittest.main()
