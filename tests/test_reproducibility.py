"""Fail-closed restoration and exact reconstruction controls (also run with -O)."""
import gzip
import hashlib
import importlib.util
import io
import json
from pathlib import Path
from fractions import Fraction
from math import gcd
import sys
import tarfile
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / 'part_XIII')]
import restore_evidence as restore
from reconstruction import rational_reconstruction


def digest(data):
    return hashlib.sha256(data).hexdigest()


class Restoration(unittest.TestCase):
    def fixture(self, root, path='evidence.bin', missing=False):
        data = b'independent evidence fixture\x00\xff'
        output = io.BytesIO()
        with tarfile.open(fileobj=output, mode='w:gz') as archive:
            if not missing:
                member = tarfile.TarInfo(digest(data));member.size = len(data)
                archive.addfile(member, io.BytesIO(data))
        packed = output.getvalue();(root/'part').write_bytes(packed)
        (root/'PAYLOAD_INDEX.json').write_text(json.dumps(dict(parts=['part'],sha256=digest(packed),files=[dict(path=path,sha256=digest(data),bytes=len(data))])))
        (root/'CHUNKED_FILES.json').write_text('[]')
        (root/'IMPORT_MANIFEST.json').write_text('{"files": []}')
        return data

    def test_valid_and_tampered_payload(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);data=self.fixture(root);restore.restore(root)
            self.assertEqual((root/'evidence.bin').read_bytes(),data)
            (root/'part').write_bytes(b'corruption')
            with self.assertRaisesRegex(ValueError,'SHA-256 mismatch'):
                restore.restore(root)

    def test_missing_member_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);self.fixture(root,missing=True)
            with self.assertRaisesRegex(ValueError,'Missing payload'):
                restore.restore(root)

    def test_path_escape_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);self.fixture(root,path='../escape')
            with self.assertRaisesRegex(ValueError,'Unsafe evidence path'):
                restore.restore(root)


class Reconstruction(unittest.TestCase):
    def test_exhaustive_small_fractions(self):
        # 2*22^2 < 1009 gives uniqueness within the stated numerator/denominator bounds.
        modulus=1009
        for denominator in range(1,23):
            for numerator in range(-22,23):
                if gcd(numerator,denominator)==1:
                    residue=numerator*pow(denominator,-1,modulus)%modulus
                    self.assertEqual(rational_reconstruction(residue,modulus),Fraction(numerator,denominator))

    def test_invalid_modulus(self):
        with self.assertRaises(ValueError):
            rational_reconstruction(1,1)


if __name__=='__main__':
    unittest.main()
