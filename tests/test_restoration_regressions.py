"""Independent adversarial coverage for every evidence restoration stage."""
import gzip
import hashlib
import io
import json
from pathlib import Path
import sys
import tarfile
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import restore_evidence as subject


def sha(data):
    return hashlib.sha256(data).hexdigest()


class RestorationRegressions(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name) / 'repo'
        self.root.mkdir()
        self.data = b'fresh evidence\x00\xff'
        self.payload()
        self.save('CHUNKED_FILES.json', [])
        self.save('IMPORT_MANIFEST.json', {'files': []})

    def save(self, name, value):
        (self.root / name).write_text(json.dumps(value))

    def payload(self, content=None, member_type=tarfile.REGTYPE, missing=False, size=None):
        output = io.BytesIO()
        with tarfile.open(fileobj=output, mode='w:gz') as archive:
            if not missing:
                item = tarfile.TarInfo(sha(self.data))
                item.type = member_type
                item.linkname = 'elsewhere'
                body = self.data if content is None else content
                item.size = len(body) if member_type == tarfile.REGTYPE else 0
                archive.addfile(item, io.BytesIO(body))
        packed = output.getvalue()
        (self.root / 'payload').write_bytes(packed)
        self.index = dict(parts=['payload'], sha256=sha(packed), files=[
            dict(path='out/payload', sha256=sha(self.data), bytes=len(self.data) if size is None else size),
            dict(path='out/alias', sha256=sha(self.data), bytes=len(self.data))])
        self.save('PAYLOAD_INDEX.json', self.index)

    def chunk(self, content=None, size=None):
        (self.root / 'chunk').write_bytes(self.data if content is None else content)
        self.save('CHUNKED_FILES.json', [dict(path='out/chunk', parts=['chunk'],
            sha256=sha(self.data), bytes=len(self.data) if size is None else size)])

    def compressed(self, content=None, size=None):
        (self.root / 'compressed').write_bytes(gzip.compress(self.data if content is None else content))
        self.save('IMPORT_MANIFEST.json', {'files': [dict(path='out/gzip', storage='gzip',
            stored_path='compressed', sha256=sha(self.data), bytes=len(self.data) if size is None else size)]})

    def test_all_stages_and_duplicate_hash_paths(self):
        self.chunk(); self.compressed()
        subject.restore(self.root)
        for name in ('payload', 'alias', 'chunk', 'gzip'):
            self.assertEqual((self.root / 'out' / name).read_bytes(), self.data)
        subject.restore(self.root)  # idempotent replay

    def test_archive_digest(self):
        (self.root / 'payload').write_bytes(b'corrupt')
        with self.assertRaisesRegex(ValueError, 'SHA-256 mismatch: payload archive'):
            subject.restore(self.root)
        self.assertFalse((self.root / 'out/payload').exists())

    def test_member_digest(self):
        self.payload(content=b'corrupt')  # valid archive digest, wrong member bytes
        with self.assertRaisesRegex(ValueError, 'SHA-256 mismatch'):
            subject.restore(self.root)
        self.assertFalse((self.root / 'out/payload').exists())

    def test_missing_member(self):
        self.payload(missing=True)
        with self.assertRaisesRegex(ValueError, 'Missing payload members'):
            subject.restore(self.root)

    def test_non_regular_member(self):
        for kind in (tarfile.SYMTYPE, tarfile.LNKTYPE, tarfile.DIRTYPE):
            with self.subTest(kind=kind):
                self.payload(member_type=kind)
                with self.assertRaisesRegex(ValueError, 'not a file'):
                    subject.restore(self.root)

    def test_payload_size(self):
        self.payload(size=999)
        with self.assertRaisesRegex(ValueError, 'size mismatch'):
            subject.restore(self.root)

    def test_chunk_corruption_and_size(self):
        for kwargs in ({'content': b'corrupt'}, {'size': 999}):
            with self.subTest(kwargs=kwargs):
                self.chunk(**kwargs)
                with self.assertRaisesRegex(ValueError, 'mismatch'):
                    subject.restore(self.root)
                self.assertFalse((self.root / 'out/chunk').exists())

    def test_gzip_corruption_and_size(self):
        for kwargs in ({'content': b'corrupt'}, {'size': 999}):
            with self.subTest(kwargs=kwargs):
                self.compressed(**kwargs)
                with self.assertRaisesRegex(ValueError, 'mismatch'):
                    subject.restore(self.root)
                self.assertFalse((self.root / 'out/gzip').exists())

    def test_malformed_gzip(self):
        self.compressed()
        (self.root / 'compressed').write_bytes(b'not gzip')
        with self.assertRaises(gzip.BadGzipFile):
            subject.restore(self.root)

    def test_unsafe_output_and_input_paths(self):
        for value in ('../escape', '/tmp/escape', 'out/../../escape', 'out\\escape'):
            for field in ('output', 'input'):
                with self.subTest(value=value, field=field):
                    self.payload()
                    if field == 'output': self.index['files'][0]['path'] = value
                    else: self.index['parts'] = [value]
                    self.save('PAYLOAD_INDEX.json', self.index)
                    with self.assertRaisesRegex(ValueError, 'Unsafe evidence path'):
                        subject.restore(self.root)

    def test_symlink_escape_input_and_output(self):
        outside = Path(self.temp.name) / 'outside'
        outside.mkdir()
        (outside / 'sentinel').write_bytes(b'unchanged')
        (self.root / 'link').symlink_to(outside, target_is_directory=True)
        for path in ('link/sentinel', 'link/new-file'):
            with self.subTest(path=path):
                with self.assertRaisesRegex(ValueError, 'escapes root'):
                    subject.write(self.root, path, b'overwrite')
        self.index['parts'] = ['link/sentinel']
        self.save('PAYLOAD_INDEX.json', self.index)
        with self.assertRaisesRegex(ValueError, 'escapes root'):
            subject.restore(self.root)
        self.assertEqual((outside / 'sentinel').read_bytes(), b'unchanged')
        self.assertFalse((outside / 'new-file').exists())


if __name__ == '__main__':
    unittest.main()
