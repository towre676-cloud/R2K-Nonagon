"""Restore Git-hosted evidence; hash checks remain active under Python -O.

Run from any working directory. No external archive is downloaded.
"""
import gzip
import hashlib
import io
import json
from pathlib import Path, PurePosixPath
import tarfile


def checked(data, expected, label):
    if hashlib.sha256(data).hexdigest() != expected:
        raise ValueError('SHA-256 mismatch: ' + label)
    return data


def target(root, relative):
    p = PurePosixPath(relative)
    if p.is_absolute() or '..' in p.parts or '\\' in relative:
        raise ValueError('Unsafe evidence path: ' + relative)
    result = root / relative
    if not result.resolve().is_relative_to(root.resolve()):
        raise ValueError('Evidence path escapes root: ' + relative)
    return result


def write(root, relative, data):
    p = target(root, relative)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_bytes(data)


def restore(root):
    root = Path(root)
    idx = json.loads((root / 'PAYLOAD_INDEX.json').read_text())
    archive = checked(b''.join(target(root, p).read_bytes() for p in idx['parts']),
                      idx['sha256'], 'payload archive')
    wanted = {}
    for row in idx['files']:
        target(root, row['path'])
        wanted.setdefault(row['sha256'], []).append(row)
    # Visit each member once; random extraction repeatedly seeks through gzip.
    with tarfile.open(fileobj=io.BytesIO(archive), mode='r|gz') as stream:
        for member in stream:
            if member.name not in wanted:
                continue
            if not member.isfile():
                raise ValueError('Evidence member is not a file: ' + member.name)
            data = checked(stream.extractfile(member).read(), member.name, member.name)
            for row in wanted.pop(member.name):
                if len(data) != row['bytes']:
                    raise ValueError('Evidence size mismatch: ' + row['path'])
                write(root, row['path'], data)
    if wanted:
        raise ValueError('Missing payload members: ' + ', '.join(sorted(wanted)))
    chunks = json.loads((root / 'CHUNKED_FILES.json').read_text())
    for row in chunks:
        data = checked(b''.join(target(root, p).read_bytes() for p in row['parts']),
                       row['sha256'], row['path'])
        if len(data) != row['bytes']:
            raise ValueError('Chunked evidence size mismatch: ' + row['path'])
        write(root, row['path'], data)
    files = json.loads((root / 'IMPORT_MANIFEST.json').read_text())['files']
    for row in files:
        if row.get('storage') == 'gzip':
            data = checked(gzip.decompress(target(root, row['stored_path']).read_bytes()),
                           row['sha256'], row['path'])
            if len(data) != row['bytes']:
                raise ValueError('Decompressed evidence size mismatch: ' + row['path'])
            write(root, row['path'], data)
    print('Restored and SHA-256-verified all Git-hosted compressed evidence: '
          f'{len(idx["files"])} payload paths, {len(chunks)} chunked objects, '
          f'{sum(r.get("storage") == "gzip" for r in files)} gzip paths.')


if __name__ == '__main__':
    restore(Path(__file__).resolve().parent)
