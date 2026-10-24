"""Content hashing and atomic file replacement.

Every report and state file the pipeline writes is published by renaming a
completed temporary file over the destination, so an interrupted write
never leaves a partial artifact for a later resume to adopt.
"""

import csv
import hashlib
import json
import os
from contextlib import contextmanager
from pathlib import Path

_BLOCK_SIZE = 1024 * 1024


def sha256(path):
    """Stream a file's SHA-256 digest."""
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(_BLOCK_SIZE), b''):
            digest.update(block)
    return digest.hexdigest()


@contextmanager
def atomic_path(path):
    """Yield a temporary sibling path, published only on a clean exit."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + '.tmp')
    try:
        yield temporary
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise
    temporary.replace(path)


@contextmanager
def atomic_text(path, *, newline=None):
    """Yield a text handle on a temporary sibling, flushed and fsynced."""
    with atomic_path(path) as temporary:
        with temporary.open('w', encoding='utf-8', newline=newline) as handle:
            yield handle
            handle.flush()
            os.fsync(handle.fileno())


def atomic_json(path, value):
    with atomic_text(path) as handle:
        json.dump(value, handle, indent=2)
        handle.write('\n')


def atomic_csv(path, fields, rows, *, extrasaction='ignore'):
    with atomic_text(path, newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields,
                                extrasaction=extrasaction)
        writer.writeheader()
        writer.writerows(rows)
