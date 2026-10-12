#!/usr/bin/env python3
"""Inspect or finish OpenMM relaxation in an interrupted PepGLAD site directory.

Run with the PepGLAD environment. The default is read-only; --apply moves
malformed generated PDBs aside and relaxes the remaining unfinished designs.
"""

import argparse
import json
import math
import os
import subprocess
import sys
from pathlib import Path


def pdb_coordinate_error(path):
    with path.open(encoding='utf-8') as handle:
        for number, line in enumerate(handle, 1):
            if not line.startswith(('ATOM  ', 'HETATM')):
                continue
            try:
                coordinates = [float(line[i:i + 8]) for i in (30, 38, 46)]
            except ValueError:
                return f'Invalid PDB coordinate field on line {number}'
            if not all(math.isfinite(value) and abs(value) < 1000.0
                       for value in coordinates):
                return f'Out-of-range PDB coordinate on line {number}'
    return None


def inspect_site(site_dir):
    summary_path = site_dir / 'summary.jsonl'
    records = [json.loads(line) for line in summary_path.read_text().splitlines()
               if line.strip()]
    ids = [record['id'] for record in records]
    if len(ids) != len(set(ids)):
        raise ValueError(f'Duplicate IDs in {summary_path}')
    if {path.stem for path in site_dir.glob('*.pdb')} != {
            record['id'] for record in records if record.get('status') != 'invalid'}:
        raise ValueError(f'PDB files do not match {summary_path}')
    relaxed, pending, invalid = [], [], []
    for record in records:
        if record.get('status') == 'invalid':
            continue
        path = site_dir / (record['id'] + '.pdb')
        with path.open(encoding='utf-8') as handle:
            if handle.readline().startswith('REMARK   1 CREATED WITH OPENMM'):
                relaxed.append(path)
                continue
        error = pdb_coordinate_error(path)
        (invalid if error else pending).append((path, error))
    return summary_path, records, relaxed, pending, invalid


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('site_dir', type=Path)
    parser.add_argument('--apply', action='store_true',
                        help='Quarantine malformed PDBs and relax pending designs')
    parser.add_argument('--timeout-seconds', type=int, default=180,
                        help='Maximum relaxation time per design (default: 180)')
    args = parser.parse_args()
    if args.timeout_seconds < 1:
        parser.error('--timeout-seconds must be positive')
    site_dir = args.site_dir.resolve()
    summary_path, records, relaxed, pending, invalid = inspect_site(site_dir)
    print(f'{site_dir}: {len(relaxed)} relaxed, {len(pending)} pending, '
          f'{len(invalid)} malformed')
    for path, error in invalid:
        print(f'  {path.name}: {error}')
    if not args.apply:
        return

    if invalid:
        quarantine = site_dir / 'invalid_designs'
        quarantine.mkdir(exist_ok=True)
        for path, error in invalid:
            destination = quarantine / path.name
            if destination.exists():
                raise FileExistsError(destination)
            path.rename(destination)
            for record in records:
                if record['id'] == path.stem:
                    record['status'] = 'invalid'
                    record['error'] = error
                    break
        temporary = summary_path.with_suffix('.jsonl.tmp')
        with temporary.open('w', encoding='utf-8') as handle:
            for record in records:
                handle.write(json.dumps(record) + '\n')
        os.replace(temporary, summary_path)

    if pending:
        pepglad_dir = Path(__file__).resolve().parent / 'PepGLAD'
        for index, (path, _) in enumerate(pending, 1):
            print(f'Relaxing {path.name} ({index}/{len(pending)})', flush=True)
            subprocess.run([
                sys.executable, '-m', 'evaluation.dG.openmm_relaxer',
                str(path), str(path),
            ], cwd=pepglad_dir, check=True, timeout=args.timeout_seconds)
    print('Relaxation recovery complete')


if __name__ == '__main__':
    main()
