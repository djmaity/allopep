"""Drive OpenStructure validation and match its reports to the designs.

The validator itself runs under the ``ost`` interpreter, so it lives in
:mod:`allopep.validation` and is invoked by path. This module owns the
driving side: the manifest it reads, the subprocess, and the invariant
that every design gets exactly one report row.
"""

import csv
import json
import shutil
import subprocess
from pathlib import Path

from .layout import PACKAGE_DIR
from .process_control import run_logged

STATUSES = ('valid', 'invalid', 'error')


def write_manifest(designs, path, input_pdb):
    """Describe each complex for the OpenStructure worker."""
    entries = [{
        'path': design.resolved_path, 'site': design.site,
        'sample': design.sample,
        'input_pdb': str(Path(input_pdb).resolve()),
        'pep_seq': design.sequence, 'pep_chain': design.peptide_chain,
        'rec_chains': list(design.receptor_chains),
    } for design in designs]
    Path(path).write_text(json.dumps(entries, indent=2) + '\n', encoding='utf-8')
    return entries


def run_ost(manifest, csv_file, jsonl_file, peptide_chain, *,
            deadline=None, log_file=None):
    """Run the validator under ost, with a deadline when one is given."""
    ost = shutil.which('ost')
    if ost is None:
        raise FileNotFoundError(
            'OpenStructure ost executable not found in the AlloPep environment.'
        )
    command = [ost, str(PACKAGE_DIR / 'validation.py'),
               '--manifest', str(Path(manifest).resolve()),
               '--csv', str(Path(csv_file).resolve()),
               '--jsonl', str(Path(jsonl_file).resolve())]
    if peptide_chain:
        command.extend(['--peptide-chain', peptide_chain])
    if deadline is None:
        subprocess.run(command, check=True)
    else:
        run_logged(command, log_file=log_file, deadline=deadline)


def read_rows(csv_file, designs):
    """Pair every report row with its design, requiring an exact match.

    Returns the CSV field names, the rows in file order, and the
    resolved-path to design mapping the callers index with.
    """
    csv_file = Path(csv_file)
    with csv_file.open(newline='', encoding='utf-8') as handle:
        reader = csv.DictReader(handle)
        fields, rows = reader.fieldnames, list(reader)
    by_path = {design.resolved_path: design for design in designs}
    if (not fields or len(rows) != len(designs) or len(by_path) != len(designs)
            or {row['path'] for row in rows} != set(by_path)):
        raise ValueError(
            f'Validation reports differ from generated designs: {csv_file}')
    if any(row['status'] not in STATUSES for row in rows):
        raise ValueError(f'Invalid validation status; see {csv_file}')
    return fields, rows, by_path
