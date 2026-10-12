"""Score validated PepGLAD complexes with Rosetta InterfaceAnalyzer."""

import csv
import hashlib
import math
import os
import shutil
import subprocess
from pathlib import Path


WEIGHTS_PATH = Path(__file__).with_name('weights') / 'optimized_rosetta_26_terms.wts'
WEIGHTS_SHA256 = 'e193a7acc0cfe0918a15662cac9d021f232797da81d57aa89e6a1c6005983044'
RESULT_FIELDS = (
    'id', 'site', 'sample', 'path', 'receptor_chains', 'peptide_chain',
    'interface', 'status', 'error', 'rosetta_binary', 'weights_sha256',
    'scorefile', 'log_file',
)


def resolve_rosetta_binary(value=None):
    """Find the externally installed InterfaceAnalyzer before design starts."""
    candidate = value if value is not None else shutil.which('InterfaceAnalyzer')
    if candidate is None:
        raise FileNotFoundError(
            'InterfaceAnalyzer is not on PATH; install Rosetta separately or '
            'pass --rosetta-bin PATH')
    binary = Path(candidate).expanduser().resolve()
    if not binary.is_file() or not os.access(binary, os.X_OK):
        raise FileNotFoundError(f'Rosetta executable is missing or not executable: {binary}')
    return binary


def verified_weights():
    """Return the bundled fitted weights only if they match the selected fit."""
    digest = hashlib.sha256(WEIGHTS_PATH.read_bytes()).hexdigest()
    if digest != WEIGHTS_SHA256:
        raise RuntimeError(f'Bundled Rosetta weights checksum differs: {WEIGHTS_PATH}')
    return WEIGHTS_PATH.resolve()


def pdb_chain_ids(path):
    """Read the coordinate chain IDs that Rosetta will see in a complex."""
    chains = set()
    with Path(path).open(encoding='utf-8', errors='replace') as handle:
        for line in handle:
            if line.startswith(('ATOM  ', 'HETATM')) and len(line) >= 22:
                chains.add(line[21])
    return chains


def parse_scorefile(path):
    """Read the final complete Rosetta score row and its interface energy."""
    header = None
    values = None
    with Path(path).open(encoding='utf-8', errors='replace') as handle:
        for line in handle:
            fields = line.split()
            if not fields or fields[0] != 'SCORE:':
                continue
            fields = fields[1:]
            if 'description' in fields:
                header = fields
                values = None
            elif header is not None:
                values = fields
    if header is None or values is None or len(header) != len(values):
        raise ValueError(f'Rosetta scorefile has no complete score row: {path}')
    metrics = {f'rosetta_{name}': value for name, value in zip(header, values)
               if name != 'description'}
    try:
        score = float(metrics['rosetta_dG_separated'])
    except (KeyError, ValueError) as exc:
        raise ValueError(f'Missing numeric dG_separated in {path}') from exc
    if not math.isfinite(score):
        raise ValueError(f'Non-finite dG_separated in {path}')
    metrics['rosetta_dG_separated'] = score
    return metrics


def _error_excerpt(output):
    lines = [line.strip() for line in output.splitlines() if line.strip()]
    return ' | '.join(lines[-8:])[-2000:]


def _score_design(design, score_dir, binary, weights, timeout):
    """Run one isolated fixed-pose calculation and retain its diagnostics."""
    interface = f'{"".join(design.receptor_chains)}_{design.peptide_chain}'
    run_dir = score_dir / 'rosetta_runs' / design.job_id
    run_dir.mkdir(parents=True, exist_ok=True)
    scorefile = run_dir / 'interface.sc'
    logfile = run_dir / 'interface.log'
    scorefile.unlink(missing_ok=True)
    logfile.unlink(missing_ok=True)
    record = {
        'id': design.job_id,
        'site': design.site,
        'sample': design.sample,
        'path': design.resolved_path,
        'receptor_chains': ''.join(design.receptor_chains),
        'peptide_chain': design.peptide_chain,
        'interface': interface,
        'status': 'error',
        'error': '',
        'rosetta_binary': str(binary),
        'weights_sha256': WEIGHTS_SHA256,
        'scorefile': str(scorefile.resolve()),
        'log_file': str(logfile.resolve()),
    }
    try:
        expected = set(design.receptor_chains) | {design.peptide_chain}
        actual = pdb_chain_ids(design.path)
        if actual != expected:
            raise ValueError(f'PDB chains {sorted(actual)} differ from '
                             f'PepGLAD chains {sorted(expected)}')
        command = [
            str(binary), '-s', design.resolved_path,
            '-interface', interface, '-score:weights', str(weights),
            '-compute_packstat', 'true', '-compute_interface_sc', 'true',
            '-pack_input', 'false', '-pack_separated', 'false',
            '-run:constant_seed', 'true', '-run:jran', '42',
            '-out:file:score_only', str(scorefile.resolve()),
            '-out:no_nstruct_label', 'true',
        ]
        completed = subprocess.run(
            command, cwd=run_dir, capture_output=True, text=True,
            timeout=timeout, check=False)
        output = completed.stdout + completed.stderr
        logfile.write_text(output, encoding='utf-8')
        if completed.returncode:
            raise RuntimeError(f'InterfaceAnalyzer exited with code '
                               f'{completed.returncode}: {_error_excerpt(output)}')
        if not scorefile.is_file():
            raise ValueError(f'Rosetta scorefile missing after success: {scorefile}')
        record.update(parse_scorefile(scorefile))
        record['status'] = 'ok'
    except subprocess.TimeoutExpired as exc:
        output = exc.stdout or b''
        if isinstance(output, bytes):
            output = output.decode('utf-8', errors='replace')
        record['error'] = f'InterfaceAnalyzer timed out after {timeout} seconds'
        logfile.write_text(output + '\n' + record['error'] + '\n',
                           encoding='utf-8')
    except (OSError, RuntimeError, ValueError) as exc:
        record['error'] = str(exc)
        with logfile.open('a', encoding='utf-8') as handle:
            handle.write('\n' + record['error'] + '\n')
    return record


def _write_csv(path, fieldnames, rows):
    """Replace one report only after the complete CSV has been written."""
    temporary = path.with_name(path.name + '.tmp')
    with temporary.open('w', newline='', encoding='utf-8') as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def rank_pepglad_outputs(valid_designs, score_dir, prefix, rosetta_bin,
                         timeout=1800):
    """Score valid designs and rank finite fitted interface energies."""
    if not valid_designs:
        raise ValueError('No validated PepGLAD structures to score')
    score_dir = Path(score_dir)
    score_dir.mkdir(parents=True, exist_ok=True)
    weights = verified_weights()
    binary = Path(rosetta_bin).resolve()
    manifest_file = score_dir / f'{prefix}_rosetta_manifest.csv'
    result_file = score_dir / f'{prefix}_rosetta_results.csv'
    ranked_file = score_dir / f'{prefix}_rosetta_ranked.csv'
    ranked_file.unlink(missing_ok=True)
    result_file.unlink(missing_ok=True)

    identifiers = [design.job_id for design in valid_designs]
    if len(set(identifiers)) != len(identifiers):
        raise ValueError('Duplicate PepGLAD job ID')
    manifest = [
        {'id': design.job_id, 'site': design.site, 'sample': design.sample,
         'path': design.resolved_path,
         'receptor_chains': ''.join(design.receptor_chains),
         'peptide_chain': design.peptide_chain,
         'interface': f'{"".join(design.receptor_chains)}_{design.peptide_chain}',
         'rosetta_binary': str(binary),
         'weights_sha256': WEIGHTS_SHA256}
        for design in valid_designs
    ]
    _write_csv(manifest_file, (
        'id', 'site', 'sample', 'path', 'receptor_chains', 'peptide_chain',
        'interface', 'rosetta_binary', 'weights_sha256'), manifest)

    results = [_score_design(design, score_dir, binary, weights, timeout)
               for design in valid_designs]
    metric_fields = sorted({'rosetta_dG_separated'} | {
        field for row in results for field in row if field.startswith('rosetta_')})
    _write_csv(result_file, (*RESULT_FIELDS, *metric_fields), results)
    ranked = [row.copy() for row in results if row['status'] == 'ok']
    if not ranked:
        raise RuntimeError(f'All {len(results)} Rosetta scoring jobs failed; '
                           f'see {result_file}')
    ranked.sort(key=lambda row: (
        row['rosetta_dG_separated'], row['site'], row['sample']))
    for rank, row in enumerate(ranked, 1):
        row['rank'] = rank
    _write_csv(ranked_file, ('rank', *RESULT_FIELDS, *metric_fields), ranked)
    return ranked
