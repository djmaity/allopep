"""Bounded, resumable rejection sampling of physically valid PepGLAD designs."""

import csv
import fcntl
import hashlib
import json
import os
import time
from dataclasses import dataclass
from contextlib import contextmanager
from pathlib import Path

from pepglad_outputs import load_round_outputs, removed_structures, site_rounds


ROOT = Path(__file__).resolve().parent


@contextmanager
def sampling_lock(out_dir):
    directory = Path(out_dir)
    directory.mkdir(parents=True, exist_ok=True)
    with (directory / '.sampling.lock').open('a') as handle:
        try:
            fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise RuntimeError(f'Another sampling process is using {directory}') from exc
        try:
            yield
        finally:
            fcntl.flock(handle, fcntl.LOCK_UN)


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def atomic_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + '.tmp')
    with temporary.open('w', encoding='utf-8') as handle:
        json.dump(value, handle, indent=2)
        handle.write('\n')
        handle.flush()
        os.fsync(handle.fileno())
    temporary.replace(path)


def next_round_size(target, attempted, valid, previous, cap):
    """Estimate demand independently per site, limiting growth and overshoot."""
    if valid >= target or attempted >= cap:
        return 0
    if not attempted:
        wanted = 8
    elif not valid:
        wanted = previous * 2
    else:
        # Integer arithmetic avoids rounding errors in ceil(remaining / p / 8).
        numerator = (target - valid) * attempted
        wanted = 8 * ((numerator + 8 * valid - 1) // (8 * valid))
        wanted = min(wanted, previous * 2)
    return min(cap - attempted, 64, max(8, wanted))


def require_resume(out_dir, prefix, resume):
    outputs = Path(out_dir) / f'{prefix}_PepGLAD_outputs'
    state = Path(out_dir) / f'{prefix}_sampling_summary.json'
    if not resume and (state.exists() or (outputs.exists() and any(outputs.iterdir()))):
        raise ValueError('Existing PepGLAD outputs require --resume; '
                         'use a separate input/output prefix for a fresh run')


def sampling_signature(args, input_pdb):
    return {
        'input_sha256': sha256(input_pdb),
        'length_min': args.length_min, 'length_max': args.length_max,
        'peptide_chain': args.peptide_chain,
        'validator_sha256': sha256(ROOT / 'pepglad_validation.py'),
        'pepglad_patch_sha256': sha256(ROOT / 'patches/PepGLAD.patch'),
        'pepglad_api_sha256': sha256(ROOT / 'PepGLAD/api/run.py'),
        'pepglad_model_sha256': sha256(ROOT / 'PepGLAD/models/LDM/ldm.py'),
        'checkpoint_sha256': sha256(ROOT / 'PepGLAD/checkpoints/codesign.ckpt'),
    }


@dataclass
class SamplingResult:
    selected: list
    shortfall: bool
    summary_path: Path


class SamplingController:
    def __init__(self, args, input_pdb, out_dir, prefix, contexts,
                 generate_round, validate_round):
        self.args = args
        self.input_pdb = input_pdb
        self.out_dir = Path(out_dir).resolve()
        self.prefix = prefix
        self.target = args.num_samples
        self.cap = args.max_num_peptides
        self.timeout = args.round_timeout_seconds
        self.generate_round = generate_round
        self.validate_round = validate_round
        self.summary_path = self.out_dir / f'{prefix}_sampling_summary.json'
        self.score_dir = self.out_dir / f'{prefix}_rosetta_score'
        self.sites_dir = self.out_dir / f'{prefix}_allosteric_sites'
        self.sites_dir.mkdir(parents=True, exist_ok=True)
        self.score_dir.mkdir(parents=True, exist_ok=True)
        self.accepted = {}
        signature = sampling_signature(args, input_pdb)
        require_resume(self.out_dir, prefix, args.resume)
        if self.summary_path.exists():
            self.state = json.loads(self.summary_path.read_text())
            if self.state.get('version') != 1 or self.state.get('signature') != signature:
                raise ValueError('Saved sampling input, generation settings, or '
                                 'validation policy differs; use a fresh output prefix')
            if len(self.state['sites']) != len(contexts):
                raise ValueError('Saved sampling site count differs from selected GAPS sites')
        else:
            self.state = {'version': 1, 'signature': signature, 'sites': []}
        if not contexts:
            raise ValueError('No GAPS sites selected for peptide generation')
        for rank, context in enumerate(contexts, 1):
            pocket = context.as_pepglad_pocket()
            pocket_file = self.sites_dir / f'{prefix}_allosteric_site_{rank}.json'
            site_dir = (self.out_dir / f'{prefix}_PepGLAD_outputs'
                        / f'{prefix}_allosteric_site_{rank}_codesign')
            if site_dir.exists() or len(self.state['sites']) >= rank:
                if not pocket_file.is_file() or json.loads(pocket_file.read_text()) != pocket:
                    raise ValueError(f'Saved PepGLAD pocket differs for site {rank}: {pocket_file}')
            else:
                atomic_json(pocket_file, pocket)
            if len(self.state['sites']) < rank:
                rounds = []
                for entry in site_rounds(site_dir, prefix):
                    rounds.append({
                        'directory': str(entry['directory'].relative_to(self.out_dir)),
                        'start': entry['start'], 'count': entry['count'],
                        'status': 'reserved', 'legacy': entry['legacy'],
                    })
                self.state['sites'].append({'site': rank, 'pocket': pocket, 'rounds': rounds})
            site = self.state['sites'][rank - 1]
            if site['site'] != rank or site['pocket'] != pocket:
                raise ValueError(f'Saved PepGLAD pocket differs for site {rank}')
            end = 0
            for entry in site['rounds']:
                if (type(entry['start']) is not int or type(entry['count']) is not int
                        or entry['start'] != end or entry['count'] < 1):
                    raise ValueError(f'Invalid reserved sample ranges for site {rank}')
                end += entry['count']
            if end > self.cap:
                raise ValueError(f'Site {rank} already charged {end} attempts, exceeding '
                                 f'--max_num_peptides {self.cap}')
            self.accepted[rank] = []
        self.state.pop('error', None)
        self.state['status'] = 'running'

    def _directory(self, entry):
        directory = (self.out_dir / entry['directory']).resolve()
        if not directory.is_relative_to(self.out_dir):
            raise ValueError('Sampling round directory is outside the output directory')
        return directory

    def _save(self):
        self.state['target_per_site'] = self.target
        self.state['max_num_peptides_per_site'] = self.cap
        for site in self.state['sites']:
            valid = sorted(self.accepted[site['site']], key=lambda design: design.sample)
            site['attempted'] = sum(entry['count'] for entry in site['rounds'])
            site['valid_count'] = len(valid)
            site['selected'] = [{'sample': design.sample, 'path': design.resolved_path}
                                for design in valid[:self.target]]
            site['failures'] = {
                key: sum(entry.get('failures', {}).get(key, 0) for entry in site['rounds'])
                for key in ('generation_invalid', 'relaxation_error',
                            'validation_invalid', 'validation_error', 'unreported')}
            site['completion_reason'] = (
                'target_reached' if len(valid) >= self.target else
                'max_num_peptides' if site['attempted'] >= self.cap else 'pending')
        atomic_json(self.summary_path, self.state)

    def _report_paths(self, directory):
        return tuple(directory / f'{self.prefix}_pepglad_validation{suffix}'
                     for suffix in ('.csv', '.jsonl', '_manifest.json'))

    def _read_validation(self, directory, designs):
        csv_path, jsonl_path, _ = self._report_paths(directory)
        with csv_path.open(newline='', encoding='utf-8') as handle:
            reader = csv.DictReader(handle)
            fields, rows = reader.fieldnames, list(reader)
        details = [json.loads(line) for line in jsonl_path.read_text().splitlines()
                   if line.strip()]
        by_path = {design.resolved_path: design for design in designs}
        if (not fields or len(rows) != len(by_path)
                or {row['path'] for row in rows} != set(by_path)
                or len(details) != len(rows)
                or {entry['result']['path'] for entry in details} != set(by_path)):
            raise ValueError(f'Validation reports differ from generated designs: {directory}')
        accepted = []
        removed = removed_structures(directory)
        for row in rows:
            design = by_path[row['path']]
            if (row['status'] not in ('valid', 'invalid', 'error')
                    or int(row['site']) != design.site or int(row['sample']) != design.sample):
                raise ValueError(f'Invalid validation status or sample metadata: {csv_path}')
            if row['status'] != 'error':
                digest = (sha256(design.path) if design.path.is_file()
                          else removed.get(design.path.name, {}).get('sha256'))
                if (row['peptide_chain'] != design.peptide_chain
                        or row['sha256'] != digest
                        or row['input_sha256'] != self.state['signature']['input_sha256']):
                    raise ValueError(f'Validation hashes or chains differ: {csv_path}')
            if row['status'] == 'valid':
                if design.path.name in removed or not design.path.is_file():
                    raise ValueError(f'Valid structure was removed: {design.path}')
                accepted.append(design)
        return accepted, fields, rows, details

    def _finish_round(self, site, entry, *, deadline, fresh=False):
        directory = self._directory(entry)
        if not entry.get('legacy'):
            metadata = {'site': site['site'], 'start': entry['start'], 'count': entry['count']}
            metadata_path = directory / 'round.json'
            if metadata_path.exists():
                if json.loads(metadata_path.read_text()) != metadata:
                    raise ValueError(f'Saved round metadata differs: {metadata_path}')
            elif entry['status'] != 'validated':
                atomic_json(metadata_path, metadata)
            else:
                raise ValueError(f'Saved round metadata missing: {metadata_path}')
        designs, records = load_round_outputs(
            directory, site['site'], self.prefix, entry['start'], entry['count'],
            partial=not fresh, include_removed=True)
        removed = removed_structures(directory, records)
        for design in designs:
            if not design.path.is_file() and design.path.name in removed:
                continue
            with design.path.open(encoding='utf-8') as handle:
                if not handle.readline().startswith('REMARK   1 CREATED WITH OPENMM'):
                    raise RuntimeError(f'PepGLAD site {site["site"]} has unrelaxed samples; '
                                       f'run recover_pepglad_relaxation.py on {directory}')
        if entry['status'] == 'validated':
            for relative, expected in entry['hashes'].items():
                path = directory / relative
                if (not path.exists() and relative in removed
                        and expected == removed[relative]['sha256']):
                    continue
                if not path.is_file() or sha256(path) != expected:
                    raise ValueError(f'Saved sampling artifact changed: {path}')
        else:
            directory.mkdir(parents=True, exist_ok=True)
            self.validate_round(designs, directory, self.prefix,
                                self.args.peptide_chain, self.input_pdb,
                                deadline=deadline, allow_errors=True)
        accepted, fields, rows, details = self._read_validation(directory, designs)
        if entry['status'] != 'validated':
            paths = [*self._report_paths(directory), *(design.path for design in designs)]
            summary = directory / 'summary.jsonl'
            if summary.exists():
                paths.append(summary)
            if not entry.get('legacy'):
                paths.append(directory / 'round.json')
            entry['hashes'] = {str(path.relative_to(directory)): sha256(path) for path in paths}
            entry['status'] = 'validated'
            entry['failures'] = {
                'generation_invalid': sum(r.get('status') == 'invalid' for r in records.values()),
                'relaxation_error': sum(r.get('status') == 'error' for r in records.values()),
                'validation_invalid': sum(row['status'] == 'invalid' for row in rows),
                'validation_error': sum(row['status'] == 'error' for row in rows),
                'unreported': entry['count'] - len(records),
            }
        self.accepted[site['site']].extend(accepted)
        self._save()
        self._remove_rejected_structures(entry, records, rows)
        return fields, rows, details

    def _remove_rejected_structures(self, entry, records, rows):
        """Journal deletion intent before unlinking, so interrupted cleanup can resume."""
        directory = self._directory(entry)
        removed = removed_structures(directory, records)
        candidates = [(Path(row['path']), row['status']) for row in rows
                      if row['status'] in ('invalid', 'error')]
        candidates.extend((directory / 'invalid_designs' / f'{identifier}.pdb', record['status'])
                          for identifier, record in records.items()
                          if record.get('status') in ('invalid', 'error'))
        for path, status in candidates:
            if path.is_file():
                if not path.resolve().is_relative_to(directory.resolve()):
                    raise ValueError(f'Rejected structure is outside its round: {path}')
                relative = str(path.relative_to(directory))
                item = {'sha256': sha256(path), 'status': status}
                if relative in removed and removed[relative] != item:
                    raise ValueError(f'Rejected structure changed before cleanup: {path}')
                removed[relative] = item
        if not removed:
            return
        journal = directory / 'removed_structures.json'
        atomic_json(journal, removed)
        # Retain the original PDB hashes alongside the journal's hash. Missing
        # files are accepted only when the rejection and original hash agree.
        entry['hashes'][journal.name] = sha256(journal)
        self._save()
        for relative, item in removed.items():
            path = directory / relative
            if path.exists() and sha256(path) != item['sha256']:
                raise ValueError(f'Rejected structure changed before cleanup: {path}')
        for relative in removed:
            (directory / relative).unlink(missing_ok=True)
        quarantine = directory / 'invalid_designs'
        if quarantine.is_dir() and not any(quarantine.iterdir()):
            quarantine.rmdir()

    def _append_reports(self, reports, *, reset=False):
        csv_path, jsonl_path, _ = self._report_paths(self.score_dir)
        fields, rows, details = reports
        mode = 'w' if reset else 'a'
        with csv_path.open(mode, newline='', encoding='utf-8') as handle:
            writer = csv.DictWriter(handle, fieldnames=fields)
            if reset:
                writer.writeheader()
            writer.writerows(rows)
        with jsonl_path.open(mode, encoding='utf-8') as handle:
            for detail in details:
                handle.write(json.dumps(detail) + '\n')

    def _write_manifest(self):
        # Reconstruct on resume/completion, avoiding a growing JSON rewrite each round.
        _, _, destination = self._report_paths(self.score_dir)
        temporary = destination.with_name(destination.name + '.tmp')
        with temporary.open('w', encoding='utf-8') as handle:
            handle.write('[\n')
            separator = ''
            for site in self.state['sites']:
                for entry in site['rounds']:
                    if entry['status'] != 'validated':
                        continue
                    manifest = self._report_paths(self._directory(entry))[2]
                    for row in json.loads(manifest.read_text()):
                        handle.write(separator + json.dumps(row))
                        separator = ',\n'
            handle.write('\n]\n')
        temporary.replace(destination)

    def run(self):
        first_report = True
        try:
            # Validate adopted/interrupted rounds once; reuse verified cached reports.
            for site in self.state['sites']:
                for entry in site['rounds']:
                    reports = self._finish_round(
                        site, entry, deadline=time.monotonic() + self.timeout)
                    self._append_reports(reports, reset=first_report)
                    first_report = False
            self._save()
            while True:
                scheduled = False
                for site in self.state['sites']:
                    previous = site['rounds'][-1]['count'] if site['rounds'] else 8
                    count = next_round_size(self.target, site['attempted'],
                                            site['valid_count'], previous, self.cap)
                    if not count:
                        continue
                    scheduled = True
                    site_dir = (self.out_dir / f'{self.prefix}_PepGLAD_outputs'
                                / f'{self.prefix}_allosteric_site_{site["site"]}_codesign')
                    directory = site_dir / f'round_{len(site["rounds"]) + 1:06d}'
                    if directory.exists():
                        raise FileExistsError(f'Refusing to overwrite sampling round: {directory}')
                    entry = {'directory': str(directory.relative_to(self.out_dir)),
                             'start': site['attempted'], 'count': count, 'status': 'reserved'}
                    site['rounds'].append(entry)
                    self._save()  # Charge the entire range before starting any subprocess.
                    directory.mkdir(parents=True)
                    atomic_json(directory / 'round.json',
                                {'site': site['site'], 'start': entry['start'], 'count': count})
                    print(f'Site {site["site"]}: {site["valid_count"]}/{self.target} valid; '
                          f'generating {count} candidates '
                          f'({site["attempted"]}/{self.cap} attempts charged)', flush=True)
                    deadline = time.monotonic() + self.timeout
                    pocket_file = (self.sites_dir
                                   / f'{self.prefix}_allosteric_site_{site["site"]}.json')
                    self.generate_round(self.args, self.input_pdb, pocket_file,
                                        directory, count, deadline=deadline)
                    entry['status'] = 'generated'
                    self._save()
                    reports = self._finish_round(site, entry, deadline=deadline, fresh=True)
                    self._append_reports(reports, reset=first_report)
                    first_report = False
                    print(f'Site {site["site"]}: {site["valid_count"]}/{self.target} valid '
                          f'after {site["attempted"]} charged attempts '
                          f'({site["completion_reason"]})', flush=True)
                if not scheduled:
                    break
            shortfall = any(site['valid_count'] < self.target for site in self.state['sites'])
            self.state['status'] = 'shortfall' if shortfall else 'complete'
            self._save()
            self._write_manifest()
            selected = [design for site in self.state['sites']
                        for design in sorted(self.accepted[site['site']],
                                             key=lambda design: design.sample)[:self.target]]
            return SamplingResult(selected, shortfall, self.summary_path)
        except BaseException as exc:
            self.state['status'] = 'error'
            self.state['error'] = f'{type(exc).__name__}: {exc}'
            self._save()
            self._write_manifest()
            raise
