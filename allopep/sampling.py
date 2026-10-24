"""Bounded, resumable rejection sampling of physically valid PepGLAD designs."""

import csv
import fcntl
import json
import re
import time
from dataclasses import dataclass
from contextlib import contextmanager
from itertools import groupby, islice
from pathlib import Path

from .fileio import atomic_json, atomic_path, sha256
from .layout import PACKAGE_DIR, REPO_ROOT
from .outputs import is_relaxed, load_round_outputs, site_rounds
from .validation_reports import read_rows


@contextmanager
def sampling_lock(layout):
    layout.root.mkdir(parents=True, exist_ok=True)
    with layout.lock_file.open('a') as handle:
        try:
            fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise RuntimeError(
                f'Another sampling process is using {layout.root}') from exc
        try:
            yield
        finally:
            fcntl.flock(handle, fcntl.LOCK_UN)


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


def require_resume(layout, resume):
    outputs = layout.designs_dir
    state = layout.sampling_summary
    if not resume and (state.exists() or (outputs.exists() and any(outputs.iterdir()))):
        raise ValueError('Existing PepGLAD outputs require --resume; '
                         'use a separate input/output prefix for a fresh run')


@dataclass(frozen=True)
class SamplingConfig:
    """Validated generation and budget settings for one sampling run.

    Error messages name the command-line options they come from, since
    that is what the operator can change.
    """

    length_min: int
    length_max: int
    num_samples: int
    max_attempts: int
    round_timeout_seconds: int
    peptide_chain: str | None = None
    resume: bool = False
    peptides_per_site: bool = False

    def __post_init__(self):
        if self.length_min < 1 or self.length_max < self.length_min:
            raise ValueError(
                '--length_min must be positive and no greater than --length_max')
        if self.num_samples < 1:
            raise ValueError('--num_samples must be positive')
        if self.max_attempts < self.num_samples:
            raise ValueError(
                '--max_num_peptides must be positive and at least --num_samples')
        if self.round_timeout_seconds < 1:
            raise ValueError('--round-timeout-seconds must be positive')
        if (self.peptide_chain is not None
                and not re.fullmatch(r'[A-Za-z0-9]', self.peptide_chain)):
            raise ValueError('--peptide-chain must be one alphanumeric chain ID')


def sampling_signature(config, input_pdb):
    return {
        'input_sha256': sha256(input_pdb),
        'length_min': config.length_min, 'length_max': config.length_max,
        'peptide_chain': config.peptide_chain,
        'peptides_per_site': config.peptides_per_site,
        'validator_sha256': sha256(PACKAGE_DIR / 'validation.py'),
        'pepglad_patch_sha256': sha256(REPO_ROOT / 'patches/PepGLAD.patch'),
        'pepglad_api_sha256': sha256(REPO_ROOT / 'PepGLAD/api/run.py'),
        'pepglad_model_sha256': sha256(REPO_ROOT / 'PepGLAD/models/LDM/ldm.py'),
        'checkpoint_sha256': sha256(REPO_ROOT / 'PepGLAD/checkpoints/codesign.ckpt'),
    }


@dataclass
class SamplingResult:
    selected: list
    shortfall: bool
    summary_path: Path


def select_valid_designs(designs, config):
    """Limit passing complexes in stable site and sample order."""
    ordered = sorted(designs, key=lambda design: (design.site, design.sample))
    if not config.peptides_per_site:
        return ordered[:config.num_samples]
    return [design for _, group in groupby(ordered, key=lambda d: d.site)
            for design in islice(group, config.num_samples)]


class SamplingController:
    def __init__(self, config, input_pdb, layout, contexts,
                 generate_round, validate_round):
        self.config = config
        self.input_pdb = input_pdb
        self.layout = layout.resolved()
        self.out_dir = self.layout.root
        self.prefix = self.layout.prefix
        self.target = config.num_samples
        self.cap = config.max_attempts
        self.timeout = config.round_timeout_seconds
        self.generate_round = generate_round
        self.validate_round = validate_round
        self.summary_path = self.layout.sampling_summary
        self.layout.sites_dir.mkdir(parents=True, exist_ok=True)
        self.layout.score_dir.mkdir(parents=True, exist_ok=True)
        self.accepted = {}
        signature = sampling_signature(config, input_pdb)
        require_resume(self.layout, config.resume)
        if self.summary_path.exists():
            self.state = json.loads(self.summary_path.read_text())
            if (self.state.get('version') != 2
                    or self.state.get('signature') != signature):
                raise ValueError('Saved sampling input, generation settings, or '
                                 'validation policy differs; use a fresh output prefix')
            if len(self.state['sites']) != len(contexts):
                raise ValueError('Saved sampling site count differs from selected GAPS sites')
        else:
            self.state = {'version': 2, 'signature': signature, 'sites': []}
        if not contexts:
            raise ValueError('No GAPS sites selected for peptide generation')
        for rank, context in enumerate(contexts, 1):
            pocket = context.as_pepglad_pocket()
            pocket_file = self.layout.pocket_file(rank)
            site_dir = self.layout.site_dir(rank)
            if site_dir.exists() or len(self.state['sites']) >= rank:
                if not pocket_file.is_file() or json.loads(pocket_file.read_text()) != pocket:
                    raise ValueError(f'Saved PepGLAD pocket differs for site {rank}: {pocket_file}')
            else:
                atomic_json(pocket_file, pocket)
            if len(self.state['sites']) < rank:
                rounds = []
                for entry in site_rounds(site_dir, self.prefix):
                    rounds.append({
                        'directory': str(entry['directory'].relative_to(self.out_dir)),
                        'start': entry['start'], 'count': entry['count'],
                        'status': 'reserved',
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
        attempted = sum(entry['count'] for site in self.state['sites']
                        for entry in site['rounds'])
        if not config.peptides_per_site and attempted > self.cap:
            raise ValueError(f'Already charged {attempted} total attempts, '
                             f'exceeding --max_num_peptides {self.cap}')
        self.state.pop('error', None)
        self.state['status'] = 'running'

    def _directory(self, entry):
        directory = (self.out_dir / entry['directory']).resolve()
        if not directory.is_relative_to(self.out_dir):
            raise ValueError('Sampling round directory is outside the output directory')
        return directory

    def _save(self):
        self.state['peptides_per_site'] = self.config.peptides_per_site
        self.state['target'] = self.target
        self.state['max_num_peptides'] = self.cap
        self.state['attempted'] = sum(
            entry['count'] for site in self.state['sites']
            for entry in site['rounds'])
        valid_designs = [design for designs in self.accepted.values()
                         for design in designs]
        selected = select_valid_designs(valid_designs, self.config)
        selected_ids = {design.job_id for design in selected}
        self.state['valid_count'] = len(valid_designs)
        self.state['selected_count'] = len(selected)
        for site in self.state['sites']:
            valid = sorted(self.accepted[site['site']], key=lambda design: design.sample)
            site['attempted'] = sum(entry['count'] for entry in site['rounds'])
            site['valid_count'] = len(valid)
            site['selected'] = [{'sample': design.sample, 'path': design.resolved_path}
                                for design in valid
                                if design.job_id in selected_ids]
            site['failures'] = {
                key: sum(entry.get('failures', {}).get(key, 0) for entry in site['rounds'])
                for key in ('generation_invalid', 'relaxation_error',
                            'validation_invalid', 'validation_error', 'unreported')}
            count = (len(valid) if self.config.peptides_per_site
                     else self.state['valid_count'])
            attempted = (site['attempted'] if self.config.peptides_per_site
                         else self.state['attempted'])
            site['completion_reason'] = (
                'target_reached' if count >= self.target else
                'max_num_peptides' if attempted >= self.cap else 'pending')
        atomic_json(self.summary_path, self.state)

    def _report_paths(self, directory):
        return self.layout.validation_reports(directory)

    def _read_validation(self, directory, designs):
        csv_path, jsonl_path, _ = self._report_paths(directory)
        fields, rows, by_path = read_rows(csv_path, designs)
        details = [json.loads(line) for line in jsonl_path.read_text().splitlines()
                   if line.strip()]
        if (len(details) != len(rows)
                or {entry['result']['path'] for entry in details} != set(by_path)):
            raise ValueError(f'Validation details differ from generated designs: {directory}')
        accepted = []
        for row in rows:
            design = by_path[row['path']]
            if (int(row['site']) != design.site
                    or int(row['sample']) != design.sample):
                raise ValueError(f'Invalid validation status or sample metadata: {csv_path}')
            if row['status'] != 'error':
                if (row['peptide_chain'] != design.peptide_chain
                        or row['sha256'] != sha256(design.path)
                        or row['input_sha256'] != self.state['signature']['input_sha256']):
                    raise ValueError(f'Validation hashes or chains differ: {csv_path}')
            if row['status'] == 'valid':
                accepted.append(design)
        return accepted, fields, rows, details

    def _finish_round(self, site, entry, *, deadline, fresh=False):
        directory = self._directory(entry)
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
            partial=not fresh)
        for design in designs:
            if not is_relaxed(design.path):
                raise RuntimeError(f'PepGLAD site {site["site"]} has unrelaxed samples; '
                                   f'run python -m allopep.recover_relaxation on {directory}')
        if entry['status'] == 'validated':
            for relative, expected in entry['hashes'].items():
                path = directory / relative
                if not path.is_file() or sha256(path) != expected:
                    raise ValueError(f'Saved sampling artifact changed: {path}')
        else:
            directory.mkdir(parents=True, exist_ok=True)
            self.validate_round(designs, directory, self.layout,
                                self.config.peptide_chain, self.input_pdb,
                                deadline=deadline, allow_errors=True)
        accepted, fields, rows, details = self._read_validation(directory, designs)
        if entry['status'] != 'validated':
            paths = [*self._report_paths(directory), *(design.path for design in designs)]
            summary = directory / 'summary.jsonl'
            if summary.exists():
                paths.append(summary)
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
        return fields, rows, details

    def _append_reports(self, reports, *, reset=False):
        csv_path, jsonl_path, _ = self._report_paths(self.layout.score_dir)
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
        _, _, destination = self._report_paths(self.layout.score_dir)
        with atomic_path(destination) as temporary:
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
                # Resume at sites with fewer rounds before revisiting others.
                for site in sorted(self.state['sites'], key=lambda s: (
                        len(s['rounds']), s['site'])):
                    previous = site['rounds'][-1]['count'] if site['rounds'] else 8
                    target, cap = self.target, self.cap
                    if not self.config.peptides_per_site:
                        target = (site['valid_count'] + self.target
                                  - self.state['valid_count'])
                        cap = (site['attempted'] + self.cap
                               - self.state['attempted'])
                    count = next_round_size(target, site['attempted'],
                                            site['valid_count'], previous, cap)
                    if not count:
                        continue
                    scheduled = True
                    directory = self.layout.round_dir(
                        site['site'], len(site['rounds']) + 1)
                    if directory.exists():
                        raise FileExistsError(f'Refusing to overwrite sampling round: {directory}')
                    entry = {'directory': str(directory.relative_to(self.out_dir)),
                             'start': site['attempted'], 'count': count, 'status': 'reserved'}
                    site['rounds'].append(entry)
                    self._save()  # Charge the entire range before starting any subprocess.
                    directory.mkdir(parents=True)
                    atomic_json(directory / 'round.json',
                                {'site': site['site'], 'start': entry['start'], 'count': count})
                    scope = ('per site' if self.config.peptides_per_site
                             else 'across sites')
                    print(
                        f'Site {site["site"]}: target {self.target} {scope}; '
                        f'generating {count} candidates '
                        f'({self.state["attempted"]} total attempts charged)',
                        flush=True)
                    deadline = time.monotonic() + self.timeout
                    pocket_file = self.layout.pocket_file(site['site'])
                    self.generate_round(self.config, self.input_pdb, pocket_file,
                                        directory, count, deadline=deadline)
                    entry['status'] = 'generated'
                    self._save()
                    reports = self._finish_round(site, entry, deadline=deadline, fresh=True)
                    self._append_reports(reports, reset=first_report)
                    first_report = False
                    print(
                        f'Site {site["site"]}: '
                        f'{site["valid_count"]} valid here, '
                        f'{self.state["valid_count"]} across sites '
                        f'({site["completion_reason"]})', flush=True)
                if not scheduled:
                    break
            selected = select_valid_designs(
                [design for designs in self.accepted.values()
                 for design in designs], self.config)
            target = self.target * (len(self.state['sites'])
                                    if self.config.peptides_per_site else 1)
            shortfall = len(selected) < target
            self.state['status'] = 'shortfall' if shortfall else 'complete'
            self._save()
            self._write_manifest()
            return SamplingResult(selected, shortfall, self.summary_path)
        except BaseException as exc:
            self.state['status'] = 'error'
            self.state['error'] = f'{type(exc).__name__}: {exc}'
            self._save()
            self._write_manifest()
            raise
