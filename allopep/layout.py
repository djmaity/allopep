"""Repository anchors and the output directory layout.

Keep this module dependency-free apart from the standard library: it is
imported from the ``gaps`` and ``PepGLAD`` environments as well as from
``allopep``.
"""

import re
from dataclasses import dataclass
from pathlib import Path

PACKAGE_DIR = Path(__file__).resolve().parent
REPO_ROOT = PACKAGE_DIR.parent
OUTPUT_BASE = Path('output')

_VALIDATION_SUFFIXES = ('.csv', '.jsonl', '_manifest.json')


@dataclass(frozen=True)
class OutputLayout:
    """Every path AlloPep writes for one input structure.

    Names are derived here and nowhere else, so renaming an artifact is a
    single edit. Site numbers start at 1.
    """

    root: Path
    prefix: str

    @classmethod
    def for_input(cls, pdb_file, base=OUTPUT_BASE):
        prefix = Path(pdb_file).stem
        return cls(Path(base) / f'{prefix}_output', prefix)

    @classmethod
    def at(cls, root, prefix):
        """Adopt an existing output directory, as the resume paths do."""
        return cls(Path(root), prefix)

    def resolved(self):
        return OutputLayout(self.root.resolve(), self.prefix)

    # Files directly under the output directory.

    @property
    def input_pdb(self):
        return self.root / f'{self.prefix}.pdb'

    @property
    def gaps_scores(self):
        return self.root / f'{self.prefix}_GAPS_scores.json'

    @property
    def gaps_output_pdb(self):
        return self.root / f'{self.prefix}_GAPS_output.pdb'

    @property
    def sampling_summary(self):
        return self.root / f'{self.prefix}_sampling_summary.json'

    @property
    def peptides_csv(self):
        return self.root / f'{self.prefix}_peptides.csv'

    @property
    def lock_file(self):
        return self.root / '.sampling.lock'

    # Directories.

    @property
    def sites_dir(self):
        return self.root / f'{self.prefix}_allosteric_sites'

    @property
    def designs_dir(self):
        return self.root / f'{self.prefix}_PepGLAD_outputs'

    @property
    def valid_complexes_dir(self):
        return self.root / f'{self.prefix}_valid_complexes'

    @property
    def score_dir(self):
        return self.root / f'{self.prefix}_vina_score'

    # Per-site and per-round paths.

    def pocket_file(self, site):
        return self.sites_dir / f'{self.prefix}_allosteric_site_{site}.json'

    def site_dir(self, site):
        return (self.designs_dir
                / f'{self.prefix}_allosteric_site_{site}_codesign')

    def round_dir(self, site, index):
        return self.site_dir(site) / f'round_{index:06d}'

    def valid_complex_file(self, site, sample):
        return (self.valid_complexes_dir / f'site_{site}'
                / f'{self.prefix}_site_{site}_peptide_{sample}.pdb')

    def site_number(self, path):
        """Recover a site number from a design or pocket path, or None."""
        pattern = re.compile(
            rf'^{re.escape(self.prefix)}_allosteric_site_([1-9]\d*)'
            r'(?:_codesign|\.json)$')
        match = pattern.fullmatch(Path(path).name)
        return int(match.group(1)) if match else None

    # Reports. The same three validation reports are written both into
    # the score directory and into each sampling round directory.

    def validation_reports(self, directory):
        return tuple(
            Path(directory) / f'{self.prefix}_pepglad_validation{suffix}'
            for suffix in _VALIDATION_SUFFIXES)

    @property
    def vina_manifest(self):
        return self.score_dir / f'{self.prefix}_vina_manifest.csv'

    @property
    def vina_results(self):
        return self.score_dir / f'{self.prefix}_vina_results.csv'

    @property
    def vina_ranked(self):
        return self.score_dir / f'{self.prefix}_vina_ranked.csv'

    # PyMOL session stems; the session and its script share a stem.

    @property
    def site_session_stem(self):
        return f'{self.prefix}_allosteric_site'

    @property
    def aligned_session_stem(self):
        return f'{self.prefix}_pepglad_aligned'
