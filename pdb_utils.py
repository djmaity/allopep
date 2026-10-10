"""Helpers for selecting and splitting protein and peptide PDB chains."""

from pathlib import Path

from Bio.PDB import PDBIO, PDBParser, Select


def detect_peptide_chain(structure):
    """Choose the shortest nonempty standard-residue chain in the first model."""

    lengths = [
        (sum(residue.id[0] == ' ' for residue in chain), chain.id)
        for chain in structure[0]
    ]
    nonempty = [item for item in lengths if item[0] > 0]
    if not nonempty:
        raise ValueError('No peptide or protein chains found in the PDB file')
    return min(nonempty, key=lambda item: item[0])[1]


class _ChainSelection(Select):
    def __init__(self, peptide_chain, select_peptide):
        self.peptide_chain = peptide_chain
        self.select_peptide = select_peptide

    def accept_chain(self, chain):
        is_peptide = chain.id == self.peptide_chain
        return is_peptide if self.select_peptide else not is_peptide


def split_protein_peptide(pdb_file, peptide_chain=None, output_prefix=None):
    """Write separate peptide and protein PDBs and return the parsed structure."""

    pdb_path = Path(pdb_file)
    structure = PDBParser(QUIET=True).get_structure(pdb_path.stem, str(pdb_path))
    if peptide_chain is None:
        peptide_chain = detect_peptide_chain(structure)

    prefix = Path(output_prefix) if output_prefix is not None else pdb_path.with_suffix('')
    io = PDBIO()
    io.set_structure(structure)
    io.save(f'{prefix}_peptide.pdb', _ChainSelection(peptide_chain, True))
    io.save(f'{prefix}_protein.pdb', _ChainSelection(peptide_chain, False))
    return structure
