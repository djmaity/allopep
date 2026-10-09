#!/usr/bin/env python3
"""OpenStructure checks for PepGLAD receptor-peptide complexes.

Run with ``ost pepglad_validation.py --manifest ... --csv ... --jsonl ...``.
Hydrogens are excluded from stereochemistry checks: OpenMM's hydrogen bond
lengths can disagree with OpenStructure's reference tables after relaxation.
"""

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
import ost
from ost import io
from ost.mol.alg import stereochemistry


FIELDS = (
    'path', 'sha256', 'ost_version', 'site', 'sample', 'peptide_chain', 'peptide_length',
    'status', 'reasons', 'peptide_clashes', 'interface_clashes',
    'receptor_clashes', 'peptide_bad_bonds', 'peptide_bad_angles',
    'inverted_ca', 'flat_ca', 'missing_backbone_atoms', 'backbone_breaks',
    'contact_residues_5a', 'minimum_interface_distance_a',
)


def _position(atom):
    return np.array((atom.pos.x, atom.pos.y, atom.pos.z), dtype=float)


def _violation_scope(atoms, peptide_chain):
    names = {atom.residue.chain.name for atom in atoms}
    if names == {peptide_chain}:
        return 'peptide'
    if peptide_chain in names:
        return 'interface'
    return 'receptor'


def _stereochemistry_details(entity, peptide_chain):
    """Collect heavy-atom violations without modifying the structure."""
    heavy = entity.Select('ele!=H and ele!=D')
    checks = {
        'clashes': stereochemistry.GetClashes(heavy, tolerance=1.5),
        'bad_bonds': stereochemistry.GetBadBonds(heavy, tolerance=12),
        'bad_angles': stereochemistry.GetBadAngles(heavy, tolerance=12),
    }
    details = {}
    for kind, violations in checks.items():
        details[kind] = []
        for violation in violations:
            atoms = (violation.a1, violation.a2)
            if kind == 'bad_angles':
                atoms += (violation.a3,)
            details[kind].append({
                'scope': _violation_scope(atoms, peptide_chain),
                **violation.ToJSON(),
            })
    return details


def _peptide_geometry(peptide, details):
    """Check backbone continuity and L amino-acid chirality."""
    peptide_atoms = []
    residue_for_atom = []
    for residue in peptide.residues:
        residue_id = f'{peptide.name}.{residue.number.num}.{residue.name}'
        missing = [name for name in ('N', 'CA', 'C', 'O')
                   if not residue.FindAtom(name).IsValid()]
        if missing:
            details['missing_backbone_atoms'].append({
                'residue': residue_id, 'atoms': missing,
            })
        if residue.name != 'GLY':
            atoms = [residue.FindAtom(name) for name in ('CA', 'N', 'C', 'CB')]
            if all(atom.IsValid() for atom in atoms):
                ca, n, c, cb = (_position(atom) for atom in atoms)
                # L amino acids have a positive signed volume in this order.
                volume = float(np.dot(np.cross(n - ca, c - ca), cb - ca))
                if volume < -0.5:
                    details['inverted_ca'].append({
                        'residue': residue_id, 'signed_volume_a3': round(volume, 3),
                    })
                elif volume <= 0.5:
                    details['flat_ca'].append({
                        'residue': residue_id, 'signed_volume_a3': round(volume, 3),
                    })
            else:
                details['flat_ca'].append({
                    'residue': residue_id, 'reason': 'missing chiral-center atom',
                })
        for atom in residue.atoms:
            if atom.element not in ('H', 'D'):
                peptide_atoms.append(_position(atom))
                residue_for_atom.append(residue_id)
    for left, right in zip(peptide.residues, peptide.residues[1:]):
        c, n = left.FindAtom('C'), right.FindAtom('N')
        if c.IsValid() and n.IsValid():
            distance = float(np.linalg.norm(_position(c) - _position(n)))
            if not 1.1 <= distance <= 1.7:
                details['backbone_breaks'].append({
                    'from': left.number.num, 'to': right.number.num,
                    'c_n_distance_a': round(distance, 3),
                })
    return peptide_atoms, residue_for_atom


def validate_structure(path, site, sample, peptide_chain=None):
    """Return a per-complex result with counts and violation details."""
    path = Path(path).resolve()
    result = {'path': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
              'ost_version': ost.__version__, 'site': site, 'sample': sample}
    entity = io.LoadPDB(str(path))
    chains = [chain for chain in entity.chains if chain.residues]
    if not chains:
        raise ValueError('No residues parsed from PDB')
    if peptide_chain is None:
        peptide_chain = min(chains, key=lambda chain: len(chain.residues)).name
    peptides = [chain for chain in chains if chain.name == peptide_chain]
    if len(peptides) != 1 or len(chains) < 2:
        raise ValueError(f'Expected peptide chain {peptide_chain} and receptor chain(s)')
    peptide = peptides[0]
    result['peptide_chain'] = peptide_chain
    result['peptide_length'] = len(peptide.residues)

    details = _stereochemistry_details(entity, peptide_chain)
    details.update({
        'inverted_ca': [], 'flat_ca': [],
        'missing_backbone_atoms': [], 'backbone_breaks': [],
    })
    for scope in ('peptide', 'interface', 'receptor'):
        result[f'{scope}_clashes'] = sum(x['scope'] == scope for x in details['clashes'])
    result['peptide_bad_bonds'] = sum(x['scope'] != 'receptor' for x in details['bad_bonds'])
    result['peptide_bad_angles'] = sum(x['scope'] != 'receptor' for x in details['bad_angles'])

    peptide_atoms, residue_for_atom = _peptide_geometry(peptide, details)

    receptor_atoms = [_position(atom) for chain in chains if chain.name != peptide_chain
                      for residue in chain.residues for atom in residue.atoms
                      if atom.element not in ('H', 'D')]
    if not peptide_atoms or not receptor_atoms:
        raise ValueError('Missing peptide or receptor heavy atoms')
    distances = np.linalg.norm(np.asarray(peptide_atoms)[:, None, :] -
                               np.asarray(receptor_atoms)[None, :, :], axis=2)
    result['minimum_interface_distance_a'] = round(float(distances.min()), 3)
    result['contact_residues_5a'] = len({residue_for_atom[i] for i in np.flatnonzero(
        (distances <= 5.0).any(axis=1))})
    for key in ('inverted_ca', 'flat_ca', 'missing_backbone_atoms', 'backbone_breaks'):
        result[key] = len(details[key])

    reasons = [key for key in ('peptide_clashes', 'interface_clashes',
               'peptide_bad_bonds', 'peptide_bad_angles', 'inverted_ca',
               'flat_ca', 'missing_backbone_atoms', 'backbone_breaks') if result[key]]
    if result['contact_residues_5a'] == 0:
        reasons.append('no_receptor_contact_5a')
    result['status'] = 'invalid' if reasons else 'valid'
    result['reasons'] = ';'.join(reasons)
    return result, details


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest', required=True, type=Path)
    parser.add_argument('--csv', required=True, type=Path)
    parser.add_argument('--jsonl', required=True, type=Path)
    parser.add_argument('--peptide-chain')
    args = parser.parse_args()
    manifest = json.loads(args.manifest.read_text(encoding='utf-8'))
    rows = []
    with args.jsonl.open('w', encoding='utf-8') as handle:
        for entry in manifest:
            try:
                result, details = validate_structure(entry['path'], entry['site'],
                                                     entry['sample'], args.peptide_chain)
            except Exception as exc:
                result = {'path': entry['path'], 'site': entry['site'],
                          'sample': entry['sample'], 'status': 'error',
                          'reasons': f'{type(exc).__name__}: {exc}'}
                details = {}
            rows.append(result)
            handle.write(json.dumps({'result': result, 'details': details}) + '\n')
    with args.csv.open('w', newline='', encoding='utf-8') as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(rows)
    print(f'OpenStructure checked {len(rows)} complexes: '
          f'{sum(r["status"] == "valid" for r in rows)} valid, '
          f'{sum(r["status"] == "invalid" for r in rows)} invalid, '
          f'{sum(r["status"] == "error" for r in rows)} errors')


if __name__ == '__main__':
    main()
