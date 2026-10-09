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
import math
import os
import sys
from pathlib import Path

import numpy as np
import ost
from iotbx import pdb
from mmtbx.validation import ramalyze, rotalyze
from ost import conop, io
from ost.mol.alg import Accessibility, stereochemistry

monomer_library = Path(sys.prefix) / 'share' / 'monomers'
if monomer_library.is_dir():
    os.environ.setdefault('MMTBX_CCP4_MONOMER_LIB', str(monomer_library))


FIELDS = (
    'path', 'sha256', 'input_pdb', 'input_sha256', 'ost_version',
    'site', 'sample', 'peptide_chain', 'peptide_length',
    'status', 'reasons', 'peptide_clashes', 'interface_clashes',
    'receptor_clashes', 'peptide_bad_bonds', 'peptide_bad_angles',
    'inverted_ca', 'flat_ca', 'missing_backbone_atoms', 'backbone_breaks',
    'sequence_mismatches', 'missing_heavy_atoms',
    'rama_evaluated', 'rama_outliers', 'rotamer_evaluated', 'rotamer_outliers',
    'omega_evaluated', 'omega_unassessable', 'twisted_peptide_bonds',
    'cis_peptide_bonds', 'nonpro_cis_peptide_bonds',
    'max_omega_planarity_deviation_degrees',
    'contact_atom_pairs_4a', 'contact_atom_pairs_5a',
    'contact_residues_4a', 'contact_residues_5a',
    'receptor_contact_residues_4a', 'receptor_contact_residues_5a',
    'minimum_interface_distance_a', 'buried_sasa_a2',
    'receptor_ca_count', 'receptor_ca_rmsd_a', 'receptor_ca_max_displacement_a',
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


def _sequence_and_atoms(peptide, expected, details):
    observed = ''.join(residue.one_letter_code for residue in peptide.residues)
    if observed != expected:
        details['sequence_mismatches'].append({'expected': expected, 'observed': observed})
    library = conop.GetDefaultLib()
    for residue in peptide.residues:
        compound = library.FindCompound(residue.name)
        if compound is None:
            raise ValueError(f'No atom definition for peptide residue {residue.name}')
        required = {atom.name for atom in compound.atom_specs
                    if atom.element not in ('H', 'D') and atom.name != 'OXT'}
        actual = {atom.name for atom in residue.atoms if atom.element not in ('H', 'D')}
        missing = sorted(required - actual)
        if missing:
            details['missing_heavy_atoms'].append({
                'residue': f'{peptide.name}.{residue.number.num}.{residue.name}',
                'atoms': missing,
            })


def _conformation(path, peptide, details):
    hierarchy = pdb.input(file_name=str(path)).construct_hierarchy()
    for name, validator in (('rama', ramalyze.ramalyze),
                            ('rotamer', rotalyze.rotalyze)):
        checked = validator(hierarchy, quiet=True)
        results = [item for item in checked.results if item.chain_id == peptide.name]
        details[f'{name}_evaluated'] = len(results)
        details[f'{name}_outliers'] = [
            {'residue': f'{item.chain_id}.{item.resid.strip()}.{item.resname}',
             'score': round(item.score, 4),
             **({'phi': round(item.phi, 2), 'psi': round(item.psi, 2)}
                if name == 'rama' else {})}
            for item in results if item.outlier
        ]
    details['omega_evaluated'] = 0
    for left, right in zip(peptide.residues, peptide.residues[1:]):
        torsion = right.omega_torsion
        if not torsion.IsValid():
            details['omega_unassessable'].append({
                'from': left.number.num, 'to': right.number.num,
            })
            continue
        degrees = (math.degrees(torsion.angle) + 180) % 360 - 180
        deviation = min(abs(degrees), 180 - abs(degrees))
        details['omega_evaluated'] += 1
        bond = {'from': left.number.num, 'to': right.number.num,
                'omega_degrees': round(degrees, 2),
                'planarity_deviation_degrees': round(deviation, 2)}
        details['omega_bonds'].append(bond)
        if deviation > 30:
            details['twisted_peptide_bonds'].append(bond)
        if abs(degrees) < 30:
            details['cis_peptide_bonds'].append(bond)
            if right.name != 'PRO':
                details['nonpro_cis_peptide_bonds'].append(bond)


def _receptor_displacement(entity, input_pdb, receptor_names):
    original = io.LoadPDB(str(input_pdb))
    reference, model = [], []
    for name in receptor_names:
        old, new = original.FindChain(name), entity.FindChain(name)
        if not old.IsValid() or not new.IsValid():
            raise ValueError(f'Missing receptor chain {name} in input or model')
        if [r.name for r in old.residues] != [r.name for r in new.residues]:
            raise ValueError(f'Receptor sequence differs in chain {name}')
        for before, after in zip(old.residues, new.residues):
            old_ca, new_ca = before.FindAtom('CA'), after.FindAtom('CA')
            if not old_ca.IsValid() or not new_ca.IsValid():
                raise ValueError(f'Missing receptor CA in chain {name}')
            reference.append(_position(old_ca))
            model.append(_position(new_ca))
    if len(reference) < 3:
        raise ValueError('At least three receptor CA atoms are needed for alignment')
    reference, model = np.asarray(reference), np.asarray(model)
    x, y = model - model.mean(axis=0), reference - reference.mean(axis=0)
    u, _, vt = np.linalg.svd(x.T @ y)
    rotation = u @ np.diag([1, 1, np.linalg.det(u @ vt)]) @ vt
    deviations = np.linalg.norm(x @ rotation - y, axis=1)
    return len(deviations), float(np.sqrt(np.mean(deviations ** 2))), float(deviations.max())


def validate_structure(entry, peptide_chain=None):
    """Return a per-complex result with counts and violation details."""
    path = Path(entry['path']).resolve()
    input_pdb = Path(entry['input_pdb']).resolve()
    result = {'path': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
              'input_pdb': str(input_pdb),
              'input_sha256': hashlib.sha256(input_pdb.read_bytes()).hexdigest(),
              'ost_version': ost.__version__, 'site': entry['site'], 'sample': entry['sample']}
    entity = io.LoadPDB(str(path))
    chains = [chain for chain in entity.chains if chain.residues]
    if not chains:
        raise ValueError('No residues parsed from PDB')
    peptide_chain = peptide_chain or entry['pep_chain']
    if peptide_chain != entry['pep_chain']:
        raise ValueError('Peptide chain override differs from PepGLAD summary')
    peptides = [chain for chain in chains if chain.name == peptide_chain]
    if len(peptides) != 1 or len(chains) < 2:
        raise ValueError(f'Expected peptide chain {peptide_chain} and receptor chain(s)')
    receptor_names = entry['rec_chains']
    if {chain.name for chain in chains} != {*receptor_names, peptide_chain}:
        raise ValueError('Complex chains differ from PepGLAD summary')
    peptide = peptides[0]
    result['peptide_chain'] = peptide_chain
    result['peptide_length'] = len(peptide.residues)

    details = _stereochemistry_details(entity, peptide_chain)
    details.update({
        'inverted_ca': [], 'flat_ca': [],
        'missing_backbone_atoms': [], 'backbone_breaks': [],
        'sequence_mismatches': [], 'missing_heavy_atoms': [],
        'omega_unassessable': [], 'twisted_peptide_bonds': [],
        'cis_peptide_bonds': [], 'nonpro_cis_peptide_bonds': [], 'omega_bonds': [],
    })
    for scope in ('peptide', 'interface', 'receptor'):
        result[f'{scope}_clashes'] = sum(x['scope'] == scope for x in details['clashes'])
    result['peptide_bad_bonds'] = sum(x['scope'] != 'receptor' for x in details['bad_bonds'])
    result['peptide_bad_angles'] = sum(x['scope'] != 'receptor' for x in details['bad_angles'])

    peptide_atoms, residue_for_atom = _peptide_geometry(peptide, details)
    _sequence_and_atoms(peptide, entry['pep_seq'], details)
    _conformation(path, peptide, details)

    receptor_atoms, receptor_residue_for_atom = [], []
    for chain in chains:
        if chain.name == peptide_chain:
            continue
        for residue in chain.residues:
            for atom in residue.atoms:
                if atom.element not in ('H', 'D'):
                    receptor_atoms.append(_position(atom))
                    receptor_residue_for_atom.append(
                        f'{chain.name}.{residue.number.num}.{residue.name}')
    if not peptide_atoms or not receptor_atoms:
        raise ValueError('Missing peptide or receptor heavy atoms')
    distances = np.linalg.norm(np.asarray(peptide_atoms)[:, None, :] -
                               np.asarray(receptor_atoms)[None, :, :], axis=2)
    result['minimum_interface_distance_a'] = round(float(distances.min()), 3)
    for cutoff, label in ((4.0, '4a'), (5.0, '5a')):
        contacts = distances <= cutoff
        result[f'contact_atom_pairs_{label}'] = int(contacts.sum())
        result[f'contact_residues_{label}'] = len({
            residue_for_atom[i] for i in np.flatnonzero(contacts.any(axis=1))})
        result[f'receptor_contact_residues_{label}'] = len({
            receptor_residue_for_atom[i] for i in np.flatnonzero(contacts.any(axis=0))})
    rec_view = entity.Select('cname=' + ','.join(receptor_names))
    pep_view = entity.Select('cname=' + peptide_chain)
    complex_sasa = Accessibility(entity)
    receptor_sasa = Accessibility(rec_view)
    peptide_sasa = Accessibility(pep_view)
    buried = (receptor_sasa + peptide_sasa - complex_sasa) / 2
    if buried < -1:
        raise ValueError(f'Negative buried SASA: {buried:.2f} A^2')
    result['buried_sasa_a2'] = round(max(0, buried), 3)
    count, rmsd, maximum = _receptor_displacement(
        entity, input_pdb, receptor_names)
    result['receptor_ca_count'] = count
    result['receptor_ca_rmsd_a'] = round(rmsd, 3)
    result['receptor_ca_max_displacement_a'] = round(maximum, 3)
    for key in ('inverted_ca', 'flat_ca', 'missing_backbone_atoms', 'backbone_breaks',
                'sequence_mismatches', 'missing_heavy_atoms', 'rama_outliers',
                'rotamer_outliers', 'omega_unassessable', 'twisted_peptide_bonds',
                'cis_peptide_bonds', 'nonpro_cis_peptide_bonds'):
        result[key] = len(details[key])
    result['missing_heavy_atoms'] = sum(
        len(item['atoms']) for item in details['missing_heavy_atoms'])
    for key in ('rama_evaluated', 'rotamer_evaluated', 'omega_evaluated'):
        result[key] = details[key]
    result['max_omega_planarity_deviation_degrees'] = max(
        (bond['planarity_deviation_degrees'] for bond in details['omega_bonds']),
        default=None)

    reasons = [key for key in ('peptide_clashes', 'interface_clashes',
               'peptide_bad_bonds', 'peptide_bad_angles', 'inverted_ca',
               'flat_ca', 'missing_backbone_atoms', 'backbone_breaks',
               'sequence_mismatches', 'missing_heavy_atoms', 'rama_outliers',
               'rotamer_outliers', 'omega_unassessable', 'twisted_peptide_bonds',
               'nonpro_cis_peptide_bonds') if result[key]]
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
                result, details = validate_structure(entry, args.peptide_chain)
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
