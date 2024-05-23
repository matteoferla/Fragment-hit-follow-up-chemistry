"""
'Chemical isomorphism' (or 'crystallographic ambiguity')
refers to the phenomenon where elements of similar electron density,
such as carbon (C), nitrogen (N), and oxygen (O), cannot be easily distinguished in the electron density map.
This is because these atoms scatter X-rays similarly, resulting in overlapping or indistinguishable electron density.

An example is terminal amides, which makes elaborations fail if incorrectly mapped.

This script expand the list of molecules to include all of these possible configurations.
"""

from typing import List, Dict
from rdkit import Chem

def get_chemical_isomorphisms(mol: Chem.Mol, suffix='-iso') -> List[Chem.Mol]:
    """
    Given a molecule, return a list of molecules where crystallographical ambiguity is resolved.
    """
    query: Chem.Mol = create_unelemental_query(mol)
    # these include all isomorphisms
    to_idx_map = lambda indxs: dict(zip(range(mol.GetNumHeavyAtoms()), indxs))  # noqa
    idx_maps = [to_idx_map(qidxs) for qidxs in mol.GetSubstructMatches(query, uniquify=False)]
    # these are replacement isomorphism
    idx_maps = remove_elemental_isomorphisms(idx_maps, mol)
    if len(idx_maps) == 0:
        raise ValueError(f'An unexplained error: no replacement_isomorphic mappings - before filtering: {idx_maps}')
    elif len(idx_maps) == 1:
        return [mol]
    else:
        return maps_to_mols(idx_maps, mol, suffix)


# dependent functions

def maps_to_mols(idx_maps: List[Dict[int, int]], mol: Chem.Mol, suffix) -> List[Chem.Mol]:
    conformer: Chem.Conformer = mol.GetConformer()  # noqa
    copies = []
    name = mol.GetProp('_Name') if mol.HasProp('_Name') else '_'
    for mi, mapping in enumerate(idx_maps):
        if mi == 0:
            # original was spiked to the front of the mappings
            copies.append(mol)
            continue
        copy = Chem.Mol(mol)
        coconfomer = copy.GetConformer()  # noqa
        for q, t in mapping.items():
            coconfomer.SetAtomPosition(q, conformer.GetAtomPosition(t))
        copy.SetProp('_Name', f'{name}{suffix}{mi}')
        copies.append(copy)
    return copies


def remove_elemental_isomorphisms(idx_maps: List[Dict[int, int]], mol: Chem.Mol) -> List[Dict[int, int]]:
    """
    Benzene can be mapped twelve ways. I do not want this. Only maps where the elements differ
    """
    unique = []
    seen = []
    # spiking in the original for zeroth position
    for m in [dict(zip(range(mol.GetNumHeavyAtoms(), mol.GetNumHeavyAtoms())))] + idx_maps:  # noqa
        # this could be done mathematically, but this is quicker to write
        element_hash = ''.join([mol.GetAtomWithIdx(i).GetSymbol() for i in m.values()])
        if element_hash in seen:
            continue
        unique.append(m)
        seen.append(element_hash)
    return unique


def create_unelemental_query(mol: Chem.Mol) -> Chem.Mol:
    """
    A Query mol where the elements are stripped.
    Changing the atomic Zahlen to 0, will not work as the atoms have to be QueryAtoms,
    which is what `arthorian_quest.enquire` does.
    """
    from arthorian_quest import enquire

    subs = {}
    for atom in mol.GetAtoms():   # noqa
        if atom.GetSymbol() == 'H':
            subs[atom.GetIdx()] = None
        elif atom.GetSymbol() not in 'CNO':
            continue
        elif not atom.GetIsAromatic():
            subs[atom.GetIdx()] = '[C,N,O]'
        else:
            subs[atom.GetIdx()] = 'a'
    return enquire(mol, subs)

def main():
    import argparse
    import time
    import pandas as pd

    parser = argparse.ArgumentParser(description='Expand a list of molecules to include their replacement isomorphisms')
    parser.add_argument('input', help='Input SDF file')
    parser.add_argument('output', help='Output SDF file')

    args = parser.parse_args()
    with Chem.SDMolSupplier(args.input) as sdfh:
        mols = [mol for mol in sdfh]
    isomols = []
    for mol in mols:
        isomols.extend(get_chemical_isomorphisms(mol))
    with Chem.SDWriter(args.output) as sdfh:
        for mol in isomols:
            sdfh.write(mol)
