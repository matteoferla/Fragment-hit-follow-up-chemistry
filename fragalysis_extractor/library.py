__all__ = ['populate_bad_inchi', 'get_matching_rows', 'read_metadata', 'to_unstereo_inchi',
           'to_singleorder_inchi', 'to_singleorder_unstereo_inchi']

import functools, operator, itertools
from rdkit import Chem
from rdkit.Chem import AllChem, PandasTools, Draw
import pandas as pd
from typing import List, Tuple
from IPython.display import Image

from collections import defaultdict

# this is kept for reference
synonyms = defaultdict(list)

def flatgrid(mols, *args, **kwargs) -> Image:
    copies: List[Chem.Mol] = [Chem.Mol(m) for m in mols]
    *map(AllChem.Compute2DCoords, copies),   # noqa, it's in place
    if 'legends' not in kwargs:
        kwargs['legends'] = [m.GetProp('_Name') if m.HasProp('_Name') else '-' for m in mols]
    return Draw.MolsToGridImage(copies, *args, **kwargs)

def to_unstereo_inchi(mol):
    mol = Chem.Mol(mol)
    Chem.RemoveStereochemistry(mol)
    return Chem.MolToInchiKey(mol)

def to_singleorder_inchi(original_mol):
    mol = AllChem.RemoveAllHs(original_mol)
    Chem.SanitizeMol(mol)
    bond: Chem.Bond
    atom: Chem.Atom
    for atom in mol.GetAtoms():
        atom.SetIsAromatic(False)
    for bond in mol.GetBonds():
        bond.SetIsAromatic(False)
        bond.SetBondType(Chem.BondType.SINGLE)
    for atom in mol.GetAtoms():
        atom.UpdatePropertyCache()
    Chem.SanitizeMol(mol)
    return Chem.MolToInchiKey(mol)

def to_singleorder_unstereo_inchi(original_mol):
    mol = AllChem.RemoveAllHs(original_mol)
    Chem.SanitizeMol(mol)
    Chem.RemoveStereochemistry(mol)
    bond: Chem.Bond
    atom: Chem.Atom
    for atom in mol.GetAtoms():
        atom.SetIsAromatic(False)
    for bond in mol.GetBonds():
        bond.SetIsAromatic(False)
        bond.SetBondType(Chem.BondType.SINGLE)
    for atom in mol.GetAtoms():
        atom.UpdatePropertyCache()
    Chem.SanitizeMol(mol)
    Chem.RemoveStereochemistry(mol)
    return Chem.MolToInchiKey(mol)

def populate_bad_inchi(df, mol_col: str):
    mols = df[mol_col]
    df['inchi'] = mols.apply(Chem.MolToInchiKey)
    df['unstereo_inchi'] = mols.apply(to_unstereo_inchi)
    df['singleorder_inchi'] = mols.apply(to_singleorder_inchi)
    df['singleorder_unstereo_inchi'] = mols.apply(to_singleorder_unstereo_inchi)

def get_matching_rows(row: pd.Series, ref: pd.DataFrame) -> Tuple[pd.DataFrame, str, str]:
    inchi_cols = ('inchi', 'unstereo_inchi', 'singleorder_inchi', 'singleorder_unstereo_inchi')
    for inchicol1,inchicol2 in itertools.product(inchi_cols, repeat=2):
        match = ref.loc[ref[inchicol2] == row[inchicol1]]
        if not match.empty:
            return match, inchicol1,inchicol2
    return pd.DataFrame(), '', ''

def read_metadata(path, library: pd.DataFrame) -> pd.DataFrame:
    meta = pd.read_csv(path)
    for comp_code_col in ('alternate_name', 'Compound code', 'Catalog ID'):
        if comp_code_col in meta.columns:
            break
    for smiles_col in ('SMILES', 'Smiles', 'smiles', 'smile'):
        if smiles_col in meta.columns:
            break
    for xstal_col in ('RealCrystalName', 'crystal_name', 'Code'):
        if xstal_col in meta.columns:
            break
    meta['desalted_smiles'] = meta[smiles_col].apply(lambda smiles: sorted(smiles.split('.'), key=len, reverse=True)[0])
    meta['mol'] = meta.desalted_smiles.apply(Chem.MolFromSmiles)
    meta.apply(lambda row: row.mol.SetProp('_Name', str(row[comp_code_col])), 1)
    meta['xcode'] = meta[xstal_col].str.extract(r'(x\d+)', expand=True)[0]
    populate_bad_inchi(meta, 'mol')
    # true_poised = df.loc[df[comp_code_col].isin(library['Catalog ID'].to_list())]
    verdicts = {}
    new_rows = []
    for i, row in meta.iterrows():
        # vanilla
        m, i1, i2 = get_matching_rows(row, library)
        if len(m) > 1:
            raise ValueError('Multiple matches')
        elif len(m) == 0:
            neo = pd.Series({'library': 'no DSiPoised match'})
            neo.name = i
            new_rows.append(neo)
        else: # single
            neo = m.iloc[0].copy()
            neo.name = i
            neo['matching_mode'] = f'{i1} -> {i2}'
            if row[comp_code_col] == neo['Catalog ID']: # full match
                neo['library'] = 'DSiPoised'
            else:
                neo['library'] = 'DSiPoised-synonym'
                synonyms[row[comp_code_col]].append(neo['Catalog ID'])
                # the meta['Catalog ID'] will be the official name
            new_rows.append(neo)
    new_df = pd.DataFrame(new_rows).copy()
    return pd.concat([meta, new_df], axis=1)
