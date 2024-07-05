__all__ = ['convert_perceived_bonds', 'convert_protonated_pH', 'convert_perceived_bonds_pH',
           'embed_props', 'parse_mol', 'parse_mols', 'remove_duplicated']

from rdkit import Chem
from rdkit.Chem import AllChem
import json
from openbabel import openbabel as ob
from typing import List
import pandas as pd
import re
from typing import List
from rdkit import Chem, RDLogger
import contextlib


def convert_perceived_bonds(mol: Chem.Mol) -> Chem.Mol:
    """
    Use OpenBabel to percieve bonds in a RDKit molecule

    :param mol:
    :return: different Chem.Mol returned
    """
    obmol = ob.OBMol()
    conv = ob.OBConversion()
    conv.SetInFormat("pdb")
    if not mol.HasProp('_Name'):
        mol.SetProp('_Name', 'molecule')
    conv.ReadString(obmol, Chem.MolToPDBBlock(mol))
    obmol.PerceiveBondOrders()
    conv.SetOutFormat("mol")
    return Chem.MolFromMolBlock(str(conv.WriteString(obmol)))


def convert_protonated_pH(mol: Chem.Mol, pH=7.4) -> Chem.Mol:
    """
    Add hydrogens to a molecule using OpenBabel to a rdkit molecule

    :param mol:
    :param pH:
    :return: different Chem.Mol returned
    """
    conv = ob.OBConversion()
    conv.SetInAndOutFormats('mol', 'mol')
    obmol = ob.OBMol()
    if not mol.HasProp('_Name'):
        mol.SetProp('_Name', 'molecule')
    conv.ReadString(obmol, Chem.MolToMolBlock(mol))
    obmol.CorrectForPH(pH)
    neomol: Chem.Mol = Chem.MolFromMolBlock(conv.WriteString(obmol).strip())
    for k, v in mol.GetPropsAsDict().items():
        with contextlib.suppress(Exception):
            if isinstance(v, bool):
                neomol.SetBoolProp(k, v)
            elif isinstance(v, int):
                neomol.SetIntProp(k, v)
            elif isinstance(v, float):
                neomol.SetDoubleProp(k, v)
            else:
                neomol.SetProp(k, v)
    return neomol


def convert_perceived_bonds_pH(mol: Chem.Mol, pH: float=7.4) -> Chem.Mol:
    """
    Add bond order and pH to a rdkit molecule

    :param mol:
    :param pH:
    :return: different Chem.Mol returned
    """
    obmol = ob.OBMol()
    conv = ob.OBConversion()
    conv.SetInFormat("pdb")
    if not mol.HasProp('_Name'):
        mol.SetProp('_Name', 'molecule')
    conv.ReadString(obmol, Chem.MolToPDBBlock(mol))
    obmol.PerceiveBondOrders()
    obmol.CorrectForPH(pH)
    conv.SetOutFormat("mol")
    neomol: Chem.Mol = Chem.MolFromMolBlock(conv.WriteString(obmol).strip())
    for k, v in mol.GetPropsAsDict().items():
        if isinstance(v, bool):
            neomol.SetBoolProp(k, v)
        elif isinstance(v, int):
            neomol.SetIntProp(k, v)
        elif isinstance(v, float):
            neomol.SetDoubleProp(k, v)
        else:
            neomol.SetProp(k, v)
    return neomol

def embed_props(mol: Chem.Mol, detail):
    name = detail['name']
    mol.SetProp('_Name', str(name))
    mol.SetProp('SMILES', str(Chem.MolToSmiles(mol)))
    for key, value in detail.items():
        with contextlib.suppress(Exception):
            if key in ('pdb_block', 'mol'):
                continue
            elif isinstance(value, float) and str(value) != 'nan':
                mol.SetDoubleProp(key, value)
            elif isinstance(value, float) and str(value) == 'nan':
                continue
            elif isinstance(value, int):
                mol.SetIntProp(key, value)
            elif isinstance(value, bool):
                mol.SetBoolProp(key, value)
            elif isinstance(value, str):
                mol.SetProp(key, value)
            elif isinstance(value, dict) or isinstance(value, list):
                mol.SetProp(key, json.dumps(value))


def parse_mol(detail) -> Chem.Mol:
    extracted = Chem.MolFromPDBBlock(detail['pdb_block'], proximityBonding=False)
    if 'SMILES' in detail and detail['SMILES'] and str(detail['SMILES']) != 'nan':
        detail['bond_annotation_method'] = 'from SMILES'
        fixer = Chem.MolFromSmiles(detail['SMILES'])
        fixed = AllChem.AssignBondOrdersFromTemplate(fixer, extracted)
        neofixed = convert_protonated_pH(fixed)
        if neofixed is None:
            print(f'Could not fix pH {detail["name"]}')
        else:
            fixed = neofixed
    else:
        detail['bond_annotation_method'] = 'OB percieved bonds'
        fixed = convert_perceived_bonds_pH(extracted)
    if fixed is None:
        fixed = extracted
        print(f'Could not fix {detail["name"]}')
    embed_props(fixed, detail)
    return fixed

def parse_mols(details, meta: pd.DataFrame) -> List[Chem.Mol]:
    for comp_code_col in ('alternate_name', 'Compound code', 'Catalog ID'):
        if comp_code_col in meta.columns:
            break
    hits = []
    for detail in details:
        name = detail.get('crystal_name', 'error')
        if not isinstance(name, str) or not re.search(r'x\d+', name):
            continue
        detail['xcode'] = re.search(r'(x\d+)', name).group(1)  # 'crystal_name': 'x0041_B_1' or MID2A-x0041
        if detail['xcode'] not in meta['xcode'].values:
            raise ValueError(f'xcode {detail["xcode"]} not found in meta')
        meta_info: pd.Series = meta.drop_duplicates('xcode').set_index('xcode').loc[detail['xcode']]
        cols = [c for c in ['pdb_entry', 'desalted_smiles', 'library'] if c in meta_info.index]
        detail.update(meta_info[cols].to_dict())
        detail['SMILES'] = detail['desalted_smiles'] if detail['desalted_smiles'] else ''
        detail['name'] = detail['crystal_name']
        # this is the problem that some hits have a synonym relative to enamine name
        detail['catalogue_id'] = meta_info['Catalog ID'] if meta_info['Catalog ID'] else meta_info[comp_code_col]
        try:
            RDLogger.DisableLog('rdApp.*')
            hits.append(parse_mol(detail))
            RDLogger.EnableLog('rdApp.*')
        except Exception as e:
            print(f'{e.__class__.__name__} in {name}: {e}')
    return remove_duplicated(hits)


def remove_duplicated(hits: List[Chem.Mol]) -> List[Chem.Mol]:
    """
    Crude deduplication based on InChIKey and coordinates sum.

    :param hits:
    :return:
    """
    seen = []
    def is_new(hit):
        h = hash(f'{Chem.MolToInchiKey(hit)}{hit.GetConformer().GetPositions().sum():.0f}')
        if h in seen:
            return False
        seen.append(h)
        return True

    return [hit for hit in hits if is_new(hit)]
