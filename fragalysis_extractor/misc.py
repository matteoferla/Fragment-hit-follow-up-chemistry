__all__ = ['add_nitrogen_charges', 'read_charged_sdf', 'gold_sdf_to_df']

import pandas as pd
from rdkit import Chem, RDLogger

def add_nitrogen_charges(mol):
    RDLogger.DisableLog('rdApp.*')
    mol.UpdatePropertyCache(strict=False)
    ps = Chem.DetectChemistryProblems(mol)
    if not ps:
        Chem.SanitizeMol(mol)
        return mol
    for p in ps:
        if p.GetType()=='AtomValenceException':
            at = mol.GetAtomWithIdx(p.GetAtomIdx())
            if at.GetAtomicNum()==7 and at.GetFormalCharge()==0 and at.GetExplicitValence()==4:
                at.SetFormalCharge(1)
    Chem.SanitizeMol(mol)
    RDLogger.EnableLog('rdApp.*')
    return mol

def read_charged_sdf(sdf_path):
    """
    SDMolSupplier does not like nitrogens without formal charges.

    :param sdf_path:
    :return:
    """
    mols = []
    with Chem.SDMolSupplier(sdf_path,sanitize=False) as sdfh:
        for mol in sdfh:
            mols.append(add_nitrogen_charges(mol))
    return mols

def gold_sdf_to_df(sdf_path) -> pd.DataFrame:
    """
    Gold SDF might not be liked by RDKit due to nitrogens.
    This function will add formal charges to nitrogens to make them happy.

    :param sdf_path:
    :return:
    """
    _data = []
    mols = read_charged_sdf(sdf_path)
    for mol in mols:
        _data.append({'index': int(mol.GetProp('_Name').split('|')[-2]),
                     'id':    mol.GetProp('_Name').split('|')[0],
                      'mol': mol,
                     **{k: v for k, v in mol.GetPropsAsDict().items() if 'Gold.Protein' not in k}
                    })
    return pd.DataFrame(_data)
