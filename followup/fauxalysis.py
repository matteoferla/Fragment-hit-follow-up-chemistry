from typing import List
import os
from rdkit import Chem

def make_fauxalysis(hits: List[Chem.Mol], target_name: str, base_folder='.') -> None:
    """
    Given a list of hits, make a fragalysis-download-like folder structure

    :param hits:
    :param target_name:
    :param base_folder:
    :return:

    .. codeblock::python
    target_name = '👾👾👾'
    make_fauxalysis(hits, target_name, os.path.join(os.getcwd(), 'fauxalysis'))
    """
    if not os.path.exists(base_folder):
        os.mkdir(base_folder)
    os.environ['FRAGALYSIS_DATA_DIR'] = base_folder
    for hit in hits:
        hit_name: str = hit.GetProp('_Name')
        hit_path = os.path.join(base_folder, f'{target_name}', 'aligned', f'{target_name}-{hit_name}')
        os.makedirs(hit_path, exist_ok=True)
        Chem.MolToMolFile(hit, os.path.join(hit_path, f'{target_name}-{hit_name}.mol'))


