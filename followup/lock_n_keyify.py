import requests
import zipfile
import io
from rdkit import Chem
from fragmenstein import Wictor  # Victor, but RDKit only
from typing import List, Tuple, Dict, Any

def get_target_data(target: str) -> dict:
    """
    Fetch target data from Fragalysis
    """
    response: requests.Response = requests.get(f'https://fragalysis.diamond.ac.uk/api/targets/?title={target}')
    info = response.json()
    assert info['next'] is None, '?? New feature: capped?'
    assert len(info['results']) == 1, 'ambiguous?'
    return info['results'][0]


def get_apo_pdbblocks(zip_url):
    """
    Fetch the apo pdbblocks from Fragalysis
    Get url from ``get_target_data``.
    """
    apo_pdbblocks = {}
    dejavu = []
    response: requests.Response = requests.get(zip_url)
    response.raise_for_status()
    with zipfile.ZipFile(io.BytesIO(response.content)) as z:
        for filename in sorted(z.namelist()):
            if '_apo.pdb' not in filename:
                continue
            # aligned/A71EV2A-x0469_0A/A71EV2A-x0469_0A_apo.pdb
            xchem_code = filename.split('/')[-1].split('-')[-1].replace('_apo.pdb', '')
            # code x0123_0A is the same structure as x0123_1A
            if xchem_code.split('_')[0] in dejavu:
                continue
            else:
                dejavu.append(xchem_code.split('_')[0])
            with z.open(filename) as f:
                apo_pdbblocks[xchem_code] = f.read().decode('utf-8')
    return apo_pdbblocks

def remove_altloc(pdbblock: str, remove_hetatm: bool=True) -> str:
    """Removes all altlocs and actually segi duplicates.
    Test line:

    ATOM      1  N   MET A   1       0.000  0.000  0.000  1.00 60.69           N
    """
    lines: List[str] = []
    seen: List[Tuple[str, str, int]] = []
    for line in pdbblock.split('\n'):
        if 'ANISOU' in line:
            continue # skip
        if line[:4] != 'ATOM' and line[:6] != 'HETATM':
            lines.append(line)
            continue
        if remove_hetatm and line[:6] == 'HETATM':
            continue
        atom_info = line[12:16].strip(), line[21].strip(), int(line[22:26].strip())
        if atom_info not in seen:
            lines.append(f'{line[:16]} {line[17:]}')
            seen.append(atom_info)
        else: # skip
            pass
    return '\n'.join(lines)

def score(vc: Chem.Mol, template_name:str, apo_pdbblock: str, **settings) -> Dict[str, Any]:
        """
        Minimised the compound (Fragalysis upload ready),
        in a frozen neighbourhood of the provided template.
        The compounds need to be in position already.

        Done using frozen neighbourhood discussed here:
        https://www.blopig.com/blog/2023/11/the-workings-of-fragmensteins-rdkit-neighbour-aware-minimisation/
        """
        wicky = Wictor([vc], pdb_block=apo_pdbblock, **settings)
        wicky.place(vc.GetProp('original SMILES'), long_name=f'{vc.GetProp("_Name")}_on_{template_name}')
        info = wicky.summarize()
        for k, v in vc.GetPropsAsDict().items():
            if isinstance(v, float):
                wicky.minimized_mol.SetDoubleProp(k, v)
            elif isinstance(v, int):
                wicky.minimized_mol.SetIntProp(k, v)
            else:
                wicky.minimized_mol.SetProp(k, v)
        wicky.minimized_mol.SetProp('ref_pdb', template_name)
        info['mol'] = wicky.minimized_mol
        info['template'] = template_name
        info['query'] = vc.GetProp('_Name')
        return info
