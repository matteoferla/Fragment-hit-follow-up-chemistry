Usage:


```python
from pathlib import Path

fragalys_download_path = Path('fragalysis-downloads')
library_path = Path('Enamine_DSI-poised_Library_plated_860cmds_20210302.sdf')
special_chains = {'XX01ZVNS2B': 'AB', 'NSP14': 'D'}
skip = []

import pymol2
from pathlib import Path
from rdkit import Chem

import pymol2
from rdkit import Chem
import fragalysis_extractor as fragex
import pandas as pd
from typing import List
from rdkit import Chem
from rdkit.Chem import AllChem, PandasTools, SaltRemover

library_cmds: List[Chem.Mol] = misc.read_charged_sdf(library_path.as_posix())
library = pd.DataFrame([{'ref_mol': mol, **mol.GetPropsAsDict()} for mol in library_cmds])
# desalt library
sr = SaltRemover.SaltRemover()
library['ref_mol'] = library.ref_mol.apply(sr.StripMol)
library['SMILES'] = library.ref_mol.apply(Chem.MolToSmiles)
library['MW (desalted)'] = library['MW (desalted)'].astype(float)
fragex.populate_bad_inchi(library, 'ref_mol')

for path in Path(fragalys_download_path).glob('*'):
    if not path.is_dir() or '.zip' in path.name:
        continue
    if (path / 'cleaned_hits.sdf').exists():
        continue
    if not (path / 'metadata.csv').exists():
        continue
    if not (path / 'reference.pdb').exists():
        continue
        
    target = path.stem
    if target in skip:
        continue
    print('*'*10, target, '*'*10)
    print('Parsing metadata')
    
    monomer_chains=special_chains.get(target, 'A')
    meta = fragex.read_metadata(path / 'metadata.csv', library)
    if sum(meta.library.str.contains('DSiPoised')) == 0:
        continue
    
    print('Parsing hits')
    with pymol2.PyMOL() as pymol:
        try:
            details = fragex.collate_details(pymol, path, monomer_chains=monomer_chains)
        except BaseException as error:
            pymol.cmd.save('test.pse')
    hits = fragex.parse_mols(details, meta)
    assert len(hits)
    with Chem.SDWriter((path / 'cleaned_hits.sdf').as_posix()) as sdfh:
        for hit in hits:
            sdfh.write(hit)

special_chains = {'XX01ZVNS2B': 'AB', 'NSP14': 'D'}

for path in Path('fragalysis-downloads').glob('*'):
    if not path.is_dir() or '.zip' in path.name:
        continue
    if (path / 'cleaned_hits.sdf').exists():
        continue
    if not (path / 'metadata.csv').exists():
        continue
    if not (path / 'reference.pdb').exists():
        continue
        
    target = path.stem
    if target in skip:
        continue
    print('*'*10, target, '*'*10)
    print('Parsing metadata')
    
    monomer_chains=special_chains.get(target, 'A')
    meta = fragex.read_metadata(path / 'metadata.csv', library)
    if sum(meta.library.str.contains('DSiPoised')) == 0:
        continue
    
    print('Parsing hits')
    with pymol2.PyMOL() as pymol:
        try:
            details = fragex.collate_details(pymol, path, monomer_chains=monomer_chains)
        except BaseException as error:
            pymol.cmd.save('test.pse')
    hits = fragex.parse_mols(details, meta)
    assert len(hits)
    with Chem.SDWriter((path / 'cleaned_hits.sdf').as_posix()) as sdfh:
        for hit in hits:
            sdfh.write(hit)
```