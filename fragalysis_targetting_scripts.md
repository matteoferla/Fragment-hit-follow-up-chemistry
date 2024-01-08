# Fragalysis

## XChem public targets

In [:octocat:: munged-Fragalysis-targets](https://github.com/matteoferla/munged-Fragalysis-targets)
the Fragalysis targets are explored by Blast to Xref with Uniprot,
yielding a [table of useful information](https://github.com/matteoferla/munged-Fragalysis-targets/blob/main/targets.csv).

NB. Not all data in Fragalysis has the same layout.

NB. The table is not polished by species popularity: a few protein are identical in some close species,
and the wrong one is picked —nobody, not even chimpanzees, cares for NUDT4 in chimpanzees.

To get the email address of a researcher involved in a target, asking is the best way.

> :construction: Rerun with newest targets

## Folder layout

Typical folder format in Fragalysis is:

* smiles.smi — optional
* metadata.csv — optional
* reference.pdb — optional
* extra_files/* — optional
* aligned/👻👻-👾👾/👾👾_bound.pdb.pdb

The latter can be replicated with `make_fauxalysis` from [fragment_elaboration_scripts.fauxalysis](fragment_elaboration_scripts/fauxalysis.py)

The dummy Fragalysis download structure is needed for a few things.

NB. `fauxalysis` is a pun on faux and Fragalysis stolen from John Chodera's group: therein it has a different role, I think.

## Fragment extraction

The problem is as follows:

1. A PDB structure may be a dimer in the crystallographic asymmetric unit, but a monomer in vivo. A compound may bind in both or in a single chain.
2. A compound may be listed as a salt, e.g. `Na.CC(=O)[O-]`, but obviously is solubilised in DMSO.
3. A compound may be racemic but only one isomer binds. 
4. A unique compound may bind more than once.
5. A compound may have crystallogrically ambiguous conformers (think glutamine and histidine orientations), the PDBs will have partial density and alt IDs.

As of July 2023, these have these effects:

1. This is addressed by aligning relevant monomer to the reference structure (see note), yielding structures with say _0A and _0B suffixes.
2. This is still results in the loss of bond order. Do not trust `_combined.sdf` file without fixing it  (see code)
3. Historically, the entry was either edited or a new entry created for the isomer. The latter is the official way.
4. Each will be its own structure, but not in legacy data

Legacy data:

There have been a few iterations of how Fragalysis aligns the crystals.
In a first iteration (by the Fragment-5) the chains in the assymetric unit were aligned to the reference in PyMOL.
As this cause corner case issues, the code was changed (by Tyler) to aligning the asymmetric unit.
Now, the code (by Connor) aligns key residues in the binding site.

In notebooks/frag-extraction.ipynb is the old case.

For the fixing of `_combined.sdf` data, the following code can be used:

```python
target_name = '...'
keep_prefix = False   # Having target_name prefix makes names longer... but easier to split etc.

# ------------------------------------------------

import pandas as pd
import operator
from rdkit import Chem
from rdkit.Chem import AllChem, PandasTools

metadata = pd.read_csv('metadata.csv')
# the crystal_name values are like TARGET-COMPOUND_\n\w
if keep_prefix:
    metadata['code'] = metadata.crystal_name
else:
    metadata['code'] = metadata.crystal_name.str.replace(f'{target_name}-', '')
metadata = metadata.set_index('crystal_name')
# remove the salts from the SMILES
metadata['smiles'] = metadata.new_smiles.apply(lambda v: sorted(v.split('.'), key=len)[-1])
# add the mol
with Chem.SDMolSupplier(f'{target_name}_combined.sdf') as sdf_r:
    hits = list(sdf_r)
metadata['mol'] = {h.GetProp('_Name'): h for h in hits}
metadata['mol'] = metadata.apply(lambda row: AllChem.AssignBondOrdersFromTemplate(Chem.MolFromSmiles(row.smiles), row.mol), axis=1)

PandasTools.WriteSDF(df=metadata.reset_index(), 
                     out=f'{target_name}.corrected.sdf',
                     molColName='mol', idName='code',
                     properties=['crystal_name', 'RealCrystalName', 'smiles', 'alternate_name', 'site_name'])
```

There is a RDKit fuction to split compounds and remove inorganics, but there is an acetate in there somewhere.
So brutally getting the longest SMILES is easy and foolproof.

    metadata.new_smiles.apply(lambda v: sorted(v.split('.'), key=len)[-1])


## Fragalysis Download

In [fragment_elaboration_scripts.fragalysis_download](fragment_elaboration_scripts/fragalysis_download.py) there is a script to download a Fragalysis target.

```python
from fragment_elaboration_scripts,fragalysis_download import QuickDownloader
import pandas as pd
from typing import List
print(f'Default settings are: {QuickDownloader.api_data}')

# Check if the target name is right
target_names: List[str] = QuickDownloader.retrieve_target_names()
target_name='Mpro'
assert target_name in target_names, f'Target named "{target_name}" not found in the list of targets'

# Download the data
quick = QuickDownloader(target_name=target_name)
quick.write_all(directory='downloads')
hits: pd.DataFrame = quick.to_pandas(star_dummy=True)

# Not all files have the reference pdb block, so if it does not the template is returned:
reference_pdbblock: str = quick.reference_pdbblock
```

NB. For legacy targets in particular, the zip may be broken.

## Fragalysis Upload

This is just to prep.
For a Colab notebook see: [upload_prep.ipynb](https://colab.research.google.com/github/matteoferla/Fragment-hit-follow-up-chemistry/blob/main/colab/upload_prep.ipynb)
The upload form is at https://fragalysis.diamond.ac.uk/viewer/upload_cset/
> :construction: The upload code exists in Fragalysis-API old version. Hunt it down.

```python
from fragment_elaboration_scripts.prep_fragalysis import prep, generate_header
header = generate_header(method='Foo',
                         ref_url='https://www.example.com',
                         submitter_name='unknown',
                         submitter_email='a@b.c',
                         submitter_institution='Nowehere',
                         generation_date='2012-12-12',
                         smiles='CN1C=NC2=C1C(=O)N(C(=O)N2C)C',
                         extras={'∆∆G': 'W. Josiah Gibbs'})
prep(df, header, 'mol', 'id', 'foo.sdf', 'x1234', 'x1234', ['∆∆G'])
```
Alternatively
```python
from gist_import import GistImporter
from types import ModuleType

fu: ModuleType = GistImporter.from_github('https://raw.githubusercontent.com/matteoferla/Fragment-hit-follow-up-chemistry/main/fragment_elaboration_scripts/prep_fragalysis.py').to_module()

header: Chem.Mol = fu.generate_header(...)  # see below
fu.prep(header, ...)   # see below
```


The header molecule has a bunch of compulsory fields. The extras are optional, but will dictate the table.

The last step is to make an SDF file with the header and the molecules.
In truth, most often doing it manually is needed.
The molecules have to have 'ref_pdb' and 'ref_mols' and 'original SMILES' properties set,
along with whatever extra was in the header.

The major drama happens with covalents: the star dummy is not accepted.
Due to stupid reasons, Xe was adopted in the fragment network. But this is not needed.

```python
from rdkit import Chem
from rkdik.Chem import AllChem
from fragment_elaboration_scripts.prep_fragalysis import DummyMasker  # from rdkit_to_params.utils

with Chem.SDWriter('molecules-for-upload.sdf') as writer:
    writer.write(header)
    writer.SetProps(['pdb_code', 'description', 'ref_pdb', 'ref_mols', 'original SMILES'])
    for mol in mols:
        mol.SetProp('ref_pdb', 'x0404_0B')
        mol.SetProp('ref_mols', 'x0404_0B')
        with DummyMasker(mol, placekeeper_zahl=16) as mask:  # Sulfur is 16. Xenon is 54
            try:
                AllChem.SanitizeMol(mask.mol)
                mask.mol.SetProp('original SMILES', Chem.MolToSmiles(mask.mol))
                writer.write(mask.mol)
            except Exception as error:
                print(mol.GetProp('_Name'), error.__class__.__name__, str(error))
```