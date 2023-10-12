# Fragment hit follow chemistry collection

> :construction: This repository is under construction :construction:

> Herein sections with the icon :construction: are incomplete.

> In notebooks/unsorted there are notebooks that need annotating.


A collection of notebooks and scripts for the prediction of follow-up compounds 
I used for follow-up chemistry for XChem targets,
_ie._ small fragments identified crystallographically.

It is not intended as an exhaustive repository of all possible methods.
For now I am adding my code (:construction:).

Furthermore, this is narrowly focused on the hit discovery stage of drug discovery,
before any activity data can be determined and any that does exist are deep in noise
(so no talk of QSAR, activity cliffs, ADMET etc etc).

For now see: [:octocat:: EV-D68-3C-protease](https://github.com/matteoferla/EV-D68-3C-protease)

For the questionnaire for now target analysis MF uses see [questionnaire.md](questionnaire.md).

## XChem compound libraries

### Dataset generation

[combined XChem libraries](combined-XChem-libraries.csv) is a CSV file of the compounds in the XChem libraries.

The code to generate it is [notebooks/combine_libraries.ipynb](notebooks/combine_libraries.ipynb).

NB. Some compounds are in two or more libraries. Each row is a library entry. See `also_in` column for duplicates.

### Lookup

Given a list of compounds, get its source library. This is not an inchikey lookup as bond order might be wrong.
So silly amount of pain.

```python
from followup.row_getter import RowGetter

library_collection = pd.read_csv('combined-XChem-libraries.csv', index_col=0)
PandasTools.AddMoleculeColumnToFrame(library_collection, 'SMILES', 'molecule', True)
matches: pd.DataFrame = RowGetter(library_collection, 'molecule').get_matching_hits(mols)
```

> :construction: CHECK there are two versions of this script. There is one that returns a single row per hit. This is the other?

### Comparison

A very common discussion point is if this compound was a hit,
what other compounds in the library are similar?

> :construction: Hunt down

## XChem public targets

In [:octocat:: munged-Fragalysis-targets](https://github.com/matteoferla/munged-Fragalysis-targets)
the Fragalysis targets are explored by Blast to Xref with Uniprot,
yielding a [table of useful information](https://github.com/matteoferla/munged-Fragalysis-targets/blob/main/targets.csv).

NB. Not all data in Fragalysis has the same layout.

NB. The table is not polished by species popularity: a few protein are identical in some close species,
and the wrong one is picked —nobody, not even chimpanzees, cares for NUDT4 in chimpanzees.

To get the email address of a researcher involved in a target, asking is the best way.

> :construction: Rerun with newest targets

## Fragment theft

> See [hit-theft notebook](notebooks/hit-theft.ipynb)

Run a Blastp search against PDB entries in NCBI to find homologues and then steal their ligands,
by aligning in PyMOL against a reference structure.

Code originally adapted from
[colab-pyrosetta-migrate_ligands.ipynb](https://github.com/matteoferla/pyrosetta-help/blob/main/colab_notebooks/colab-pyrosetta-migrate_ligands.ipynb)
(modded to use PyMOL over PyRosetta).

For a chemical component kill list see
[relevant blog post](https://blog.matteoferla.com/2019/11/go-away-glycerol.html)

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
5. Each will be its own structure, but not in legacy data

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

## Decomposition

BRICS is good, but has two problems for my use.
It does not keep the metadata and it does not split fused/spiro rings, which is understandable,
but in this case I want to split them.

In [followup.fragment](followup/fragment.py) there is a function (`fragmént`) to fragment a molecule into frágments.

```python
from followup.fragment import fragmént, remove_duplicated

fhits: List[Chem.Mol] = []
for hit in hits:
    fhits.extend( fragmént(hit, minFragmentSize=7, fused_splitting=True) )

fhits = remove_duplicated(fhits)
```

Using the remote version:

```python
fmodule = GistImporter.from_github('https://raw.githubusercontent.com/matteoferla/Fragment-hit-follow-up-chemistry/main/followup/fragment.py')
remove_duplicated = fmodule['remove_duplicated']
fragmént = fmodule['fragmént']
```

## PLIP

PLIP extract interactions from PDB files.

NB. `intxn` is non-standard and unclear to many, please don't copy.

In [followup.plip](followup/plip.py) there is a class to run PLIP on a PDB file
and converts the results into a uniform format.

```python
plipper = SerialPLIPper(pdb_block,
                        resn=lig_resn,
                        resi=lig_resi,
                        chain=lig_chain)
# list of dictionaries with interactions spelled out.
details = plipper.summarize_interactions(atom_names)

# or dictionary of interaction counts for a given mol (assuming apo was provided)
plipper = SerialPLIPper(apo_block)
details = plipper(mol)
```

A few things to note: the atom names are lost for the non-canonical ligand,
hence why they are passed.

The bond order and protonation is guessed with no option for SMILES template.
As a result in [plip_initials notebook](notebooks/plip_initials.ipynb),
PyRosetta is used to generate the protonated form with CONECT records.
allowing the happy usage of:

```python
from plip.basic import config

config.NOHYDRO = True
```

See also https://github.com/matteoferla/PLIP-PyRosetta-hotspots-test
and https://blog.matteoferla.com/2023/07/a-note-no-plip-interactions.html

## Fragmenstein

In [followup.fragmenstein_merge_sw_place](followup/fragmenstein_merge_sw_place.py)
is a script to run Fragmenstein on a PDB file and SDF hits, merge them,
find analogues via SmallWorld and place them in the PDB.


> For more [:octocat::Fragmenstein](https://github.com/matteoferla/Fragmenstein)

> :construction: Migrate origins is being migrated to Fragmenstein repo.

## Fragment network

> :construction: Clean-up

## Arthorian quest

> See [:octocat::Arthorian quest](https://github.com/matteoferla/arthorian-quest)

> :construction: Repo under construction

## Filtering

> :construction: Clean-up

## CoPriNet

> See [:octocat::CoPriNet](https://github.com/oxpig/CoPriNet)

Predict the price for a compound.

> :construction: hunt down API code

## Fragalysis

Typical folder format in Fragalysis is:

* smiles.smi — optional
* metadata.csv — optional
* reference.pdb — optional
* extra_files/* — optional
* aligned/👻👻-👾👾/👾👾_bound.pdb.pdb

### Fauxalysis

The latter can be replicated with `make_fauxalysis` from [followup.fauxalysis](followup/fauxalysis.py)

The dummy Fragalysis download structure is needed for a few things.

### Fragalysis Download

In [followup.fragalysis_download](followup/fragalysis_download.py) there is a script to download a Fragalysis target.

```python
from followup/fragalysis_download import QuickDownloader
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

### Fragalysis Upload

This is just to prep.
The upload form is at https://fragalysis.diamond.ac.uk/viewer/upload_cset/
> :construction: The upload code exists in Fragalysis-API old version. Hunt it down.

```python
from followup.prep_fragalysis import prep, generate_header
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
from followup.prep_fragalysis import DummyMasker  # from rdkit_to_params.utils

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

## Michelanglo

Upload to Fragalysis is not always okay.

> :construction: Clean-up

## OE ROCS

This task has three parts:

1. Subsetting of Enamine BB, see below.
2. Conformer generation with OE Omega.
3. FastROCS with the conformers, or OE ROCS GUI derp around.

In [:octocat:: Functional-subsets-of-Enamine-BB](https://github.com/matteoferla/Functional-subsets-of-Enamine-BB)
Enamine BB is split into functional subsets.
Due to size, the OE Omega conformers are not included.

> :construction: Copy 2 & 3 off cluster

> :construction: Hunt down the fluorescent subset!

## MMCIF annotation

> :construction: Hunt down and combine

## Minor

### Zinc data

`ZincInformer` in [zinc_data.py](followup/zinc_data.py) is a class to get data from a series of ZINC IDs.

### Enamine REAL DB download

The download of the Enamine REAL DB was done using [this script](https://gist.github.com/matteoferla/b1eee8656079d006835f2d8dc159fbb5)
requires registration.
The other databases in OPIG bitbucket are from different sources: see README.md files therein!

### ...

> :construction: ...