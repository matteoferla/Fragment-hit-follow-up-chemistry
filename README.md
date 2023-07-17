# Fragment hit follow chemistry collection

> :construction: This repository is under construction :construction:

A collection of notebooks and scripts for the prediction of follow-up compounds 
I used for follow-up chemistry for XChem targets,
_ie._ small fragments identified crystallographically.

It is not intended as an exhaustive repository of all possible methods.

Furthermore, this is narrowly focused on the hit discovery stage of drug discovery,
before any activity data can be determined and any that does exist are deep in noise
(so no talk of QSAR, activity cliffs, ADMET etc etc).

For now see: https://github.com/matteoferla/EV-D68-3C-protease

## XChem libraries

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
RowGetter(library_collection, 'molecule').get_matching_hits(mols)
```

## Fragment theft

> See [hit-theft notebook](notebooks/hit-theft.ipynb)

Run a Blastp search against PDB entries in NCBI to find homologues and then steal their ligands.
Code originally adapted from
[colab-pyrosetta-migrate_ligands.ipynb](https://github.com/matteoferla/pyrosetta-help/blob/main/colab_notebooks/colab-pyrosetta-migrate_ligands.ipynb)

For a chemical component kill list see
[relevant blog post](https://blog.matteoferla.com/2019/11/go-away-glycerol.html)

## Fragment extraction

:hammer: To polish

**Any script here present is for fragment cleanup is not the official Fragalysis way.**

But briefly, the problem is as follows:

* A compound may be racemic but only one isomer binds. It will be split later on.
* A unique compound may bind more than once. It will have a single XChem code.
* A PDB structure may be a dimer in the crystallographic asymmetric unit, but a monomer in vivo. A compound may bind in both or in a single chain. Each chain will have a different PDB.
* A compound may have ambiguous conformers (think glutamine and histidine orientations), the PDBs will have partial density and alt IDs.

## Fragment PLIP

:hammer: To polish

## Fragmenstein

> For more [:octocat::Fragmenstein](https://github.com/matteoferla/Fragmenstein)

...
## Fragment network

:hammer: To polish

## Arthorian quest

> See [:octocat::Arthorian quest](https://github.com/matteoferla/arthorian-quest)
> (:construction: under construction :construction:)

...

## Filtering

:hammer: To polish

## CoPriNet

> See [:octocat::CoPriNet](https://github.com/oxpig/CoPriNet)

Predict the price for a compound.

...

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

### Fragalysis Upload

This is just to prep.
The upload form is at https://fragalysis.diamond.ac.uk/viewer/upload_cset/
> :warning: The upload code exists in Fragalysis-API old version. Hunt it down.

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

## Michelanglo

Upload to Fragalysis is not always okay.

## Minor

### Zinc data

`ZincInformer` in [zinc_data.py](followup/zinc_data.py) is a class to get data from a series of ZINC IDs.

