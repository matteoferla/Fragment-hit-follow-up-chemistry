# Ultility scripts
For full details see [![Read the Docs](https://img.shields.io/readthedocs/fragment-hit-follow-up-chemistry)](https://fragment-hit-follow-up-chemistry.readthedocs.io/en/latest/fragment_elaboration_scripts.html)

## Fragment theft

Requires PyMOL. See [retrieve_PDB_ligands.py](fragment_elaboration_scripts/retrieve_PDB_ligands.py) for yaml.
`conda env update -n thief --file thief.yml` or `conda install -c schrodinger pymol-open-source`


CLI: `python retrieve_PDB_ligands.py --reference template.pdb --sdf compounds.sdf`


> See [hit-theft notebook](notebooks/hit-theft.ipynb)

Run a Blastp search against PDB entries in NCBI to find homologues and then steal their ligands,
by aligning in PyMOL against a reference structure.

Code originally adapted from
[colab-pyrosetta-migrate_ligands.ipynb](https://github.com/matteoferla/pyrosetta-help/blob/main/colab_notebooks/colab-pyrosetta-migrate_ligands.ipynb)
(modded to use PyMOL over PyRosetta).

For a chemical component kill list see
[relevant blog post](https://blog.matteoferla.com/2019/11/go-away-glycerol.html)

## Query external databases

### Zinc data

`ZincInformer` in [zinc_data.py](fragment_elaboration_scripts/zinc_data.py) is a class to get data from a series of ZINC IDs.

CLI: `zinc-data ZINC00000001 ZINC00000002 ZINC00000003 -o zinc.csv`

### Enamine Store

The functions ``search`` and ``get_price`` in [enamine_store.py](fragment_elaboration_scripts/enamine_store.py) are used to search the Enamine REAL database.

CLI: `enamine-store price Z12345678 --catalogue REALDB --currency EUR`

### Enamine Catalogue Download

`DownloadEnamine` in [enamine_catalogue_download.py](fragment_elaboration_scripts/enamine_catalog_download.py)
is a class to download the Enamine catalogues on a remote server.

CLI: `enamine-catalogue-download`

## Lookup XChem screen collection

Given a list of compounds, get its source library. This is not an inchikey lookup as bond order might be wrong.
So silly amount of pain.

```python
from fragment_elaboration_scripts.row_getter import RowGetter
import pandas as pd
from rdkit import Chem
from rdkit.Chem import PandasTools
from typing import List

mols: List[Chem.Mol] = ...
library_collection = pd.read_csv('combined-XChem-libraries.csv', index_col=0)
PandasTools.AddMoleculeColumnToFrame(library_collection, 'SMILES', 'molecule', True)
matches: pd.DataFrame = RowGetter(library_collection, 'molecule').get_matching_hits(mols)
```

> :construction: CHECK there are two versions of this script. There is one that returns a single row per hit. This is the other?

## Pyrosetta minimise

Requires PyRosetta. For installation, see official documentation or 
for remote clusters and avoiding conda [pyrosetta_help](https://github.com/matteoferla/pyrosetta-help).
```bash
pip install pyrosetta-help
install_pyrosetta -u 👾👾👾 -p 👾👾👾
```

See [pyrosetta_minimise.py](fragment_elaboration_scripts/pyrosetta_minimise.py) for functions.

As a CLI script it does not make too much sense as Rosetta directly is CLI based,
and this is just fast relax with constraints.
The benefit comes from tinkering it in a notebook...
```bash
pyrosetta-minimize --target template.pdb --constraint_weight 10 --cycles 15 --output minimised.pdb
```

## Get name via Cactus

`smiles_to_iupac`

## Comparison

A very common discussion point is if this compound was a hit,
what other compounds in the library are similar?

> :construction: Hunt down the code

## MMCIF annotation

> :construction: Hunt down and combine

## Decomposition

In [fragment_elaboration_scripts.fragment](fragment_elaboration_scripts/fragment.py) there is a function (`fragmént`) to fragment a molecule into frágments.

```python
from fragment_elaboration_scripts.fragment import fragmént, remove_duplicated
from rdkit import Chem
from typing import List

hits: List[Chem.Mol] = ...
fhits: List[Chem.Mol] = [frag for hit in hits for frag in fragmént(hit, minFragmentSize=7, fused_splitting=True)]
fhits = remove_duplicated(fhits)
```

CLI: `fragment compounds.sdf -o fragments.sdf`

### PLIP

PLIP extract interactions from PDB files.

NB. `intxn` is non-standard and unclear to many, please don't copy.

In [fragment_elaboration_scripts.plip](fragment_elaboration_scripts/plip.py) there is a class to run PLIP on a PDB file
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

Note that the hydrogens need to be stripped from the PDB as there is an issue with PLIP and hydrogens:

```python
from plip.basic import config

config.NOHYDRO = True
```

See also https://github.com/matteoferla/PLIP-PyRosetta-hotspots-test
and https://blog.matteoferla.com/2023/07/a-note-no-plip-interactions.html

## Michelanglo

Upload to Fragalysis is not always okay.

> :construction: Clean-up