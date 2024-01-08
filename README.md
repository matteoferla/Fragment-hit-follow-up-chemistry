# Scripts and pipelines to elaborate fragment hits

> :construction: This repository is under construction :construction:

> Herein sections with the icon :construction: are incomplete.

> In notebooks/unsorted there are notebooks that need annotating.


A collection of notebooks and scripts to elaborate fragment hits.
I used for elaboration chemistry for XChem targets,
_ie._ small fragments identified crystallographically.

It is not intended as an exhaustive repository of all possible methods.
For now I am adding my code (:construction:).

Furthermore, this is narrowly focused on the hit discovery stage of drug discovery,
before any activity data can be determined and any that does exist are deep in noise
(so no talk of QSAR, activity cliffs, ADMET etc etc). Many scripts will be Diamond Light Source XChem specific.

For now see: [:octocat:: EV D68-3C protease](https://github.com/matteoferla/EV-D68-3C-protease) or
[:octocat:: EV A71-2A protease](https://github.com/matteoferla/EV-A71-2A-elaborations).

## Colab notebooks

Given a Uniprot ID make a sequence logo of the cleavage sites (viral polyprotein)
[![colab](https://img.shields.io/badge/Run_cleavege-logo.ipynb-f9ab00?logo=googlecolab)](https://colab.research.google.com/github/matteoferla/Fragment-hit-follow-chemistry/blob/main/colab/make_cleavage_logo.ipynb)

Given a Fragalysis valid upload set, fix it to pick the best Fragalysis template and do minor tweaks to make it fit there.
I.e. pretend an induced fit compound binds by lock and key.
[![colab](https://img.shields.io/badge/Run_cleavege-logo.ipynb-f9ab00?logo=googlecolab)](https://colab.research.google.com/github/matteoferla/Fragment-hit-follow-up-chemistry/blob/main/colab/upload_fix.ipynb).

For the preparation of an upload file
[![colab](https://img.shields.io/badge/Run_cleavege-logo.ipynb-f9ab00?logo=googlecolab)](https://colab.research.google.com/github/matteoferla/Fragment-hit-follow-up-chemistry/blob/main/colab/upload_prep.ipynb).

## Scripts
These scripts in `fragment_elaboration_scripts`.
Some contain in their docstring a description of the script and a yaml block for installation with conda if needed.

I have used somely via Python in Jupyter notebooks, I have converted them to CLI scripts for the convenience of others.
But I have not tested them all.

### Installation (without dependencies)
These can be installed with pip:

```bash
pip install fragment_elaboration_scripts
```

or individually fetched from the repo (which means it will be the latest):

```bash
pip install gist_import
```

For example:
```python
from gist_import import GistImporter
from types import ModuleType
from typing import List
from rdkit import Chem

url = 'https://raw.githubusercontent.com/matteoferla/Fragment-hit-follow-up-chemistry/main/fragment_elaboration_scripts/fragment.py'
frag: ModuleType = GistImporter.from_github(url).to_module()
mols: List[Chem.Mol] = frag.fragment(...)
mols: List[Chem.Mol] = frag.remove_duplicated(...)
```

### Fragalysis

See [fragalysis_targetting_scripts.md](fragalysis_targetting_scripts.md).

### Utilities

See [utilities_scripts.md](utilities_scripts.md).

### Elaboration pipelines

See [elaboration_scripts.md](elaboration_scripts.md).

## Other

For the questionnaire given to user to direct which target analyses MF uses see [questionnaire.md](questionnaire.md).

### Dataset generation

[combined XChem libraries](combined-XChem-libraries.csv) is a CSV file of the compounds in the XChem libraries.

The raw code to generate it is [notebooks/combine_libraries.ipynb](notebooks/combine_libraries.ipynb).

NB. Some compounds are in two or more libraries. Each row is a library entry. See `also_in` column for duplicates.