# Scripts

These scripts in `fragment_elaboration_scripts`.
Some contain in their docstring a description of the script and a yaml block for installation with conda if needed.

I have used somely via Python in Jupyter notebooks, I have converted them to CLI scripts for the convenience of others.
But I have not tested them all.

## Installation (without dependencies)
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

url = 'https://raw.githubusercontent.com/matteoferla/Fragment-hit-follow-up-chemistry/main/fragment_elaboration_scripts/fragment.py'
frag: ModuleType = GistImporter.from_github(url).to_module()
mols: List[Chem.Mol] = frag.fragment(...)
mols: List[Chem.Mol] = frag.remove_duplicated(...)
```

## Fragalysis

See [fragalysis_targetting_scripts.md](fragalysis_targetting_scripts.md).

## Utilities

See [utilities_scripts.md](utilities_scripts.md).

## Elaboration pipelines

See [elaboration_scripts.md](elaboration_scripts.md).





### ...

> :construction: ...

