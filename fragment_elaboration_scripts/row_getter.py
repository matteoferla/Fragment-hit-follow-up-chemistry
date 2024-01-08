import pandas as pd
import requests
from typing import List, Union
from rdkit import Chem
from rdkit.Chem import AllChem

options = []

class RowGetter:
    """
    Given a pandas table with a column w/ a molecule,
    find the _first_ match of a given molecule.

    This is not inchi based, which may be silly
    Furthermore, what if a molecules is from different libraries?

    .. codeblock::python
        library_collection = pd.read_csv('../XChem-targets/combined_libraries.csv', index_col=0)
        PandasTools.AddMoleculeColumnToFrame(library_collection, 'SMILES', 'molecule', True)
        RowGetter(library_collection, 'molecule').get_matching_hits(mols)
    """
    def __init__(self, df: pd.DataFrame, mol_col: str):
        self.df = df
        self.mol_col = mol_col
        self.df['disordered'] = self.df[self.mol_col].apply(self.disorder)

    def disorder(self, mol: Chem.Mol) -> Chem.Mol:
        """
        As in, lacks bond order
        """
        disorded = Chem.Mol(mol)
        for bond in disorded.GetBonds():
            bond.SetBondType(Chem.BondType.SINGLE)
        return disorded


    def __call__(self, hit: Chem.Mol) -> Union[None, pd.Series]:
        hit = AllChem.RemoveAllHs(hit)
        m = (self.df[self.mol_col].apply(lambda m: hit.GetNumHeavyAtoms() == m.GetNumHeavyAtoms())) & (self.df[self.mol_col] >= hit)
        if sum(m) == 0:
            m = (self.df[self.mol_col].apply(lambda m: hit.GetNumHeavyAtoms() == m.GetNumHeavyAtoms())) & (self.df.disordered >= hit)
        if sum(m) == 0:
            return None
        elif sum(m) == 1:
            return self.df.loc[m].reset_index().iloc[0]
        else:
            matches = self.df[m].reset_index()
            return pd.Series(dict(Id='Ambigous:' + '|'.join(set(matches.Id.to_list())),
                        SMILES=matches.iloc[0].SMILES,
                        molecule=matches.iloc[0].molecule,
                        library='Ambigous:' + '|'.join(set(matches.library.to_list())),
                        Name=matches.iloc[0].Name,
                        ))

    def get_matching_hits(self, hits: List[Chem.Mol]) -> pd.DataFrame:
        rows = []
        for hit in hits:
            row = self(hit)
            if row is None:
                rows.append(pd.Series({'3D_molecule': hit, 'in_library': False, 'query': hit, 'library': 'Unknown'}))
            else:
                row['query'] = hit
                row['3D_molecule'] = hit
                row['in_library'] = True
                rows.append(row)
        return pd.DataFrame(rows)

    @classmethod
    def smiles_to_iupac(cls, smiles: str, raise_error:bool=True) -> str:
        """
        Given a SMILES, get the IUPAC name from NIH Cactus server
        """
        rep = "iupac_name"
        CACTUS = "https://cactus.nci.nih.gov/chemical/structure/{0}/{1}"
        url = CACTUS.format(smiles, rep)
        response = requests.get(url)
        if response.ok:
            return response.text
        elif raise_error:
            response.raise_for_status()
        else:
            return ''
