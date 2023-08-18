"""
This uses the new Enamine Store which does not require authentication.
Nor does the old one work anymore.
"""

import enum, requests
from rdkit import Chem
from typing import Union
import pandas as pd

class StoreTypes(enum.Enum):
    ID = enum.auto()
    CAS = enum.auto()
    MFCD = enum.auto()
    SMARTS = enum.auto()

class StoreSSTypes(enum.Enum):
    EXACT = enum.auto()
    SUB = enum.auto()
    SIM = enum.auto()

class StoreCatalogues(enum.Enum):
    SCR = enum.auto()
    BB = enum.auto()
    REALDB = enum.auto()
    MADE = enum.auto()
    EBC = enum.auto()

class StoreCurrency(enum.Enum):
    USD = enum.auto()
    EUR = enum.auto()

def search(mol_or_smarts: Union[Chem.Mol, str],
           catalogue: StoreCatalogues = StoreCatalogues.REALDB,
           search_type: StoreTypes = StoreTypes.SMARTS,
           structural_type: StoreSSTypes = StoreSSTypes.SIM,
           size: int = 100,
           sim: float = 0.1) -> pd.DataFrame:
    """
    Search the Enamine Store.

    :param mol_or_smarts:
    :param catalogue:
    :param search_type:
    :param structural_type:
    :param size:
    :param sim:
    :return:
    """
    smarts = Chem.MolToSmiles(mol_or_smarts) if isinstance(mol_or_smarts, Chem.Mol) else str(mol_or_smarts)
    response = requests.get(f'https://new.enaminestore.com/api/v1/catalog/',
                            dict(q=smarts,
                                 cat=catalogue.name if isinstance(catalogue, StoreCatalogues) else str(catalogue),
                                 type=search_type.name if isinstance(search_type, StoreTypes) else str(search_type),
                                 sstype=structural_type.name if isinstance(structural_type, StoreSSTypes) else str(structural_type),
                                 curPage=1, # does nothing
                                 pageSize=int(size),   # does nothing
                                 sim=max(float(sim), 0.01), # does nothing
                                 )
                            )
    response.raise_for_status()
    if 'searchResults' not in response.json():
        raise ValueError(f'Incorrect query: {response.text}')
    return pd.DataFrame(response.json()['searchResults'])

def get_price(enamine_code: str,
              catalogue:StoreCatalogues = StoreCatalogues.REALDB,
              currency: StoreCurrency=StoreCurrency.EUR) -> float:
    """
    Get the price of a compound from the Enamine Store.

    :param enamine_code:
    :param catalogue:
    :param currency:
    :return:
    """
    db_name = catalogue.name if isinstance(catalogue, StoreCatalogues) else str(catalogue)
    cur_name = currency.name if isinstance(currency, StoreCurrency) else str(currency)
    enamine_code = str(enamine_code).strip()
    response = requests.get(f'https://new.enaminestore.com/api/v1/catalog/price?id={enamine_code}&cat={db_name}&cur={cur_name}')
    data = response.json()
    for datum in data['samples']:
        if datum['amount'] != 1.0 or datum['measure'] != 'mg':
            continue
        return datum['price']
    return 0.
