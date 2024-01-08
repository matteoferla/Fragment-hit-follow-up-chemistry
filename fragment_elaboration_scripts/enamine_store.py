"""
This uses the new Enamine Store which does not require authentication.
Nor does the old one work anymore.

The new one is a REST API, but the documentation is not public.
This is meant to be used by the browser. Safeguards are in place and will not work on residential IPs.
I set a sleep of 30 seconds between calls to be safe.

.. code-block:: python

    from fragment_elaboration_scripts.enamine_store import search, StoreCatalog, StoreCurrency
    get_price(enamine_code, catalogue=StoreCatalog.REALDB, currency=StoreCurrency.EUR)

.. code-block:: python

    import pandas as pd
    from rdkit import Chem
    from fragment_elaboration_scripts.enamine_store import search, StoreCatalog, StoreTypes, StoreSSTypes
    mol: Chem.Mol
    df: pd.DataFrame = search(mol, catalogue=StoreCatalog.REALDB, search_type=StoreTypes.SMARTS)

------------
CLI Usage
------------

.. code-block:: bash

    $ enamine-store price Z12345678 --catalogue REALDB --currency EUR

or

.. code-block:: bash

    $ enamine-store search Cn1cnc2c1c(=O)n(C)c(=O)n2C > caffeine_analogues.csv

-----------
Installation
-----------

No special installation requirements beyond ``pip install fragment_elaboration_scripts``
"""

import enum, requests
from rdkit import Chem
from typing import Union
import pandas as pd


class StoreTypes(enum.Enum):
    """
    EnamineStore types of searches: ID, CAS, MFCD, SMARTS
    """
    ID = enum.auto()
    CAS = enum.auto()
    MFCD = enum.auto()
    SMARTS = enum.auto()


class StoreSSTypes(enum.Enum):
    """
    EnamineStore types of structural searches: EXACT, SUB, SIM
    """
    EXACT = enum.auto()
    SUB = enum.auto()
    SIM = enum.auto()


class StoreCatalog(enum.Enum):
    """
    EnamineStore catalogues: SCR, BB, REALDB, MADE, EBC
    """
    SCR = enum.auto()
    BB = enum.auto()
    REALDB = enum.auto()
    MADE = enum.auto()
    EBC = enum.auto()


class StoreCurrency(enum.Enum):
    """
    EnamineStore currencies: USD, EUR

    No other currencies are supported.
    """
    USD = enum.auto()
    EUR = enum.auto()

def search(mol_or_smarts: Union[Chem.Mol, str],
           catalogue: StoreCatalog = StoreCatalog.REALDB,
           search_type: StoreTypes = StoreTypes.SMARTS,
           structural_type: StoreSSTypes = StoreSSTypes.SIM,
           size: int = 100,
           sim: float = 0.1) -> pd.DataFrame:
    """
    Search the Enamine Store.

    :param mol_or_smarts: rdkit Mol or SMARTS string
    :param catalogue: see StoreCatalog enum. Options: SCR, BB, REALDB, MADE, EBC
    :param search_type: see StoreTypes enum. Options: ID, CAS, MFCD, SMARTS
    :param structural_type: see StoreSSTypes enum. Options: EXACT, SUB, SIM
    :param size: number of results to return
    :param sim: similarity threshold
    :return: pandas DataFrame
    """
    if isinstance(mol_or_smarts, Chem.Mol):
        smarts = Chem.MolToSmarts(mol_or_smarts)
    else:
        smarts = str(mol_or_smarts)
        assert Chem.MolFromSmarts(smarts) is not None, f'Could not parse SMARTS {smarts}'
    # spoof as a browser
    headers = {'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/102.0.0.0 Safari/537.36'}
    response = requests.get(f'https://new.enaminestore.com/api/v1/catalog/',
                            dict(q=smarts,
                                 cat=catalogue.name if isinstance(catalogue, StoreCatalog) else str(catalogue),
                                 type=search_type.name if isinstance(search_type, StoreTypes) else str(search_type),
                                 sstype=structural_type.name if isinstance(structural_type, StoreSSTypes) else str(
                                     structural_type),
                                 curPage=1,  # does nothing
                                 pageSize=int(size),  # does nothing
                                 sim=max(float(sim), 0.01),  # does nothing
                                 ),
                            headers=headers
                            )
    response.raise_for_status()
    if response.text[0] == '<':
        raise ValueError(f'Incorrect response (not JSON): {response.text}... Did you pass a SMILES not a SMARTS?')
    if 'searchResults' not in response.json():
        raise ValueError(f'Incorrect query: {response.text}')
    return pd.DataFrame(response.json()['searchResults'])


def get_price(enamine_code: str,
              catalogue: StoreCatalog = StoreCatalog.REALDB,
              currency: StoreCurrency = StoreCurrency.EUR) -> float:
    """
    Get the price of a compound from the Enamine Store.
    
    .. code-block:: python

        price: float = get_price(enamine_code, catalogue=StoreCatalog.REALDB, currency=StoreCurrency.EUR)

    :param enamine_code: str. Enamine code, e.g. Z12345678
    :param catalogue: see StoreCatalog enum. Options: SCR, BB, REALDB, MADE, EBC
    :param currency: see StoreCurrency enum. Options: USD, EUR
    :return: price for 1 mg (fudged math for higher amounts)
    """
    db_name = catalogue.name if isinstance(catalogue, StoreCatalog) else str(catalogue)
    cur_name = currency.name if isinstance(currency, StoreCurrency) else str(currency)
    enamine_code = str(enamine_code).strip()
    headers = {
        'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/102.0.0.0 Safari/537.36'}

    response = requests.get(
        f'https://new.enaminestore.com/api/v1/catalog/price?id={enamine_code}&cat={db_name}&cur={cur_name}',
        headers=headers
    )
    data = response.json()
    for datum in data['samples']:
        if datum['amount'] != 1.0 or datum['measure'] != 'mg':
            continue
        return datum['price']
    if len(data['samples']) == 0:
        return 0.
    datum = data['samples'][0]
    if datum['measure'] == 'mg':
        return data['samples'][0]['price'] / datum['amount']
    elif datum['measure'] == 'g':
        return data['samples'][0]['price'] / datum['amount'] / 1_000
    elif datum['measure'] == 'kg':
        return data['samples'][0]['price'] / datum['amount'] / 1_000_000
    else:
        return 0.

def main():
    import argparse
    import time
    import pandas as pd

    parser = argparse.ArgumentParser(description='Search the Enamine Store')
    subparsers = parser.add_subparsers(help='search or price')

    def parse_price(args: argparse.Namespace):
        print(get_price(args.enamine_code, catalogue=args.catalogue, currency=args.currency))

    price_parser = subparsers.add_parser('price', help='get price')
    price_parser.add_argument('enamine_code', type=str, help='Enamine code')
    price_parser.add_argument('--catalogue', type=StoreCatalog, default=StoreCatalog.REALDB, help='Catalogue')
    price_parser.add_argument('--currency', type=StoreCurrency, default=StoreCurrency.EUR, help='Currency')
    price_parser.set_defaults(func=parse_price)

    def parse_search(args: argparse.Namespace):
        df = search(args.smiles, catalogue=args.catalogue, search_type=args.search_type,
                    structural_type=args.structural_type, size=args.size, sim=args.sim)
        print(df.to_csv(index=False))

    search_parser = subparsers.add_parser('search', help='search')
    search_parser.add_argument('smiles', type=str, help='SMILES string')
    search_parser.add_argument('--catalogue', type=StoreCatalog, default=StoreCatalog.REALDB, help='Catalogue')
    search_parser.add_argument('--search_type', type=StoreTypes, default=StoreTypes.SMARTS, help='Search type')
    search_parser.add_argument('--structural_type', type=StoreSSTypes, default=StoreSSTypes.SIM, help='Structural type')
    search_parser.add_argument('--size', type=int, default=100, help='Size')
    search_parser.add_argument('--sim', type=float, default=0.1, help='Similarity')
    search_parser.set_defaults(func=parse_search)

    args = parser.parse_args()
    args.func(args)
