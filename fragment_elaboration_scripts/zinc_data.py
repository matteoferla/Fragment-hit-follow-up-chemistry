"""
A very simple script to get info from Zinc.
Zinc, as far as I can tell, does not have an API,
so the HTML has to be scraped.

NB. I have not asked for permission to query ZINC programmatically

See ``ZincInformer`` for the class.

-------------
CLI Usage
-------------

.. code-block:: bash

    $ python -m zinc_data ZINC00000001 ZINC00000002 ZINC00000003 > zinc.json

or

.. code-block:: bash

    $ python -m zinc_data ZINC00000001 ZINC00000002 ZINC00000003 -o zinc.csv

-------------
Installation
-------------

No special installation requirements beyond ``pip install fragment_elaboration_scripts``
"""

import requests, collections, contextlib, json
from bs4 import BeautifulSoup


class ZincInformer(collections.abc.MutableMapping):
    """
    A simple class to get Zinc info.

    A class that stores the retieved values —in ``.data`` (``.dump`` and ``.load`` to store).
    The values can be accessed as a subscript or by calling the instance,
    the latter captures errors declared during initialisation by the argument ``suppressed_exception``.

    The instance is callable, so can be used in a ``pandas.Series.apply``:

    .. code-block::python

        zinfo: Callable = ZincInformer()
        data: pd.DataFrame = series.apply(zinfo)

    The data is stored in a dictionary, so can be dumped to a JSON file and loaded back in.
    This is useful for caching without repeating requests.

    .. code-block::python

        zinfo: Callable = ZincInformer()
        zinfo.load('zinc.json')
        data: Dict = zinfo['ZINC00000001']
        print(zinfo.data)
        zinfo.dump('zinc.json')

    The data is fetched via a call by ``get_soup`` and then parsed by ``get_zinc_info``,
    which calls ``get_dl`` and ``polísh``.
    """

    def __init__(self, suppressed_exception=Exception):
        self.data = {}
        self.suppressed_exception = suppressed_exception

    def __getitem__(self, zinc_id):
        if zinc_id not in self.data:
            soup = self.get_soup(zinc_id)
            self.data[zinc_id] = self.get_zinc_info(zinc_id, soup)
        return self.data[zinc_id]

    def __call__(self, zinc_id):
        with contextlib.suppress(self.suppressed_exception):
            return self[zinc_id]
        self[zinc_id] = {}
        return {}

    def __setitem__(self, zinc_id: str, info: dict):
        self.data[zinc_id] = info

    def __delitem__(self, zinc_id):
        del self.data[zinc_id]

    def __iter__(self):
        return iter(self.data)

    def __len__(self):
        return len(self.data)

    def dump(self, filename: str = 'zinc.json'):
        with open(filename, 'w') as fh:
            json.dump(self.data, fh)

    def load(self, filename: str = 'zinc.json'):
        with open(filename, 'r') as fh:
            self.data = json.load(fh)

    # ======== specific methods

    @classmethod
    def get_soup(self, zinc_id: str) -> BeautifulSoup:
        """parse HTML. Return a soup"""
        response = requests.get(f'https://zinc.docking.org/substances/{zinc_id}')
        response.raise_for_status()
        return BeautifulSoup(response.text, 'html.parser')

    @classmethod
    def get_dl(self, soup: BeautifulSoup) -> dict:
        """
        Data tables are organised in dl entries, dt headers and dd fields.
        """
        keys, values = [], []
        for dl in soup.findAll("dl"):
            for dt in dl.findAll("dt"):
                keys.append(dt.text.strip())
            for dd in dl.findAll("dd"):
                values.append(dd.text.strip())
        return dict(zip(keys, values))

    @classmethod
    def polísh(self, dl: dict) -> dict:
        """
        There's a bunch of words that get in the way...
        """
        return {k.replace('Bioactive', '').replace('Natural Products', '').replace('Building Blocks', '').strip(): v for
                k, v in dl.items()}

    @classmethod
    def get_zinc_info(self, zinc_id, soup):
        """
        These fields ought to always exist!
        """
        return {'query_name': zinc_id,
                'title': soup.title.text.strip(),
                'SMILES': soup.find('input', dict(id="substance-smiles-field")).attrs['value'].strip(),
                'inchi': soup.find('input', dict(id="substance-inchi-field")).attrs['value'].strip(),
                'inchikey': soup.find('input', dict(id="substance-inchikey-field")).attrs['value'].strip(),
                **self.polísh(self.get_dl(soup))
                }


# ----------------- CLI ---------------------------------------------------------------------------------------------

def main():
    import argparse
    import pandas as pd
    from typing import List, Dict

    parser = argparse.ArgumentParser(description='Get info from Zinc.')
    parser.add_argument('zinc_ids', nargs='+', help='Zinc IDs')
    parser.add_argument('-o', '--output', help='Output CSV file', default='')
    parser.add_argument('-c', '--cache', help='Saved cache file', default='')
    args = parser.parse_args()
    zinformer = ZincInformer()
    if args.cache:
        zinformer.load(args.cache)
    details: List[Dict] = []
    zinc_id: str
    for zinc_id in args.zinc_ids:
        info = {'id': zinc_id, **zinformer(zinc_id)}
        details.append(info)
    if args.cache:
        zinformer.dump(args.cache)
    if args.output:
        pd.DataFrame(details).to_csv(args.output, index=False)
    else:
        print(zinformer.data)

if __name__ == '__main__':
    main()
