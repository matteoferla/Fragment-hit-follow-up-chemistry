import requests

class DownloadEnamine:
    """
    This class is set up to download Enamine REAL database on a remote machine.
    Instatiation requires plain ``username`` and ``password``.

    .. code-block::python
        de = DownloadEnamine('foo.bar@baz.ac.uk', 'Foo123')
        de.download_all('REAL')

    Note, this is copied off the route of the web page and not the Enamine Store API.
    Plus the official documentation (emailed Word document) is for the old Store and
    no longer applies anyway (plain text username and password in GET header "Authorization").

    The URLs pointing to the download pages were copied off manually.
    """
    REAL = [
        'Enamine_REAL_HAC_6_21_420M_CXSMILES.cxsmiles.bz2',
        'Enamine_REAL_HAC_22_23_471M_CXSMILES.cxsmiles.bz2',
        'Enamine_REAL_HAC_24_394M_CXSMILES.cxsmiles.bz2',
        'Enamine_REAL_HAC_25_557M_CXSMILES.cxsmiles.bz2',
        'Enamine_REAL_HAC_26_833M_Part_1_CXSMILES.cxsmiles.bz2',
        'Enamine_REAL_HAC_26_833M_Part_2_CXSMILES.cxsmiles.bz2',
        'Enamine_REAL_HAC_27_1.1B_Part_1_CXSMILES.cxsmiles.bz2',
        'Enamine_REAL_HAC_27_1.1B_Part_2_CXSMILES.cxsmiles.bz2',
        'Enamine_REAL_HAC_28_1.2B_Part_1_CXSMILES.cxsmiles.bz2',
        'Enamine_REAL_HAC_28_1.2B_Part_2_CXSMILES.cxsmiles.bz2',
        'Enamine_REAL_HAC_29_38_988M_Part_1_CXSMILES.cxsmiles.bz2',
        'Enamine_REAL_HAC_29_38_988M_Part_2_CXSMILES.cxsmiles.bz2',
    ]

    def __init__(self, username, password):
        self.sesh = requests.Session()
        response = self.sesh.get('https://enamine.net/compound-collections/real-compounds/real-database',
                                 params={'username': username,
                                         'password': password,
                                         'Submit': 'Login',
                                         'remember': 'yes',
                                         'option': 'com_users',
                                         'task': 'user.login'})
        response.raise_for_status()

    def download_all(self, catalogue='REAL'):
        """
        The URLs of the databases files are in the class attribute of that same catalogue name (i.e. ``REAL``).
        """
        for filename in getattr(self, catalogue):
            self.download('REAL', filename)

    def check(self, catalogue='REAL'):
        for filename in getattr(self, catalogue):
            with self.sesh.get(f'https://ftp.enamine.net/download/{catalogue}/{filename}', stream=True) as r:
                r.raise_for_status()  # requests.exceptions.HTTPError
                for chunk in r.iter_content(chunk_size=8192):
                    break

    def download(self, catalogue, filename):
        """
        Downloads the ``filename`` of the given ``catalogue``
        """
        with self.sesh.get(f'https://ftp.enamine.net/download/{catalogue}/{filename}', stream=True) as r:
            r.raise_for_status()
            with open(filename, 'wb') as f:
                for chunk in r.iter_content(chunk_size=8192):
                    f.write(chunk)