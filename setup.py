from setuptools import setup, find_packages

setup(
    name='fragment_elaboration_scripts',
    version='0.1.1',
    packages=find_packages(),
    url='https://github.com/matteoferla/Fragment-hit-follow-up-chemistry',
    license='MIT',
    author='Matteo Ferla',
    author_email='matteo.ferla@stats.ox.ac.uk',
    description='',
    classifiers=[  # https://pypi.org/classifiers/
        'Development Status :: 4 - Beta',  # Development Status :: 5 - Production/Stable
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Chemistry',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
    ],
    package_data={'fragment_elaboration_scripts': ['combined-XChem-libraries.csv']},
    entry_points={
        'console_scripts': ['zinc-data=fragment_elaboration_scripts.zinc_data:main',
                            'enamine-store=fragment_elaboration_scripts.enamine_store:main',
                            'enamine-catalog-download=fragment_elaboration_scripts.enamine_catalog_download:main',
                            'retrieve-PDB-ligands=fragment_elaboration_scripts.retrieve_PDB_ligands:main',
                            'fragment=fragment_elaboration_scripts.fragment:main',
                            'pyrosetta-minimize=fragment_elaboration_scripts.pyrosetta_min:main',
                            'oeb_to_sdf=fragment_elaboration_scripts.oe_conformer_gen:oeb_to_sdf',
                            ],
    },
    install_requires=['requests', 'beautifulsoup4', 'pandas', 'rdkit', 'biopython'],
    extras_require={},
)
