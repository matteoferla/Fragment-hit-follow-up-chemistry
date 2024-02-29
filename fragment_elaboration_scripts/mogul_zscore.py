"""
Mogul (Cambridge Crystallographic Data Centre) gives a histogram of bond lengths, angles, torsions, and rings
of a molecule. This script gets the mean and max of these z-scores of these distributions.
Uncommon

## Usage:

.. code-block:: bash

        $ python mogul_zscore.py -i input.sdf -o output.json

## Installation

Requires the CCDC software and Python API to be installed.
Assuming in ``$HOME/CCDC`` and the conda environment is called ``CSD``:

.. code-block:: bash

    conda create -n CSD python=3.9
    conda activate CSD
    conda install --channel=https://conda.ccdc.cam.ac.uk csd-python-api
    export CCDC_HOME=$HOME/CCDC
    conda env config vars set CCDC_TOOLKIT_ASER_DATABASE=$CCDC_HOME/ccdc-data/csd/as544be
    conda env config vars set CCDC_TOOLKIT_SQLITE_DATABASE=$CCDC_HOME/ccdc-data/csd/as544be_ASER.sqlite
    conda env config vars set CCDC_ISOSTAR_DATA_DIRECTORY=$CCDC_HOME/ccdc-data/isostar
    conda env config vars set CCDC_MOGUL_DATA=$CCDC_HOME/ccdc-data/mogul
    conda env config vars set CCDC_MOGUL_INITIALISATION_FILE=$CONDA_PREFIX/lib/python3.9/site-packages/ccdc/parameter_files/mogul.ini
    conda env config vars set CCDC_CROSSMINER_DATABASE=$CCDC_HOME/ccdc-data/crossminer/csd_pdb_crossminer.feat
    conda env config vars set CCDC_CROSSMINER_FEATURE_DEFINITIONS=$CCDC_HOME/ccdc-data/crossminer/feature_definitions
    conda env config vars set GOLD_DIR=$CCDC_HOME/ccdc-software/gold/GOLD
    # always good:
    conda env config vars set PYTHONUSERBASE=$CONDA_PREFIX
"""

from ccdc import conformer, io, molecule
from typing import List, Tuple
from scipy.stats import norm
import numpy as np
import json
print( conformer._mogul_version() )

def get_mean_and_max_z_score(mol: molecule.Molecule) -> Tuple[float, float]:
    mol.assign_bond_types(which='unknown')
    mol.standardise_aromatic_bonds()
    mol.standardise_delocalised_bonds()
    mol.add_hydrogens()

    engine = conformer.GeometryAnalyser()
    geometry_analysed_mol: molecule.Molecule = engine.analyse_molecule(mol)

    analysed: List[conformer.GeometryAnalyser.Analysis]
    z_scores: List[float] = []
    #counts: List[int] = []
    for analysed in (geometry_analysed_mol.analysed_bonds,
                     geometry_analysed_mol.analysed_angles,
                     geometry_analysed_mol.analysed_torsions,
                     geometry_analysed_mol.analysed_rings):
        for a in analysed:
            if a.enough_hits:
                z_scores.append(a.z_score)
                #counts.append(a.nhits)
    z_scores = np.array(z_scores, dtype=np.float64)  # the None will be nan
    return np.nanmean(z_scores), np.nanmax(z_scores)

def parse(infile, outfile):
    assert '.json' in outfile, 'outfile must be a json file.
    mol_reader = io.MoleculeReader(infile)
    mol: molecule.Molecule
    data = []
    for mol in mol_reader:
        try:
            mu_z, max_z = get_mean_and_max_z_score(mol)
            data.append(dict(title=mol.identifier,
                             molecular_weight=round(mol.molecular_weight, 4),
                             mean_z_score=mu_z,
                                max_z_score=max_z)
                        )
        except Exception as e:
            print(e)
            continue
    mol_reader.close()
    with open(outfile, 'w') as fh:
        json.dump(data, fh, indent=2)

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Get mean and max z-scores of how CCDC typical bonds etc. from a molecule file')
    parser.add_argument('-i', '--infile', help='Input SDF file', required=True, type=str)
    parser.add_argument('-o', '--outfile', help='Output JSON file', required=True, type=str)
    args = parser.parse_args()
    parse(args.infile, args.outfile)


# attempt at weighted Fisher's method
# p_values = np.nan_to_num(2 * (1 - norm.cdf(np.abs(z_scores))), nan=1)
# test_statistics = [norm.ppf(1 - p_value/2) for p_value in p_values]
# weighted_test_statistics = [test_statistic * np.sqrt(frequency) for test_statistic, frequency in zip(test_statistics, counts)]
# sum_squared_weighted_test_statistics = np.sum(np.square(weighted_test_statistics))
# chi_squared = sum_squared_weighted_test_statistics
# p_value_chi2 = 1 - norm.cdf(chi_squared)
# combined_z_score = norm.ppf(1 - p_value_chi2/2)
# print(combined_z_score)