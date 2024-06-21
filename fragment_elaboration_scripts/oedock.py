# heavily adapted from the following source:
# https://github.com/volkamerlab/DocKin/blob/master/dockin/oe_docking.py

input_sdf = '👾👾.sdf'
output_sdf = '👾👾_docked.sdf'
protein_filename = '👾👾pdb'
# fix box ref mol too!


# if this code block fails,
# it is because the openeye license is not set, expired or cover OEDocking.
import os, re
import numpy as np
from pathlib import Path
oelicence = Path(os.environ['OE_LICENSE'])
assert oelicence.exists()
print('OE Licence', re.search('EXP_DATE.*', oelicence.read_text()).group() )

from openeye import oechem, oedocking, oequacpac, oeomega
import pandas as pd
from openbabel import openbabel as ob
from rdkit.Chem import AllChem
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem


def calculate_centroid(molecule):
    conf = molecule.GetConformer()
    positions = conf.GetPositions()
    centroid = positions.mean(axis=0)
    return centroid

def translate_to_coordinates(molecule, target_coords):
    centroid = calculate_centroid(molecule)
    translation_vector = target_coords - centroid
    conf = molecule.GetConformer()
    for i in range(molecule.GetNumAtoms()):
        pos = np.array(conf.GetAtomPosition(i))
        new_pos = pos + translation_vector
        conf.SetAtomPosition(i, new_pos)
    return molecule

def smiles_to_mol_bonds_pH(smiles, name: str, pH: float=7.4) -> Chem.Mol:
    obmol = ob.OBMol()
    conv = ob.OBConversion()
    conv.SetInFormat("smi")
    conv.ReadString(obmol, smiles)
    obmol.PerceiveBondOrders()
    obmol.CorrectForPH(pH)
    conv.SetOutFormat("mol")
    obmol.SetTitle( name )
    mol: Chem.Mol = Chem.MolFromMolBlock(conv.WriteString(obmol).strip())
    return mol

from typing import List

def enumerate_tautomers(mol: oechem.OEMol) -> List[oechem.OEMol]:
    tautomer_options = oequacpac.OETautomerOptions()
    tautomer_options.SetMaxTautomersGenerated(4096)
    tautomer_options.SetMaxTautomersToReturn(16)
    tautomer_options.SetCarbonHybridization(True)
    tautomer_options.SetMaxZoneSize(50)
    tautomer_options.SetApplyWarts(True)
    pKa_norm = True
    tautomers = [oechem.OEMol(tautomer) for tautomer in
                 oequacpac.OEGetReasonableTautomers(mol, tautomer_options, pKa_norm)]
    return tautomers
# -----------------------------------------
# ## Load mols

mols = []
with oechem.oemolistream(input_sdf) as ifs:
    for mol in ifs.GetOEMols():
        mols.append(oechem.OEGraphMol(mol))
print(len(mols), 'ligands')

# -----------------------------------------
# ## Load protein

with oechem.oemolistream() as ifs:
    ifs.SetFlavor(oechem.OEFormat_PDB,
            oechem.OEIFlavor_PDB_Default | oechem.OEIFlavor_PDB_DATA | oechem.OEIFlavor_PDB_ALTLOC)
    ifs.open(protein_filename)
    protein = oechem.OEGraphMol()
    assert oechem.OEReadMolecule(ifs, protein), 'failed reading'

c = np.array(list(protein.GetCoords().values()))
print(np.sqrt(np.square(c.max(axis=0) - c.min(axis=0)).sum()))

# -----------------------------------------
# ## Prepare protein

receptor = oechem.OEGraphMol()

# ### choice A: use the ligand
# oedocking.OEMakeReceptor(receptor, protein, molecule)

# ### choice B: use a box
ref_mol: Chem.Mol = ... # rdkit Chem.Mol
pos = ref_mol.GetConformer().GetPositions()
box_dims = np.min(pos, axis=0) - 5, np.max(pos, axis=0) + 5
ref_cen = calculate_centroid(ref_mol)
box = oedocking.OEBox(*box_dims[0], *box_dims[1])
oedocking.OEMakeReceptor(receptor, protein, box)

# -----------------------------------------
# ## Dock initialise

dock = oedocking.OEDock()
dock_resolution = oedocking.OESearchResolution_High
dock_method = oedocking.OEDockMethod_Hybrid2
#dock = oedocking.OEDock(dock_method, dock_resolution)
dock = oedocking.OEDock()
dock.Initialize(receptor)

# -----------------------------------------
# ## Dock

num_poses = 100

def score(molecule, dock=dock):
    """Return the docking score."""
    value = oechem.OEGetSDData(molecule, dock.GetName())
    return float(value)

omega_options = oeomega.OEOmegaOptions()
omega_options.SetMaxSearchTime(60.0)  # time out
omega = oeomega.OEOmega(omega_options)
omega.SetStrictStereo(False)  # enumerate stereochemistry if uncertain

results = {}
oemol: oechem.OEMol
oetautomol: oechem.OEMol
for oemol in mols:
    print('='*10, oemol.GetTitle(), '='*10)
    docked_tautomers = list()
    tautomers = enumerate_tautomers(oemol)
    for oetautomol in tautomers:
        docked_mol = oechem.OEMol()
        # expand conformers
        omega.Build(oetautomol)

        # dock molecule
        return_code = dock.DockMultiConformerMolecule(docked_mol, oetautomol, num_poses)
        if return_code != oedocking.OEDockingReturnCode_Success:
            print(f'Docking failed for molecule with title {oetautomol.GetTitle()} with error code '
                  f'{oedocking.OEDockingReturnCodeGetName(return_code)}.')
            continue

        # store docking data
        oedocking.OESetSDScore(docked_mol, dock, dock.GetName())

        # expand conformations
        for conformation in docked_mol.GetConfs():
            docked_tautomers.append(oechem.OEGraphMol(conformation))

    # sort all conformations of all tautomers by score
    docked_tautomers.sort(key=score)

    # keep number of conformations as specified by num_poses
    docked = docked_tautomers[:num_poses]
    results[oemol.GetTitle()] = docked

# -----------------------------------------
# ## Write results

with oechem.oemolostream(output_sdf) as ofs:
    for name, docked in results.items():
        for mol in docked:
            mol.SetData('Id', name)
            oechem.OEWriteMolecule(ofs, mol)