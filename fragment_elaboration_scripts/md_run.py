"""
Given the enviroment variable $VC of a mol file, and a $TEMPLATE receptor/template/protein pdb file,
run $NSTEPS 40 fs steps of OpenMM simulation in explicit solvent
and save a pickled DataFrame of the RMSD of the ligand to the template.

At 40 fs steps, 500_000 steps is 20 ns.
1_250_000 steps is 50 ns.
2_500_000 steps is 100 ns.


Bash script to run this:

.. code-block:: bash

    #!/bin/bash
    source $CONDA_PREFIX/etc/profile.d/conda.sh
    conda install -n openmm -c conda-forge -c omnia python=3.10 rdkit openff-toolkit openforcefield pdbfixer torch
    conda activate openmm
    # don't forget to make sure CUDA words:
    # I install torch and use that cudatoolkit as I am too lazy to disinstall the wrong cuda-toolkit
    # that openmm installs.
    conda env config set var LD_LIBRARY_PATH = $CONDA_PREFIX/lib/python3.10/site-packages/triton/third_party/cuda
    conda env config set var CUDA_HOME = $CONDA_PREFIX/lib/python3.10/site-packages/triton/third_party/cuda

    cd 👾👾
    conda activate openmm
    export SLACK_WEBHOOK='https://hooks.slack.com/services/👾👾/👾👾'
    export TEMPLATE=👾👾/👾👾.pdb
    export VC=👾👾/👾👾.mol;
    export NSTEPS=500_000
    python 👾👾/👾👾/openmm_validation.py;
    curl -X POST -H 'Content-type: application/json' --data '{"text":"OpenMM '$VC' complete"}' $SLACK_WEBHOOK;
"""

import os, sys, io
from pathlib import Path
from typing import Sequence, Union, Optional

from rdkit import Chem
from rdkit.Chem import AllChem

import openmm as mm
import openmm.app as mma
import openmm.unit as mmu
import openff.toolkit.topology as fft
from openmmforcefields.generators import SMIRNOFFTemplateGenerator

from fragmenstein import Fritz
from fragmenstein.demo import MPro

from pathlib import Path
import os

implicit_solvent = False
n_equlibration_steps = 1_000
n_save_steps = 1_000
n_steps = os.environ.get('NSTEPS', 500_000)

# apo_block = MPro.get_template()
# mol = MPro.get_mol('x2646')

template_filename = os.environ['TEMPLATE']
apo_block = Path(template_filename).read_text()
mol_name = os.environ['VC']
mol = Chem.MolFromMolBlock(Path(f'picks/{mol_name}.mol').read_text())
assert mol is not None, 'Failed to read Mol'


# print(mm.Platform.getPluginLoadFailures())
# assert mm.Platform.getPlatformByName('CUDA').supportsKernels('CUDA')


def pdbblock_to_PDB(pdb_block: str) -> mma.PDBFile:
    """
    Read a PDB block (string) into a PDBFile.

    Copied from Fritz.pdbblock_to_PDB https://github.com/matteoferla/Fragmenstein/blob/e4001e7781e50a552d27a79176db0846693ab494/fragmenstein/openmm/fritz.py#L170C5-L179C19
    """
    assert isinstance(pdb_block, str), 'pdb_block must be a string'
    assert len(pdb_block) > 0, 'pdb_block must be a non-empty string'
    iostr = io.StringIO(pdb_block)
    iostr.seek(0)
    pdb = mma.PDBFile(iostr)
    return pdb


def create_forcefield(mol: Chem.Mol,
                      forcefield_names: Sequence[str]) -> mma.ForceField:
    """set up forcefield with ligand

    Copied from https://github.com/matteoferla/Fragmenstein/blob/e4001e7781e50a552d27a79176db0846693ab494/fragmenstein/openmm/fritz.py#L181
    """
    ideal_mol: Chem.Mol = AllChem.AddHs(mol, addCoords=True)
    ideal_mol.RemoveAllConformers()
    AllChem.EmbedMolecule(ideal_mol, enforceChirality=True)
    molecule = fft.Molecule.from_rdkit(ideal_mol, allow_undefined_stereo=True)
    smirnoff = SMIRNOFFTemplateGenerator(molecules=molecule)
    forcefield = mma.ForceField(*forcefield_names)
    forcefield.registerTemplateGenerator(smirnoff.generator)
    return forcefield


def plonk(apo: Union[str, mma.PDBFile, mma.Modeller], mol: Optional[Chem.Mol] = None) -> mma.Modeller:
    """
    Plonk the ligand into the apo structure.
    """
    if mol is None:
        mol = self.prepped_mol
    if isinstance(apo, str) and '.pdb' in apo and apo.find('\n') == -1:
        raise ValueError('you passed a filename, not a block')
    elif isinstance(apo, str):  # PDB block
        apo: mma.PDBFile = pdbblock_to_PDB(apo)
    # make a copy:
    holo: mma.Modeller = mma.Modeller(apo.topology, apo.positions)
    lig: mma.Modeller = rdkit_to_openMM(mol)
    holo.add(lig.topology, lig.positions)  # noqa mmu.Quantity is okay
    return holo


def rdkit_to_openMM(mol: Chem.Mol) -> mma.Modeller:
    # rdkit AssignStereochemistryFrom3D previously applied
    lig = fft.Molecule.from_rdkit(mol, allow_undefined_stereo=True)
    # minor corrections:
    # there is no need to fix via
    # lig_topo._chains[0]._residues[0].name = 'LIG'
    for r_a, m_a in zip(mol.GetAtoms(), lig.atoms):
        assert r_a.GetSymbol() == m_a.symbol, 'Discrepancy'
        r_name = r_a.GetPDBResidueInfo().GetName()
        m_a.name = r_name
    # convert and merge
    lig_topo: mma.Topology = fft.Topology.from_molecules([lig]).to_openmm()
    lig_pos: mmu.Quantity = lig.conformers[0].to_openmm()
    return mma.Modeller(lig_topo, lig_pos)  # noqa mmu.Quantity is okay


# --------------------------


mol = AllChem.AddHs(mol, addCoords=True)
AllChem.SanitizeMol(mol)
Fritz.correct_pdbinfo(mol=mol, resn='LIG', resi=1, chain='X')
AllChem.SanitizeMol(mol)

pdb: mma.PDBFile = pdbblock_to_PDB(apo_block)
apo = mma.Modeller(pdb.topology, pdb.positions)
holo = plonk(apo, mol)

forcefield_names = ('amber14-all.xml', 'implicit/gbn2.xml' if implicit_solvent else 'amber14/tip3pfb.xml')
forcefield: mma.ForceField = create_forcefield(mol, forcefield_names)

holo.addHydrogens(forcefield, pH=7.0)
if not implicit_solvent:
    holo.addSolvent(forcefield, model="tip3p", padding=12.0 * mmu.angstroms)

system: mm.System = forcefield.createSystem(holo.topology,
                                            nonbondedMethod=mma.NoCutoff,
                                            nonbondedCutoff=1 * mmu.nanometer,
                                            constraints=mma.HBonds)
integrator = mm.LangevinMiddleIntegrator(300 * mmu.kelvin, 1 / mmu.picosecond, 0.004 * mmu.picoseconds)
simulation = mma.Simulation(holo.topology, system, integrator)
simulation.context.setPositions(holo.positions)
assert simulation.context.getPlatform().getName() == 'CUDA'
simulation.minimizeEnergy()
positions: mmu.Quantity = simulation.context.getState(getPositions=True).getPositions()

state: mm.State = simulation.context.getState(getEnergy=True)
print(
    f'The system has {state.getPotentialEnergy().value_in_unit(mmu.kilocalorie_per_mole):.1f} kcal/mol in potential energy')
print(
    f'The system has {state.getKineticEnergy().value_in_unit(mmu.kilocalorie_per_mole):.1f} kcal/mol in kinitic energy')

prepped_filename = 'openmm_output/' + template_filename.replace('.pdb', f'-{mol_name}.pdb')
mma.PDBFile.writeFile(simulation.topology, positions, prepped_filename)

# Minimize the energy
simulation.minimizeEnergy()
min_filename = 'openmm_output/' + template_filename.replace('.pdb', f'-{mol_name}.min.pdb')
mma.PDBFile.writeFile(simulation.topology, positions, min_filename)

# Equilibrate
simulation.step(n_equlibration_steps)

# Run a short trajectory
traj_filename = 'openmm_output/' + template_filename.replace('.pdb', f'-{mol_name}.dcd')
simulation.reporters.append(mma.DCDReporter(traj_filename, n_save_steps))
simulation.reporters.append(
    mma.StateDataReporter(sys.stdout, n_save_steps, step=True, potentialEnergy=True, kineticEnergy=True,
                          temperature=True))
simulation.step(n_steps)
end_filename = 'openmm_output/' + template_filename.replace('.pdb', f'-{mol_name}.end.pdb')
mma.PDBFile.writeFile(simulation.topology, positions, end_filename)

# analyse RMSD

import MDAnalysis as mda
from MDAnalysis.analysis import align
from MDAnalysis.analysis import rms
import plotly.express as px
import plotly.io as pio

pio.renderers.default = 'iframe'
import pandas as pd

# Load the topology and trajectory
u = mda.Universe(prepped_filename, traj_filename)

# superpose frames
aligner = align.AlignTraj(u, u, select="backbone", in_memory=True)
aligner.run()

# Set up the RMSD calculation
selection = "resname LIG"
R = rms.RMSD(u.select_atoms(selection), u.select_atoms(selection), select=selection, ref_frame=0)
R.run()

# Create a DataFrame for the RMSD values
rmsd_df = pd.DataFrame({
    'Frame': R.results.rmsd[:, 1],
    'RMSD': R.results.rmsd[:, 2]
})

rmsd_df.to_pickle(f'openmm_output/{mol_name}.pkl.gz')
