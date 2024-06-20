"""
Minimize a PDB file using PyRosetta

Usage:

.. code-block:: bash

    $ pyrosetta-minimize -t 1uqb.pdb -o 1uqb_min.pdb

"""

from pathlib import Path
import pyrosetta
import pyrosetta_help as ph
from types import ModuleType
from IPython.display import display, HTML

from collections import Counter
prc: ModuleType = pyrosetta.rosetta.core
prp: ModuleType = pyrosetta.rosetta.protocols
prcc: ModuleType = pyrosetta.rosetta.core.conformation
pr_scoring: ModuleType = pyrosetta.rosetta.core.scoring

def relax(original: pyrosetta.Pose,
          constraint_weight: float=5,
          cycles: int=15,
          relax_to_start_coords:bool=True) -> pyrosetta.Pose:
    pose: pyrosetta.Pose = original.clone()
    # configure constraints
    scorefxn: pr_scoring.ScoreFunction = pyrosetta.get_fa_scorefxn()
    scorefxn.set_weight(pr_scoring.ScoreType.coordinate_constraint, constraint_weight)
    scorefxn.set_weight(pr_scoring.ScoreType.angle_constraint, constraint_weight)
    scorefxn.set_weight(pr_scoring.ScoreType.atom_pair_constraint, constraint_weight)
    pyrosetta.rosetta.basic.options.set_boolean_option('relax:constrain_relax_to_start_coords', relax_to_start_coords)
    pyrosetta.rosetta.basic.options.set_boolean_option('relax:coord_constrain_sidechains', relax_to_start_coords)
    # set up the relax sampler
    pyrosetta.rosetta.protocols.relax.FastRelax.register_options()
    relax = pyrosetta.rosetta.protocols.relax.FastRelax(scorefxn, cycles)
    relax.constrain_relax_to_start_coords(relax_to_start_coords)
    relax.apply(pose)
    return pose

def init_pyrosetta(ignore_unrecognized_res=True, load_PDB_components=False):
    logger = ph.configure_logger()
    pyrosetta.distributed.maybe_init(extra_options=ph.make_option_string(no_optH=False,
                                                                         ex1=None,
                                                                         ex2=None,
                                                                         # mute='all',
                                                                         ignore_unrecognized_res=ignore_unrecognized_res,
                                                                         load_PDB_components=load_PDB_components,
                                                                         ignore_waters=True)
                                     )
    pyrosetta.rosetta.basic.options.set_boolean_option('run:ignore_zero_occupancy', False)
    pyrosetta.rosetta.basic.options.set_boolean_option('in:auto_setup_metals', True)

def main():
    import argparse

    parser = argparse.ArgumentParser(description='Minimise target without electron density')
    parser.add_argument('-t', '--target', help='Reference PDB file', required=True, type=str)
    parser.add_argument('-o', '--outfile', help='Outfile', default='minimised.pdb', type=str)
    parser.add_argument('-w', '--constraint_weight', help='constraint_weight', default=5, type=float)
    parser.add_argument('-c', '--cycles', help='cycles', default=15, type=int)
    parser.add_argument('-s', '--single_chain', help='keep only first chain', default=False, type=bool)
    parser.add_argument('--ignore_unrecognized_res',
                        help='If PyRosetta does not recognise a residue (chemical component) it will ignore it',
                        default=True, type=bool)
    parser.add_argument('--load_PDB_components',
                        help='Rosetta has a cache of residue (chemical component) as they appear in the PDB (mixed accuracy)',
                        default=False, type=bool)
    parser.add_argument('--relax_to_start_coords',
                        help='relax_to_start_coords same as setting `-w 0`',
                        default=True, type=bool)
    # ## Get values
    args: argparse.Namespace = parser.parse_args()
    target: str = args.target
    outfile: str = args.outfile
    constraint_weight: float = args.constraint_weight
    cycles:int = args.cycles
    single_chain:bool = args.single_chain
    # checks
    assert Path(target).exists(), f'{target} does not exist'
    assert not Path(outfile).exists(), f'{outfile} already exists'
    # init
    init_pyrosetta(ignore_unrecognized_res=args.ignore_unrecognized_res,
                   load_PDB_components=args.load_PDB_components)
    # load
    original: pyrosetta.Pose = pyrosetta.pose_from_file(target)
    con_tally = Counter(
        [type(con).__name__ for con in original.constraint_set().get_all_constraints()]).most_common()
    print('chains', prc.pose.conf2pdb_chain(original))
    print('sequence', original.sequence())
    print('Constraints present', con_tally)
    if single_chain:
        original: pyrosetta.Pose = original.split_by_chain(1)
    # relax
    relaxed = relax(original=original,
                    constraint_weight=constraint_weight,
                    cycles=cycles,
                    relax_to_start_coords=args.relax_to_start_coords)
    prp.toolbox.CA_superimpose(original, relaxed)
    print(prc.simple_metrics.metrics.RMSDMetric(original).calculate(relaxed))
    relaxed.dump_pdb(outfile)

if __name__ == '__main__':
    main()
