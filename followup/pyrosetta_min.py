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

def main(target: str,
             outfile: str,
             constraint_weight: float=5,
             cycles: int=15,
             single_chain: bool = True):
    # checks
    assert Path(target).exists(), f'{target} does not exist'
    assert not Path(outfile).exists(), f'{outfile} already exists'
    # init
    init_pyrosetta()
    # load
    original: pyrosetta.Pose = pyrosetta.pose_from_file(target)
    con_tally: dict = Counter([type(con).__name__ for con in original.constraint_set().get_all_constraints()]).most_common()
    print('chains', prc.pose.conf2pdb_chain(original))
    print('sequence', original.sequence())
    print('Constraints present', con_tally)
    if single_chain:
        original: pyrosetta.Pose = original.split_by_chain(1)
    # relax
    relaxed = relax(original, constraint_weight,cycles)
    prp.toolbox.CA_superimpose(original, relaxed)
    print(prc.simple_metrics.metrics.RMSDMetric(original).calculate(relaxed))
    relaxed.dump_pdb(outfile)

def relax(original: pyrosetta.Pose, constraint_weight: float=5, cycles: int=15) -> pyrosetta.Pose:
    pose: pyrosetta.Pose = original.clone()
    scorefxn: pr_scoring.ScoreFunction = pyrosetta.get_fa_scorefxn()
    scorefxn.set_weight(pr_scoring.ScoreType.coordinate_constraint, constraint_weight)
    scorefxn.set_weight(pr_scoring.ScoreType.angle_constraint, constraint_weight)
    scorefxn.set_weight(pr_scoring.ScoreType.atom_pair_constraint, constraint_weight)
    pyrosetta.rosetta.basic.options.set_boolean_option('relax:constrain_relax_to_start_coords', True)
    pyrosetta.rosetta.basic.options.set_boolean_option('relax:coord_constrain_sidechains', True)
    pyrosetta.rosetta.protocols.relax.FastRelax.register_options()
    relax = pyrosetta.rosetta.protocols.relax.FastRelax(scorefxn, cycles)
    relax.constrain_relax_to_start_coords(True)
    relax.apply(pose)
    return pose

def init_pyrosetta():
    logger = ph.configure_logger()
    pyrosetta.distributed.maybe_init(extra_options=ph.make_option_string(no_optH=False,
                                                                         ex1=None,
                                                                         ex2=None,
                                                                         # mute='all',
                                                                         ignore_unrecognized_res=True,
                                                                         load_PDB_components=False,
                                                                         ignore_waters=True)
                                     )
    pyrosetta.rosetta.basic.options.set_boolean_option('run:ignore_zero_occupancy', False)
    pyrosetta.rosetta.basic.options.set_boolean_option('in:auto_setup_metals', True)


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Minimise target without electron density')
    parser.add_argument('-t', '--target', help='Reference PDB file', required=True)
    parser.add_argument('-o', '--outfile', help='Outfile', default='blast_hits.csv')
    parser.add_argument('-w', '--constraint_weight', help='constraint_weight', default=5, type=float)
    parser.add_argument('-c', '--cycles', help='cycles', default=15, type=int)
    args = parser.parse_args()

    main(target=args.target,
         outfile=args.outfile,
         constraint_weight=args.constraint_weight,
         cycles=args.cycles
         )
