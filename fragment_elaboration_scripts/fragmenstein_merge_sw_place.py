"""

Example weights file:

    {"N_rotatable_bonds": 3, "\\u2206\\u2206G": 3, "interaction_uniqueness_metric": -20, "N_unconstrained_atoms": 0.5, "N_constrained_atoms": -0.2, "N_interactions": -5, "N_interactions_lost": 10, "max_hit_Tanimoto": -1}

Example of multiple dbs

    --sw_databases Enamine-SC-Stock-Mar2022.smi.anon Enamine-BB-Stock-Mar2022.smi.anon REAL-Database-22Q1.smi.anon
"""

from rdkit.ML.Cluster import Butina
from rdkit.Chem import rdMolDescriptors as rdmd
from rdkit.Chem import Descriptors

from typing import List, Dict, Any, Optional
import operator, os, re, logging, random, time, argparse, string, itertools, json, contextlib, requests
from warnings import warn
from pathlib import Path

import pyrosetta
import pyrosetta_help as ph

import pandas as pd
import pandera.typing as pdt
from pandarallel import pandarallel
from smallworld_api import SmallWorld


from rdkit import Chem, rdBase, DataStructs
from rdkit.Chem import AllChem, Draw, PandasTools, BRICS
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import rdFingerprintGenerator as rdfpg
from rdkit.Chem.rdfiltercatalog import FilterCatalogParams, FilterCatalog, FilterCatalogEntry


Chem.SetDefaultPickleProperties(Chem.PropertyPickleOptions.AllProps)

from fragmenstein import Igor, Victor, Laboratory, Monster
from fragmenstein.laboratory.validator import place_input_validator

pandarallel.initialize()

SUPRESSED_EXCEPTION = Exception  # debug purposes...
DEFAULT_WEIGHTS = {"N_rotatable_bonds": 3, "\u2206\u2206G": 3, "interaction_uniqueness_metric": -20, "N_unconstrained_atoms": 0.5, "N_constrained_atoms": -0.2, "N_interactions": -5, "N_interactions_lost": 10, "max_hit_Tanimoto": -1, "N_PAINS": 20, "strain_per_HA": 1}
# ------------------------------------------------------

logger = ph.configure_logger()
if logger.handlers:
    logger.handlers[0].setLevel(logging.ERROR)  # logging.WARNING = 30

Igor.init_pyrosetta()

try:
    sws = SmallWorld()
    chemical_databases: pd.DataFrame = sws.retrieve_databases()
except requests.HTTPError as error:
    warn(f'sw.docking.org is down. {error}')

# ------------------------------------------------------

def set_up(output, cutoff, quick, suffix, **kwargs):
    os.makedirs(output, exist_ok=True)
    Victor.work_path = output
    Victor.monster_throw_on_discard = True  # stop this merger if a fragment cannot be used.
    Victor.monster_joining_cutoff = cutoff  # Å
    Victor.quick_reanimation = quick  # for the impatient
    Victor.error_to_catch = Exception  # stop the whole laboratory otherwise
    Victor.enable_stdout(logging.CRITICAL)
    Victor.enable_logfile(os.path.join(output, f'{suffix}.log'), logging.ERROR)


def config_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', '--template', help='Template PDB file', required=True)
    parser.add_argument('-i', '--hits', help='Hits SDF file', required=True)
    parser.add_argument('-o', '--output', help='Output folder', default='output')
    parser.add_argument('-r', '--ranking', help='Ranking method', default='∆∆G')
    parser.add_argument('-c', '--cutoff', help='Joining cutoff', default=5)
    parser.add_argument('-q', '--quick', help='Quick reanimation', default=False)
    parser.add_argument('-d', '--sw_dist', help='SmallWorld distance', default=25, type=int)
    parser.add_argument('-l', '--sw_length', help='SmallWorld length', default=50, type=int)
    parser.add_argument('-b', '--sw_databases', help='SmallWorld databases. Accepts multiple',
                        nargs='+', default=sws.REAL_dataset)
    parser.add_argument('-s', '--suffix', help='Suffix for output files', default='')
    parser.add_argument('-n', '--n_cores', help='Number of cores', default=55, type=int)
    parser.add_argument('-m', '--combination_size', help='Number of hits to combine in one step', default=2, type=int)
    parser.add_argument('-k', '--top_mergers', help='Max number of mergers to followup up on', default=500, type=int)
    parser.add_argument('-e', '--timeout', help='Timeout for each merger', default=240, type=int)
    parser.add_argument('-x', '--max_tasks', help='Max number of combinations to try', default=0, type=int)
    parser.add_argument('-z', '--blacklist', help='Blacklist file', default='')
    parser.add_argument('-j', '--weights', help='JSON weights file', default='')
    return parser


# ------------------------------------------------------

def replace_hits(pdbblock, hits, n_cores, timeout, suffix, **settings):
    # place themselves for a ∆∆G score
    lab = Laboratory(pdbblock=pdbblock, covalent_resi=None, run_plip=True)
    selfies = pd.DataFrame([dict(name=hit.GetProp('_Name'),
                                 hits=[hit],
                                 smiles=Chem.MolToSmiles(hit)
                                 ) for hit in hits])
    replacements: pd.DataFrame = lab.place(place_input_validator(selfies), n_cores=n_cores, timeout=timeout)
    fix_intxns(replacements)
    replacements['bleached_name'] = replacements['name']
    replacements['name'] = replacements.hit_mols.apply(lambda ms: ms[0].GetProp('_Name'))
    replacements.to_pickle(f'fragmenstein_hit_replacements{suffix}.pkl.gz')
    # replacements.to_csv(f'fragmenstein_hit_replacements{suffix}.csv')
    return replacements

def UFF_Gibbs(mol):
    # free energy cost of bound conformer
    if not isinstance(mol, Chem.Mol) or mol.GetNumHeavyAtoms() == 0:
        return float('nan')
    with contextlib.suppress(SUPRESSED_EXCEPTION):
        AllChem.SanitizeMol(mol)
        # this is actually UFF
        copy = Chem.Mol(mol)
        return Monster.MMFF_score(None, mol, True)
    return float('nan')


params = FilterCatalogParams()
params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
catalog = FilterCatalog(params)

def get_pains(mol) -> List[str]:
    with contextlib.suppress(Exception):
        entry: FilterCatalogEntry
        if not isinstance(mol, Chem.Mol) or mol.GetNumHeavyAtoms() == 0:
            return []
        AllChem.SanitizeMol(mol)
        return [entry.GetDescription() for entry in catalog.GetMatches(mol)]

def merge(hits, pdbblock, suffix, n_cores, combination_size, timeout, max_tasks, blacklist, **settings) -> pd.DataFrame:
    lab = Laboratory(pdbblock=pdbblock, covalent_resi=None)
    lab.blacklist = blacklist
    tick = time.time()
    combinations: pd.DataFrame = lab.combine(hits,
                                             n_cores=n_cores,
                                             timeout=timeout,
                                             combination_size=combination_size,
                                             max_tasks=max_tasks)
    with rdBase.BlockLogs():
        combinations['simple_smiles'] = combinations.unminimized_mol.apply(Victor.to_simple_smiles)
    combinations.to_pickle(f'fragmenstein_mergers{suffix}.pkl.gz')
    combinations.to_csv(f'fragmenstein_mergers{suffix}.csv')
    print(tick - time.time())
    return combinations


# ------------------------------------------------------
def get_custom_map(row: pd.Series) -> Dict[str, Dict[int, int]]:
    """
    SmallWorld returns a mapping of the indices.
    This functions maps the indices back to the original hit molecules.

    :param row: SW pd.Series (or dict)
    :return:
    """
    temp = Victor(row.hits, pdb_block=pdbblock)
    temp.monster.positioned_mol = row.unminimized_merger
    temp.minimized_mol = row.minimized_merger
    return temp.migrate_sw_origins(row)


def search(combinations: pd.DataFrame, suffix: str,
           sw_dist:int, sw_length:int, top_mergers: int,
           ranking: str, sw_db:str,ranking_ascending: Optional[bool]=None, **setting) -> pd.DataFrame:
    if ranking_ascending is None:
        ranking_ascending = False if ranking in ('LE', 'N_interactions') else True
    queries = combinations.sort_values(ranking, ascending=ranking_ascending) \
        .loc[(combinations.outcome == 'acceptable')] \
        .drop_duplicates('simple_smiles') \
        .reset_index() \
        .head(top_mergers)
    similars = sws.search_many(queries.simple_smiles.to_list(),
                               dist=sw_dist,
                               length=sw_length,
                               db=sw_db,
                               tolerated_exceptions=Exception)
    print(f'Found {len(similars)} analogues')
    # query_index was added clientside to keep track!
    similars['catalogue'] = sw_db
    similars['query_name'] = similars.query_index.map(queries.name.to_dict())
    similars['hits'] = similars.query_index.map(queries.hit_mols.to_dict())
    similars['hit_names'] = similars.hits.apply(lambda m: [mm.GetProp('_Name') for mm in m])
    similars['minimized_merger'] = similars.query_index.map(queries.minimized_mol.to_dict())
    similars['unminimized_merger'] = similars.query_index.map(queries.unminimized_mol.to_dict())
    similars['name'] = similars['id'] + ':' + similars['query_name']
    similars['smiles'] = similars.hitSmiles.str.split(expand=True)[0]
    similars['custom_map'] = similars.apply(get_custom_map, axis=1)
    similars.to_pickle(f'fragmenstein_similars{suffix}.{sw_db}.pkl.gz')
    return similars


# ------------------------------------------------------

def place(similars, pdbblock, n_cores, timeout, suffix, **settings) -> pd.DataFrame:
    lab = Laboratory(pdbblock=pdbblock, covalent_resi=None, run_plip=True)
    placements: pd.DataFrame = lab.place(place_input_validator(similars), n_cores=n_cores, timeout=timeout)
    placements.loc[(placements['∆∆G'] > -1) & (placements.outcome == 'acceptable'), 'outcome'] = 'weak'
    placements['unminimized_mol'] = placements.unminimized_mol.fillna(Chem.Mol())
    fix_intxns(placements)
    placements.to_pickle(f'fragmenstein_placed{suffix}.pkl.gz')
    placements.to_csv(f'fragmenstein_placed{suffix}.csv')
    print(placements.outcome.value_counts())
    return placements


def fix_intxns(df):
    intxn_names = [c for c in df.columns if isinstance(c, tuple)]
    for intxn_name in intxn_names:
        df[intxn_name] = df[intxn_name].fillna(0).astype(int)


def write(placements, suffix):
    valids = placements.loc[placements.outcome == 'acceptable'].sort_values('∆∆G')
    valids['combined_name'] = valids['name']
    valids['name'] = valids['enamine_name']
    valids.to_pickle(f'fragmenstein_acceptables{suffix}.pkl.gz')
    valids.path.apply(Path.exists).value_counts()


def correct_weaklings(hit_replacements, df):
    df['hit_names'] = df.hit_mols \
        .apply(lambda v: v if isinstance(v, list) else []) \
        .apply(lambda v: [m.GetProp('_Name') for m in v])
    dG_mapping = hit_replacements.set_index('name')['∆∆G'].to_dict()
    dG_mapping.update({k.replace('-', '_'): v for k, v in dG_mapping.items()})
    get_lowest = lambda names: min([0] + [dG_mapping.get(name, 0) for name in names])
    df['lowest_hit_∆∆G'] = df.hit_names.apply(get_lowest)
    worseness_mask = (df['∆∆G'] > df['lowest_hit_∆∆G'] * 0.8) & (
            df.outcome == 'acceptable')
    df.loc[worseness_mask, 'outcome'] = 'weaker'


def core_ops(hit_replacements, **settings):
    combinations: pd.DataFrame = merge(**settings)
    correct_weaklings(hit_replacements, combinations)
    uncat_similars: List[pd.DataFrame] = []
    for sw_db in settings['sw_databases']:
        s = search(combinations, sw_db=sw_db, **settings)
        if len(s):
            uncat_similars.append(s)
    assert uncat_similars, 'No analogues were found!'
    similars: pd.DataFrame = pd.concat(uncat_similars, ignore_index=True).drop_duplicates('smiles')
    similars.to_pickle(f'fragmenstein_similars{settings["suffix"]}.pkl.gz')
    placements: pd.DataFrame = place(similars, **settings)
    return placements


# ----- Scoring -----------------------------------
class GetRowSimilarity:
    def __init__(self, hits):
        self.fpgen = rdfpg.GetRDKitFPGenerator()
        self.hit2fp = {h.GetProp('_Name'): self.fpgen.GetFingerprint(h) for h in hits}

    def __call__(self, row: pd.Series):
        with contextlib.suppress(SUPRESSED_EXCEPTION):
            if not isinstance(row.minimized_mol, Chem.Mol):
                return float('nan')
            elif isinstance(row.hit_names, str):
                hit_names = row.hit_names.split(',')
            elif isinstance(row.hit_names, list):
                hit_names = row.hit_names
            else:
                return float('nan')
            fp = self.fpgen.GetFingerprint(AllChem.RemoveHs(row.minimized_mol))
            return max([DataStructs.TanimotoSimilarity(fp, self.hit2fp[name]) for name in hit_names])
        return float('nan')


class HitIntxnTallier:

    def __init__(self, hit_replacements):
        self.slim_hits = self.slim_down(hit_replacements)

    def slim_down(self, hit_replacements):
        # the bleaching was fixed cf bleached_name
        # hit_replacements['new_name'] = hit_replacements.name.str.replace('-', '_')
        # undoing bleaching
        hit_replacements['new_name'] = hit_replacements.hit_mols.apply(lambda ms: ms[0].GetProp('_Name'))
        columns = [c for c in hit_replacements.columns if isinstance(c, tuple)]
        return hit_replacements.set_index('new_name')[columns].fillna(0).copy()

    def __call__(self, row: pd.Series):
        with contextlib.suppress(SUPRESSED_EXCEPTION):
            if not isinstance(row.minimized_mol, Chem.Mol) or isinstance(row.hit_names, float):
                return float('nan'), float('nan')
            present_tally = 0
            absent_tally = 0
            for hit_name in list(row.hit_names):
                if hit_name not in self.slim_hits.index:
                    raise Exception('Name' + hit_name)
                hit_row = self.slim_hits.loc[hit_name]
                for intxn_name, hit_value in hit_row.items():
                    if not hit_value:
                        continue
                    elif intxn_name not in row.index:
                        absent_tally += 1 if intxn_name[0] != 'hydroph_interaction' else 0.5
                    elif row[intxn_name]:
                        absent_tally += 1 if intxn_name[0] != 'hydroph_interaction' else 0.5
                    else:
                        present_tally += 1 if intxn_name[0] != 'hydroph_interaction' else 0.5
            return present_tally, absent_tally
        return float('nan'), float('nan')


class UniquenessMeter:
    def __init__(self, tallies, intxn_names, k=0.5):
        self.tallies = tallies
        self.intxn_names = intxn_names
        self.k = k

    def __call__(self, row):
        with contextlib.suppress(SUPRESSED_EXCEPTION):
            return sum([(row[name] / self.tallies[name]) ** self.k for name in self.intxn_names if
                        row[name] and self.tallies[name]])
        return float('nan')

    def tally_interactions(self, row):
        return sum([row[c] if self.intxn_names[0] != 'hydroph_interaction' else row[c] * 0.5 for c in self.intxn_names])


class PenaltyMeter:
    def __init__(self, weights):
        self.weights = weights

    def __call__(self, row):
        with contextlib.suppress(SUPRESSED_EXCEPTION):
            penalty = 0
            if row.outcome != 'acceptable':
                return float('inf')
            for col, w in self.weights.items():
                if col not in row.index:
                    warn(f'{col} column is missing from df')
                    continue
                penalty += row[col] * w
            return penalty
        return float('nan')



def butina_cluster(mol_list, cutoff=0.35):
    # https://github.com/PatWalters/workshop/blob/master/clustering/taylor_butina.ipynb
    fp_list = [rdmd.GetMorganFingerprintAsBitVect(AllChem.RemoveAllHs(m), 3, nBits=2048) for m in mol_list]
    dists = []
    nfps = len(fp_list)
    for i in range(1, nfps):
        sims = DataStructs.BulkTanimotoSimilarity(fp_list[i], fp_list[:i])
        dists.extend([1 - x for x in sims])
    mol_clusters = Butina.ClusterData(dists, nfps, cutoff, isDistData=True)
    cluster_id_list = [0] * nfps
    for idx, cluster in enumerate(mol_clusters, 1):
        for member in cluster:
            cluster_id_list[member] = idx
    return cluster_id_list


def score(placements, hit_replacements, suffix, hits, weights, **settings):
    # tanimoto
    get_similarity = GetRowSimilarity(hits)
    placements['max_hit_Tanimoto'] = placements.apply(get_similarity, axis=1)
    # properties
    m = placements.minimized_mol.apply(lambda m: m if isinstance(m, Chem.Mol) else Chem.Mol())
    # macrocyclics... yuck.
    placements['largest_ring'] = m.apply(lambda mol: max([0] + list(map(len, mol.GetRingInfo().AtomRings()))))
    # interactions
    intxn_names = [c for c in placements.columns if isinstance(c, tuple)]
    for intxn_name in intxn_names:
        placements[intxn_name] = placements[intxn_name].fillna(0).astype(int)
    tally_hit_intxns = HitIntxnTallier(hit_replacements)
    hit_checks = placements.apply(tally_hit_intxns, axis=1)
    placements['N_interactions_kept'] = hit_checks.apply(operator.itemgetter(0))  # .fillna(0).astype(int)
    placements['N_interactions_lost'] = hit_checks.apply(operator.itemgetter(1))  # .fillna(99).astype(int)
    tallies = placements[intxn_names].sum()
    ratioed = UniquenessMeter(tallies, intxn_names, k=0.5)
    placements['interaction_uniqueness_metric'] = placements.apply(ratioed, axis=1)
    placements['N_HA'] = m.apply(Chem.Mol.GetNumHeavyAtoms)
    placements['N_rotatable_bonds'] = m.apply(Chem.rdMolDescriptors.CalcNumRotatableBonds)
    placements['N_interactions'] = placements.apply(ratioed.tally_interactions, axis=1)
    placements['PAINSes'] = placements.minimized_mol.apply(get_pains)
    placements['N_PAINS'] = placements.PAINSes.apply(len)
    placements['UFF_Gibbs'] = placements.minimized_mol.apply(UFF_Gibbs)
    placements['strain_per_HA'] = placements.UFF_Gibbs / min(placements.N_HA, 0.0001)
    penalize = PenaltyMeter(weights)
    placements['ad_hoc_penalty'] = placements.apply(penalize, axis=1)
    # placements['N_interactions'] = placements[[c for c in placements.columns if isinstance(c, tuple)]].fillna(0).sum(axis=1)
    with contextlib.suppress(SUPRESSED_EXCEPTION):
        placements['cluster'] = butina_cluster(m.to_list())


if __name__ == '__main__':
    parser = config_parser()
    # load
    settings: Dict[str, Any] = vars(parser.parse_args())
    set_up(**settings)
    with Chem.SDMolSupplier(settings['hits'].strip()) as sd:
        # hitdex: Dict[str, Chem.Mol] = {mol.GetProp('_Name'): mol for mol in sd}
        hits: List[Chem.Mol] = list(sd)
    settings['hits'] = hits
    if settings['blacklist']:
        with open(settings['blacklist'].strip()) as fh:
            settings['blacklist'] = [line.strip() for line in fh.readlines()]
    else:
        settings['blacklist'] = []
    print(f'N hits: {len(hits)}')
    with open(settings['template'].strip()) as fh:
        pdbblock = fh.read()
    settings['pdbblock'] = pdbblock
    if settings['weights']:
        with open(settings['weights'].strip()) as fh:
            settings['weights'] = json.load(fh)
    else:
        settings['weights'] = DEFAULT_WEIGHTS
    # self
    hit_replacements: pd.DataFrame = replace_hits(**settings)
    # run
    max_tasks = settings['max_tasks']
    hitnames = [h.GetProp('_Name') for h in hits]
    all_names = list(map('-'.join, itertools.permutations(hitnames, settings['combination_size'])))
    if max_tasks == 0 or max_tasks > len(all_names):
        core_ops(hit_replacements, **settings)
        exit()
    base_suffix = settings['suffix']
    all_placements = pd.DataFrame()
    letters = iter(string.ascii_uppercase)
    for i in range(0, len(all_names) + max_tasks, max_tasks):
        settings['suffix'] = base_suffix + next(letters)
        placements: pd.DataFrame = core_ops(hit_replacements, **settings)
        settings['blacklist'] += all_names[i:i + max_tasks]
        all_placements = pd.concat([all_placements, placements], ignore_index=True)
    correct_weaklings(hit_replacements, all_placements)
    settings['suffix'] = base_suffix
    all_placements.to_pickle(f'fragmenstein_placed{base_suffix}.pkl.gz')
    score(all_placements, hit_replacements, **settings)
    all_placements.to_pickle(f'fragmenstein_placed{base_suffix}.pkl.gz')
    all_placements.to_csv(f'fragmenstein_placed{base_suffix}.csv')
