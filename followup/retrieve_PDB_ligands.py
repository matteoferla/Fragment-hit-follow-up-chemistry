"""
This script formerly called ``hit_theft``, finds sequence homologues in the PDB and steals their ligands.

It does the following step:

* Get the sequence of the longest chain in the PDB file
* Blast it against the PDB dataset in NCBI's servers


# --- YAML file for installation -------
# conda env update -n thief --file thief.yml

name: thief
channels:
  - conda-forge
  - schrodinger
  - defaults
dependencies:
  - python==3.12
  # pymol-bundle
  - pymol-open-source
  - pip
  - pip:
    - requests
    # pyrosetta-help
    - rdkit
    - biopython

"""
import pandas as pd
import pymol2
import re, logging
from typing import Dict, List, Tuple
import chempy
from typing import List, Set, Dict
from collections import defaultdict
from IPython.display import Image as IPyImage
from PIL.PngImagePlugin import PngImageFile
import json
import pickle
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
import contextlib


def get_log(level: int = logging.INFO) -> logging.Logger:
    log = logging.getLogger()
    if len(log.handlers) > 0:
        log.setLevel(level)
        return log
    handler = logging.StreamHandler()
    log.addHandler(handler)
    log.setLevel(level)
    handler.setLevel(level)
    return log


def load_reference(pymol, filename_or_block: str) -> None:
    if filename_or_block.count('\n') == 0:
        pymol.cmd.load(filename_or_block, 'reference')
    else:
        pymol.cmd.read_pdbstr(filename_or_block, 'reference')


def get_sequences(reference: str) -> List[Tuple[str, str]]:
    """
    Get the sequences of the chains in the PDB file
    :param reference: can be a filename or a PDB block (loaded as 'reference' `by load_reference``)
    :return: list of chain to sequence sorted by length (longest first)
    """
    with pymol2.PyMOL() as pymol:
        load_reference(pymol, reference)
        pymol.cmd.remove('not polymer')
        fasta: str = pymol.cmd.get_fastastr('reference')
        seqs: Dict[str, str] = {re.match('.*_(\w)\n', seq).group(1): ''.join(seq.split('\n')[1:]) for seq in
                                fasta.strip('>').split('>')}
        return sorted(seqs.items(), key=lambda cs: len(cs[1]), reverse=True)


def extract_pdbblocks(reference: str,
                      reference_chain: str,
                    df: pd.DataFrame,
                    fluff_marker: str = '§') \
                    -> Tuple[Dict[str, str], Dict[str, Dict], Dict[str, int]]:
    """
    :param reference:
    :param df:
    :param fluff_marker: name + fluff + index, say MPRO-x0102_0B§1
    :return: pdb_blocks, details, chemcomp_tally
    """
    # safety checks for custom pd.DataFrame
    for key in ('chain', 'pdb_code', 'ligand_codes'):
        assert key in df.columns, f'Column {key} missing in dataframe. This was not outputted by ``ph.LigandHunter``'
    # fill these:
    pdb_blocks: Dict[str, str] = {}
    details: Dict[str, Dict] = {}
    chemcomp_tally = defaultdict(int)
    with pymol2.PyMOL() as pymol, contextlib.suppress(Exception):
        pymol.cmd.set('pdb_conect_all', 'on')
        load_reference(pymol, reference)
        for i, row in df.iterrows():
            pymol.cmd.fetch(row.pdb_code)
            pymol_name = row.pdb_code
            xstal_name = row.pdb_code
            pymol.cmd.align(f'%{row.pdb_code} and chain {row.chain}', f'reference and chain {reference_chain}')
            pymol.cmd.remove(f'%{row.pdb_code} and not (byres reference around 5)')
            for lig_resn in set(row.ligand_codes):
                # ## Determine what is unique ligand residue
                # get atoms of ligand
                lig_atoms: List[chempy.Atom] = pymol.cmd.get_model(f'%{pymol_name} and resn {lig_resn}').atom
                # get set of tuple of resi chain segi alt
                selectors: Set[tuple] = {(atom.resi, atom.chain, atom.segi, atom.alt) for atom in lig_atoms}
                # iterate for all unique
                for resi, chain, segi, alt in selectors:
                    sele = f'%{pymol_name} and resn {lig_resn} and resi {resi}'
                    if chain:
                        sele += f' and chain {chain}'
                    if segi:
                        sele += f' and segi {segi}'
                    if alt:
                        sele += f' and alt {alt}'
                    # MOD: no asymmetric chain business here: too dangerous
                    xsele = f'(not ({sele})) and bound_to ({sele})'
                    covalency = False
                    if pymol.cmd.count_atoms(xsele):
                        # crosslink
                        covalency = True
                        pymol.cmd.alter(xsele, 'elem="T"')  # Chem.MolFromPDBBlock does not like *, R, X
                        pymol.cmd.alter(xsele, 'name = " X  "')
                        pymol.cmd.sort()
                    pymol.cmd.create('copied', f'bound_to ({sele})')
                    pymol.cmd.alter('%copied', 'alt=""')
                    pymol.cmd.alter('%copied', 'segi=""')
                    pymol.cmd.alter('%copied', 'chain="X"')
                    pymol.cmd.alter('%copied', 'resi="1"')
                    pymol.cmd.alter('%copied', f'resn="LIG"')
                    pymol.cmd.sort()
                    # MOD
                    # neoname = f'{xstal_name}{fluff_marker}{i}'
                    chemcomp_tally[lig_resn] += 1
                    neoname = f'{lig_resn}{fluff_marker}{chemcomp_tally[lig_resn]}'
                    pdb_blocks[neoname] = pymol.cmd.get_pdbstr('%copied').replace('ATOM  ', 'HETATM')
                    pymol.cmd.delete('copied')
                    details[neoname] = {**row.to_dict(),  # mod
                                        'covalent': covalency,
                                        **dict(base_name=xstal_name,
                                               lig_resn=lig_resn,
                                               lig_resi=resi,
                                               lig_chain=chain,
                                               lig_segi=segi,
                                               lig_alt=alt,
                                               polymer_chain='A', )  # mod
                                        }
    return pdb_blocks, details, chemcomp_tally

def munge_to_rdkit(pdb_blocks: Dict[str, str],
                   details: Dict[str, Dict],
                   smiles_col_name='single_smiles') -> Tuple[List[Chem.Mol], List[Dict], List[str], List[str]]:
    mols: List[Chem.Mol] = []
    baddies: List[Dict] = []
    inorganics: List[str] = []
    broken: List[str] = []
    for name, block in pdb_blocks.items():
        if block.strip().count('\n') == 0:
            broken.append(name)
            continue
        detail: dict = details[name]
        mol = Chem.MolFromPDBBlock(block, proximityBonding='CONECT' not in block)
        if mol is None:
            # print(f'Issue with sanitisation, trying without for {name}')
            mol = Chem.MolFromPDBBlock(block, proximityBonding='CONECT' not in block, sanitize=False)
        if len(mol.GetAtomsMatchingQuery(AllChem.AtomNumEqualsQueryAtom(6))) == 0:
            inorganics.append(name)
            continue
        tritium = mol.GetAtomsMatchingQuery(AllChem.IsotopeEqualsQueryAtom(3))
        for atom in tritium:  #: Chem.Atom
            atom.SetIsotope(0)
            atom.SetAtomicNum(0)  # dummy atom is zahl zero
        assert mol, f'{name} failed to load'
        # MOD
        # assert mol.GetNumAtoms() > 0 and mol.GetNumBonds() > 0, f'{name} failed to load'
        assert mol.GetNumAtoms() > 0, f'{name} failed to load'
        smiles = detail[smiles_col_name]
        if smiles is None:
            smiles = ''
        ref = Chem.MolFromSmiles(smiles)
        try:
            mol = AllChem.AssignBondOrdersFromTemplate(ref, mol)
        except Exception as e:
            baddies.append(dict(name=name, block=block, mol=mol, ref=ref, detail=detail, exception=e))
            # continue # MOD: These are PDB depositions they have bond order.
        mol.SetProp('XChem_code', name.split('_')[0])
        mol.SetProp('Occupancy', json.dumps([a.GetPDBResidueInfo().GetOccupancy() for a in mol.GetAtoms()]))
        mol.SetProp('TempFactor', json.dumps([a.GetPDBResidueInfo().GetTempFactor() for a in mol.GetAtoms()]))
        unprefixed = name.replace('D68EV3CPROA-', '')
        mol.SetProp('_Name', unprefixed)
        for k, v in detail.items():
            if isinstance(v, dict):
                v = json.dumps(v)
            mol.SetProp(k, str(v))
        for atom in mol.GetAtoms():
            name = atom.GetPDBResidueInfo().GetName()
            atom.SetProp('molFileAlias', name)
        mols.append(mol)
    return mols,  baddies, inorganics, broken


def polymer():
    raise NotImplementedError('this code has not been ported yet')
    # HACK FOR POLYMER

    xstal_name: str = '👾👾👾'

    from collections import defaultdict

    with pymol2.PyMOL() as pymol:
        pymol.cmd.load(reference, 'reference')
        xstal_name = '2b0f'
        pymol.cmd.fetch(xstal_name)
        pymol.cmd.align(f'%{xstal_name} and chain A', 'reference')
        pymol.cmd.delete('reference')
        pymol.cmd.remove('chain A')
        atom_tally = defaultdict(int)
        for atom in pymol.cmd.get_model('chain B').atom:
            atom_tally[atom.symbol] += 1
            pymol.cmd.alter(f'name {atom.name}', f'name="{atom.symbol: >2}{atom_tally[atom.symbol]: <2}"')
        pymol.cmd.alter('*', 'alt=""')
        pymol.cmd.alter('*', 'segi=""')
        pymol.cmd.alter('*', 'chain="X"')
        pymol.cmd.alter('*', 'resi="1"')
        pymol.cmd.alter('*', 'resn="LIG"')
        pymol.cmd.sort()
        neoname = 'polymer'
        pymol.cmd.set('pdb_conect_all', 'on')
        pdb_blocks[neoname] = pymol.cmd.get_pdbstr('*').replace('ATOM  ', 'HETATM')
        details[neoname] = dict(base_name=xstal_name,
                                lig_resn=None,
                                lig_resi=1,
                                lig_chain='B',
                                lig_segi='',
                                lig_alt='',
                                polymer_chain='A', )

def get_homolog_chemcomps(reference: str,
                          blast_filename: str = 'blast_hits.csv',
                          image_filename: str = 'compounds.png',
                          sdf_filename: str = 'compounds.sdf',
                          error_sdf_filename: str = 'errors.sdf',
                          log_level: int = logging.INFO,
                          fluff_marker='§') -> None:
    """
    Master function.

    :param reference: can be a filename or a PDB block (loaded as 'reference' `by load_reference``)
    :return:
    """
    log = get_log(log_level)
    log.debug('Exacting sequences...')
    seqs: List[Tuple[str, str]] = get_sequences(reference)
    chain, seq = seqs[0]
    log.info(f'Structure contained {len(seqs)}, using chain {chain} which is {len(seq)} residues long')
    # ----------------------------------
    # blast
    hunter = LigandHunter(seq)
    log.debug('Most common ligands:')
    for lig, c in hunter.get_most_common_ligands()[:10]:
        log.debug(lig, c, hunter.ligand_data[lig][0]['name'])
    structures: pd.DataFrame = hunter.to_dataframe()
    structures.to_csv(blast_filename)
    # ----------------------------------
    log.info(f'Extracting from {len(structures)} structures')
    pdb_blocks: Dict[str, str]  # accession: pdbblock
    details: Dict[str, Dict]    # accession: details
    chemcomp_tally: Dict[str, int]  # ligand (chem comp name): tally
    pdb_blocks, details, chemcomp_tally = extract_pdbblocks(reference=reference,
                                                            reference_chain=chain,
                                                            df=structures,
                                                            fluff_marker=fluff_marker)
    smilesdex = {detail['lig_resn']: get_smiles(detail['lig_resn']) for detail in details.values() if
                 detail['lig_resn'] != 'LIG' and detail['lig_resn'] != 'UNL' and detail['lig_resn'] is not None}
    log.info(f'Found {sum(chemcomp_tally.values())} ligands ({len(chemcomp_tally)} unique)')
    for detail in details.values():
        detail['single_smiles'] = smilesdex.get(detail['lig_resn'], None)
    # ----------------------------------
    mols: List[Chem.Mol]
    baddies: List[Dict]
    inorganics: List[str]
    broken: List[str]
    mols, baddies, inorganics, broken = munge_to_rdkit(pdb_blocks, details)
    log.info(f'Inorganics: {inorganics}')
    log.info(f'Broken: {broken}')
    for bad in baddies:
        log.warning(f'Issue with bond orders for {bad["name"]} - ' +
                    f'{bad["exception"].__class__.__name__}{bad["exception"]}')
    assert mols, 'No mols!'
    image = Draw.MolsToGridImage(mols)
    if isinstance(image, IPyImage):  # we are in a notebook
        with open(image_filename, 'wb') as fh:
            fh.write(image.data)
    elif isinstance(image, PngImageFile):
        image.save(image_filename)
    else:
        log.critical('Could not save image - no idea where we are!')
    with Chem.SDWriter(sdf_filename) as sdf:
        for mol in sorted(mols, key=Chem.Mol.GetNumHeavyAtoms, reverse=True):
            if mol is None:
                continue
            sdf.write(mol)
    with Chem.SDWriter(error_sdf_filename) as sdf:
        for bad in baddies:
            if bad.mol is None:
                continue
            sdf.write(bad.mol)
    log.debug('Completed.')
    return

# ---------------------------------------------------------
# copied straight out of pyrosetta-help
from collections import Counter
from io import StringIO
from typing import (Tuple, Dict, List)

import pandas as pd
import requests
from Bio import SearchIO
from Bio.Blast.NCBIWWW import qblast

def get_smiles(ligand_code: str) -> str:
    """
    Get the smiles of a ligand.
    Remember that PDBe smiles need to charged to pH 7.
    """
    ligand_code = ligand_code.upper()
    ligand_data = requests.get(f'https://www.ebi.ac.uk/pdbe/api/pdb/compound/summary/{ligand_code}').json()
    return ligand_data[ligand_code][0]['smiles'][0]['name']
class LigandHunter:
    """
    Given a sequence find homologues (Blast) and find if they have ligands (PDBe API).
    The definition of what ligand is a cofactor comes from PDBe and does not count ions or
    triphospho-nucleotides as cofactors, but does count NADH adducts.

    * ``.data`` is a list of dictionaries.
    * ``.to_dataframe()`` converts into a pandas dataframe for further analysis.
    * ``.candidate_ligands`` list of ligand residue 3-letter codes
    """
    _cofactor_reference = {}

    @property
    def cofactor_reference(self) -> Dict[str, List[dict]]:
        if not self._cofactor_reference: #cached
            self._cofactor_reference = requests.get('https://www.ebi.ac.uk/pdbe/api/pdb/compound/cofactors').json()
        return self._cofactor_reference

    missing_codes = ['ATP', 'GTP', 'CA', 'MG', 'W']

    @property
    def cofactor_codes(self):
        _grouped_cofactor_codes = {name: value[0]['cofactors'] for name, value in self.cofactor_reference.items()}
        return [code for codes in _grouped_cofactor_codes.values() for code in codes] + self.missing_codes

    def __init__(self, sequence: str):
        """

        :param sequence: is assumed clean protein sequence
        """
        self.sequence = sequence
        self._blast_result = self._retrieve_homologues()  #:StringIO
        self.data = self._parse_blast_result()
        self._retrieve_ligands()  # inplace.
        self._ligand_data = None

    def to_dataframe(self):
        """
        Converts ``.data`` to a pandas dataframe
        """
        return pd.DataFrame(self.data).transpose()

    @property
    def candidate_ligands(self) -> List[str]:
        return list(set([code for datum in self.data.values() for code in datum['ligand_codes']]))

    @property  # cached the old way
    def ligand_data(self):
        """
        Data for the ligands.

        Note that LigandNicker has a get smiles method, that works like this, but is unrelated.
        """
        if self._ligand_data is None:
            self._ligand_data = requests.post('https://www.ebi.ac.uk/pdbe/api/pdb/compound/summary/',
                                              data=','.join(self.candidate_ligands)).json()
        return self._ligand_data

    def get_most_common_ligands(self) -> List[Tuple[str, int]]:
        # Uses collections.Counter not typing.Counter:
        c = Counter([code for datum in self.data.values() for code in datum['ligand_codes']])  # noqa
        return c.most_common()

    def get_pdb_entry_by_ligand(self, ligand_code: str) -> dict:
        """
        get pdb **entry** by ligand.
        Returns the first, which should be lowest e-value
        """
        for datum in self.data.values():
            if ligand_code in datum['ligand_codes']:
                return datum
        else:
            raise ValueError(f'{ligand_code} not found in any of the {len(self.data)} hits')

    # ------------ initialisation methods ------------------------

    def _retrieve_homologues(self) -> StringIO:
        return qblast('blastp', 'pdb', self.sequence)

    def _parse_blast_result(self) -> Dict[str, dict]:
        results = {}
        self._blast_result.seek(0)
        for query in SearchIO.parse(self._blast_result, 'blast-xml'):
            for hit in query:
                datum = dict(accession=hit.accession,
                             description=hit.description,
                             evalue=hit.hsps[0].evalue)
                datum['pdb_code'], datum['chain'] = hit.accession.split('_')
                results[datum['pdb_code'].upper()] = datum
        return results

    def _retrieve_ligands(self) -> None:
        query = ','.join(self.data.keys())
        reply = requests.post(url='https://www.ebi.ac.uk/pdbe/api/pdb/entry/ligand_monomers/',
                              data=query)
        for code, datum in reply.json().items():
            entry = self.data[code.upper()]
            entry['ligands'] = datum
            entry['ligand_codes'] = [inner['chem_comp_id'] for inner in datum]
            entry['cofactor_codes'] = [code for code in entry['ligand_codes'] if code in self.cofactor_codes]
            entry['has_cofactor'] = bool(entry['cofactor_codes'])
        # TODO find the part of the code that is somehow truncated between commits!!


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Get chemical components from homologs in the PDB')
    parser.add_argument('-r', '--reference', help='Reference PDB file or block', required=True)
    parser.add_argument('-b', '--blast', help='Blast output file', default='blast_hits.csv')
    parser.add_argument('-i', '--image', help='Image output file', default='compounds.png')
    parser.add_argument('-s', '--sdf', help='SDF output file', default='compounds.sdf')
    parser.add_argument('-e', '--error_sdf', help='Error SDF output file', default='errors.sdf')
    parser.add_argument('-l', '--log_level', help='Log level', default=logging.INFO, type=int)
    parser.add_argument('-f', '--fluff_marker', help='Fluff marker', default='§')
    args = parser.parse_args()

    get_homolog_chemcomps(reference=args.reference,
                          blast_filename=args.blast,
                          image_filename=args.image,
                          sdf_filename=args.sdf,
                          error_sdf_filename=args.error_sdf,
                            log_level=args.log_level,
                            fluff_marker=args.fluff_marker)
