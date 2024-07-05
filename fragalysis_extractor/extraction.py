__all__ = ['Extractor']

import re
from pathlib import Path
import pymol2, chempy  # noqa
from typing import List, Set, Dict, Optional, Tuple
from collections import defaultdict
import logging

logger = logging.getLogger(__name__)

class Extractor:

    def __init__(self, folder: Path, ref_path: Optional[Path] = None,
                 monomer_chains: str = 'A',
                 lig_resn='LIG', fluff_marker='_'
                 ):
        self.folder = Path(folder)
        self.paths, self.reference_path = self.parse_paths(folder, ref_path)
        self.pymol = None
        self.monomer_chains = monomer_chains
        self.lig_resn = lig_resn
        self.fluff_marker = fluff_marker

    def __call__(self, pymol=None):
        if pymol is None:
            with pymol2.PyMOL() as pymol:
                self.pymol = pymol
                return self.collate_details()
        else:
            return self.collate_details()

    def parse_paths(self, folder: Optional[Path]=None, ref_path: Optional[Path] = None) -> Tuple[List[Path], Path]:
        if folder is None:
            folder = self.folder
        for subfolder in  ('crystallographic', 'combined'):
            if (folder / subfolder).exists():
                logger.debug('Using crystallographic folder')
                paths = sorted(Path(folder / subfolder).glob('*_bound.pdb'))
                if len(paths) == 0:
                    paths = sorted(Path(folder / subfolder).glob('*.pdb'))
                break
        else:
            if (folder / 'aligned').exists():
                paths = sorted(Path(folder / 'aligned').glob('*/*_bound.pdb'))
                # if path.stem[:2] != '._' and 'min' not in path.stem
            else:
                raise Exception('NEW CASE')
        assert len(paths), 'Could not find PDB files'
        if ref_path is not None:
            return paths, ref_path
        elif hasattr(self, 'reference_path') and self.reference_path.exists():
            return paths, self.reference_path
        else:
            ref_path = folder / 'reference.pdb'
            assert ref_path.exists(), ref_path.as_posix()
            return paths, ref_path


    def get_detail(self, pymol_name, resn, resi, chain, segi, alt):
        detail = {}
        detail['crystal_name']: str = pymol_name
        sele = f'%{pymol_name} and resn {resn} and resi {resi}'
        if chain:
            sele += f' and chain {chain}'
        if segi:
            sele += f' and segi {segi}'
        if alt:
            sele += f' and alt {alt}'
        # logger.debug(pymol_name, sele)
        self.pymol.cmd.create('copied', sele)
        self.pymol.cmd.alter('%copied', 'alt=""')
        self.pymol.cmd.alter('%copied', 'segi=""')
        self.pymol.cmd.alter('%copied', 'chain="X"')
        self.pymol.cmd.alter('%copied', 'resi="1"')
        self.pymol.cmd.sort()
        detail['pdb_block'] = self.pymol.cmd.get_pdbstr('%copied')
        detail['names'] = [a.name for a in self.pymol.cmd.get_model('%copied').atom]
        detail['occupancies'] = [a.q for a in self.pymol.cmd.get_model('%copied').atom]
        detail['b-factors'] = [a.b for a in self.pymol.cmd.get_model('%copied').atom]
        self.pymol.cmd.delete('%copied')
        detail.update(dict(
            lig_resn=resn,
            lig_resi=resi,
            lig_chain=chain,
            lig_segi=segi,
            lig_alt=alt)
        )
        return detail


    def get_objects(self):
        return [pymol_name for pymol_name in self.pymol.cmd.get_names('objects') if pymol_name != 'reference']


    def get_details(self, pymol_name):
        details = []
        lig_atoms: List[chempy.Atom] = self.pymol.cmd.get_model(f'resn {self.lig_resn} and {pymol_name}').atom
        selectors: Set[tuple] = {(atom.resn, atom.resi, atom.chain, atom.segi, atom.alt) for atom in lig_atoms}
        for resn, resi, chain, segi, alt in sorted(selectors, key=lambda s: s[2] + s[0]):
            detail = self.get_detail(pymol_name, resn, resi, chain, segi, alt)
            details.append(detail)
        return details


    def collate_details(self) -> List[Dict]:
        """
        Use PyMOL across a Fragalsis folder

        monomer_chains as in heteromer
        """
        # ## Figure out structure
        logger.debug(len(self.paths), 'found')
        pathdex = {re.search(r'(x\d+)', path.stem).group(1): path for path in self.paths}
        # ## Parse
        details = []
        count = defaultdict(int)
        mono_sele = '(chain ' + "+".join(self.monomer_chains) + ' and polymer)'
        # ## load
        logger.debug('loading')
        self.pymol.cmd.load(str(self.reference_path), 'reference')
        self.pymol.cmd.remove(f'reference and not chain {"+".join(self.monomer_chains)}')
        self.pymol.cmd.remove(f'reference and not polymer')
        assert self.pymol.cmd.count_atoms(f'reference'), 'reference empty'
        logger.debug(f'loaded reference')
        path: Path
        for pymol_name, path in pathdex.items():
            if pymol_name in self.get_objects():
                continue
            self.pymol.cmd.load(path.as_posix(), pymol_name)
            self.pymol.cmd.align(f'{pymol_name}', f'reference and {mono_sele}')
            logger.debug(f'loaded {pymol_name}')
            assert self.pymol.cmd.count_atoms(pymol_name), f'{pymol_name} empty'
            self.pymol.cmd.remove('inorganic or solvent')
            other_chains = self.pymol.cmd.get_chains(f'{pymol_name} and polymer')
            pymol_names = [pymol_name]
            if set(other_chains).difference(self.monomer_chains) == set():
                # there are no extra chains...
                pass
            elif len(self.monomer_chains) > 1:
                raise Exception(f"Too complicated: expected {self.monomer_chains} but got also {other_chains}")
            else:
                for other in other_chains:
                    self.pymol.cmd.create(f'{pymol_name}_{other}',
                                     f'({pymol_name} and chain {other}) or (resn {self.lig_resn} and byres ({pymol_name} and chain {other}) around 4)')
                    pymol_names.append(f'{pymol_name}_{other}')
                    self.pymol.cmd.align(f'{pymol_name}_{other} and chain {other}', f'reference and {mono_sele}')
                    self.pymol.cmd.remove(
                        f'{pymol_name}_{other} and not (byres (reference and {mono_sele}) around 3.5)')  # kill ligands from the other chain
                    self.pymol.cmd.remove(f'{pymol_name}_{other} and polymer and not reference')  # I dont need protein
                    details.extend(self.get_details(f'{pymol_name}_{other}'))
                    self.pymol.cmd.delete(f'{pymol_name}_{other}')
            self.pymol.cmd.remove(
                f'{pymol_name} and not (byres (reference and {mono_sele}) around 4)')  # kill ligands from the other chain
            self.pymol.cmd.remove(f'{pymol_name} and polymer and not reference')  # I dont need protein
            details.extend(self.get_details(pymol_name))
            self.pymol.cmd.delete(pymol_name)
        # annotate number
        for detail in details:
            count[detail['crystal_name']] += 1
            detail['crystal_name'] = detail['crystal_name'] + self.fluff_marker + str(count[detail['crystal_name']])
        return details
