from rdkit import Chem, RDLogger
from rdkit.Chem import BRICS, AllChem
from typing import List, Sequence
from collections import Counter
from molecular_rectifier import Rectifier

def _get_fused(mol):
    ri = mol.GetRingInfo()
    rings_idxs = ri.AtomRings()
    if not rings_idxs:
        return [], ()
    commons = [i for i, c in Counter([i for r in rings_idxs for i in r]).items() if c > 1]
    return commons, tuple([r for r in rings_idxs if set(commons).union(r)])

def _get_ring(mol, commons: List[int], rings_idxs: Sequence[Sequence[int]], keep_ring: int) -> Chem.Mol:
    """
    Required by split_fused

    :param mol:
    :param commons: see _get_fused
    :param rings_idxs: see _get_fused
    :param keep_ring: index of ring to keep
    :return:
    """
    m = Chem.RWMol(mol)
    m.BeginBatchEdit()
    aromaticities = [all([m.GetAtomWithIdx(i).GetIsAromatic() for i in ring_idxs if i not in commons]) for ring_idxs in rings_idxs]
    for rn, ring_idxs  in enumerate(rings_idxs):
        for i in ring_idxs:
            if i in commons:
                continue
            if rn != keep_ring:
                m.RemoveAtom(i)
            else:
                m.GetAtomWithIdx(i).SetBoolProp('_nice_ring', True)
    if not aromaticities[keep_ring] and any(aromaticities):
        for c in commons:
            m.GetAtomWithIdx(c).SetIsAromatic(False)
        if len(commons) == 2: # Really? I mean a spiro could not have been aromatic...
            m.GetBondBetweenAtoms(*commons).SetBondType(Chem.BondType.DOUBLE)
    m.CommitBatchEdit()
    mol = m.GetMol()
    try:
        ms = Chem.GetMolFrags(mol, asMols=True)
    except Exception as error:
        mol = Rectifier(mol).fix().mol
        ms = Chem.GetMolFrags(mol, asMols=True)
    return [m for m in ms if any(a.HasProp('_nice_ring') for a in m.GetAtoms())][0]

def split_fused(mol) -> List[Chem.Mol]:
    """
    Return the possible fused/spiro ring splits if present.
    (If there are no fused rings it returns an empty list)

    :param mol:
    :return:
    """
    commons, rings_idxs = _get_fused(mol)
    if not commons:
        return []
    return [_get_ring(mol, commons, rings_idxs, k) for k in range(len(rings_idxs))]

def fragment(mol: Chem.Mol, minFragmentSize=4,
             fused_splitting=True,
             remove_new_dummies=True,
             temp_zahl=85) -> List[Chem.Mol]:
    """
    Wrapper for BRICSDecompose, which keeps metadata.

    BRICS is good, but has two problems for my use.
    It does not keep the metadata and it does not split fused/spiro rings, which is understandable,
    but in this case I want to split them.

    :param mol:
    :param minFragmentSize:
    :param fused_splitting: Split fused/spiro rings? See ``split_fused``
    :param remove_new_dummies:
    :param temp_zahl: 85. Prevent pre-existant dummies from being removed. 85 is astatine. No reason.
    :return:

    NB. renamed from fragmént (fragmént = verb, frágment = noun).
    which does prevent overriding... but requires deadkeys for typing.
    Also, I am not sure if stress changes in US pronunciation from noun to verb.
    """
    if mol is None:
        return []
    RDLogger.DisableLog('rdApp.*')
    name = mol.GetProp('_Name')
    for atom in mol.GetAtoms():
        atom.SetProp('brics_ori_name', name)
        atom.SetIntProp('brics_ori_i', atom.GetIdx())
    for dummy in mol.GetAtomsMatchingQuery(AllChem.AtomNumEqualsQueryAtom(0)):
        dummy.SetAtomicNum(temp_zahl)
        dummy.SetBoolProp('dummy_atom', True)
    i = 0
    raw_fragments = []
    fragments = []
    # BRICS
    for fragment in BRICS.BRICSDecompose(mol, keepNonLeafNodes=True, returnMols=True, minFragmentSize=minFragmentSize):
        dummies = list(fragment.GetAtomsMatchingQuery(AllChem.AtomNumEqualsQueryAtom(0)))
        if len(dummies) == 0:
            fragment.SetProp('_Name', name)
        else:
            i+=1
            fragment.SetProp('_Name', f'{name}S{i}')
        for dummy in fragment.GetAtomsMatchingQuery(AllChem.AtomNumEqualsQueryAtom(0)):
            dummy.SetAtomicNum(1)  # hydrogen!
            dummy.SetIsotope(0)
        raw_fragments.append(AllChem.RemoveAllHs(fragment))
    # Fused split
    if fused_splitting:
        for fragment in split_fused(mol):
            i+=1
            fragment.SetProp('_Name', f'{name}S{i}')
            if fragment.GetNumAtoms() < minFragmentSize:
                continue
            raw_fragments.append(AllChem.RemoveAllHs(fragment))
    # fix metadata
    for fragment in raw_fragments:
        for k, v in mol.GetPropsAsDict().items():
            if isinstance(v, str):
                fragment.SetProp(k, v)
            elif isinstance(v, int):
                fragment.SetIntProp(k, v)
            elif isinstance(v, float):
                fragment.SetDoubleProp(k, v)
            elif isinstance(v, bool):
                fragment.SetBoolProp(k, v)
            else:
                fragment.SetProp(k, str(v))
        for atom in fragment.GetAtomsMatchingQuery(AllChem.AtomNumEqualsQueryAtom(temp_zahl)):
            atom.SetAtomicNum(0)
            atom.SetIsotope(0)
        fragments.append(fragment)
    RDLogger.EnableLog('rdApp.*')
    return fragments


def remove_duplicated(hits: List[Chem.Mol]) -> List[Chem.Mol]:
    """
    Crude deduplication based on InChIKey and coordinates sum.

    :param hits:
    :return:
    """
    seen = []
    def is_new(hit):
        h = hash(f'{Chem.MolToInchiKey(hit)}{hit.GetConformer().GetPositions().sum():.0f}')
        if h in seen:
            return False
        seen.append(h)
        return True

    return [hit for hit in hits if is_new(hit)]

fragmént = fragment

def main():
    import argparse
    from rdkit import Chem

    parser = argparse.ArgumentParser(description='Fragment a molecule')
    parser.add_argument('sdf', type=str, help='sdf file')
    parser.add_argument('output', type=str, help='sdf outfile')
    parser.add_argument('--minFragmentSize', type=int, default=7, help='Minimum fragment size')

    args = parser.parse_args()

    with Chem.SDMolSupplier(args.sdf) as suppl:
        hits: List[Chem.Mol] = [mol for mol in suppl if mol is not None]
    fhits: List[Chem.Mol] = []
    for hit in hits:
        fhits.extend(fragmént(hit, minFragmentSize=args.minFragmentSize, fused_splitting=True))
    fhits: List[Chem.Mol] = remove_duplicated(fhits)
    with Chem.SDWriter(args.output) as writer:
        for hit in fhits:
            writer.write(hit)

if __name__ == '__main__':
    main()