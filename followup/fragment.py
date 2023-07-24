from rdkit import Chem
from rdkit.Chem import BRICS, AllChem



def fragment(mol: Chem.Mol):
    fragments = []
    if mol is None:
        return []
    name = mol.GetProp('_Name')
    for atom in mol.GetAtoms():
        atom.SetProp('brics_ori_name', name)
        atom.SetIntProp('brics_ori_i', atom.GetIdx())
    for dummy in mol.GetAtomsMatchingQuery(AllChem.AtomNumEqualsQueryAtom(0)):
        # tritium
        dummy.SetAtomicNum(85)  # astatine just becase
        dummy.SetBoolProp('dummy_atom', True)
    i = 0
    for fragment in BRICS.BRICSDecompose(mol, keepNonLeafNodes=True, returnMols=True, minFragmentSize=4):
        dummies = list(fragment.GetAtomsMatchingQuery(AllChem.AtomNumEqualsQueryAtom(0)))
        if len(dummies) == 0:
            fragment.SetProp('_Name', name)
        else:
            i+=1
            fragment.SetProp('_Name', f'{name}S{i}')
        for dummy in dummies:
            dummy.SetAtomicNum(1)  # hydrogen!
            dummy.SetIsotope(0)
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
        for atom in fragment.GetAtomsMatchingQuery(AllChem.AtomNumEqualsQueryAtom(85)):
            atom.SetAtomicNum(0)
            atom.SetIsotope(0)
        fragments.append(AllChem.RemoveHs(fragment))
    return fragments