"""
Genereate a conformer library from a set of molecules for use in OpenEye ROCS
"""

import os
from pathlib import Path

assert Path(os.environ['OE_LICENSE']).exists(), f'$OE_LICENSE={os.environ["OE_LICENSE"]} does not exist'
from openeye import oechem
from openeye import oeomega


def init_omega(dense=True):
    """
    Initialise omega

    :param dense: use ``oeomega.OEOmegaSampling_Dense``?
    :return:
    """
    omegaOpts = oeomega.OEOmegaOptions(oeomega.OEOmegaSampling_Dense) if dense else oeomega.OEOmegaOptions()
    opts = oechem.OESimpleAppOptions(omegaOpts, "Omega", oechem.OEFileStringType_Mol, oechem.OEFileStringType_Mol3D)
    # print(oechem.OEConfigureOpts(opts, [], False))
    omegaOpts.UpdateValues(opts)
    omega = oeomega.OEOmega(omegaOpts)

    flipperOpts = oeomega.OEFlipperOptions()
    opts = oechem.OESimpleAppOptions(flipperOpts, "stereo_and_torsion", oechem.OEFileStringType_Mol,
                                     oechem.OEFileStringType_Mol)
    flipperOpts.UpdateValues(opts)
    return omega, flipperOpts


def sdf_to_oeb(in_filename: str, out_filename: str, dense=True, tally_printable=1_000):
    """
    Convert sdf to oeb

    :param in_filename:
    :param out_filename:
    :param dense: use ``oeomega.OEOmegaSampling_Dense``?
    :param tally_printable:
    :return:
    """
    omega, flipperOpts = init_omega(dense)
    # error checking: tally of return codes
    ret_codes = {getattr(oeomega, v): v.lstrip('OEOmegaReturnCode_') for v in dir(oeomega) if 'OEOmegaReturnCode' in v}
    ret_code_tallies = {c: 0 for c in ret_codes.values()}

    # inner
    def gen_conf_inner(mol: oechem.OEMol) -> None:
        ret_code = omega.Build(mol)
        ret_code_tallies[ret_codes[ret_code]] += 1
        if ret_code == oeomega.OEOmegaReturnCode_UnspecifiedStereo:
            for enantiomer in oeomega.OEFlipper(mol.GetActive(), flipperOpts):
                gen_conf_inner(mol)
        elif ret_code == oeomega.OEOmegaReturnCode_Success:
            oechem.OEWriteMolecule(ofs, mol)
        else:
            pass

    # I/O
    ifs = oechem.oemolistream()
    assert ifs.open(in_filename)
    ofs = oechem.oemolostream()
    assert ofs.open(out_filename)
    # iterate
    mol: oechem.OEMol
    c = 0
    for i, mol in enumerate(ifs.GetOEMols()):
        # tally
        if i % tally_printable == 0:
            print(i, flush=True)
        # skip overly linear
        if oechem.OECount(mol, oechem.OEIsRotor()) > 5:
            continue
        # fix enamine title
        identifier = oechem.OEGetSDData(mol, 'id')
        mol.SetTitle(identifier)
        # conf gen
        mol.ClearCoords()
        oechem.OEAddExplicitHydrogens(mol)
        gen_conf_inner(mol)
        # the `gen_conf_inner` fun calls `oechem.OEWriteMolecule`
        c += mol.GetMaxConfIdx()
    # done
    print(f'{i} molecules converted into {c} conformer. The following return codes were seen: {ret_code_tallies}')
    ofs.close()


def oeb_to_sdf(in_filename: str, out_filename: str, tally_printable=1_000):
    """
    Convert oeb to sdf
    """
    # I/O
    ifs = oechem.oemolistream()
    assert ifs.open(in_filename)
    ofs = oechem.oemolostream()
    ofs.SetFormat(oechem.OEFormat_SDF)
    assert ofs.open(out_filename)
    # iterate
    mol: oechem.OEMol
    for i, mol in enumerate(ifs.GetOEMols()):
        # tally
        if i % tally_printable == 0:
            print(i, flush=True)
        oechem.OEWriteMolecule(mol)
    # done
    print(f'{i} molecules saved')
    ofs.close()

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Convert sdf to oeb')
    parser.add_argument('in_filename', type=str, help='input filename')
    parser.add_argument('out_filename', type=str, help='output filename')
    args = parser.parse_args()
    oeb_to_sdf(in_filename=args.in_filename, out_filename=args.out_filename)
