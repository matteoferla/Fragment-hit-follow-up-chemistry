import requests
def smiles_to_iupac(smiles: str, raise_error: bool = True) -> str:
    """
    Given a SMILES, get the IUPAC name from NIH Cactus server
    """
    rep = "iupac_name"
    CACTUS = "https://cactus.nci.nih.gov/chemical/structure/{0}/{1}"
    url = CACTUS.format(smiles, rep)
    response = requests.get(url)
    if response.ok:
        return response.text
    elif raise_error:
        response.raise_for_status()
    else:
        return ''