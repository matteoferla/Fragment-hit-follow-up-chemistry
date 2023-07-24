def display_source(function):
    """
    Display the source code formatted
    """
    import inspect
    from IPython.display import HTML, display
    from pygments import highlight
    from pygments.lexers import PythonLexer
    from pygments.formatters import HtmlFormatter
    code:str = inspect.getsource(function)
    html:str = highlight(code, PythonLexer(), HtmlFormatter(style='colorful'))
    stylesheet:str = f"<style>{HtmlFormatter().get_style_defs('.highlight')}</style>"
    display(HTML(f"{stylesheet}{html}"))

# --------------------------------------------------------------------------------------------
from rdkit import Chem
from rdkit.Chem import Draw
from typing import List
from IPython.display import Image

def flatgrid(mols, *args, **kwargs) -> Image:
    # from https://www.blopig.com/blog/2023/06/customising-mcs-mapping-in-rdkit/
    copies: List[Chem.Mol] = [Chem.Mol(m) for m in mols]
    *map(AllChem.Compute2DCoords, copies),   # noqa, it's in place
    if 'legends' not in kwargs:
        kwargs['legends'] = [m.GetProp('_Name') if m.HasProp('_Name') else '-' for m in mols]
    return Draw.MolsToGridImage(copies, *args, **kwargs)


# --------------------------------------------------------------------------------------------

import importlib.util
from types import ModuleType

def import_path(module_path, module_name='custom_module') -> ModuleType:
    """
    Import a module from a path.

    :param module_path:
    :param module_name:
    :return:
    """
    spec = importlib.util.spec_from_file_location(module_name, module_path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod
