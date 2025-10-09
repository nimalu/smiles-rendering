from dataclasses import dataclass
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D


@dataclass
class SVGOptions:
    width: int = 500
    height: int = 500


def create_svg(mol: Chem.Mol, options: SVGOptions) -> str:
    opts = rdMolDraw2D.MolDrawOptions()
    opts.useBWAtomPalette()
    d = rdMolDraw2D.MolDraw2DSVG(options.width, options.height)
    d.SetDrawOptions(opts)
    rdMolDraw2D.PrepareAndDrawMolecule(d, mol)
    d.FinishDrawing()
    return d.GetDrawingText()
