from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D
from . import Options


def create_mol(smiles: str) -> Chem.Mol:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES string: {smiles}")
    mol = rdMolDraw2D.PrepareMolForDrawing(mol)
    return mol


def create_svg(mol: Chem.Mol, options: Options) -> str:
    """
    Generate an SVG representation of a molecule.

    Creates an SVG drawing of the molecule using RDKit's drawing capabilities,
    applying the specified options for appearance and layout.

    Args:
        mol (Chem.Mol): RDKit molecule object to render
        options (Options): Configuration options for rendering

    Returns:
        str: SVG representation as a string
    """
    opts = rdMolDraw2D.MolDrawOptions()
    opts.useBWAtomPalette()
    opts.setBackgroundColour(options.backgroundColour)

    for option_name, option_value in vars(options).items():
        if option_value is None:
            continue
        try:
            setattr(opts, option_name, option_value)
        except AttributeError as e:
            if str(e).startswith("Cannot set unknown attribute"):
                # Skip options that don't exist in rdMolDraw2D.MolDrawOptions
                continue
            raise e

    d = rdMolDraw2D.MolDraw2DSVG(options.width, options.height)
    d.SetDrawOptions(opts)
    rdMolDraw2D.PrepareAndDrawMolecule(d, mol)
    d.FinishDrawing()
    return d.GetDrawingText()


def svg_to_pil(svg: str):
    """
    Convert an SVG string to a PIL Image object.

    Uses CairoSVG to convert the SVG to PNG data, then loads it as a PIL Image.
    This is useful for generating raster image outputs from vector molecular diagrams.

    Args:
        svg (str): SVG content as a string

    Returns:
        PIL.Image.Image: PIL Image object containing the rasterized molecular diagram

    Raises:
        ValueError: If SVG conversion fails
        ImportError: If required dependencies (cairosvg, PIL) are not available
    """
    from cairosvg import svg2png
    from PIL import Image
    from io import BytesIO

    png_data = svg2png(bytestring=svg.encode("utf-8"))
    if png_data is None:
        raise ValueError("Failed to convert SVG to PNG")
    return Image.open(BytesIO(png_data))
