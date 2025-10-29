"""
Render SMILES (Simplified Molecular Input Line Entry System) strings to images.

This module provides functionality to convert SMILES strings into visual molecular representations
with coloring based on molecular fragments.
"""

import argparse
from collections import namedtuple
from pathlib import Path
from . import renderer
from . import postprocessing
from . import Options


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input", help="File containing a SMILES string")
    group = parser.add_mutually_exclusive_group(required=False)
    group.add_argument(
        "--svg", dest="svg", action="store_true", default=True, help="Output SVG files (default)"
    )
    group.add_argument(
        "--png", dest="png", action="store_true", default=False, help="Output PNG files"
    )
    args = parser.parse_args()

    input_path = Path(args.input)
    smiles = input_path.read_text().strip()
    sample = create_sample(smiles)

    if args.svg:
        write_SVGs(sample, input_path)
    if args.png:
        write_PNGs(sample, input_path)


Sample = namedtuple("Sample", ["smiles", "mol", "svg", "svg_colored"])
"""
A data structure containing molecular representation data.

Fields:
    smiles (str): The original SMILES string
    mol (Chem.Mol): RDKit molecule object
    svg (str): Standard SVG representation of the molecule
    svg_colored (str): Colored SVG representation with unique colors for different molecular fragments
"""


def create_sample(smiles: str, options: Options | None = None) -> Sample:
    """
    Create a Sample object containing molecular data and SVG representations.

    Args:
        smiles (str): A valid SMILES string representing a molecule
        options (Options, optional): Rendering options for molecular diagram generation.
            If None, default Options will be used.

    Returns:
        Sample: A namedtuple containing the original SMILES string, RDKit molecule object,
                standard SVG representation, and colored SVG representation
    """
    if options is None:
        options = Options()
    mol = renderer.create_mol(smiles)
    svg = renderer.create_svg(mol, options)
    svg = postprocessing.post_process_svg(svg, mol, options)
    svg_colored = postprocessing.color_instances(svg)
    return Sample(smiles, mol, svg, svg_colored)


def write_SVGs(sample: Sample, input_path: Path):
    """
    Write SVG representations to files.

    Creates two SVG files:
    - A standard black and white molecular diagram
    - A colored version with unique colors for different molecular fragments

    Args:
        sample (Sample): Sample object containing the SVG data
        input_path (Path): Original input file path, used to determine output filenames
    """
    with open(input_path.with_suffix(".svg"), "w") as f:
        f.write(sample.svg)
    with open(input_path.parent / (input_path.stem + "_colored.svg"), "w") as f:
        f.write(sample.svg_colored)


def write_PNGs(sample: Sample, input_path: Path):
    """
    Write PNG representations to files.

    Creates two PNG files by converting the SVG representations:
    - A standard black and white molecular diagram
    - A colored version with unique colors for different molecular fragments

    Args:
        sample (Sample): Sample object containing the SVG data to convert
        input_path (Path): Original input file path, used to determine output filenames
    """
    image = renderer.svg_to_pil(sample.svg)
    with open(input_path.with_suffix(".png"), "wb") as f:
        image.save(f, format="PNG")

    image_colored = renderer.svg_to_pil(sample.svg_colored)
    with open(input_path.parent / (input_path.stem + "_colored.png"), "wb") as f:
        image_colored.save(f, format="PNG")


if __name__ == "__main__":
    main()
