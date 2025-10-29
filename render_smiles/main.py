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
        output_svg(sample, input_path)
    if args.png:
        output_png(sample, input_path)


Sample = namedtuple("Sample", ["smiles", "mol", "svg", "svg_colored"])


def create_sample(smiles: str):
    options = Options()
    mol = renderer.create_mol(smiles)
    svg = renderer.create_svg(mol, options)
    svg = postprocessing.post_process_svg(svg, mol, options)
    svg_colored = postprocessing.color_instances(svg)
    return Sample(smiles, mol, svg, svg_colored)


def output_svg(sample: Sample, input_path: Path):
    with open(input_path.with_suffix(".svg"), "w") as f:
        f.write(sample.svg)
    with open(input_path.parent / (input_path.stem + "_colored.svg"), "w") as f:
        f.write(sample.svg_colored)


def output_png(sample: Sample, input_path: Path):
    image = renderer.svg_to_pil(sample.svg)
    with open(input_path.with_suffix(".png"), "wb") as f:
        image.save(f, format="PNG")

    image_colored = renderer.svg_to_pil(sample.svg_colored)
    with open(input_path.parent / (input_path.stem + "_colored.png"), "wb") as f:
        image_colored.save(f, format="PNG")


if __name__ == "__main__":
    main()
