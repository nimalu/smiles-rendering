import argparse
from pathlib import Path
import renderer
from rdkit import Chem
import postprocessing


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("smiles", help="SMILES string to process")
    parser.add_argument("destination", help="Destination file path")
    args = parser.parse_args()
    create_sample(args.smiles, Path(args.destination).with_suffix(""))


def create_sample(smiles: str, destination: Path):
    mol = Chem.MolFromSmiles(smiles)

    options = renderer.SVGOptions()
    svg = renderer.create_svg(mol, options)

    svg = postprocessing.post_process_svg(svg, mol)
    with open(destination.with_suffix(".svg"), "w") as f:
        f.write(svg)

    svg = postprocessing.color_instances(svg)
    with open(destination.parent / (destination.stem + "-colored.svg"), "w") as f:
        f.write(svg)


if __name__ == "__main__":
    main()
