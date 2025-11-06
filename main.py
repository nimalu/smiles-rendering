import argparse
import json
from pathlib import Path
from smiles_segmentation import Options
from smiles_segmentation.images import svg_to_pil
from smiles_segmentation.postprocessing import (
    annotate_svg_with_instances,
    color_instances,
    flatten_paths_to_polygons,
)
from smiles_segmentation.renderer import create_mol, create_svg
from smiles_segmentation.yolo import svg_to_yolo_format, build_class_to_id


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input", help="File containing a SMILES string")
    args = parser.parse_args()

    input_path = Path(args.input)
    smiles = input_path.read_text().strip()

    options = Options()
    mol = create_mol(smiles)
    svg = create_svg(mol, options)
    svg = annotate_svg_with_instances(svg, mol, options)
    svg = flatten_paths_to_polygons(svg)
    svg_colored = color_instances(svg)
    class_to_id = build_class_to_id(svg)
    yolo_str = svg_to_yolo_format(svg, class_to_id)

    with open(input_path.with_suffix(".txt"), "w") as f:
        f.write(yolo_str)
    with open(input_path.parent / (input_path.stem + "_classes.json"), "w") as f:
        json.dump(class_to_id, f)

    with open(input_path.with_suffix(".svg"), "w") as f:
        f.write(svg)
    with open(input_path.parent / (input_path.stem + "_colored.svg"), "w") as f:
        f.write(svg_colored)
    image = svg_to_pil(svg)
    with open(input_path.with_suffix(".png"), "wb") as f:
        image.save(f, format="PNG")
    image_colored = svg_to_pil(svg_colored)
    with open(input_path.parent / (input_path.stem + "_colored.png"), "wb") as f:
        image_colored.save(f, format="PNG")


if __name__ == "__main__":
    main()
