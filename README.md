
# SMILES Segmentation

A Python package for rendering SMILES (Simplified Molecular Input Line Entry System) strings into visual molecular diagrams with fragment coloring.

## Features

- **SMILES to Molecule Conversion**: Convert SMILES strings to RDKit molecule objects
- **SVG Generation**: Create scalable vector graphics of molecular structures
- **PNG Export**: Convert SVG diagrams to raster images
- **Fragment Coloring**: Automatically assign unique colors to different molecular fragments
- **Customizable Rendering**: Extensive options for controlling diagram appearance
- **Command-line Interface**: Easy batch processing of SMILES files


## Development Setup

```bash
uv sync
```

## Python API

The package provides a modular API for converting SMILES strings to visual molecular diagrams:

```python
from smiles_segmentation import Options
from smiles_segmentation.renderer import create_mol, create_svg
from smiles_segmentation.postprocessing import (
    annotate_svg_with_instances,
    color_instances,
    flatten_paths_to_polygons,
)
from smiles_segmentation.images import svg_to_pil
from smiles_segmentation.yolo import svg_to_yolo_format, build_class_to_id

# Basic usage
smiles = "CCc(c1)ccc2[n+]1ccc3c2[nH]c4c3cccc4"

# Create molecule and generate SVG
options = Options()
mol = create_mol(smiles)
svg = create_svg(mol, options)

# Add instance annotations (atom and bond identification)
svg = annotate_svg_with_instances(svg, mol, options)

# Convert paths to polygons for segmentation
svg = flatten_paths_to_polygons(svg)

# Generate colored version
svg_colored = color_instances(svg)

# Convert to PIL Image for PNG export
image = svg_to_pil(svg)
image.save("output.png", format="PNG")

# Generate YOLO format annotations
class_to_id = build_class_to_id(svg)
yolo_annotations = svg_to_yolo_format(svg, class_to_id)

# Customized rendering options
options = Options(
    width=800,
    height=600,
    backgroundColour=(0.95, 0.95, 0.95),  # Light gray background
    bondLineWidth=3,
    baseFontSize=0.8
)
```

## Example Usage

See `main.py` for a complete example that demonstrates the full pipeline:
- Reading SMILES from a file
- Generating SVG and PNG outputs
- Creating both standard and colored versions
- Exporting YOLO format annotations

Run the example:
```bash
python main.py smiles.smi
```

This will generate:
- `smiles.svg` - Standard molecular diagram
- `smiles_colored.svg` - Colored by molecular fragments
- `smiles.png` and `smiles_colored.png` - PNG versions
- `smiles.txt` - YOLO format segmentation annotations
- `smiles_classes.json` - Class-to-ID mapping for annotations


