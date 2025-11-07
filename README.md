
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

## Dataset Generation for YOLO Instance Segmentation

The `generate_dataset.py` script creates a complete dataset ready for training YOLO instance segmentation models.

### Overview

The dataset generation script:

- Generates molecular structure images from SMILES strings
- Creates corresponding YOLO format annotations
- Randomizes rotation for sample diversity
- Organizes files in the standard YOLO dataset structure
- Generates a dataset YAML configuration file

### Dataset Structure

The generated dataset follows the Ultralytics YOLO format:

```
dataset/
├── dataset.yaml          # Dataset configuration
├── images/
│   ├── train/           # Training images (PNG)
│   └── val/             # Validation images (PNG)
└── labels/
    ├── train/           # Training annotations (TXT)
    └── val/             # Validation annotations (TXT)
```

### Usage

This project uses [uv](https://docs.astral.sh/uv/) for Python package management. Make sure you have `uv` installed before running the script.

#### Configuration

The script uses constants defined at the top of `generate_dataset.py`. Edit these values to customize the dataset generation:

```python
# Dataset configuration constants
OUTPUT_DIR = Path("./dataset")
NUM_TRAIN_SAMPLES = 100
NUM_VAL_SAMPLES = 20

# Randomization ranges
ROTATION_RANGE = (0, 360)
WIDTH_RANGE = (512, 768)
HEIGHT_RANGE = (512, 768)
BASE_FONT_SIZE_RANGE = (0.5, 0.7)
ATOM_LABEL_PADDING_RANGE = (0.05, 0.15)
BOND_LINE_WIDTH_RANGE = (1.5, 2.5)
MULTIPLE_BOND_OFFSET_RANGE = (0.1, 0.2)

# Random seed for reproducibility (None for random)
RANDOM_SEED = None
```

#### Running the Script

Generate the dataset with the configured settings:

```bash
uv run generate_dataset.py
```

### Dataset Content

The dataset includes ~40 common molecules:

- **Simple organic molecules**: acetic acid, ethanol, benzene, acetone, etc.
- **Amino acids**: alanine, phenylalanine, leucine, etc.
- **Sugars**: glucose, fructose
- **Nucleotides**: guanine, thymidine, cytosine
- **Pharmaceuticals**: aspirin, ibuprofen, caffeine
- **Aromatic compounds**: anthracene, phenanthrene
- **Heterocycles**: pyridine, imidazole, indole
- And more...

Each sample randomly selects a molecule and applies a random rotation within the specified range.

### YOLO Format

Each annotation file contains one line per instance:

```
<class_id> <x1> <y1> <x2> <y2> ... <xn> <yn>
```

Where:
- `class_id`: Integer ID of the class (atom type or bond type)
- `x1 y1 ... xn yn`: Polygon vertices in absolute pixel coordinates

### Training with Ultralytics

Once the dataset is generated, you can train a YOLO segmentation model. First, install Ultralytics:

```bash
uv pip install ultralytics
```

Then train the model:

```python
from ultralytics import YOLO

# Load a pretrained model
model = YOLO('yolov8n-seg.pt')

# Train the model
results = model.train(
    data='dataset/dataset.yaml',
    epochs=100,
    imgsz=640,
    batch=16
)
```

### Notes

- The dataset directory is excluded from git (see `.gitignore`)
- Images are generated in black and white for clear atom/bond distinction
- Rotation randomization ensures model learns rotation-invariant features
- Class IDs are consistent across train/val splits


