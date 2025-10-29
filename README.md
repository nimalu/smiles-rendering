
# Render SMILES

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

## Command Line Usage

The package provides a command-line tool `render-smiles` for processing SMILES files:

```bash
# Render to SVG (default)
uv run render-smiles smiles.smi

# Render to PNG
uv run render-smiles --png smiles.smi

# This will generate:
# - smiles.svg (standard molecular diagram)
# - smiles_colored.svg (colored by molecular fragments)
# - smiles.png and smiles_colored.png (if --png is used)
```

### Python API

```python
from render_smiles import Options, create_sample

# Basic usage
smiles = "CCc(c1)ccc2[n+]1ccc3c2[nH]c4c3cccc4"
sample = create_sample(smiles)

# Access the generated content
print(sample.smiles)      # Original SMILES string
print(sample.mol)         # RDKit molecule object
print(sample.svg)         # Standard SVG representation
print(sample.svg_colored) # Colored SVG representation

# Customized rendering
options = Options(
    width=800,
    height=600,
    backgroundColour=(0.95, 0.95, 0.95),  # Light gray background
    bondLineWidth=3,
    baseFontSize=0.8
)

sample = create_sample(smiles, options)
```


