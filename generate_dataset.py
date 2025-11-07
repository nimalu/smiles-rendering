"""
Generate a YOLO format instance segmentation dataset from SMILES strings.

This script creates a dataset with molecular structure images and corresponding
YOLO format annotations for instance segmentation training.
"""

import random
from pathlib import Path
from typing import List

from smiles_segmentation import Options
from smiles_segmentation.images import svg_to_pil
from smiles_segmentation.postprocessing import (
    annotate_svg_with_instances,
    flatten_paths_to_polygons,
)
from smiles_segmentation.renderer import create_mol, create_svg
from smiles_segmentation.yolo import svg_to_yolo_format, build_class_to_id


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


# Common molecules for dataset diversity
COMMON_MOLECULES = [
    # Simple organic molecules
    ("CC(=O)O", "acetic_acid"),
    ("CCO", "ethanol"),
    ("CC(C)O", "isopropanol"),
    ("C1=CC=CC=C1", "benzene"),
    ("CC(=O)C", "acetone"),
    # Amino acids
    ("NC(C)C(=O)O", "alanine"),
    ("NC(CC(=O)O)C(=O)O", "aspartic_acid"),
    ("NC(CCC(=O)O)C(=O)O", "glutamic_acid"),
    ("NC(CC1=CC=CC=C1)C(=O)O", "phenylalanine"),
    ("NC(CC(C)C)C(=O)O", "leucine"),
    # Sugars
    ("OCC1OC(O)C(O)C(O)C1O", "glucose"),
    ("OCC1OC(O)C(O)C(O)C1O", "fructose"),
    # Nucleotides
    ("C1=NC2=C(N1)C(=O)NC(=N2)N", "guanine"),
    ("CC1=CN(C(=O)NC1=O)C2CC(C(O2)CO)O", "thymidine"),
    ("C1=CN(C=O)C(=O)NC1=N", "cytosine"),
    # Pharmaceuticals
    ("CC(=O)Oc1ccccc1C(=O)O", "aspirin"),
    ("CC(C)Cc1ccc(cc1)C(C)C(=O)O", "ibuprofen"),
    ("CN1C=NC2=C1C(=O)N(C(=O)N2C)C", "caffeine"),
    # Aromatic compounds
    ("c1ccc2c(c1)ccc3c2cccc3", "anthracene"),
    ("c1ccc2c(c1)ccc1c2cccc1", "phenanthrene"),
    ("CC1=CC=C(C=C1)C", "para-xylene"),
    # Alcohols and ethers
    ("CCCCO", "butanol"),
    ("CCOC(=O)C", "ethyl_acetate"),
    ("CCOCC", "diethyl_ether"),
    # Aldehydes and ketones
    ("CC=O", "acetaldehyde"),
    ("O=Cc1ccccc1", "benzaldehyde"),
    ("CC(=O)CC(C)C", "methyl_isobutyl_ketone"),
    # Amines
    ("CCN", "ethylamine"),
    ("c1ccc(cc1)N", "aniline"),
    ("CN(C)C", "trimethylamine"),
    # Carboxylic acids
    ("CCCC(=O)O", "butyric_acid"),
    ("c1ccc(cc1)C(=O)O", "benzoic_acid"),
    # Heterocycles
    ("c1cccnc1", "pyridine"),
    ("c1cnc[nH]1", "imidazole"),
    ("c1ccc2[nH]ccc2c1", "indole"),
    # Additional diversity
    ("CC(C)(C)O", "tert-butanol"),
    ("c1ccc(cc1)O", "phenol"),
    ("CC(C)C(=O)C(C)C", "diisobutyl_ketone"),
    ("CCCCCC", "hexane"),
    ("C1CCCCC1", "cyclohexane"),
    ("CC1=CC=C(C=C1)O", "para-cresol"),
]


def generate_sample(
    smiles: str,
    output_dir: Path,
    sample_id: int,
    rotation_range: tuple[int, int] = (0, 360),
    width_range: tuple[int, int] = (512, 768),
    height_range: tuple[int, int] = (512, 768),
    base_font_size_range: tuple[float, float] = (0.5, 0.7),
    atom_label_padding_range: tuple[float, float] = (0.05, 0.15),
    bond_line_width_range: tuple[float, float] = (1.5, 2.5),
    multiple_bond_offset_range: tuple[float, float] = (0.1, 0.2),
) -> dict[str, int]:
    """
    Generate a single sample with randomized rendering options.

    Args:
        smiles: SMILES string of the molecule
        output_dir: Directory to save the output files
        sample_id: Unique ID for this sample
        rotation_range: Min and max rotation angles in degrees
        width_range: Min and max image width in pixels
        height_range: Min and max image height in pixels
        base_font_size_range: Min and max base font size
        atom_label_padding_range: Min and max additional atom label padding
        bond_line_width_range: Min and max bond line width
        multiple_bond_offset_range: Min and max multiple bond offset

    Returns:
        dict mapping class names to IDs for this sample
    """
    # Randomize all parameters
    rotation = random.randint(*rotation_range)
    width = random.randint(*width_range)
    height = random.randint(*height_range)
    base_font_size = random.uniform(*base_font_size_range)
    atom_label_padding = random.uniform(*atom_label_padding_range)
    bond_line_width = random.uniform(*bond_line_width_range)
    multiple_bond_offset = random.uniform(*multiple_bond_offset_range)

    # Create options with randomized parameters
    options = Options(
        width=width,
        height=height,
        rotate=rotation,
        baseFontSize=base_font_size,
        additionalAtomLabelPadding=atom_label_padding,
        bondLineWidth=bond_line_width,
        multipleBondOffset=multiple_bond_offset,
        backgroundColour=(1, 1, 1),  # White background
    )

    # Generate molecular structure
    mol = create_mol(smiles)
    svg = create_svg(mol, options)
    svg = annotate_svg_with_instances(svg, mol, options)
    svg = flatten_paths_to_polygons(svg)

    # Build class mapping
    class_to_id = build_class_to_id(svg)

    # Generate YOLO format annotation
    yolo_str = svg_to_yolo_format(svg, class_to_id)

    # Save annotation file
    annotation_path = output_dir / f"{sample_id:06d}.txt"
    with open(annotation_path, "w") as f:
        f.write(yolo_str)

    # Save black and white PNG
    image = svg_to_pil(svg)
    image_path = output_dir / f"{sample_id:06d}.png"
    with open(image_path, "wb") as f:
        image.save(f, format="PNG")

    return class_to_id


def merge_class_mappings(class_mappings: List[dict[str, int]]) -> dict[str, int]:
    """
    Merge multiple class mappings into a single consistent mapping.

    Args:
        class_mappings: List of class-to-id dictionaries

    Returns:
        Unified class-to-id dictionary
    """
    all_classes = set()
    for mapping in class_mappings:
        all_classes.update(mapping.keys())

    # Sort for consistency
    return {cls: idx for idx, cls in enumerate(sorted(all_classes))}


def generate_dataset(
    output_base: Path,
    num_train: int = 100,
    num_val: int = 20,
    rotation_range: tuple[int, int] = (0, 360),
    width_range: tuple[int, int] = (512, 768),
    height_range: tuple[int, int] = (512, 768),
    base_font_size_range: tuple[float, float] = (0.5, 0.7),
    atom_label_padding_range: tuple[float, float] = (0.05, 0.15),
    bond_line_width_range: tuple[float, float] = (1.5, 2.5),
    multiple_bond_offset_range: tuple[float, float] = (0.1, 0.2),
    seed: int | None = None,
):
    """
    Generate a complete YOLO format dataset.

    Args:
        output_base: Base directory for the dataset
        num_train: Number of training samples
        num_val: Number of validation samples
        rotation_range: Min and max rotation angles in degrees
        width_range: Min and max image width in pixels
        height_range: Min and max image height in pixels
        base_font_size_range: Min and max base font size
        atom_label_padding_range: Min and max additional atom label padding
        bond_line_width_range: Min and max bond line width
        multiple_bond_offset_range: Min and max multiple bond offset
        seed: Random seed for reproducibility
    """
    if seed is not None:
        random.seed(seed)

    # Create directory structure
    train_images_dir = output_base / "images" / "train"
    train_labels_dir = output_base / "labels" / "train"
    val_images_dir = output_base / "images" / "val"
    val_labels_dir = output_base / "labels" / "val"

    for dir_path in [train_images_dir, train_labels_dir, val_images_dir, val_labels_dir]:
        dir_path.mkdir(parents=True, exist_ok=True)

    print(f"Generating dataset in {output_base}")
    print(f"Train samples: {num_train}, Validation samples: {num_val}")
    print(f"Rotation range: {rotation_range[0]}° to {rotation_range[1]}°")
    print(
        f"Image size range: {width_range[0]}-{width_range[1]}x{height_range[0]}-{height_range[1]}"
    )
    print(f"Font size range: {base_font_size_range[0]:.2f}-{base_font_size_range[1]:.2f}")
    print(f"Bond line width range: {bond_line_width_range[0]:.2f}-{bond_line_width_range[1]:.2f}")

    # Generate training samples
    print("\nGenerating training samples...")
    train_class_mappings = []
    for i in range(num_train):
        smiles, name = random.choice(COMMON_MOLECULES)
        try:
            # Save to images directory, annotations go to labels directory
            image_output_dir = train_images_dir
            label_output_dir = train_labels_dir

            # Generate image
            class_mapping = generate_sample(
                smiles=smiles,
                output_dir=image_output_dir,
                sample_id=i,
                rotation_range=rotation_range,
                width_range=width_range,
                height_range=height_range,
                base_font_size_range=base_font_size_range,
                atom_label_padding_range=atom_label_padding_range,
                bond_line_width_range=bond_line_width_range,
                multiple_bond_offset_range=multiple_bond_offset_range,
            )

            # Move annotation to labels directory
            (image_output_dir / f"{i:06d}.txt").rename(label_output_dir / f"{i:06d}.txt")

            train_class_mappings.append(class_mapping)

            if (i + 1) % 10 == 0:
                print(f"  Generated {i + 1}/{num_train} samples")
        except Exception as e:
            print(f"  Error generating sample {i} ({name}): {e}")

    # Generate validation samples
    print("\nGenerating validation samples...")
    val_class_mappings = []
    for i in range(num_val):
        smiles, name = random.choice(COMMON_MOLECULES)
        try:
            # Save to images directory, annotations go to labels directory
            image_output_dir = val_images_dir
            label_output_dir = val_labels_dir

            # Generate image
            class_mapping = generate_sample(
                smiles=smiles,
                output_dir=image_output_dir,
                sample_id=i,
                rotation_range=rotation_range,
                width_range=width_range,
                height_range=height_range,
                base_font_size_range=base_font_size_range,
                atom_label_padding_range=atom_label_padding_range,
                bond_line_width_range=bond_line_width_range,
                multiple_bond_offset_range=multiple_bond_offset_range,
            )

            # Move annotation to labels directory
            (image_output_dir / f"{i:06d}.txt").rename(label_output_dir / f"{i:06d}.txt")

            val_class_mappings.append(class_mapping)

            if (i + 1) % 10 == 0:
                print(f"  Generated {i + 1}/{num_val} samples")
        except Exception as e:
            print(f"  Error generating sample {i} ({name}): {e}")

    # Merge all class mappings
    all_class_mappings = train_class_mappings + val_class_mappings
    unified_class_mapping = merge_class_mappings(all_class_mappings)

    # Create dataset YAML file
    yaml_content = f"""# SMILES Molecular Instance Segmentation Dataset
# Generated for use with Ultralytics YOLO
# Documentation: https://docs.ultralytics.com/datasets/segment/

# Dataset root directory
path: {output_base.name}

# Train/val/test sets
train: images/train  # training images
val: images/val      # validation images
test:                # test images (optional)

# Classes
names:
"""

    for class_name, class_id in sorted(unified_class_mapping.items(), key=lambda x: x[1]):
        yaml_content += f"  {class_id}: {class_name}\n"

    yaml_path = output_base / "dataset.yaml"
    with open(yaml_path, "w") as f:
        f.write(yaml_content)

    print("\nDataset generation complete!")
    print(f"Total classes: {len(unified_class_mapping)}")
    print(f"Dataset YAML saved to: {yaml_path}")
    print("\nClass mapping:")
    for class_name, class_id in sorted(unified_class_mapping.items(), key=lambda x: x[1]):
        print(f"  {class_id}: {class_name}")


def main():
    generate_dataset(
        output_base=OUTPUT_DIR,
        num_train=NUM_TRAIN_SAMPLES,
        num_val=NUM_VAL_SAMPLES,
        rotation_range=ROTATION_RANGE,
        width_range=WIDTH_RANGE,
        height_range=HEIGHT_RANGE,
        base_font_size_range=BASE_FONT_SIZE_RANGE,
        atom_label_padding_range=ATOM_LABEL_PADDING_RANGE,
        bond_line_width_range=BOND_LINE_WIDTH_RANGE,
        multiple_bond_offset_range=MULTIPLE_BOND_OFFSET_RANGE,
        seed=RANDOM_SEED,
    )


if __name__ == "__main__":
    main()
