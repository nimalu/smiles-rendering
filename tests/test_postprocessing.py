from smiles_segmentation.postprocessing import _unique_colors


def test_unique_colors_are_different():
    """Test that unique_colors generates different colors."""
    color_generator = _unique_colors()

    # Generate a reasonable number of colors to test
    num_colors_to_test = 2500
    generated_colors = []

    for _ in range(num_colors_to_test):
        generated_colors.append(next(color_generator))

    # Verify all colors are unique
    unique_color_set = set(generated_colors)
    assert len(unique_color_set) == len(generated_colors), (
        f"Expected {len(generated_colors)} unique colors, but got {len(unique_color_set)}. "
        f"Duplicates found: {len(generated_colors) - len(unique_color_set)}"
    )
