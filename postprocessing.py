from xml.etree import ElementTree
import re
from rdkit import Chem


def color_instances(svg: str):
    """
    Color instances in the SVG output. Each unique instance-id is assigned a unique color.
    """
    matches = re.finditer(r'instance-id="(\d+)"', svg)
    instance_ids = {int(m.group(1)) for m in matches}
    svg_root = ElementTree.fromstring(svg)
    for color, instance_id in zip(unique_colors(), sorted(instance_ids)):
        # find all elements with the current instance-id
        instance_elements = []
        for element in list(svg_root.iter()):
            if element.get("instance-id") == str(instance_id):
                instance_elements.append(element)

        # set the color of these elements
        for element in instance_elements:
            # adjust style
            style = element.get("style", "")
            style_dict = dict(s.split(":") for s in style.split(";") if s)
            if "stroke" in style_dict:
                style_dict["stroke"] = color
            style = ";".join(f"{k}:{v}" for k, v in style_dict.items())
            element.set("style", style)

            # adjust fill
            if "fill" in element.attrib:
                element.set("fill", color)

    return ElementTree.tostring(svg_root, encoding="unicode")


def post_process_svg(svg: str, mol: Chem.Mol) -> str:
    """
    Post-process the SVG output by adding instance IDs and classes.
    Each atom and bond is assigned a unique instance ID and a class based on
    its type (element symbol for atoms, bond type for bonds).

    Non-atom/bond paths are removed from the SVG.
    """
    bond_types = [bond.GetBondType() for bond in mol.GetBonds()]
    atom_symbols = [atom.GetSymbol() for atom in mol.GetAtoms()]
    instances = []

    svg_root = ElementTree.fromstring(svg)
    for element in list(svg_root.iter()):
        # only process path elements
        if element.tag != "{http://www.w3.org/2000/svg}path":
            continue

        bond_index = get_bond_index(element)
        if bond_index is not None:
            instance_name = f"bond-{bond_index}"
            if instance_name not in instances:
                instances.append(instance_name)
            instance_id = instances.index(instance_name)

            attrib = {
                "instance-class": str(bond_types[bond_index]),
                "instance-id": str(instance_id),
            }
            attrib.update(element.attrib)
            element.attrib = attrib
            continue

        atom_index = get_atom_index(element)
        if atom_index is not None:
            instance_name = f"atom-{atom_index}"
            if instance_name not in instances:
                instances.append(instance_name)
            instance_id = instances.index(instance_name)

            attrib = {
                "instance-class": str(atom_symbols[atom_index]),
                "instance-id": str(instance_id),
            }
            attrib.update(element.attrib)
            element.attrib = attrib
            continue

        # remove non-atom/bond paths
        svg_root.remove(element)

    # remove class attributes
    for element in list(svg_root.iter()):
        if "class" in element.attrib:
            del element.attrib["class"]

    return ElementTree.tostring(svg_root, encoding="unicode")


def get_atom_index(element: ElementTree.Element):
    """
    Extract the atom index from the class attribute of an SVG element.
    """
    class_attr = element.get("class", "")
    match = re.search(r"atom-(\d+)", class_attr)
    if match:
        return int(match.group(1))
    return None


def get_bond_index(element: ElementTree.Element):
    """
    Extract the bond index from the class attribute of an SVG element.
    """
    class_attr = element.get("class", "")
    match = re.search(r"bond-(\d+)", class_attr)
    if match:
        return int(match.group(1))
    return None


def unique_colors():
    def hsv_to_rgb(h, s, v) -> tuple[int, int, int]:
        i = int(h * 6)
        f = h * 6 - i
        p = v * (1 - s)
        q = v * (1 - f * s)
        t = v * (1 - (1 - f) * s)

        i %= 6
        if i == 0:
            return (v, t, p)
        if i == 1:
            return (q, v, p)
        if i == 2:
            return (p, v, t)
        if i == 3:
            return (p, q, v)
        if i == 4:
            return (t, p, v)
        else:
            return (v, p, q)

    def generate_colors():
        hue = 0
        golden_ratio = 0.618033988749895
        while True:
            hue = (hue + golden_ratio) % 1
            rgb = hsv_to_rgb(hue, 0.8, 0.95)
            yield f"#{int(rgb[0] * 255):02x}{int(rgb[1] * 255):02x}{int(rgb[2] * 255):02x}"

    return generate_colors()
