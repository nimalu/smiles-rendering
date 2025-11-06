from xml.etree import ElementTree


def build_class_to_id(svg: str):
    svg_root = ElementTree.fromstring(svg)

    classes = set()
    for element in svg_root.iter():
        if element.get("instance-class") is not None:
            classes.add(element.get("instance-class"))

    return {cls: idx for idx, cls in enumerate(sorted(classes))}


def svg_to_yolo_format(svg: str, class_to_id: dict[str, int] | None = None) -> str:
    # parse svg
    svg_root = ElementTree.fromstring(svg)

    if class_to_id is None:
        class_to_id = build_class_to_id(svg)

    yolo_lines = []
    for element in svg_root.iter():
        if element.tag != "{http://www.w3.org/2000/svg}polygon":
            continue
        instance_class = element.get("instance-class")
        if instance_class is None:
            raise ValueError("SVG element missing 'instance-class' attribute")

        class_id = class_to_id[instance_class]

        points = element.get("points", "")
        point_list = []
        for point_str in points.strip().split(" "):
            x_str, y_str = point_str.split(",")
            point_list.append((float(x_str), float(y_str)))

        yolo_line = f"{class_id} " + " ".join(f"{x} {y}" for x, y in point_list)
        yolo_lines.append(yolo_line)

    return "\n".join(yolo_lines)
