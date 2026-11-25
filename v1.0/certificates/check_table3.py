"""Check row 26 against a locally saved copy of arXiv:2606.00290v3 HTML."""

import argparse
import json
import math
import re
from html.parser import HTMLParser
from pathlib import Path


class TableRow(HTMLParser):
    def __init__(self):
        super().__init__()
        self.active = False
        self.cells = []
        self.paths = []

    def handle_starttag(self, tag, attrs):
        fields = dict(attrs)
        if tag == "tr":
            self.active = fields.get("id") == "S5.T3.5.7"
        if self.active and tag == "td":
            self.cells.append("")
        if self.active and len(self.cells) == 2 and tag == "path":
            self.paths.append(fields["d"])

    def handle_endtag(self, tag):
        if tag == "tr":
            self.active = False

    def handle_data(self, data):
        if self.active and self.cells:
            self.cells[-1] += data


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("html", type=Path)
    row = TableRow()
    row.feed(parser.parse_args().html.read_text(encoding="utf-8"))
    cells = [value.strip() for value in row.cells]
    if len(cells) != 7 or [cells[i] for i in (0, 4, 5, 6)] != [
        "26",
        "0.0293174175",
        "0.0432004682",
        "8",
    ]:
        raise ValueError(
            "Row number, bounds or sample order differ from Table 3, version 3."
        )
    # In the arXiv SVG each circular vertex ends with a move to its centre.
    # Each straight edge is shortened slightly to meet the vertex boundary.
    dots, lines = [], []
    number = r"(-?\d+(?:\.\d+)?)"
    for path in row.paths:
        dot = re.search(r"\bZ M " + number + " " + number + r"$", path)
        edge = re.fullmatch(
            "M " + number + " " + number + " L " + number + " " + number, path
        )
        if dot:
            dots.append(tuple(map(float, dot.groups())))
        elif edge:
            x1, y1, x2, y2 = map(float, edge.groups())
            lines.append(((x1, y1), (x2, y2)))
        else:
            raise ValueError("Unrecognized drawing command.")
    if len(dots) != 6 or len(lines) != 5:
        raise ValueError("Unexpected graph order or size.")
    cx, cy = (sum(p[k] for p in dots) / 6 for k in (0, 1))
    dots.sort(key=lambda p: math.atan2(p[0] - cx, p[1] - cy) % math.tau)
    # A tiny rounding error in cx can move the top vertex across the angle cut.
    top = max(range(6), key=lambda i: dots[i][1])
    dots = dots[top:] + dots[:top]

    def vertex(point):
        matches = [i for i, p in enumerate(dots) if math.dist(point, p) < 3]
        if len(matches) != 1:
            raise ValueError("Ambiguous line endpoint.")
        return matches[0]

    edges = sorted(tuple(sorted(map(vertex, edge))) for edge in lines)
    if edges != [(0, 1), (0, 2), (0, 5), (1, 3), (4, 5)]:
        raise ValueError("Different graph in the F_i column.")
    print(
        json.dumps(
            {
                "version": "2606.00290v3",
                "table": 3,
                "row": 26,
                "edges": edges,
                "displayed_lower": cells[4],
                "numerical_upper": cells[5],
                "N": int(cells[6]),
            }
        )
    )


if __name__ == "__main__":
    main()
