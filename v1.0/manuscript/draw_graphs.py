"""Draw the two labelled graphs in Figure 1. Requires ReportLab."""

from itertools import combinations, permutations
from math import cos, pi, sin
from pathlib import Path

import reportlab
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.pdfgen import canvas


def draw_graph(pdf, points, labels, edges):
    pdf.setLineWidth(0.8)
    pdf.setStrokeColorRGB(0.20, 0.21, 0.22)
    for i, j in edges:
        pdf.line(*points[i], *points[j])
    pdf.setFillColorRGB(1, 1, 1)
    for x, y in points:
        pdf.circle(x, y, 6.8, stroke=1, fill=1)
    pdf.setFillColorRGB(0.07, 0.08, 0.09)
    pdf.setFont("FigureLabels", 7.5)
    for (x, y), label in zip(points, labels):
        pdf.drawCentredString(x, y - 2.6, label)


def main():
    path = Path(__file__).with_name("graphs.pdf")
    pdfmetrics.registerFont(
        TTFont("FigureLabels", str(Path(reportlab.__file__).parent / "fonts/Vera.ttf"))
    )
    pdf = canvas.Canvas(str(path), pagesize=(344, 114), invariant=True)
    pdf.setTitle("The graph F26 and the Petersen graph")
    pdf.setAuthor("Oleg Kyzlasov")
    tree = [(80, 60), (50, 42), (80, 94), (20, 24), (140, 24), (110, 42)]
    draw_graph(pdf, tree, list("012345"), [(0, 1), (0, 2), (0, 5), (1, 3), (4, 5)])
    pairs = list(combinations(range(5), 2))

    def adjacent(i, j):
        return set(pairs[i]).isdisjoint(pairs[j])

    outer = next(
        p
        for p in permutations(range(10), 5)
        if all(adjacent(p[i], p[(i + 1) % 5]) for i in range(5))
    )
    inner = tuple(
        next(v for v in range(10) if v not in outer and adjacent(v, i)) for i in outer
    )
    order = outer + inner
    assert len(set(order)) == 10
    edges = [
        (i, j) for i, j in combinations(range(10), 2) if adjacent(order[i], order[j])
    ]
    assert len(edges) == 15
    points = [
        (
            260 + radius * cos(pi / 2 - 2 * pi * i / 5),
            56 + radius * sin(pi / 2 - 2 * pi * i / 5),
        )
        for radius in (48, 22.5)
        for i in range(5)
    ]
    labels = ["".join(map(str, pairs[i])) for i in order]
    draw_graph(pdf, points, labels, edges)
    pdf.showPage()
    pdf.save()


if __name__ == "__main__":
    main()
