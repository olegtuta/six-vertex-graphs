"""Regenerate Figure 2 and check all adjacencies in the displayed copy."""

from itertools import combinations
from pathlib import Path

import reportlab
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.pdfgen.canvas import Canvas

from draw_graphs import draw_graph


def main():
    edges = {(0, 1), (0, 3), (0, 4), (1, 2), (1, 3), (3, 5), (4, 5)}
    cells = [(0, 0), (1, 1), (0, 2), (2, 2), (1, 2), (0, 1)]
    assert len(set(cells)) == 6
    for i, j in combinations(range(6), 2):
        a, b = cells[i]
        c, d = cells[j]
        assert ((i, j) in edges) == (a != c and b != d)
    pdfmetrics.registerFont(TTFont("FigureLabels", str(Path(reportlab.__file__).parent / "fonts/Vera.ttf")))
    pdf = Canvas(str(Path(__file__).with_name("rook-example.pdf")), pagesize=(254, 108), invariant=True)
    pdf.setTitle("An induced copy of F74 in the complement of the rook graph")
    pdf.setAuthor("Oleg Kyzlasov")
    draw_graph(pdf, [(17, 47), (42, 70), (42, 96), (67, 47), (17, 13), (67, 13)],
               list("012345"), sorted(edges))
    side = 30
    left, bottom = 148, 9
    pdf.setLineWidth(.7)
    for k in range(4):
        pdf.line(left + k*side, bottom, left + k*side, bottom + 3*side)
        pdf.line(left, bottom + k*side, left + 3*side, bottom + k*side)
    for i, (row, col) in enumerate(cells):
        x, y = left + (col + .5)*side, bottom + (2.5 - row)*side
        pdf.circle(x, y, 6.8, stroke=1, fill=0)
        pdf.drawCentredString(x, y - 2.6, str(i))
    pdf.showPage()
    pdf.save()


if __name__ == "__main__":
    main()
