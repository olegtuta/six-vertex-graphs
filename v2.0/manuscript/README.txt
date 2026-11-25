Substitution bounds for the inducibility of six-vertex graphs
Oleg Kyzlasov | Version 2.0 | 22 November 2026

six-vertex.tex is the main manuscript and substitutions.tex contains
the sections on F74, F29, F22, F34 and F62. The vector figures are
graphs.pdf (Figure 1) and rook-example.pdf (Figure 2). Build with

  tectonic six-vertex.tex

or run pdflatex six-vertex.tex twice. A first installation of the TeX
toolchain may need to download standard packages. All figure assets
are present. To reproduce the release build with Tectonic 0.17.0,
set SOURCE_DATE_EPOCH to 1795305600, then run

  tectonic -Z deterministic-mode six-vertex.tex

The optional scripts draw_graphs.py and draw_rook.py regenerate the two
figures with ReportLab and its bundled Bitstream Vera fonts. The latter
script also checks all fifteen pairs of vertices in the displayed copy.
The font notice is included. Proof data and checkers are in
six-vertex-certificates.zip.
