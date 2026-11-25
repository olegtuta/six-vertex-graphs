Petersen substitutions and the inducibility of the six-vertex tree F26
Oleg Kyzlasov | Manuscript source | v1.0, 18 November 2026

f26.tex is the complete LaTeX source. graphs.pdf is the vector artwork
for Figure 1. From the extracted directory, build with either:

  tectonic f26.tex

or a standard LaTeX installation:

  pdflatex f26.tex
  pdflatex f26.tex

The release PDF was built with Tectonic 0.17.0. The figure is supplied,
so the build requires no drawing software. To regenerate it, run
python draw_graphs.py in an environment with ReportLab and its bundled
Bitstream Vera fonts. The drawing script constructs the Petersen edges
from disjoint pairs. The font notice is included in this directory.

For the release build, set the environment variable SOURCE_DATE_EPOCH to
1794960000 and run: tectonic -Z deterministic-mode f26.tex

The computational data and proof checkers are in the separate archive
f26-certificates.zip deposited with the preprint.
