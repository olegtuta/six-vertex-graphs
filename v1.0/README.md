# Version 1.0 · Petersen substitutions and the inducibility of the six-vertex tree F26

Oleg Kyzlasov, 18 November 2026. [Preprint on Zenodo](https://zenodo.org/records/24341542) · [DOI 10.5281/zenodo.24341542](https://doi.org/10.5281/zenodo.24341542)

Balanced substitutions of the Petersen graph give limiting induced density
480/11111 for the six-vertex tree F26. This value is sharp among graphs
with no induced bull or house. A rational certificate gives a separate
upper bound for unrestricted graphs.

- [Manuscript source and figures](manuscript/)
- [Verification code and certificate data](certificates/)
- [Certificate documentation](certificates/README.txt)

## Verify the certificates

From the repository root:

~~~sh
cd v1.0/certificates
python -B verify.py --part elementary
python -B verify.py --part restricted
python -B verify.py
~~~

The first command after `cd` checks the Petersen copy count and finite
examples. The restricted replay also checks the coloured census, the rational
inequality and its equality calculation. Both use only Python 3.10 or later.
The full replay also checks the unrestricted upper bound. It needs a
little-endian system and a C++20 compiler with signed 128-bit arithmetic.
Use `--compiler g++` or `--compiler clang++` to select a compiler. On
Windows, run from the x64 Native Tools Command Prompt for Visual Studio
2022 and pass `--compiler cl`.

The deposited [reference results](certificates/reference-results.json)
and [SHA-256 manifest](certificates/SHA256SUMS.txt) are included.

## Build the manuscript

Run `tectonic f26.tex` in [the manuscript directory](manuscript/), or run
`pdflatex f26.tex` twice. The [source notes](manuscript/README.txt)
document the release build and figure generation.

Graph numbering follows [Bodnár et al., Table 3, version 3](https://arxiv.org/html/2606.00290v3#S5.T3).
