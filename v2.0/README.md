# Version 2.0 · Substitution bounds for the inducibility of six-vertex graphs

Oleg Kyzlasov, 22 November 2026. [Preprint on Zenodo](https://zenodo.org/records/24437581) · [DOI 10.5281/zenodo.24437581](https://doi.org/10.5281/zenodo.24437581)

The paper retains the F26 results and adds five cases:

- **F74.** Balanced substitutions of the complement of the 3-by-3 rook
  graph give limiting density 360/7381. The paper proves a sharp counting
  inequality in the rook substitution class and supplies a rational
  upper bound for unrestricted graphs.
- **F22, F29, F34 and F62.** Explicit constructions exceed the respective
  conjectured values in version 3 of the reference paper.

The [manuscript source](manuscript/) includes both figures.
The [computational supplement](certificates/) includes exact data,
verification programs and the full F26 replay from version 1.0.

## Verify the certificates

From the repository root:

~~~sh
cd v2.0/certificates
python -B verify.py
~~~

The full replay uses Python 3.10 or later and a C++20 compiler with signed
128-bit arithmetic on a little-endian system. It runs offline and takes
several minutes. Use `--compiler g++` or `--compiler clang++` to select a
compiler. On Windows, use the x64 Native Tools Command Prompt for Visual
Studio 2022 and pass `--compiler cl`. A successful run ends with
`Verified: all`.

| Option | Checks |
| --- | --- |
| `--part cases` | The five added constructions and the weighted F74 inequality. This needs only Python. |
| `--part f74-upper` | The unrestricted eight-vertex upper certificate for F74. |
| `--part f26` | Both F26 bounds and the restricted equality calculation. |

See the [certificate documentation](certificates/README.txt) for data
formats, provenance and the [SHA-256 manifest](certificates/SHA256SUMS.txt).

## Build the manuscript

Run `tectonic six-vertex.tex` in [the manuscript directory](manuscript/),
or run `pdflatex six-vertex.tex` twice. The [source notes](manuscript/README.txt)
document the release build and figure generation.

Graph and conjecture numbering follow [Bodnár et al., version 3](https://arxiv.org/abs/2606.00290v3).
