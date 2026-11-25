Substitution bounds for the inducibility of six-vertex graphs
Oleg Kyzlasov | Version 2.0 | Computational supplement

Extract this archive, then run from the directory containing this file:

  python -B verify.py

The full check runs offline on a little-endian system. Use Python 3.10+
with assertions enabled and a C++20 compiler with signed 128-bit arithmetic.
Use the x64 Native Tools command prompt on Windows (Visual Studio 2022):

  python -B verify.py --compiler cl

GCC or Clang can be selected with --compiler g++ or --compiler clang++.
All compilation and intermediate files go into a temporary directory.
Allow several minutes. Success ends with "Verified: all". A missing
file, failed inequality, wrong expected result or compiler failure gives
a nonzero exit code. The checks require no solver or third-party Python
modules. The manifest SHA256SUMS.txt covers all the other deposited files.

PARTIAL REPLAYS

  python -B verify.py --part cases
     Reconstruct the five additional lower witnesses. Includes the sharp
     weighted inequality on the nine rook cells. Python alone is sufficient.

  python -B verify.py --part f74-upper
     Recount the unrestricted eight-vertex upper certificate for F74.

  python -B verify.py --part f26
     Reproduce both F26 bounds and the restricted equality calculation.
     The contents of f26/ are the first edition's certificate archive.
     Its README documents its coloured catalogues and binary tables.

CASES AND NORMALIZATION

cases.json specifies labelled targets and templates. The characters ij
denote the edge {i,j}. A pair between distinct classes is complete if
listed and empty otherwise. "cliques" gives the constant classes with
diagonal bit one; the other constant classes have diagonal bit zero.
"recursive" replaces a class by another copy of the entire pattern.
The two-level test finishes every remaining recursive class at bit zero.

verify_cases.py counts both the surviving maps on the template and the
maps on the flattened matrix with the rational terminal weights. For F29
it visits all 1,947,792 six-subsets of the 36-vertex graph and checks the
weight formula for two disjoint components. For F74 it reconstructs all
84 six-subsets of the board, keeps the 36 targets and expands every
coefficient on every chamber modulo row, column and transpose symmetries.
Polynomial coefficients are unsigned integers. Disjoint guarded 32-bit
digits permit all coefficient comparisons in an ordering to be made
with Python integers. The bounds that prevent borrows and carries are
checked before the loop. expected-cases.json is a reference output.
The enumerations do not use that file to compute any count or sign.

F74 UPPER CERTIFICATE

f74/factors.txt consists entirely of integers. Its first line gives the
number of blocks, number of hosts and scale D=2^40. For each block the
header gives flag order, root order, root mask, dimension and extension
code width. Next come the extension classifier and the square array of
integer factors. Each host row then specifies a 28-bit simple graph mask
and its induced copy count on six-subsets. Bit positions enumerate
unordered vertex pairs in lexicographic order. The verifier recomputes
the count and uses A A^T / D^2 for each matrix A in the file.

f26/unrestricted/grammar.py checks invariance under permutations of
the unlabelled extension vertices. f74/eight.cpp checks the Gram entries,
graph orbits and the maximum over 12,346 host classes. The raw maximum is
2377428879170370923861228790 over 40320 * 2^80. The resulting upper
bound is 11321089900811290113624899/232113757366008801543585792.
This is a global upper estimate. The sharp result in the manuscript is
for rook substitutions. The code here does not test a structural theorem
for arbitrary extremal graphons.

PROVENANCE

Graph numbering, conjecture numbers and table entries refer throughout
to Bodnar, Gao, Leon, Liu, Pikhurko and Sun, arXiv:2606.00290v3,
https://arxiv.org/abs/2606.00290v3

The F74 integer factors were obtained from the authors' numerical
matrices by spectral decomposition, truncating negative eigenvalues to
zero and rounding the scaled square-root factors to integers. The proof
uses only the deposited integers. The source for the numerical proposal:
https://github.com/bodnalev/supplementary_files/blob/9f356bb3dc10dcacd4d26a4eb1d556c2650e4e1b/graph_inducibility_6vtx/certs/inexact_74.pickle
SHA-256 of that source:
c10963465f06ccae8f482bbb9ecc093bfbb5ec7e59e6be68d350b8b4b5b11c68

The first version, including the F26 source and certificates, is at
https://doi.org/10.5281/zenodo.24341542
