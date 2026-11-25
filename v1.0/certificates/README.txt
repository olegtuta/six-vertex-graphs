Petersen substitutions and the inducibility of the six-vertex tree F26
Oleg Kyzlasov | Preprint v1.0 | Computational supplement

REPRODUCING THE RESULTS

Extract the ZIP to a directory. Run the commands below from the directory
containing verify.py. Use ordinary Python 3.10 or later with assertions
enabled. The complete run needs a little-endian system and a C++20 compiler
with signed 128-bit arithmetic. The computation runs offline. No graph
library, numerical solver, Python package installation or precompiled
executable is needed. The programs use the system temporary directory.

  python -B verify.py --part elementary
      Recount the Petersen copies by isomorphism and check the replacement
      examples in Table 1 of the preprint. The reported difference and
      relative gain are exact arithmetic on the stated endpoints.

  python -B verify.py --part restricted
      Also verify the coloured catalogue, full flag census, positive
      definiteness, all host inequalities, and the equality calculation.
      These two commands require only the Python standard library.

  python -B verify.py --output fresh-results.json
      Run all of the above and the unrestricted upper certificate. This
      includes compilation from source. Allow several minutes. The exit
      code is zero on success. The optional output file can be compared
      directly with reference-results.json. The latter is a saved output;
      the proof calculations do not use its reported values.

On Linux with GCC or Clang, or macOS with the command line tools installed,
the driver can select c++ automatically. To choose a compiler explicitly:

  python -B verify.py --compiler g++ --output fresh-results.json

On Windows use the x64 Native Tools Command Prompt for Visual Studio 2022,
with the C++ workload installed:

  python -B verify.py --compiler cl --output fresh-results.json

The driver compiles eight.cpp, census.cpp and nine.cpp. The GCC/Clang
arguments are -O3 -std=c++20; the MSVC arguments are /O2 /EHsc /std:c++20.
The MSVC build uses the signed 128-bit integer supplied with the VS 2022
C++ standard library. The distributed tables are little-endian. To run
only the unrestricted replay, add --part unrestricted. The published
package was tested with 64-bit CPython 3.12 and MSVC on Windows.

MATHEMATICAL CHECKS

1. The target has edges 01, 02, 05, 13, 45. In KG(5,2), 60 of the 210
   six-subsets induce the target. Each centre lies at the centre of six
   copies. Recursive substitution gives 480/11111 in the limit.

2. restricted/data/colored-triangle-free-catalog-461 contains the host
   catalogues at orders 1 through 8. For a pair of positions i<j, take its
   bit index in lexicographic order, starting at zero. A set red bit means
   an edge between distinct quotient vertices. A set blue bit means the
   positions denote the same quotient vertex. All other pairs are white.
   The catalogue checker enumerates labelled triangle-free graphs and
   checks class-size orbits and all Petersen population weights. The counts
   at orders 1 through 8 are 1, 3, 6, 17, 41, 125, 388, 1491.

3. restricted/data/colored-triangle-free-products-462 describes the flag
   types and flag representatives. The sparse count matrix has one row per
   host. moments.ptr.u64 holds the row endpoints, moments.indices.u32 the
   packed column indices, and moments.counts.u32 their multiplicities.
   moments.P-moments.u64 is the separate Petersen population cross-check.
   In each symmetric block column (i,j), i<=j, use j(j+1)/2+i plus the
   block offset. The diagonal stores M_ii and the off-diagonal stores
   2M_ij. The checker independently reconstructs all 534006 nonzero
   coefficients and checks the zero entries by comparing whole sparse rows.

4. restricted/data/colored-exact-certificate-480/certificate.json contains
   the matrices of the coloured inequality. kernel_columns gives B by
   sparse integer columns. Divide base_Gram_numerators by the global
   base_Gram_denominator, then add the listed rational symmetric
   corrections to obtain Q. operator_divisor is d in the manuscript.
   The lower triangular factor and its denominator give C; the verifier
   checks strict diagonal dominance of Q-C C^T in exact arithmetic.
   It then evaluates every original host inequality, with 276 tight rows.
   Direct computations in the six-root block verify the equality argument.
   Its index is 142 (zero-based), with 28 flags and rank 27 for B.

5. unrestricted/eight.txt begins with the counts 169, 12346 and the
   denominator 2^42. Each type is given by a flag order, root order, root
   mask, number of flags and extension width, followed by a lookup row
   and the square integer factor. The tail gives graph masks and their
   six-subset counts, which the checker recounts. unrestricted/one.txt
   holds the additional root-one classifier and its 90-by-90 factor.
   unrestricted/hosts/nine.u64 contains the nine-vertex adjacency masks.
   nine.target.u8 contains their target counts, and nine.sub8.u16 the
   nine deletion indices per graph. The census verifier independently
   checks every deletion and every orbit. A Gram factor has entries of
   absolute value below 2^47 and at most 272 columns. Thus an unscaled
   Gram entry has absolute value below 2^103. Summing at most 3*8! entries
   and the target leaves the eight-vertex numerator below 2^121 in
   absolute value. The nine lifted scores and the 630 root-one terms
   (with multiplier 576) together remain below 2^126. The native
   checkers use signed 128-bit arithmetic for these calculations.

   The nine-vertex maximum is
     233973195106737927117969229 / 5415987671873538702683668480.
   Subtracting 480/11111 gives
     2088331666530919595233019 / 60177039022186888525518240481280.

6. The deduction for all (bull, house)-free graphs uses primeness of F
   and Theorem 1.3 of Belavadi and Hoang. The preprint contains the
   induction, the equality proof and the sampling arguments connecting
   both finite certificates to arbitrary graph orders.

PROVENANCE AND INTEGRITY

SHA256SUMS.txt lists the distributed proof inputs and source files with
SHA-256 digests. The driver validates these before replay. Archive input
files are read-only during verification. The reference report is included
for comparison and is not a source of certificate coefficients.

Row 26, Table 3, Section 5 of the reference article:
https://arxiv.org/html/2606.00290v3#S5.T3
The graph can be read clockwise from the top in its F_i column. To repeat
the source check, save the HTML for that version locally and run
  python -B check_table3.py article-v3.html
It reads both displayed bounds, N, and the graph's five drawn edges directly
from the specified row. This optional source check is separate from the
offline mathematical verification.

The numerical starting data for the unrestricted order-eight calculation
come from Bodnar, Gao, Leon, Liu, Pikhurko and Sun's supplement:
https://github.com/bodnalev/supplementary_files/blob/9f356bb3dc10dcacd4d26a4eb1d556c2650e4e1b/graph_inducibility_6vtx/certs/inexact_26.pickle
SHA-256 d3847ea2f168027593f997a5d69851782d5bf1e4d1547ab2acf0c200a75ad605
The proof uses the integer factors in this archive. The original numerical
matrices and their floating-point diagnostics are unnecessary for replay.
The external decomposition theorem is available at
https://arxiv.org/html/2604.27594v3
