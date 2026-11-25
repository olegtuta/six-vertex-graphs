"""Replay every original colored flag coefficient using standard Python.

Root subsets and explicit root isomorphisms replace the C++ generator's
ordered-tuple enumeration. No numerical or graph libraries are used.
"""

from pathlib import Path
from itertools import combinations, permutations
from collections import defaultdict
from array import array
import argparse, json, sys, time

ROOT = Path(__file__).resolve().parent
DEFAULT_DATA = ROOT / "data"


def read_array(path, kind):
    values = array(kind)
    with path.open("rb") as f:
        values.fromfile(f, path.stat().st_size // values.itemsize)
    if sys.byteorder != "little":
        values.byteswap()
    return values


def matrix(n, red, blue):
    A = [[0] * n for _ in range(n)]
    assert (
        not (red & blue)
        and red >= 0
        and blue >= 0
        and (red | blue) < 1 << (n * (n - 1) // 2)
    )
    for z, (i, j) in enumerate(combinations(range(n), 2)):
        A[i][j] = A[j][i] = 2 if blue >> z & 1 else 1 if red >> z & 1 else 0
    return A


def key(A, S):
    r = b = 0
    for z, (i, j) in enumerate(combinations(S, 2)):
        if A[i][j] == 1:
            r |= 1 << z
        elif A[i][j] == 2:
            b |= 1 << z
    return r | (b << (len(S) * (len(S) - 1) // 2))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-root", default=str(DEFAULT_DATA))
    ap.add_argument("--name", default="colored-portable-full-census-485")
    args = ap.parse_args()
    data = Path(args.data_root)
    start = time.monotonic()
    p = data / "colored-triangle-free-products-462/moments"
    schema = json.loads(Path(str(p) + ".schema.json").read_text())
    blocks = schema["blocks"]
    n = schema["host_count"]
    columns = schema["packed_columns"]
    cat = json.loads(
        (data / "colored-triangle-free-catalog-461/catalog-8.json").read_text()
    )
    assert n == len(cat) == 1491
    ptr = read_array(Path(str(p) + ".ptr.u64"), "Q")
    ix = read_array(Path(str(p) + ".indices.u32"), "I")
    count = read_array(Path(str(p) + ".counts.u32"), "I")
    assert len(ptr) == n + 1 and ptr[0] == 0 and ptr[-1] == len(ix) == len(count)
    types = {}
    flags = []
    offset = 0
    for t, block in enumerate(blocks):
        k = block["k"]
        nf = block["n"]
        d = block["dimension"]
        assert k in (0, 2, 4, 6) and nf == (8 + k) // 2
        assert block["offset"] == offset and d == len(block["flag_representatives"])
        offset += d * (d + 1) // 2
        typecode = block["type_red"] | (block["type_blue"] << (k * (k - 1) // 2))
        assert (k, typecode) not in types
        types[k, typecode] = t
        lookup = {}
        bits = nf * (nf - 1) // 2
        for i, code in enumerate(block["flag_representatives"]):
            A = matrix(nf, code & ((1 << bits) - 1), code >> bits)
            assert key(A, tuple(range(k))) == typecode
            for ext in permutations(range(k, nf)):
                c = key(A, tuple(range(k)) + ext)
                assert c not in lookup or lookup[c] == i
                lookup[c] = i
        flags.append(lookup)
    assert offset == columns and len(blocks) == 146
    expected_types = {(0, 0)}
    for k in (2, 4, 6):
        for r in json.loads(
            (data / f"colored-triangle-free-catalog-461/catalog-{k}.json").read_text()
        ):
            expected_types.add(
                (k, r["red_mask"] | (r["blue_mask"] << (k * (k - 1) // 2)))
            )
    assert set(types) == expected_types
    root_lookup = {k: {} for k in (0, 2, 4, 6)}
    for (k, c), t in types.items():
        bits = k * (k - 1) // 2
        A = matrix(k, c & ((1 << bits) - 1), c >> bits)
        for q in permutations(range(k)):
            raw = key(A, q)
            inv = tuple(q.index(i) for i in range(k))
            if raw not in root_lookup[k]:
                root_lookup[k][raw] = (t, [])
            assert root_lookup[k][raw][0] == t
            root_lookup[k][raw][1].append(inv)
    pop = json.loads(
        (
            data / "colored-triangle-free-catalog-461/Petersen-population-8.json"
        ).read_text()
    )
    weights = pop["integer_weights"]
    assert len(weights) == n
    moments = [0] * columns
    nnz = 0
    for h, entry in enumerate(cat):
        A = matrix(8, entry["red_mask"], entry["blue_mask"])
        expected = defaultdict(int)
        for k in (0, 2, 4, 6):
            u = (8 - k) // 2
            mult = {0: 576, 2: 36, 4: 4, 6: 1}[k]
            for S in combinations(range(8), k):
                match = root_lookup[k].get(key(A, S))
                if match is None:
                    continue
                t, orders = match
                remaining = tuple(v for v in range(8) if v not in S)
                off = blocks[t]["offset"]
                for p0 in orders:
                    roots = tuple(S[v] for v in p0)
                    for E in combinations(remaining, u):
                        F = tuple(v for v in remaining if v not in E)
                        i = flags[t][key(A, roots + E)]
                        j = flags[t][key(A, roots + F)]
                        if i > j:
                            i, j = j, i
                        expected[off + j * (j + 1) // 2 + i] += mult
        actual = {}
        for z in range(ptr[h], ptr[h + 1]):
            c = ix[z]
            assert c < columns and c not in actual and count[z] > 0
            actual[c] = count[z]
            moments[c] += weights[h] * count[z]
        assert dict(expected) == actual, (h, "original flag census differs")
        nnz += len(actual)
        if (h + 1) % 250 == 0:
            print(
                json.dumps(
                    dict(
                        status="portable_all_host_census",
                        hosts=h + 1,
                        seconds=time.monotonic() - start,
                    )
                ),
                flush=True,
            )
    native = read_array(Path(str(p) + ".P-moments.u64"), "Q")
    assert list(native) == moments
    out = dict(
        status="complete_original_flag_census_independently_replayed_with_standard_Python",
        colored_hosts=n,
        root_types=len(types),
        packed_columns=columns,
        exact_nonzero_coefficients=nnz,
        every_original_flag_product_coefficient_checked=True,
        every_P_moment_entry_checked=True,
        no_external_libraries_used=True,
        normalization=40320,
        seconds=time.monotonic() - start,
    )
    (Path.cwd() / (args.name + ".json")).write_text(json.dumps(out, indent=2))
    print(json.dumps(out), flush=True)


if __name__ == "__main__":
    main()
