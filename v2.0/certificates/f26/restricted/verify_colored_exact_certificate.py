"""Standard-Python independent PSD and ALL-host replay in ORIGINAL flags.

The recovery program used projected Kronecker operators and FLINT.
This verifier uses neither: it forms BQB^T with Python integers, verifies
Cholesky-error diagonal dominance, then applies the original flag census.
Full census/catalog completeness audits are separate required evidence.
"""

from pathlib import Path
from fractions import Fraction as Fr
from math import lcm
from array import array
from itertools import combinations
import argparse, json, time, sys, hashlib

ROOT = Path(__file__).resolve().parent
BIN = ROOT / "data"
if hasattr(sys, "set_int_max_str_digits"):
    sys.set_int_max_str_digits(0)


def read_array(path, kind):
    a = array(kind)
    with path.open("rb") as f:
        a.fromfile(f, path.stat().st_size // a.itemsize)
    if sys.byteorder != "little":
        a.byteswap()
    return a


def main():
    global BIN
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-root", default=str(BIN))
    ap.add_argument("--certificate", default=None)
    ap.add_argument("--name", default="colored-exact-independent-481")
    args = ap.parse_args()
    BIN = Path(args.data_root)
    start = time.monotonic()
    path = (
        Path(args.certificate)
        if args.certificate
        else BIN / "colored-exact-certificate-480/certificate.json"
    )
    cert = json.loads(path.read_text())
    assert cert["schema"] == "colored_TF_order8_exact_PSD_certificate_v1"
    source = BIN / cert["source_flag_products"]
    p = source / "moments"
    schema = json.loads(Path(str(p) + ".schema.json").read_text())
    assert sorted(record["original_block"] for record in cert["blocks"]) == list(
        range(len(schema["blocks"]))
    )
    catalog = json.loads(
        (BIN / cert["source_colored_catalog"] / "catalog-8.json").read_text()
    )
    n = len(catalog)
    D = cert["base_Gram_denominator"]
    assert D > 0
    ptr = read_array(Path(str(p) + ".ptr.u64"), "Q")
    ix = read_array(Path(str(p) + ".indices.u32"), "I")
    count = read_array(Path(str(p) + ".counts.u32"), "I")
    assert len(ptr) == n + 1 and ptr[-1] == len(ix) == len(count)
    columns = [[] for _ in range(schema["packed_columns"])]
    for h in range(n):
        for z in range(ptr[h], ptr[h + 1]):
            columns[ix[z]].append((h, int(count[z])))
    # Recount the objective from adjacency, independently of the generator's
    # isomorphism forms and its stored six-vertex deck.
    target = []
    for entry in catalog:
        red = [[False] * 8 for _ in range(8)]
        blue = [[False] * 8 for _ in range(8)]
        for z, (i, j) in enumerate(combinations(range(8), 2)):
            red[i][j] = red[j][i] = bool(entry["red_mask"] >> z & 1)
            blue[i][j] = blue[j][i] = bool(entry["blue_mask"] >> z & 1)
            assert not (red[i][j] and blue[i][j])
        f = b = 0
        for S in combinations(range(8), 6):
            if all(blue[i][j] for i, j in combinations(S, 2)):
                b += 1
            degrees = {i: sum(red[i][j] for j in S) for i in S}
            if sorted(degrees.values()) != [1, 1, 1, 2, 2, 3]:
                continue
            center = next(i for i in S if degrees[i] == 3)
            if sorted(degrees[j] for j in S if red[center][j]) != [1, 2, 2]:
                continue
            reached = {center}
            todo = [center]
            while todo:
                v = todo.pop()
                for w in S:
                    if red[v][w] and w not in reached:
                        reached.add(w)
                        todo.append(w)
            if len(reached) == 6:
                f += 1
        numerator = 11111 * f + 480 * b - 13440
        assert numerator == entry["target_integer_numerator"]
        assert f == entry["F_six_subsets"] and b == entry["blue_six_subsets"]
        target.append(numerator)
    scores = [Fr(v, 13440) for v in target]
    margins = []
    Gram_dimensions = 0
    terms = 0
    for b, record in enumerate(cert["blocks"]):
        t = record["original_block"]
        s = schema["blocks"][t]
        k = len(record["base_Gram_numerators"])
        assert k == len(record["kernel_columns"])
        B = [[] for _ in range(s["dimension"])]
        for a, col in enumerate(record["kernel_columns"]):
            for i, v in col:
                assert 0 <= i < len(B) and isinstance(v, int)
                B[i].append((a, v))
        Q0 = record["base_Gram_numerators"]
        assert len(Q0) == k and all(len(row) == k for row in Q0)
        assert all(
            Q0[i][j] == Q0[j][i] and isinstance(Q0[i][j], int)
            for i in range(k)
            for j in range(k)
        )
        den = D
        for c in record["corrections"]:
            assert int(c["denominator"]) > 0
            den = lcm(den, int(c["denominator"]))
        Q = [[v * (den // D) for v in row] for row in Q0]
        for c in record["corrections"]:
            i = c["i"]
            j = c["j"]
            v = int(c["numerator"]) * (den // int(c["denominator"]))
            assert 0 <= i <= j < k
            Q[i][j] += v
            if i != j:
                Q[j][i] += v
        F = record["PSD_factor_lower_numerators"]
        FD = record["PSD_factor_denominator"]
        assert (
            FD > 0 and len(F) == k and all(len(row) == i + 1 for i, row in enumerate(F))
        )
        common = lcm(den, FD**2)
        a = common // den
        bb = common // (FD**2)
        offsum = [0] * k
        diag = [0] * k
        for j in range(k):
            for i in range(j + 1):
                gram = sum(F[i][z] * F[j][z] for z in range(i + 1))
                res = Q[i][j] * a - gram * bb
                if i == j:
                    diag[i] = res
                else:
                    offsum[i] += abs(res)
                    offsum[j] += abs(res)
        margin = min(diag[i] - offsum[i] for i in range(k))
        assert margin > 0
        assert Fr(margin, common) == Fr(
            record["residual_strict_diagonal_dominance_margin"]
        )
        margins.append(Fr(margin, common))
        Gram_dimensions += k
        # Expand the actual flag Gram matrix, without projected operators.
        off = s["offset"]
        div = record["operator_divisor"]
        assert isinstance(div, int) and div > 0
        packed = 0
        for j in range(s["dimension"]):
            for i in range(j + 1):
                coefficient = sum(v * w * Q[u][z] for u, v in B[i] for z, w in B[j])
                if coefficient:
                    for h, m in columns[off + packed]:
                        scores[h] += Fr(m * coefficient, den * div)
                        terms += 1
                packed += 1
        assert packed == s["dimension"] * (s["dimension"] + 1) // 2
        if b % 30 == 0:
            print(
                json.dumps(
                    dict(
                        status="independent_exact_PSD_and_original_flag_replay",
                        blocks_done=b + 1,
                    )
                ),
                flush=True,
            )
    assert len(cert["blocks"]) == 146 and all(v <= 0 for v in scores)
    pop = json.loads(
        (
            BIN / cert["source_colored_catalog"] / "Petersen-population-8.json"
        ).read_text()
    )
    assert all(scores[h] == 0 for h, w in enumerate(pop["integer_weights"]) if w)
    out = dict(
        status="all_PSD_and_original_flag_host_inequalities_independently_verified",
        colored_hosts=n,
        Gram_blocks=146,
        total_kernel_Gram_dimensions=Gram_dimensions,
        exact_nonpositive_host_inequalities=n,
        exact_zero_host_inequalities=sum(v == 0 for v in scores),
        original_flag_coefficient_contributions_replayed=terms,
        minimum_exact_PSD_DD_margin=str(min(margins)),
        no_projected_operator_data_used=True,
        no_external_arithmetic_libraries_used=True,
        objective_independently_recounted_from_all_actual_red_blue_adjacencies=True,
        certificate_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        seconds=time.monotonic() - start,
    )
    (Path.cwd() / (args.name + ".json")).write_text(json.dumps(out, indent=2))
    print(json.dumps(out), flush=True)


if __name__ == "__main__":
    main()
