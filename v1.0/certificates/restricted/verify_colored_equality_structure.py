"""Exact equality structure extracted from the positive-definite F-root block.

This checks the finite ingredients of the triangle-free equality lemma.
The PSD premise is supplied by verify_colored_exact_certificate.py.
"""

from pathlib import Path
from itertools import combinations, permutations
from collections import Counter
import argparse, json

ROOT = Path(__file__).resolve().parent


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-root", default=str(ROOT / "data"))
    ap.add_argument("--name", default="colored-equality-structure-487")
    args = ap.parse_args()
    data = Path(args.data_root)
    schema = json.loads(
        (data / "colored-triangle-free-products-462/moments.schema.json").read_text()
    )
    cert = json.loads(
        (data / "colored-exact-certificate-480/certificate.json").read_text()
    )

    def matrix(n, mask):
        A = [[0] * n for _ in range(n)]
        for bit, (i, j) in enumerate(combinations(range(n), 2)):
            A[i][j] = A[j][i] = int(mask >> bit & 1)
        return A

    def isF(A, S):
        d = {i: sum(A[i][j] for j in S) for i in S}
        if sorted(d.values()) != [1, 1, 1, 2, 2, 3]:
            return False
        root = next(i for i in S if d[i] == 3)
        if sorted(d[j] for j in S if A[root][j]) != [1, 2, 2]:
            return False
        seen = {root}
        todo = [root]
        while todo:
            v = todo.pop()
            for w in S:
                if A[v][w] and w not in seen:
                    seen.add(w)
                    todo.append(w)
        return len(seen) == 6

    candidates = [
        (t, b)
        for t, b in enumerate(schema["blocks"])
        if b["k"] == 6
        and b["type_blue"] == 0
        and isF(matrix(6, b["type_red"]), tuple(range(6)))
    ]
    assert len(candidates) == 1
    t, b = candidates[0]
    record = next(r for r in cert["blocks"] if r["original_block"] == t)
    d = b["dimension"]
    assert d == 28 and len(record["kernel_columns"]) == d - 1
    F = matrix(6, b["type_red"])
    pairs = list(combinations(range(5), 2))
    P = [[int(set(a).isdisjoint(c)) for c in pairs] for a in pairs]
    roots = next(
        p
        for p in permutations(range(10), 6)
        if all(F[i][j] == P[p[i]][p[j]] for i, j in combinations(range(6), 2))
    )
    flags = b["flag_representatives"]
    index = {c: i for i, c in enumerate(flags)}
    extension_indices = []
    for v in range(10):
        labels = roots + (v,)
        red = blue = 0
        for bit, (i, j) in enumerate(combinations(range(7), 2)):
            if labels[i] == labels[j]:
                blue |= 1 << bit
            elif P[labels[i]][labels[j]]:
                red |= 1 << bit
        extension_indices.append(index[red | (blue << 21)])
    v = [extension_indices.count(i) for i in range(d)]
    assert sorted(v) == [0] * 18 + [1] * 10
    B = [[0] * (d - 1) for _ in range(d)]
    for j, col in enumerate(record["kernel_columns"]):
        for i, a in col:
            B[i][j] += a
        assert sum(v[i] * a for i, a in col) == 0
    # A nonzero pivot minor modulo a prime proves full rank over the reals.
    modulus = 1000000007
    A = [[a % modulus for a in row] for row in B]
    rank = 0
    pivots = []
    for j in range(d - 1):
        pivot = next((i for i in range(rank, d) if A[i][j]), None)
        assert pivot is not None
        A[rank], A[pivot] = A[pivot], A[rank]
        inv = pow(A[rank][j], -1, modulus)
        A[rank] = [(a * inv) % modulus for a in A[rank]]
        for i in range(rank + 1, d):
            a = A[i][j]
            if a:
                A[i] = [(x - a * y) % modulus for x, y in zip(A[i], A[rank])]
        pivots.append(pivot)
        rank += 1
    assert rank == 27
    # Six blue flags force each root atom to mass1/10. Four remaining normal
    # profiles have mass1/10; every atom in such a profile extends five roots
    # to an F, forcing that atom itself to mass1/10 by the same F-root block.
    root_flags = extension_indices[:]
    assert len(set(extension_indices)) == 10
    witnesses = []
    for v0 in range(10):
        if v0 in roots:
            continue
        S = next(S for S in combinations(roots, 5) if isF(P, tuple(S) + (v0,)))
        witnesses.append(dict(outside_P_vertex=v0, five_root_vertices=list(S)))
    assert len(witnesses) == 4
    order = roots + tuple(v0 for v0 in range(10) if v0 not in roots)
    T = [[P[order[i]][order[j]] for j in range(10)] for i in range(10)]
    unknown = list(combinations(range(6, 10), 2))
    hist = Counter()
    maxcopies = -1
    winners = []
    for mask in range(64):
        G = [row[:] for row in T]
        for bit, (i, j) in enumerate(unknown):
            G[i][j] = G[j][i] = int(mask >> bit & 1)
        if any(
            G[i][j] and G[i][k] and G[j][k] for i, j, k in combinations(range(10), 3)
        ):
            continue
        f = sum(isF(G, S) for S in combinations(range(10), 6))
        hist[f] += 1
        if f > maxcopies:
            maxcopies = f
            winners = [mask]
        elif f == maxcopies:
            winners.append(mask)
    actual = sum(T[i][j] << bit for bit, (i, j) in enumerate(unknown))
    assert maxcopies == 60 and winners == [actual]
    out = dict(
        status="weighted_triangle_free_equality_structure_exactly_extracted",
        F_root_type=t,
        F_root_type_red_mask=b["type_red"],
        flag_dimension=d,
        kernel_rank_exact=rank,
        kernel_rank_proof_modulus=modulus,
        kernel_orthogonal_complement="span of ten normal flag indicators",
        normal_flag_indices=extension_indices,
        rooted_embedding_into_Petersen=list(roots),
        outside_profile_F_extension_witnesses=witnesses,
        triangle_free_completion_F_count_histogram=dict(sorted(hist.items())),
        maximum_F_six_subsets=maxcopies,
        unique_maximizing_completion_mask=actual,
        equality_theorem="After removing zero weights, equality with positive F density holds exactly for Petersen with all ten weights1/10.",
        strict_positive_definiteness_premise="Verified in verify_colored_exact_certificate.py.",
    )
    (Path.cwd() / (args.name + ".json")).write_text(json.dumps(out, indent=2))
    print(json.dumps(out))


if __name__ == "__main__":
    main()
