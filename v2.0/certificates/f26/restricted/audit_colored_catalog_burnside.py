"""Independent Burnside counts and embedding/stabilizer P population law.

This uses graph-isomorphism backtracking, not the catalog generator's
degree-ordered canonical search or its occupancy enumeration.
"""

from pathlib import Path
from itertools import combinations
from math import factorial, prod
from collections import defaultdict
import argparse, json, time

ROOT = Path(__file__).resolve().parent
BIN = ROOT / "data"


def main():
    global BIN
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-root", default=str(BIN))
    ap.add_argument("--verify-labelled-counts", action="store_true")
    ap.add_argument("--name", default="colored-catalog-independent-476")
    args = ap.parse_args()
    BIN = Path(args.data_root)
    start = time.monotonic()
    src = BIN / "colored-triangle-free-catalog-461"
    cats = {n: json.loads((src / f"catalog-{n}.json").read_text()) for n in range(1, 9)}

    def matrix(n, mask):
        A = [[0] * n for _ in range(n)]
        for k, (i, j) in enumerate(combinations(range(n), 2)):
            A[i][j] = A[j][i] = (mask >> k) & 1
        return A

    def embeddings(A, B, onto=False):
        n = len(A)
        m = len(B)
        order = sorted(range(n), key=lambda i: -sum(A[i]))
        mapping = [-1] * n
        result = []

        def extend(k, used):
            if k == n:
                result.append(tuple(mapping))
                return
            i = order[k]
            for j in range(m):
                if used >> j & 1 or (
                    sum(A[i]) != sum(B[j]) if onto else sum(A[i]) > sum(B[j])
                ):
                    continue
                if any(A[i][u] != B[j][mapping[u]] for u in order[:k]):
                    continue
                mapping[i] = j
                extend(k + 1, used | (1 << j))
            mapping[i] = -1

        extend(0, 0)
        return result

    labels = list(combinations(range(5), 2))
    P = [[int(set(a).isdisjoint(b)) for b in labels] for a in labels]
    quotients = {}
    groups = {}
    Pemb = {}
    for n in range(1, 9):
        q = [r for r in cats[n] if r["quotient_order"] == n]
        assert len(q) == [0, 1, 2, 3, 7, 14, 38, 107, 410][n]
        buckets = defaultdict(list)
        for r in q:
            mask = r["quotient_mask"]
            A = matrix(n, mask)
            assert not any(
                A[i][j] and A[i][k] and A[j][k] for i, j, k in combinations(range(n), 3)
            )
            key = tuple(sorted(map(sum, A)))
            for B in buckets[key]:
                assert not embeddings(A, B, onto=True)
            buckets[key].append(A)
            quotients[n, mask] = A
            groups[n, mask] = embeddings(A, A, onto=True)
            Pemb[n, mask] = len(embeddings(A, P))
    # Burnside: a fixed composition is constant on each permutation cycle.
    predicted = {n: 0 for n in range(1, 9)}
    for (k, mask), G in groups.items():
        for total in range(k, 9):
            fixedsum = 0
            for p in G:
                seen = set()
                cycles = []
                for v in range(k):
                    if v in seen:
                        continue
                    u = v
                    length = 0
                    while u not in seen:
                        seen.add(u)
                        length += 1
                        u = p[u]
                    cycles.append(length)
                dp = [0] * (total - k + 1)
                dp[0] = 1
                for length in cycles:
                    for d in range(length, len(dp)):
                        dp[d] += dp[d - length]
                fixedsum += dp[-1]
            assert fixedsum % len(G) == 0
            predicted[total] += fixedsum // len(G)
    assert predicted == {n: len(cats[n]) for n in cats}
    for n, entries in cats.items():
        seen = set()
        for r in entries:
            k = r["quotient_order"]
            mask = r["quotient_mask"]
            m = r["class_sizes"]
            assert len(m) == k and sum(m) == n and all(x > 0 for x in m)
            T = quotients[k, mask]
            labels = [i for i, v in enumerate(m) for _ in range(v)]
            red = blue = 0
            for bit, (i, j) in enumerate(combinations(range(n), 2)):
                if labels[i] == labels[j]:
                    blue |= 1 << bit
                elif T[labels[i]][labels[j]]:
                    red |= 1 << bit
            assert (red, blue) == (r["red_mask"], r["blue_mask"])
            canonical_m = min(tuple(m[p[i]] for i in range(k)) for p in groups[k, mask])
            key = k, mask, canonical_m
            assert key not in seen
            seen.add(key)
    orbit_sums = {
        k: sum(factorial(k) // len(G) for (n, mask), G in groups.items() if n == k)
        for k in range(1, 9)
    }
    labelled = {}
    if args.verify_labelled_counts:
        counts = [0] * 9
        adj = [0] * 8

        def enumerate_graphs(k):
            counts[k] += 1
            valid = [True] + [False] * ((1 << k) - 1)
            for s in range(1, 1 << k):
                bit = s & -s
                v = bit.bit_length() - 1
                rest = s ^ bit
                valid[s] = valid[rest] and not (adj[v] & rest)
            if k == 7:
                counts[8] += sum(valid)
                return
            for s, ok in enumerate(valid):
                if not ok:
                    continue
                adj[k] = s
                bits = s
                while bits:
                    bit = bits & -bits
                    v = bit.bit_length() - 1
                    adj[v] |= 1 << k
                    bits ^= bit
                enumerate_graphs(k + 1)
                bits = s
                while bits:
                    bit = bits & -bits
                    v = bit.bit_length() - 1
                    adj[v] ^= 1 << k
                    bits ^= bit

        enumerate_graphs(0)
        labelled = {k: counts[k] for k in range(1, 9)}
        assert labelled == orbit_sums
    pop = json.loads((src / "Petersen-population-8.json").read_text())
    replayed = []
    for r in cats[8]:
        k = r["quotient_order"]
        mask = r["quotient_mask"]
        m = r["class_sizes"]
        G = groups[k, mask]
        stabilizer = sum(all(m[i] == m[p[i]] for i in range(k)) for p in G)
        aut = stabilizer * prod(factorial(a) for a in m)
        assert factorial(8) % aut == 0
        replayed.append((factorial(8) // aut) * Pemb[k, mask])
    assert replayed == pop["integer_weights"] and sum(replayed) == 10**8
    out = dict(
        status="complete_colored_catalog_counts_and_full_P_population_independently_replayed",
        colored_host_counts=predicted,
        all_triangle_free_quotient_representatives_pairwise_nonisomorphic=True,
        every_catalog_count_checked_by_Burnside_cycle_generating_functions=True,
        P_population_all1491_weights_checked_by_induced_embeddings_and_size_stabilizers=True,
        P_supported_hosts=sum(v > 0 for v in replayed),
        P_population_denominator=10**8,
        seconds=time.monotonic() - start,
    )
    out.update(
        actual_red_blue_masks_and_weight_orbit_uniqueness_checked_for_all_catalogs=True,
        quotient_labelled_orbit_sums=orbit_sums,
        independently_enumerated_labelled_triangle_free_graph_counts=labelled,
        quotient_completeness_self_contained_labelled_count_audit=args.verify_labelled_counts,
    )
    (Path.cwd() / (args.name + ".json")).write_text(json.dumps(out, indent=2))
    print(json.dumps(out), flush=True)


if __name__ == "__main__":
    main()
