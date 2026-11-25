"""Elementary exact prerequisites and sharp Petersen example; standard Python."""

from itertools import combinations
from fractions import Fraction as F
from pathlib import Path
import json


def main():
    def graph(n, edges):
        A = [[0] * n for _ in range(n)]
        for i, j in edges:
            A[i][j] = A[j][i] = 1
        return A

    def prime(A):
        n = len(A)
        for k in range(2, n):
            for S in combinations(range(n), k):
                if all(len({A[v][w] for w in S}) == 1 for v in range(n) if v not in S):
                    return False
        return True

    T = graph(6, [(0, 1), (0, 2), (0, 5), (1, 3), (4, 5)])
    bull = graph(5, [(0, 1), (1, 2), (0, 2), (0, 3), (1, 4)])
    house = graph(5, [(i, j) for i, j in combinations(range(5), 2) if j != i + 1])
    assert all(prime(A) for A in (T, bull, house))
    alpha = max(
        len(S)
        for k in range(7)
        for S in combinations(range(6), k)
        if all(not T[i][j] for i, j in combinations(S, 2))
    )
    assert alpha == 3
    pairs = list(combinations(range(5), 2))
    P = [[int(set(a).isdisjoint(b)) for b in pairs] for a in pairs]
    assert all(sum(row) == 3 for row in P) and not any(
        P[i][j] and P[i][k] and P[j][k] for i, j, k in combinations(range(10), 3)
    )
    # Classify connected six-vertex trees by their degree profile and center arms.
    copies = 0
    for S in combinations(range(10), 6):
        degrees = {i: sum(P[i][j] for j in S) for i in S}
        if sorted(degrees.values()) != [1, 1, 1, 2, 2, 3]:
            continue
        center = next(i for i in S if degrees[i] == 3)
        if sorted(degrees[j] for j in S if P[center][j]) != [1, 2, 2]:
            continue
        reached = {center}
        todo = [center]
        while todo:
            v = todo.pop()
            for w in S:
                if P[v][w] and w not in reached:
                    reached.add(w)
                    todo.append(w)
        if len(reached) == 6:
            copies += 1
    assert copies == 60
    q = F(720 * copies, 10**6)
    b = F(10, 10**6)
    L = F(480, 11111)
    assert q == F(27, 625) and q / (1 - b) == L
    qt = F(0)
    for t in range(1, 9):
        qt = q + b * qt
        n = 10**t
        assert qt == L * (1 - F(1, n**5))
        count = qt * n**6 / 720
        assert count == F(2 * (n**6 - n), 33333) and count.denominator == 1
    out = dict(
        status="F26_primality_alpha_and_Petersen_sharpness_exactly_verified",
        F26_bull_house_prime=True,
        F26_independence_number=alpha,
        Petersen_F_six_subsets=copies,
        sharp_corrected_capacity=str(L),
        finite_equality_orders_checked=8,
        general_problem_closed=False,
    )
    Path("elementary-verification.json").write_text(json.dumps(out, indent=2))
    print(json.dumps(out))


if __name__ == "__main__":
    main()
