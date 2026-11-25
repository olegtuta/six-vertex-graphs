"""Reproduce the counting and coefficient certificates in the second edition."""
from collections import Counter
from fractions import Fraction as Q
from itertools import combinations, combinations_with_replacement, permutations, product
from math import comb, factorial, prod
from pathlib import Path
import json


def edges(words):
    return {tuple(map(int, s)) for s in words.split()}


def matrix(n, pairs, loops=()):
    return [[int((min(i,j), max(i,j)) in pairs or i == j and i in loops)
             for j in range(n)] for i in range(n)]


def orbit(target):
    index = {e: k for k, e in enumerate(combinations(range(6), 2))}
    return {sum(1 << index[tuple(sorted((p[a],p[b])))] for a,b in target)
            for p in permutations(range(6))}


def induced_code(adj, vertices):
    return sum(adj[i][j] << k for k,(i,j) in enumerate(combinations(vertices, 2)))


def compatible_maps(target, adj, wild=()):
    """Generate assignments of one labelling, allowing both bits within wild parts."""
    def visit(prefix):
        v = len(prefix)
        if v == 6:
            yield tuple(prefix)
            return
        for i in range(len(adj)):
            for u,j in enumerate(prefix):
                if i == j and i in wild:
                    continue
                if adj[j][i] != int((u,v) in target):
                    break
            else:
                yield from visit(prefix + [i])
    yield from visit([])


def perturbation(case):
    target = edges(case["edges"])
    rec = set(case["recursive"])
    r = len(rec)
    adj = matrix(case["parts"], edges(case["between"]), case["cliques"])
    n = len(adj)
    assert n == 2*r
    labels = len(orbit(target))
    x = Q(case["x"])
    y = 1/Q(r) - x
    coeff = 0
    for part in compatible_maps(target, adj, rec):
        if len(set(part)) == 1:
            continue
        inside = [i for i in part if i in rec]
        assert len(inside) == len(set(inside)) == 2
        coeff += labels
    assert coeff == {2:180,3:2160}[r]
    density = lambda t: coeff*t**2*(1/Q(r)-t)**4 / (1-r*t**6)
    old = Q(case["proposed"])
    assert density(1/Q(3*r)) == old
    # Flatten two levels, retain clique diagonals, freeze remaining recursion at zero.
    parts = [(i,j) for i in range(n) for j in (range(n) if i in rec else [None])]
    weight = [x if i in rec else y for i in range(n)]
    mass = [weight[i]*(weight[j] if j is not None else 1) for i,j in parts]
    assert sum(mass) == 1
    flat = [[adj[i][k] if i != k or i not in rec else adj[j][l]
             for k,l in parts] for i,j in parts]
    finite = labels*sum(prod(mass[i] for i in p) for p in compatible_maps(target, flat))
    assert finite == density(x)*(1-(r*x**6)**2) > old
    return dict(cross_coefficient=coeff, classes=len(parts), x=str(x),
                finite=str(finite), excess_over_proposal=str(finite-old), recursive=str(density(x)))


def cycle_and_isolate(case):
    target = edges(case["edges"])
    base = matrix(6,target)
    labels = orbit(target)
    parts = list(product(range(6),repeat=2))
    adj = [[base[i][j] if i != j else base[a][b] for j,b in parts] for i,a in parts]
    count = sum(induced_code(adj, s) in labels for s in combinations(range(36),6))
    assert count == 46782
    # Check every pair has an outside distinguisher, including adjacent pairs.
    assert all(any(base[i][k] != base[j][k] for k in range(6) if k not in (i,j))
               for i,j in combinations(range(6),2))
    density = Q(720*count,36**6)
    c = Q(120,5**5-5)
    d = 60*c/(5**6-5)
    assert (c,d) == (Q(1,26),Q(3,20306))
    a = Q(167216,10**6)
    b = 1-a
    improved = d*(a**6+b**6)+6*c*(a**5*b+a*b**5)
    assert improved > density == Q(12995,839808) > Q(case["proposed"])
    return dict(six_sets_checked=comb(36,6), copies=count, blowup=str(density),
                excess=str(density-Q(case["proposed"])), two_components=str(improved))


def packed(values):
    # Little-endian base 2^32 digits. All individual coefficients are nonnegative.
    return int.from_bytes(b''.join(v.to_bytes(4,'little') for v in values),'little')


def rook(case):
    motif = matrix(6, edges(case["edges"]))
    for size in range(2,6):
        for group in combinations(range(6), size):
            assert any(len({motif[v][i] for i in group}) == 2
                       for v in range(6) if v not in group)
    cells = list(product(range(3),repeat=2))
    adj = [[int(a!=c and b!=d) for c,d in cells] for a,b in cells]
    codes = orbit(edges(case["edges"]))
    sixes = list(combinations(range(9),6))
    copies = [s for s in sixes if induced_code(adj,s) in codes]
    value = Q(case["value"])
    assert len(copies) == 36 and Q(720*len(copies),9**6-9) == value

    # The column for ranks S is the polynomial product of the six prefix sums.
    terms = list(combinations_with_replacement(range(9),6))
    columns = {}
    for ranks in sixes:
        poly = Counter({():1})
        for rank in ranks:
            extended = Counter()
            for monomial, coeff in poly.items():
                for j in range(rank+1):
                    extended[tuple(sorted((*monomial,j)))] += coeff
            poly = extended
        columns[ranks] = packed(poly[t] for t in terms)
    # Coefficients of (sum x)^6 - sum x^6 in the same ordered coordinates.
    def coefficient(term):
        multi = 720 // prod(factorial(v) for v in Counter(term).values())
        return multi*(prod(9-j for j in term)-(9-max(term)))
    totals = [coefficient(t) for t in terms]
    guard = 1 << 31
    flags = packed(guard for t in terms)
    initial = packed(guard+v for v in totals)
    assert max(totals) < guard and 14762*len(copies)*720 < guard
    # Each column coefficient is <= 6!. The guard digit prevents borrows.
    # Reading its sign bit therefore tests each coefficient separately.
    positive = [terms.index((i,)*6) for i in range(1,8)]
    positive.append(terms.index((0,0,0,0,8,8)))
    perms = list(permutations(range(3)))
    group = {tuple(3*p[a]+q[b] for a,b in cells) for p,q in product(perms,repeat=2)}
    group |= {tuple(3*p[b]+q[a] for a,b in cells) for p,q in product(perms,repeat=2)}
    assert len(group) == 72
    assert all(adj[i][j] == adj[p[i]][p[j]] for p in group for i,j in combinations(range(9),2))
    stabilizer = [p for p in group if p[0] == 0]
    seen = set()
    checked = 0
    for tail in permutations(range(1,9)):
        order = (0,*tail)
        if order in seen:
            continue
        seen.update(tuple(p[i] for i in order) for p in stabilizer)
        rank = {v:i for i,v in enumerate(order)}
        subtotal = sum(columns[tuple(sorted(rank[v] for v in s))] for s in copies)
        digits = initial - 14762*subtotal
        assert digits & flags == flags
        assert all(((digits >> (32*j)) & (guard-1)) > 0 for j in positive)
        checked += 1
    assert len(seen) == factorial(8) and checked*len(group) == factorial(9)
    return dict(copies=len(copies), value=str(value), orders=checked,
                coefficients_per_order=len(terms), total_coefficients=checked*len(terms))


def main():
    if not __debug__:
        raise SystemExit("Use Python with assertions enabled.")
    data = json.loads(Path(__file__).with_name('cases.json').read_text(encoding='utf-8'))
    results = {k:perturbation(data[k]) for k in ('F22','F34','F62')}
    results['F29'] = cycle_and_isolate(data['F29'])
    results['F74'] = rook(data['F74'])
    print(json.dumps(results, indent=2,sort_keys=True))


if __name__ == '__main__':
    main()
