"""Verify the exact certificate's flag grammar without trusting author data.

A classifier of extensions is a legitimate flag function if it is invariant
under relabelling unlabelled vertices. It may merge classes; positivity remains.
This validates that invariant directly on every possible extension bit pattern.
The C++ verifier checks Gram factors, all host inequalities, target and coverage.
"""

import hashlib
from itertools import combinations, permutations
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent


def check(path):
    blocks = 0
    patterns = 0
    with Path(path).open(encoding="ascii") as stream:
        count, hosts, denominator = map(int, next(stream).split())
        if (count, hosts) != (169, 12346) or not 0 < denominator <= 2**42:
            raise ValueError("invalid certificate dimensions")
        seen = set()
        for _ in range(count):
            n, k, root, size, width = map(int, next(stream).split())
            u = n - k
            if (
                k not in (2, 4, 6)
                or 2 * n - k != 8
                or not 0 <= root < 2 ** (k * (k - 1) // 2)
                or not 1 <= size <= 272
            ):
                raise ValueError("invalid rooted graph grammar")
            if width != k * u + u * (u - 1) // 2 or (k, root) in seen:
                raise ValueError("invalid or duplicate type")
            seen.add((k, root))
            codes = list(map(int, next(stream).split()))
            if len(codes) != 2**width or any(not 0 <= i < size for i in codes):
                raise ValueError("invalid extension classifier")
            edges = [(a, b) for a, b in combinations(range(n), 2) if b >= k]
            edge_index = {edge: i for i, edge in enumerate(edges)}
            for extra in permutations(range(k, n)):
                order = tuple(range(k)) + extra
                remap = [
                    edge_index[tuple(sorted((order[a], order[b])))] for a, b in edges
                ]
                for value in range(2**width):
                    permuted = sum(
                        ((value >> old) & 1) << new for new, old in enumerate(remap)
                    )
                    if codes[value] != codes[permuted]:
                        raise ValueError(
                            "classifier depends on an unlabelled vertex name"
                        )
            patterns += len(codes)
            blocks += 1
            for _ in range(size):
                row = list(map(int, next(stream).split()))
                if len(row) != size or any(abs(v) >= 2**47 for v in row):
                    raise ValueError("invalid integer Gram factor")
        for _ in range(hosts):
            graph, target = map(int, next(stream).split())
            if not 0 <= graph < 2**28 or not 0 <= target <= 28:
                raise ValueError("invalid finite host")
        if stream.read().strip():
            raise ValueError("trailing certificate data")
    return {
        "status": "verified_exact_flag_grammar",
        "typed_blocks": blocks,
        "extension_patterns_checked": patterns,
        "integer_gram_denominator": denominator,
        "certificate_sha256": hashlib.sha256(Path(path).read_bytes()).hexdigest(),
    }
