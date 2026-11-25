"""Replay the mathematical certificates in a fresh temporary directory."""

import argparse
import hashlib
import importlib.util
import json
import os
import shlex
import shutil
import subprocess
import sys
import tempfile
from fractions import Fraction
from itertools import combinations, permutations
from pathlib import Path

HERE = Path(__file__).resolve().parent
LOWER = Fraction(480, 11111)
UPPER = Fraction(233973195106737927117969229, 5415987671873538702683668480)


def check(condition, message):
    if not condition:
        raise ValueError(message)


def read(path):
    return json.loads(path.read_text(encoding="utf-8"))


def execute(command, cwd):
    subprocess.run([str(part) for part in command], cwd=cwd, check=True)


def verify_files():
    for line in (HERE / "SHA256SUMS.txt").read_text(encoding="utf-8").splitlines():
        digest, name = line.split("  ", 1)
        check(
            hashlib.sha256((HERE / name).read_bytes()).hexdigest() == digest,
            "SHA-256 mismatch for " + name,
        )


def elementary():
    # Isomorphism by explicit permutation, independent of the degree/arm
    # classifier used by the coloured certificate.
    f_edges = ((0, 1), (0, 2), (0, 5), (1, 3), (4, 5))
    orbit = {
        frozenset(tuple(sorted((p[i], p[j]))) for i, j in f_edges)
        for p in permutations(range(6))
    }
    pairs = list(combinations(range(5), 2))
    centers = [0] * 10
    for vertices in combinations(range(10), 6):
        edges = frozenset(
            (i, j)
            for i, j in combinations(range(6), 2)
            if set(pairs[vertices[i]]).isdisjoint(pairs[vertices[j]])
        )
        if edges in orbit:
            (center,) = [v for v in range(6) if sum(v in e for e in edges) == 3]
            centers[vertices[center]] += 1
    check(centers == [6] * 10, "Petersen copy count")
    check(Fraction(720 * sum(centers), 10**6 - 10) == LOWER, "construction limit")
    # Table 1 in the manuscript, using two-digit Kneser labels. Check
    # both the profile and the claimed replacement by graph isomorphism.
    roots = ("01", "02", "12", "03", "13", "24")
    profile_rows = [
        ("04", (2, 4), 3),
        ("14", (1, 3), 4),
        ("23", (0,), 2),
        ("34", (0, 1, 2), 5),
    ]
    for label, neighbors, omit in profile_rows:
        actual = tuple(i for i, r in enumerate(roots) if set(label).isdisjoint(r))
        check(actual == neighbors, "root neighbourhood for " + label)
        alt = tuple(label if i == omit else r for i, r in enumerate(roots))
        edges = frozenset(
            (i, j)
            for i, j in combinations(range(6), 2)
            if set(alt[i]).isdisjoint(alt[j])
        )
        check(edges in orbit, "replacement witness for " + label)
    gap = UPPER - LOWER
    table_lower = Fraction("0.0293174175")
    table_upper = Fraction("0.0432004682")
    check(0 < gap < Fraction(34704, 10**12) and UPPER < table_upper, "order of bounds")
    return dict(
        petersen_copies=sum(centers),
        copies_by_center=centers,
        construction_limit=str(LOWER),
        upper_endpoint_used_for_comparison=str(UPPER),
        interval_width=str(gap),
        replacement_examples_checked=len(profile_rows),
        table_lower_relative_gain_percent=str(100 * (LOWER / table_lower - 1)),
        table_numerical_upper_minus_lower=str(table_upper - LOWER),
    )


def restricted(work):
    src = HERE / "restricted"
    data = src / "data"

    def run(script, *arguments):
        execute([sys.executable, "-I", "-B", src / script, *arguments], work)

    run("verify_F26_elementary_portable.py")
    run(
        "audit_colored_catalog_burnside.py",
        "--data-root",
        data,
        "--verify-labelled-counts",
        "--name",
        "catalog",
    )
    run(
        "verify_colored_full_census_portable.py",
        "--data-root",
        data,
        "--name",
        "census",
    )
    run("verify_colored_exact_certificate.py", "--data-root", data, "--name", "gram")
    run(
        "verify_colored_equality_structure.py",
        "--data-root",
        data,
        "--name",
        "equality",
    )
    cat, census, gram, eq = [
        read(work / (name + ".json"))
        for name in ("catalog", "census", "gram", "equality")
    ]
    check(cat["quotient_completeness_self_contained_labelled_count_audit"], "catalogue")
    check(census["every_original_flag_product_coefficient_checked"], "product census")
    check(
        gram["exact_nonpositive_host_inequalities"] == census["colored_hosts"] == 1491,
        "host coverage",
    )
    check(eq["F_root_type"] == 142 and eq["kernel_rank_exact"] == 27, "equality block")
    check(
        eq["triangle_free_completion_F_count_histogram"]
        == {"26": 1, "34": 1, "38": 2, "44": 2, "60": 1},
        "completion census",
    )
    return dict(
        classes_by_order=cat["colored_host_counts"],
        gram_blocks=gram["Gram_blocks"],
        total_gram_dimension=gram["total_kernel_Gram_dimensions"],
        coefficients=census["exact_nonzero_coefficients"],
        packed_columns=census["packed_columns"],
        tight_host_rows=gram["exact_zero_host_inequalities"],
        petersen_supported_rows=cat["P_supported_hosts"],
        equality_f_root=eq["F_root_type"],
        equality_flags=eq["flag_dimension"],
        kernel_column_rank=eq["kernel_rank_exact"],
        completions=eq["triangle_free_completion_F_count_histogram"],
        remaining_edge_mask=eq["unique_maximizing_completion_mask"],
    )


def compile_programs(src, work, chosen):
    compiler = chosen or os.environ.get("CXX")
    if not compiler:
        compiler = next(
            (name for name in ("c++", "g++", "cl", "clang++") if shutil.which(name)), ""
        )
    command = (
        [compiler]
        if Path(compiler).is_file()
        else shlex.split(compiler, posix=os.name != "nt")
    )
    if not command or not shutil.which(command[0]):
        raise SystemExit(
            "A C++20 compiler is required for the unrestricted bound. On Windows open the x64 Native Tools Command Prompt for VS 2022. See README.txt."
        )
    is_msvc = Path(command[0]).stem.lower() == "cl"
    programs = {}
    for name in ("eight", "census", "nine"):
        binary = work / (name + (".exe" if os.name == "nt" else ""))
        arguments = (
            [
                "/nologo",
                "/O2",
                "/EHsc",
                "/std:c++20",
                src / (name + ".cpp"),
                "/Fe:" + str(binary),
            ]
            if is_msvc
            else ["-O3", "-std=c++20", src / (name + ".cpp"), "-o", binary]
        )
        execute(command + arguments, work)
        programs[name] = binary
    return programs


def unrestricted(work, compiler):
    check(
        sys.byteorder == "little",
        "The order-nine binary tables require a little-endian system.",
    )
    src = HERE / "unrestricted"
    programs = compile_programs(src, work, compiler)
    spec = importlib.util.spec_from_file_location("flag_grammar", src / "grammar.py")
    grammar = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(grammar)
    geometry = grammar.check(src / "eight.txt")
    with (src / "one.txt").open(encoding="ascii") as stream:
        dim, denominator = map(int, next(stream).split())
        classes = list(map(int, next(stream).split()))
        factors = [int(x) for row in stream for x in row.split()]
    check(
        dim == 90 and denominator == geometry["integer_gram_denominator"] == 2**42,
        "Gram dimensions",
    )
    check(len(factors) == dim**2 and max(map(abs, factors)) < 2**47, "root-one factor")
    check(
        len(classes) == 1024 and min(classes) >= 0 and max(classes) < dim,
        "root-one classifier",
    )
    pairs = list(combinations(range(5), 2))
    for perm in permutations(range(1, 5)):
        p = (0, *perm)
        indices = [pairs.index(tuple(sorted((p[i], p[j])))) for i, j in pairs]
        for m in range(1024):
            image = sum(((m >> i) & 1) << j for j, i in enumerate(indices))
            check(classes[m] == classes[image], "root-one flag relabelling")
    execute(
        [
            programs["eight"],
            src / "eight.txt",
            "33701886732571024057023867254",
            work / "eight.json",
            work / "eight_progress",
            work / "eight_rows.json",
        ],
        work,
    )
    rows = read(work / "eight_rows.json")
    check(
        len(rows["host_numerators"]) == 12346
        and int(rows["denominator"]) == 40320 * denominator**2,
        "eight-host scores",
    )
    (work / "values").write_text(
        "12346 "
        + rows["denominator"]
        + "\n"
        + "\n".join(rows["host_numerators"])
        + "\n",
        encoding="ascii",
    )
    with (src / "eight.txt").open(encoding="ascii") as stream:
        nb, nh, _ = map(int, next(stream).split())
        for _ in range(nb):
            dimension = int(next(stream).split()[3])
            for _ in range(dimension + 1):
                next(stream)
        masks = [next(stream).split()[0] for _ in range(nh)]
    (work / "masks").write_text("\n".join(masks) + "\n", encoding="ascii")
    prefix = src / "hosts/nine"
    for suffix, size in ((".u64", 8), (".sub8.u16", 18), (".target.u8", 1)):
        check(
            Path(str(prefix) + suffix).stat().st_size == 274668 * size,
            "order-nine table length",
        )
    execute([programs["census"], prefix, work / "masks", work / "orbits.json"], work)
    orbits = read(work / "orbits.json")
    check(orbits["labelled_graphs_covered"] == 1 << 36, "unlabelled orbit coverage")
    execute(
        [
            programs["nine"],
            src / "one.txt",
            work / "values",
            prefix,
            "independent",
            work / "nine.json",
        ],
        work,
    )
    nine = read(work / "nine.json")
    value = Fraction(int(nine["numerator"]), int(nine["denominator"]))
    check(
        value == UPPER and int(nine["denominator"]) == 362880 * denominator**2,
        "global maximum",
    )
    return dict(
        eight_graphs=12346,
        nine_graphs=orbits["host_classes"],
        labelled_eight=read(work / "eight.json")["labelled_graphs_covered"],
        labelled_nine=orbits["labelled_graphs_covered"],
        even_blocks=geometry["typed_blocks"],
        odd_block_dimension=dim,
        denominator_of_factors=denominator,
        bound=str(value),
    )


def main():
    if sys.flags.optimize:
        raise SystemExit("Use ordinary Python. Optimization disables proof assertions.")
    sys.dont_write_bytecode = True
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--part",
        choices=("all", "restricted", "unrestricted", "elementary"),
        default="all",
    )
    parser.add_argument(
        "--compiler", help="C++20 compiler executable for the unrestricted certificate"
    )
    parser.add_argument(
        "--output", type=Path, help="write the resulting JSON receipt to this file"
    )
    args = parser.parse_args()
    verify_files()
    report = {"elementary": elementary()}
    with tempfile.TemporaryDirectory(prefix="f26-proof-") as tmp:
        work = Path(tmp)
        if args.part in ("all", "restricted"):
            report["restricted"] = restricted(work)
        if args.part in ("all", "unrestricted"):
            report["unrestricted"] = unrestricted(work, args.compiler)
    report["verified"] = args.part
    output = json.dumps(report, indent=2, ensure_ascii=True) + "\n"
    if args.output:
        args.output.write_text(output, encoding="utf-8")
    print(output)


if __name__ == "__main__":
    main()
