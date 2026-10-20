"""Classify finite three-house deductions; Python 3, standard library only.

This is independently executable research code. It neither calls the production
solver nor claims that an obstruction is a historically new Sudoku technique.
Bit digit-1 records a candidate. Domains, including any earlier eliminations,
are premises. No uniqueness premise is used.
"""

import argparse
from collections import Counter
from functools import lru_cache
import hashlib
import itertools
import json
from pathlib import Path

if not __debug__:
    raise SystemExit("Run without -O: this research classifier uses assertions.")

FULL = 511
ROWS = [tuple(9 * r + c for c in range(9)) for r in range(9)]
COLUMNS = [tuple(9 * r + c for r in range(9)) for c in range(9)]
BOXES = [
    tuple(9 * (3 * br + r) + 3 * bc + c for r in range(3) for c in range(3))
    for br in range(3) for bc in range(3)
]
HOUSES = ROWS + COLUMNS + BOXES
HERE = Path(__file__).resolve().parent


def bits(mask):
    while mask:
        bit = mask & -mask
        yield bit
        mask -= bit


def count(mask):
    return bin(mask).count("1")


def label(cell):
    return "r%dc%d" % (cell // 9 + 1, cell % 9 + 1)


def digits(mask):
    return [d for d in range(1, 10) if mask & (1 << (d - 1))]


@lru_cache(maxsize=32768)
def matching_support(domains, allowed=FULL):
    """Return every matching-supported edge, or None if no covering matching.

    The recursion enumerates a subset DP; the second traversal unions edges of
    all accepting paths. It does not decide the complete Sudoku.
    """
    order = sorted(range(len(domains)), key=lambda i: count(domains[i] & allowed))
    ordered = tuple(domains[i] & allowed for i in order)

    @lru_cache(None)
    def good(i, used):
        if i == len(ordered):
            return True
        return any(good(i + 1, used | bit) for bit in bits(ordered[i] & ~used))

    if not good(0, 0):
        return None
    support = [0] * len(domains)
    visited = set()

    def visit(i, used):
        if i == len(ordered) or (i, used) in visited:
            return
        visited.add((i, used))
        for bit in bits(ordered[i] & ~used):
            if good(i + 1, used | bit):
                support[order[i]] |= bit
                visit(i + 1, used | bit)

    visit(0, 0)
    return tuple(support)


def cluster(center):
    r, c = divmod(center, 9)
    row, column, box = set(ROWS[r]), set(COLUMNS[c]), set(BOXES[3 * (r // 3) + c // 3])
    atoms = [
        (center,), tuple(sorted((row & box) - {center})),
        tuple(sorted((column & box) - {center})), tuple(sorted(row - box)),
        tuple(sorted(column - box)), tuple(sorted(box - row - column)),
    ]
    assert list(map(len, atoms)) == [1, 2, 2, 6, 6, 4]
    return atoms, [tuple(sorted(row)), tuple(sorted(column)), tuple(sorted(box))]


def boundary_signatures():
    """All 3,780 assignments of digit SETS to the three boundary atoms."""
    for center in bits(FULL):
        for row_pair in itertools.combinations(tuple(bits(FULL ^ center)), 2):
            row = row_pair[0] | row_pair[1]
            for column_pair in itertools.combinations(tuple(bits(FULL ^ (center | row))), 2):
                column = column_pair[0] | column_pair[1]
                yield center, row, column


SIGNATURES = tuple(boundary_signatures())
assert len(SIGNATURES) == 3780 and len(set(SIGNATURES)) == 3780


def atom_scan(domains, center=0, first_only=False):
    atoms, _ = cluster(center)
    atom_domains = [tuple(domains[c] for c in atom) for atom in atoms]
    supported = {c: 0 for atom in atoms for c in atom}
    signature_count = 0
    for digit, row, column in SIGNATURES:
        allowed = (digit, row, column, FULL ^ (digit | row),
                   FULL ^ (digit | column), FULL ^ (digit | row | column))
        matchings = []
        for ds, mask in zip(atom_domains, allowed):
            support = matching_support(ds, mask)
            if support is None:
                break
            matchings.append(support)
        if len(matchings) != 6:
            continue
        if first_only:
            return True
        signature_count += 1
        for atom, support in zip(atoms, matchings):
            for cell, mask in zip(atom, support):
                supported[cell] |= mask
    if first_only:
        return False
    return supported, signature_count


def pair_feasible(domains, first, second):
    shared = tuple(sorted(set(first) & set(second)))
    private_a = tuple(sorted(set(first) - set(second)))
    private_b = tuple(sorted(set(second) - set(first)))
    for choice in itertools.combinations(tuple(bits(FULL)), len(shared)):
        mask = sum(choice)
        if matching_support(tuple(domains[c] for c in shared), mask) is None:
            continue
        if matching_support(tuple(domains[c] for c in private_a), FULL ^ mask) is None:
            continue
        if matching_support(tuple(domains[c] for c in private_b), FULL ^ mask) is not None:
            return True
    return False


def assume(domains, target):
    result = list(domains)
    result[target["cell"]] = 1 << (target["digit"] - 1)
    return result


def conditional_hall_closure(domains, houses, locked=False):
    """GAC of the selected houses, optionally plus intersection propagation.

    Failure means ordinary Hall propagation under the assumed target suffices;
    it does not imply that a static ALS/AIC detector recognizes the deduction.
    """
    domains = list(domains)
    rounds = 0
    while True:
        changed = False
        rounds += 1
        for house in houses:
            support = matching_support(tuple(domains[c] for c in house))
            if support is None:
                return False, rounds
            for cell, mask in zip(house, support):
                if domains[cell] != mask:
                    domains[cell] = mask
                    changed = True
        if locked:
            for first in houses:
                for second in houses:
                    if first == second:
                        continue
                    intersection = set(first) & set(second)
                    if not intersection:
                        continue
                    outside_first = set(first) - intersection
                    outside_second = set(second) - intersection
                    for digit in bits(FULL):
                        if (any(domains[c] & digit for c in intersection)
                                and not any(domains[c] & digit for c in outside_first)):
                            for cell in outside_second:
                                if domains[cell] & digit:
                                    domains[cell] &= ~digit
                                    changed = True
        if not changed:
            return True, rounds


def classify(domains, target, center=0):
    assert len(domains) == 81 and all(0 < mask <= FULL for mask in domains)
    assert domains[target["cell"]] & (1 << (target["digit"] - 1))
    forced = assume(domains, target)
    _, houses = cluster(center)
    assert not atom_scan(forced, center, first_only=True), "Not a three-house deduction"
    single_refutations = [name for name, house in zip(("R", "C", "B"), houses)
                          if matching_support(tuple(forced[c] for c in house)) is None]
    pair_refutations = [name for name, a, b in (("RB", 0, 2), ("CB", 1, 2), ("RC", 0, 1))
                        if not pair_feasible(forced, houses[a], houses[b])]
    hall_ok, hall_rounds = conditional_hall_closure(forced, houses)
    locked_ok, locked_rounds = conditional_hall_closure(forced, houses, locked=True)
    category = ("single_house" if single_refutations else
                "pair_reducible" if pair_refutations else "three_house_irreducible")
    return {
        "category": category, "single_refutations": single_refutations,
        "pair_refutations": pair_refutations,
        "conditional_hall_refutes": not hall_ok, "conditional_hall_rounds": hall_rounds,
        "conditional_hall_locked_refutes": not locked_ok,
        "conditional_hall_locked_rounds": locked_rounds,
    }


def all_native_als(domains):
    """Enumerate every 1..8-cell ALS in a native house, without engine bounds."""
    found = {}
    for house in HOUSES:
        for size in range(1, 9):
            for cells in itertools.combinations(house, size):
                union = 0
                for cell in cells:
                    union |= domains[cell]
                if count(union) == size + 1:
                    found[tuple(sorted(cells))] = union
    return sorted(found.items())


def sees(first, second):
    if first == second:
        return False
    ar, ac = divmod(first, 9)
    br, bc = divmod(second, 9)
    return ar == br or ac == bc or (ar // 3 == br // 3 and ac // 3 == bc // 3)


def als_xz_explanation(domains, target, als_sets):
    """Find the smallest maximum-ALS-size standard, disjoint ALS-XZ witness.

    This is a complete search of this ONE known rule over native-house ALSes;
    it is not a complete search of ALS chains, overlapping ALSes, or rank logic.
    """
    z = 1 << (target["digit"] - 1)
    eligible = [(cells, union) for cells, union in als_sets
                if union & z and target["cell"] not in cells
                and all(sees(target["cell"], c) for c in cells if domains[c] & z)]
    best = None
    for (first, first_union), (second, second_union) in itertools.combinations(eligible, 2):
        if set(first) & set(second):
            continue
        for x in bits((first_union & second_union) & ~z):
            if all(sees(a, b) for a in first if domains[a] & x
                   for b in second if domains[b] & x):
                key = (max(len(first), len(second)), len(first) + len(second),
                       first, second, x)
                if best is None or key < best:
                    best = key
    if best is None:
        return None
    maximum, total, first, second, x = best
    # Explicitly recheck the theorem premises independently of search filters.
    for cells in (first, second):
        assert any(set(cells) <= set(house) for house in HOUSES)
        assert len(set(cells)) == len(cells)
        assert count(sum_union(domains[c] for c in cells)) == len(cells) + 1
        assert all(sees(target["cell"], c) for c in cells if domains[c] & z)
    assert not (set(first) & set(second))
    assert all(sees(a, b) for a in first if domains[a] & x for b in second if domains[b] & x)
    return {"first": list(first), "second": list(second), "rcc_digit": digits(x)[0],
            "elimination_digit": target["digit"], "maximum_als_size": maximum,
            "total_cells": total}


def sum_union(masks):
    union = 0
    for mask in masks:
        union |= mask
    return union


def canonical_fingerprint(domains, target, center=0):
    """Exact local equivalence, not equivalence of entire Sudoku boards.

    Quotient digit renaming, row/column transposition, and permutations of cells
    WITHIN a membership atom. The target cell and target digit are distinguished.
    Arbitrary within-atom permutations need not extend to a native-grid symmetry.
    """
    atoms, _ = cluster(center)
    target_atom = next(i for i, atom in enumerate(atoms) if target["cell"] in atom)
    ordinary = [tuple(domains[c] for c in atom if c != target["cell"]) for atom in atoms]
    target_mask = domains[target["cell"]]
    distinct = sorted(set([target_mask] + [m for atom in ordinary for m in atom]))
    original_bits = [1 << d for d in range(9) if d != target["digit"] - 1]
    fixed_bit = 1 << (target["digit"] - 1)
    best = None
    for destinations in itertools.permutations([1 << d for d in range(1, 9)]):
        mapping = dict(zip(original_bits, destinations))
        mapping[fixed_bit] = 1
        renamed = {m: sum(mapping[bit] for bit in bits(m)) for m in distinct}
        groups = [tuple(sorted(renamed[m] for m in atom)) for atom in ordinary]
        for atom_order in ((0, 1, 2, 3, 4, 5), (0, 2, 1, 4, 3, 5)):
            form = (atom_order.index(target_atom), renamed[target_mask],
                    tuple(groups[i] for i in atom_order))
            if best is None or form < best:
                best = form
    encoded = json.dumps(best, separators=(",", ":"))
    return hashlib.sha256(encoded.encode("ascii")).hexdigest()


def minimize_premises(domains, target, center=0, reverse=False):
    """Greedy 1-minimal zero-edge premises, NOT globally minimum cardinality."""
    atoms, _ = cluster(center)
    region = sorted({c for atom in atoms for c in atom}, reverse=reverse)
    result = [FULL] * 81
    for cell in region:
        result[cell] = domains[cell]
    assert not atom_scan(assume(result, target), center, first_only=True)
    checks = 0
    for cell in region:
        old = result[cell]
        result[cell] = FULL
        checks += 1
        if atom_scan(assume(result, target), center, first_only=True):
            result[cell] = old
    for cell in region:
        for bit in sorted(bits(FULL ^ result[cell]), reverse=reverse):
            result[cell] |= bit
            checks += 1
            if atom_scan(assume(result, target), center, first_only=True):
                result[cell] ^= bit
    assert not atom_scan(assume(result, target), center, first_only=True)
    # Monotonicity makes the single pass 1-minimal; independently check it.
    for cell in region:
        for bit in bits(FULL ^ result[cell]):
            relaxed = result[:]
            relaxed[cell] |= bit
            assert atom_scan(assume(relaxed, target), center, first_only=True)
    return {
        "order": "descending" if reverse else "ascending", "checks": checks,
        "restricted_cells": sum(m != FULL for m in result),
        "excluded_candidates": sum(9 - count(m) for m in result),
        "domains": result, "classification": classify(result, target, center),
    }


def hall_witness(domains, cells, allowed):
    for size in range(1, len(cells) + 1):
        for indices in itertools.combinations(range(len(cells)), size):
            union = 0
            for i in indices:
                union |= domains[cells[i]] & allowed
            if count(union) < size:
                return tuple(cells[i] for i in indices), union
    raise AssertionError("Failed matching had no Hall witness")


def rejection_certificate_statistics(domains, target, center=0):
    """Independently recompute a Hall rejection for every boundary signature."""
    forced = assume(domains, target)
    atoms, _ = cluster(center)
    reasons = Counter()
    witnesses = set()
    for digit, row, column in SIGNATURES:
        allowed = (digit, row, column, FULL ^ (digit | row),
                   FULL ^ (digit | column), FULL ^ (digit | row | column))
        for index, (cells, mask) in enumerate(zip(atoms, allowed)):
            if matching_support(tuple(forced[c] for c in cells), mask) is None:
                subset, neighbors = hall_witness(forced, cells, mask)
                assert count(neighbors) < len(subset)
                reasons[str(index)] += 1
                witnesses.add((index, subset, neighbors))
                break
        else:
            raise AssertionError("Certificate had an unrefuted signature")
    assert sum(reasons.values()) == 3780
    return {"covered_signatures": 3780, "first_failing_atom_counts": dict(reasons),
            "distinct_hall_witnesses": len(witnesses)}


def compact_subset_cover(domains, target, center=0):
    """Greedy cover using atom/cell-subset witnesses without fixed digit sets."""
    forced = assume(domains, target)
    atoms, _ = cluster(center)
    signatures = [(d, r, c, FULL ^ (d | r), FULL ^ (d | c), FULL ^ (d | r | c))
                  for d, r, c in SIGNATURES if forced[center] & d]
    candidates = []
    for atom in range(1, 6):
        for local in range(1, 1 << len(atoms[atom])):
            cells = [c for i, c in enumerate(atoms[atom]) if local & (1 << i)]
            neighbors = sum_union(forced[c] for c in cells)
            covered = sum(1 << i for i, signature in enumerate(signatures)
                          if count(neighbors & signature[atom]) < len(cells))
            if covered:
                candidates.append((atom, local, neighbors, covered))
    uncovered = (1 << len(signatures)) - 1
    chosen = []
    while uncovered:
        best = min(candidates, key=lambda x: (-count(x[3] & uncovered), x[0], x[1]))
        assert best[3] & uncovered, "An exact Hall subset cover was incomplete"
        atom, local, neighbors, covered = best
        chosen.append({"atom": atom, "cells": local, "neighbor_union": neighbors,
                       "newly_covered": count(covered & uncovered)})
        uncovered &= ~covered
    return {"signatures": len(signatures), "potential_nonempty_subsets": 147,
            "greedy_blockers": chosen, "minimum_cardinality_established": False}


def permutation_oracle(domains, target, center=0):
    """Independent factorial enumeration joined on actual ordered cell values."""
    forced = assume(domains, target)
    _, houses = cluster(center)
    valid = []
    for house in houses:
        valid.append([p for p in itertools.permutations(range(1, 10))
                      if all(forced[c] & (1 << (digit - 1)) for c, digit in zip(house, p))])
    rb = sorted(set(houses[0]) & set(houses[2]))
    cb = sorted(set(houses[1]) & set(houses[2]))
    indices = [{c: i for i, c in enumerate(house)} for house in houses]

    def project(assignment, house, shared):
        return tuple(assignment[indices[house][c]] for c in shared)

    row = {project(p, 0, rb) for p in valid[0]}
    column = {project(p, 1, cb) for p in valid[1]}
    row_box = [p for p in valid[2] if project(p, 2, rb) in row]
    column_box = [p for p in valid[2] if project(p, 2, cb) in column]
    joint = [p for p in row_box if project(p, 2, cb) in column]
    row_column = ({p[indices[0][center]] for p in valid[0]}
                  & {p[indices[1][center]] for p in valid[1]})
    assert not joint
    category = classify(domains, target, center)
    actual_pairs = [name for name, feasible in (("RB", bool(row_box)),
                                              ("CB", bool(column_box)),
                                              ("RC", bool(row_column))) if not feasible]
    assert actual_pairs == category["pair_refutations"]
    return {"individual_assignments": list(map(len, valid)), "row_box": len(row_box),
            "column_box": len(column_box), "row_column_digits": sorted(row_column),
            "triple_join": len(joint)}


def check_completion(domains, completion):
    board = list(map(int, completion))
    assert len(board) == 81
    assert all(domains[c] & (1 << (d - 1)) for c, d in enumerate(board))
    assert all({board[c] for c in house} == set(range(1, 10)) for house in HOUSES)


def classify_corpus(corpus, canonical=True):
    rows = []
    for state in corpus["states"]:
        domains = state["domains"]
        check_completion(domains, corpus["completion"])
        supported, signatures = atom_scan(domains, state["center"])
        assert signatures == state["compatible_signatures"]
        actual = sorted((cell, digit) for cell, support in supported.items()
                        for digit in digits(domains[cell] & ~support))
        expected = sorted((t["cell"], t["digit"]) for t in state["eliminations"])
        assert actual == expected
        als_sets = all_native_als(domains)
        for target in state["eliminations"]:
            result = {"trial": state["trial"], "cell": target["cell"],
                      "digit": target["digit"], **classify(domains, target, state["center"])}
            result["known_als_xz"] = als_xz_explanation(domains, target, als_sets)
            if canonical and result["category"] == "three_house_irreducible":
                result["canonical_fingerprint"] = canonical_fingerprint(domains, target, state["center"])
            rows.append(result)
    return {
        "schema_version": 1, "states": len(corpus["states"]), "eliminations": len(rows),
        "categories": dict(Counter(row["category"] for row in rows)),
        "conditional_hall_refutable": sum(row["conditional_hall_refutes"] for row in rows),
        "conditional_hall_locked_refutable": sum(row["conditional_hall_locked_refutes"] for row in rows),
        "known_als_xz_explanations": sum(row["known_als_xz"] is not None for row in rows),
        "known_als_xz_minimum_maximum_sizes": dict(Counter(
            str(row["known_als_xz"]["maximum_als_size"]) for row in rows
            if row["known_als_xz"] is not None)),
        "canonical_irreducible_classes": (len({row["canonical_fingerprint"] for row in rows
                                               if "canonical_fingerprint" in row}) if canonical else None),
        "cases": rows,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--no-canonical", action="store_true", help="Skip the exact 8! digit-renaming quotient")
    parser.add_argument("--minimize", action="store_true", help="Also derive two 1-minimal premise sets for trial 847")
    parser.add_argument("--permutation-check", action="store_true", help="Independently enumerate all 9! permutations for the seven irreducible cases")
    parser.add_argument("--output", type=Path, help="Write the deterministic result as JSON")
    options = parser.parse_args()
    corpus = json.loads((HERE / "three_house_classification_corpus.json").read_text())
    result = classify_corpus(corpus, canonical=not options.no_canonical)
    for name, expected in corpus["expected_classification"].items():
        if name == "canonical_irreducible_classes" and options.no_canonical:
            continue
        assert result[name] == expected, (name, result[name], expected)
    fixture = json.loads((HERE / "three-house-847.json").read_text())
    check_completion(fixture["domains"], fixture["completion"])
    result["fixture_847_certificate_statistics"] = rejection_certificate_statistics(fixture["domains"], fixture["target"])
    result["fixture_847_known_als_xz"] = als_xz_explanation(fixture["domains"], fixture["target"], all_native_als(fixture["domains"]))
    assert result["fixture_847_known_als_xz"] == {
        "first": [0, 9, 18, 36, 45, 54, 72], "second": [1, 2], "rcc_digit": 2,
        "elimination_digit": 1, "maximum_als_size": 7, "total_cells": 9,
    }
    result["fixture_847_compact_subset_cover"] = compact_subset_cover(fixture["domains"], fixture["target"])
    if options.permutation_check:
        result["independent_permutation_checks"] = []
        for case in result["cases"]:
            if case["category"] != "three_house_irreducible":
                continue
            state = next(s for s in corpus["states"] if s["trial"] == case["trial"])
            target = {"cell": case["cell"], "digit": case["digit"]}
            result["independent_permutation_checks"].append({"trial": case["trial"], **target,
                **permutation_oracle(state["domains"], target, state["center"])})
    if options.minimize:
        result["fixture_847_minimized"] = [minimize_premises(fixture["domains"], fixture["target"], reverse=reverse)
                                          for reverse in (False, True)]
    encoded = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if options.output:
        options.output.write_text(encoded)
    print(json.dumps({k: v for k, v in result.items() if k not in ("cases", "fixture_847_minimized")}, sort_keys=True))
    if options.minimize:
        print("Greedy minima:", [(r["order"], r["restricted_cells"], r["excluded_candidates"])
                                for r in result["fixture_847_minimized"]])


if __name__ == "__main__":
    main()
