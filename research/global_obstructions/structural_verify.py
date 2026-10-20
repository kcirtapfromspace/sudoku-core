"""Independently check the synthetic 13-cell, four-digit global obstruction.

No graph edges, independence bounds, floating-point results, or colors supplied
by a searcher are trusted. This script uses only Python's standard library and
the earlier independent finite matching baseline in this repository.
"""

import argparse
from fractions import Fraction
import hashlib
import itertools
import json
from pathlib import Path
import sys

if not __debug__:
    raise SystemExit("Run without -O: assertions are part of the research checks.")

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "nonmodular"))
from three_house_classification import HOUSES, FULL, atom_scan, matching_support


def shares_house(first, second):
    return first != second and any(first in h and second in h for h in HOUSES)


def bit_values(mask):
    while mask:
        bit = mask & -mask
        yield bit
        mask -= bit


def check_full_completion(board, domains):
    assert len(board) == 81
    assert all(domains[c] & (1 << (digit - 1)) for c, digit in enumerate(board))
    assert all({board[c] for c in h} == set(range(1, 10)) for h in HOUSES)


def independence_check(vertices):
    counts = {}
    independent_triple = None
    for size in range(1, 5):
        count = 0
        for subset in itertools.combinations(vertices, size):
            if not any(shares_house(a, b) for a, b in itertools.combinations(subset, 2)):
                count += 1
                if size == 3:
                    independent_triple = list(subset)
        counts[str(size)] = count
    assert counts == {"1": 13, "2": 43, "3": 32, "4": 0}
    assert independent_triple is not None
    # Every same-digit set is independent in the reconstructed native graph.
    # Four available digits therefore cover at most 4*3=12 of these 13 cells.
    assert 4 * 3 < len(vertices)
    return {"independence_number": 3, "independent_sets_by_size": counts,
            "four_subsets_checked": 715, "independent_triple": independent_triple}


def check_hand_capacity_proof(vertices):
    """Check the short four-box argument separately from exhaustive subsets."""
    groups = [set(vertices) & set(HOUSES[i]) for i in (18, 20, 21, 23)]
    assert groups == [{9}, {15, 16, 17, 25}, {27, 36, 46, 47}, {34, 43, 51, 53}]
    assert set().union(*groups) == set(vertices)
    assert sum(map(len, groups)) == len(vertices)
    # Four independent cells would need one per box. The singleton first box
    # forces cell 9, then cell 25, then two mutually conflicting row-six sets.
    box_three = {c for c in groups[1] if not shares_house(9, c)}
    box_four = {c for c in groups[2] if not shares_house(9, c)}
    assert box_three == {25} and box_four == {46, 47}
    box_six = {c for c in groups[3] if not shares_house(25, c)}
    assert box_six == {51, 53}
    assert all(shares_house(a, b) for a in box_four for b in box_six)
    return {"covering_boxes": [1, 3, 4, 6], "forced_cells": [9, 25],
            "last_box_choices": [sorted(box_four), sorted(box_six)],
            "last_choices_conflict_in_row": 6}


def exact_linear_witness(domains, palette, roles):
    """Check all 81 cell + 243 house/digit equations using exact Fractions."""
    assert len(roles) == 81 and set(roles) == {0, 1}
    assert all(sum(roles[c] for c in h) == 4 for h in HOUSES)
    palette = set(palette)
    values = [[Fraction(1, 4) if roles[c] and digit in palette else
               Fraction(1, 5) if not roles[c] and digit not in palette else Fraction(0)
               for digit in range(1, 10)] for c in range(81)]
    for c, row in enumerate(values):
        assert sum(row) == 1
        for digit, value in enumerate(row, 1):
            assert 0 <= value <= 1
            if value:
                assert domains[c] & (1 << (digit - 1))
    for h in HOUSES:
        for d in range(9):
            assert sum(values[c][d] for c in h) == 1
    return {"native_equalities": 324, "nonzero_variables": sum(bool(v) for row in values for v in row),
            "nonzero_denominators": sorted({v.denominator for row in values for v in row if v})}


def check_hall_and_crossing_fixed_point(domains):
    for house in HOUSES:
        supported = matching_support(tuple(domains[c] for c in house))
        assert supported == tuple(domains[c] for c in house)
    for center in range(81):
        supported, count = atom_scan(domains, center)
        assert count > 0
        assert all(domains[c] == mask for c, mask in supported.items())


def selected_house_triple_witness(house_ids, domains, palette, roles):
    """Construct a coloring using ONLY the three selected all-different rules.

    Other native peer edges on the same union of cells are deliberately absent:
    this baseline is the conjunction of three selected native constraints.
    """
    selected = [set(HOUSES[i]) for i in house_ids]
    union = set().union(*selected)
    adjacency = {c: set().union(*(h for h in selected if c in h)) - {c} for c in union}
    coloring = {}
    for role, colors in ((1, set(palette)), (0, set(range(1, 10)) - set(palette))):
        vertices = {c for c in union if roles[c] == role}
        remaining = vertices.copy()
        ordering = []
        # Reconstruct and verify a perfect elimination ordering, not a trusted
        # assertion that the selected graph happens to be chordal.
        while remaining:
            simplicial = None
            for c in sorted(remaining):
                neighbors = adjacency[c] & remaining
                if all(b in adjacency[a] for a, b in itertools.combinations(neighbors, 2)):
                    simplicial = c
                    break
            assert simplicial is not None, (house_ids, role, "not chordal")
            ordering.append(simplicial)
            remaining.remove(simplicial)
        for c in reversed(ordering):
            used = {coloring[n] for n in adjacency[c] if n in coloring}
            available = colors - used
            assert available, (house_ids, role, "palette insufficient")
            digit = min(available)
            assert domains[c] & (1 << (digit - 1))
            coloring[c] = digit
    for house in selected:
        assert {coloring[c] for c in house} == set(range(1, 10))
    return tuple((c, coloring[c]) for c in sorted(coloring))


def check_all_native_house_triples(domains, palette, roles):
    digest = hashlib.sha256()
    count = 0
    for ids in itertools.combinations(range(27), 3):
        witness = selected_house_triple_witness(ids, domains, palette, roles)
        digest.update(json.dumps([ids, witness], separators=(",", ":")).encode("ascii"))
        count += 1
    assert count == 2925
    return {"selected_house_triples": count, "constructive_witnesses_sha256": digest.hexdigest(),
            "scope": "Only the three selected native all-different constraints, not every native edge on their union"}


def full_branch_coloring(vertices, palette, node_limit=1000000):
    """Independent full-branch MRV search: no symmetry skipping or memoization.

    Count empty-domain leaves as real nodes. Return a candidate-respecting
    coloring when possible. This is research validation, not a production proof.
    """
    vertices = tuple(sorted(vertices))
    adjacency = [tuple(j for j, other in enumerate(vertices) if shares_house(cell, other))
                 for cell in vertices]
    full = sum(1 << (digit - 1) for digit in palette)
    calls = 0
    empty_leaves = 0

    def search(domains, assigned):
        nonlocal calls, empty_leaves
        calls += 1
        if calls > node_limit:
            raise RuntimeError("Research coloring node budget exhausted")
        if any(mask == 0 for mask in domains):
            empty_leaves += 1
            return None
        if len(assigned) == len(vertices):
            return domains
        cell = min((i for i in range(len(vertices)) if i not in assigned),
                   key=lambda i: (bin(domains[i]).count("1"),
                                  -sum(n not in assigned for n in adjacency[i]), i))
        for bit in bit_values(domains[cell]):
            updated = list(domains)
            updated[cell] = bit
            for neighbor in adjacency[cell]:
                if neighbor not in assigned:
                    updated[neighbor] &= ~bit
            result = search(updated, assigned | {cell})
            if result is not None:
                return result
        return None

    witness = search([full] * len(vertices), set())
    if witness is not None:
        assignment = {c: mask.bit_length() for c, mask in zip(vertices, witness)}
        assert all(assignment[a] != assignment[b] for a, b in itertools.combinations(vertices, 2)
                   if shares_house(a, b))
    else:
        assignment = None
    return assignment, {"nodes": calls, "empty_domain_leaves": empty_leaves,
                        "symmetry_pruning": False, "memoization": False}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--skip-crossing-matching", action="store_true",
                        help="Skip the independent 81-crossing matching fixed-point checks")
    args = parser.parse_args()
    state = json.loads((HERE / "structural_four_digit_13.json").read_text())
    base, assumed = state["base_domains"], state["assumed_domains"]
    assert len(base) == len(assumed) == 81
    assert all(0 < a <= FULL and a & b == a for a, b in zip(assumed, base))
    target, vertices, palette = state["target"], state["support"], state["palette"]
    assert target == {"cell": 16, "digit": 1, "value": False, "inference": "placement"}
    assert vertices == sorted(set(vertices)) and len(vertices) == 13
    mask = sum(1 << (d - 1) for d in palette)
    raw_assumption = base.copy()
    raw_assumption[target["cell"]] &= ~(1 << (target["digit"] - 1))
    assert all(raw_assumption[c] == mask for c in vertices)
    assert all(assumed[c] == mask for c in vertices)
    assert base[target["cell"]] & 1 and not assumed[target["cell"]] & 1
    board = list(map(int, state["completion"]))
    check_full_completion(board, base)
    assert board[target["cell"]] == target["digit"]
    result = {"schema_version": 1, "support_size": 13,
              "maximum_house_support": max(len(set(vertices) & set(h)) for h in HOUSES),
              **independence_check(vertices)}
    result["hand_capacity_proof"] = check_hand_capacity_proof(vertices)
    assert result["maximum_house_support"] == 4
    # The complete respecting solution uses exactly five colors on the support;
    # the independent-set bound excludes four, so chi is exactly five.
    assert len({board[c] for c in vertices}) == 5
    result["chromatic_number"] = 5
    result["linear_relaxation"] = exact_linear_witness(assumed, palette, state["role_witness"])
    result["all_three_house_relaxations"] = check_all_native_house_triples(
        assumed, palette, state["role_witness"])
    if not args.skip_crossing_matching:
        check_hall_and_crossing_fixed_point(base)
        check_hall_and_crossing_fixed_point(assumed)
        result["hall_and_81_crossing_fixed_points"] = ["base", "assumed"]
    witness, result["full_branch_refutation"] = full_branch_coloring(vertices, palette)
    assert witness is None
    result["vertex_deletion_colorings"] = {}
    for deleted in vertices:
        witness, _ = full_branch_coloring([c for c in vertices if c != deleted], palette)
        assert witness is not None
        result["vertex_deletion_colorings"][str(deleted)] = witness
    encoded = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(encoded)
    print(json.dumps({k: v for k, v in result.items() if k != "vertex_deletion_colorings"}, sort_keys=True))


if __name__ == "__main__":
    main()
