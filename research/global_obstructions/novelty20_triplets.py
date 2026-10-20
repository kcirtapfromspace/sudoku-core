"""Exact witnesses excluding every direct three-digit list contradiction.

For every palette of three labels, color ALL native cells whose effective
opposite-state lists are contained in it. Restricting a proper list coloring to
any subset remains a proper list coloring. Thus no impossible three-digit
pattern, catalogued or not, can occur directly on such lists in this snapshot.
Guarded patterns, OR relations, conditional branches and changed snapshots are
not excluded. Generation is optional; acceptance checks explicit witnesses.
"""
import argparse
import hashlib
import itertools
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
FIXTURE = HERE / "different_lead.json"
OUTPUT = HERE / "novelty20_triplets.json"


def adjacent(a, b):
    return a != b and (a//9 == b//9 or a%9 == b%9 or
                       (a//27, a%9//3) == (b//27, b%9//3))


def generate(case):
    domains = case["assumed_domains"]
    records = []
    for palette in itertools.combinations(range(1, 10), 3):
        mask = sum(1 << (d-1) for d in palette)
        cells = [c for c, m in enumerate(domains) if m & ~mask == 0]
        neighbors = {c: [q for q in cells if adjacent(c, q)] for c in cells}
        assigned = {}

        def search():
            available = {c: domains[c] & ~sum({1 << (assigned[q]-1)
                         for q in neighbors[c] if q in assigned})
                         for c in cells if c not in assigned}
            if not available:
                return True
            c = min(available, key=lambda q: (bin(available[q]).count("1"),
                                               -len(neighbors[q]), q))
            for d in palette:
                if available[c] & (1 << (d-1)):
                    assigned[c] = d
                    if search():
                        return True
                    del assigned[c]
            return False

        if not search():
            raise ValueError("Direct contradiction found for palette " + str(palette))
        records.append({"palette": list(palette), "cells": cells,
                        "coloring": [assigned[c] for c in cells]})
    return {"schema_version": 1,
            "fixture_sha256": hashlib.sha256(FIXTURE.read_bytes()).hexdigest(),
            "target_assumption": case["target"], "snapshot": "assumed_domains",
            "records": records,
            "scope": "All direct three-digit list obstructions excluded in this snapshot; guarded OR relations and their consequences are outside this certificate."}


def verify(case, proof):
    if (proof["schema_version"] != 1 or proof["snapshot"] != "assumed_domains"
            or proof["fixture_sha256"] != hashlib.sha256(FIXTURE.read_bytes()).hexdigest()
            or proof["target_assumption"] != case["target"]):
        raise ValueError("Source binding mismatch")
    domains = case["assumed_domains"]
    expected = list(itertools.combinations(range(1, 10), 3))
    if [tuple(row["palette"]) for row in proof["records"]] != expected:
        raise ValueError("Every three-digit palette must occur exactly once")
    assignments = 0
    for row in proof["records"]:
        mask = sum(1 << (d-1) for d in row["palette"])
        cells = [c for c, m in enumerate(domains) if m & ~mask == 0]
        if cells != row["cells"] or len(cells) != len(row["coloring"]):
            raise ValueError("Witness omits an eligible native cell")
        for i, (c, d) in enumerate(zip(cells, row["coloring"])):
            if type(d) is not int or d not in row["palette"] or not domains[c] & (1 << (d-1)):
                raise ValueError("Color violates its actual list")
            for q, other in zip(cells[:i], row["coloring"][:i]):
                r, col = divmod(c, 9)
                s, other_col = divmod(q, 9)
                linked = r == s or col == other_col or (r//3, col//3) == (s//3, other_col//3)
                if linked and d == other:
                    raise ValueError("Native peer conflict")
        assignments += len(cells)
    return {"all_palettes": len(expected), "checked_cell_assignments": assignments,
            "direct_three_digit_list_obstruction": False,
            "scope": proof["scope"]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--generate", action="store_true")
    args = parser.parse_args()
    case = json.loads(FIXTURE.read_text())
    if args.generate:
        result = generate(case)
        verify(case, result)
        OUTPUT.write_text(json.dumps(result, indent=2) + "\n")
    result = json.loads(OUTPUT.read_text())
    print(json.dumps(verify(case, result)))
    bad = json.loads(json.dumps(result))
    bad["records"].pop()
    try:
        verify(case, bad)
    except ValueError:
        pass
    else:
        raise ValueError("Missing palette was accepted")
    bad = json.loads(json.dumps(result))
    bad["records"][0]["coloring"][0] = 0
    try:
        verify(case, bad)
    except ValueError:
        pass
    else:
        raise ValueError("Invalid color was accepted")


if __name__ == "__main__":
    main()
