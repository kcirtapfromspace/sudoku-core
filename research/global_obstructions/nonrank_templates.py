"""Discover exact primal witnesses for the full single-digit template LP.

There are 46,656 legal placements of one digit in classic Sudoku. A variable is
the probability of one domain-compatible template for one digit. Nine digit
mass equations and 81 cell coverage equations couple the nine distributions.
This contains every valid inequality of each individual digit's placement
polytope. It is not the integer hull of the multi-digit Sudoku problem.

SciPy proposes a basic point; acceptance below is exact rational arithmetic.
The separate nonrank_audit_* checker does not import this search implementation.
"""
import os
for _name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_name] = "1"

import argparse
from fractions import Fraction
from functools import lru_cache
import hashlib
import itertools
import json
from pathlib import Path
import time

from linear_witness import rational_basic_solution, source_state

RELAXATION = "native-single-digit-template-v1"


def premises(fixture):
    values, _, domains, target = source_state(fixture)
    payload = {"schema": "native-template-lp-source-v1", "values": values,
               "domains": domains, "target": target}
    digest = hashlib.sha256(json.dumps(payload, sort_keys=True,
                                      separators=(",", ":")).encode()).hexdigest()
    return values, domains, target, digest


@lru_cache(None)
def all_templates():
    templates = []
    for columns in itertools.permutations(range(9)):
        if any(len({columns[r] // 3 for r in range(start, start + 3)}) != 3
               for start in (0, 3, 6)):
            continue
        cells = tuple(9 * row + column for row, column in enumerate(columns))
        templates.append((sum(1 << cell for cell in cells), cells))
    if len(templates) != 46656:
        raise ValueError("Template enumeration count changed")
    return tuple(templates)


def compatible_templates(domains):
    answer = []
    for digit in range(1, 10):
        bit = 1 << (digit - 1)
        allowed = sum(1 << cell for cell, mask in enumerate(domains) if mask & bit)
        forced = sum(1 << cell for cell, mask in enumerate(domains) if mask == bit)
        answer.append([(mask, cells) for mask, cells in all_templates()
                       if mask & ~allowed == 0 and mask & forced == forced])
    return answer


def verify(fixture, witness):
    values, domains, target, digest = premises(fixture)
    if (witness.get("schema_version") != 1 or witness.get("relaxation") != RELAXATION
            or witness.get("state_hash") != digest or witness.get("target") != target):
        raise ValueError("Template witness binding mismatch")
    totals = [Fraction(0)] * 90
    indicators = [[Fraction(0)] * 9 for _ in range(81)]
    seen = set()
    for entry in witness["terms"]:
        if set(entry) != {"digit", "cells", "numerator", "denominator"}:
            raise ValueError("Malformed template term")
        digit, cells = entry["digit"], entry["cells"]
        numerator, denominator = entry["numerator"], entry["denominator"]
        if (any(type(value) is not int for value in (digit, numerator, denominator))
                or not 1 <= digit <= 9 or numerator <= 0 or denominator <= 0):
            raise ValueError("Malformed template digit or rational")
        if (len(cells) != 9 or any(type(c) is not int or not 0 <= c < 81 for c in cells)
                or cells != sorted(set(cells))):
            raise ValueError("Malformed template cells")
        if (len({c // 9 for c in cells}) != 9 or len({c % 9 for c in cells}) != 9
                or len({3 * (c // 27) + c % 9 // 3 for c in cells}) != 9):
            raise ValueError("Template is not one cell per native house")
        bit = 1 << (digit - 1)
        if any(not domains[c] & bit for c in cells):
            raise ValueError("Template violates effective domains")
        if any(value == digit and c not in cells for c, value in enumerate(values)):
            raise ValueError("Template omits an already placed digit")
        key = digit, tuple(cells)
        if key in seen:
            raise ValueError("Duplicate template")
        seen.add(key)
        weight = Fraction(numerator, denominator)
        totals[digit - 1] += weight
        for cell in cells:
            totals[9 + cell] += weight
            indicators[cell][digit - 1] += weight
    if any(total != 1 for total in totals):
        raise ValueError("A digit mass or cell coverage equation failed")
    if indicators[target["cell"]][target["digit"] - 1] != int(target["assumed"]):
        raise ValueError("Target probability differs from assumption")
    return {"template_terms": len(seen), "equalities": 90,
            "implied_native_equalities": 324,
            "max_denominator": max(term["denominator"] for term in witness["terms"])}


def search(fixture, time_limit=15, max_columns=200000):
    import numpy as np
    from scipy.optimize import linprog
    from scipy.sparse import csc_matrix
    started = time.monotonic()
    _, domains, target, digest = premises(fixture)
    compatible = compatible_templates(domains)
    counts = list(map(len, compatible))
    status = {"relaxation": RELAXATION, "compatible_templates_by_digit": counts,
              "template_universe_count": 46656, "columns": sum(counts)}
    if not all(counts):
        return {**status, "status": "a digit has no compatible template", "witness": None}
    if sum(counts) > max_columns:
        return {**status, "status": "column budget cutoff", "witness": None}
    columns = [(digit, cells) for digit, terms in enumerate(compatible, 1) for _, cells in terms]
    rows, column_indices = [], []
    for index, (digit, cells) in enumerate(columns):
        rows.extend([digit - 1] + [9 + cell for cell in cells])
        column_indices.extend([index] * 10)
    matrix = csc_matrix((np.ones(len(rows)), (rows, column_indices)),
                        shape=(90, len(columns)))
    result = linprog(np.zeros(len(columns)), A_eq=matrix, b_eq=np.ones(90),
                     bounds=(0, None), method="highs-ds", options={"time_limit": time_limit})
    if not result.success:
        return {**status, "status": "floating optimizer did not return a primal",
                "solver_status": int(result.status), "solver_message": result.message,
                "seconds": time.monotonic() - started, "witness": None}
    active = [i for i, value in enumerate(result.x) if value > 1e-9]
    rational = [Fraction(float(result.x[i])).limit_denominator(100000) for i in active]

    def make_witness(weights):
        return {"schema_version": 1, "relaxation": RELAXATION,
                "state_hash": digest, "target": target,
                "terms": [{"digit": columns[i][0], "cells": list(columns[i][1]),
                           "numerator": weight.numerator, "denominator": weight.denominator}
                          for i, weight in zip(active, weights) if weight]}

    witness = make_witness(rational)
    try:
        checked = verify(fixture, witness)
    except ValueError:
        rational = rational_basic_solution(matrix[:, active].toarray(), result.x[active])
        witness = make_witness(rational)
        checked = verify(fixture, witness)
    return {**status, "status": "exact rational primal verified", "check": checked,
            "seconds": time.monotonic() - started, "witness": witness}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("fixture", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--time-limit", type=float, default=15)
    parser.add_argument("--max-columns", type=int, default=200000)
    args = parser.parse_args()
    result = search(json.loads(args.fixture.read_text()), args.time_limit, args.max_columns)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({k: value for k, value in result.items() if k != "witness"}))


if __name__ == "__main__":
    main()
