"""Search for, then EXACTLY check, full native Sudoku LP feasibility witnesses.

Verification uses Python's standard library only. Search additionally uses numpy
and scipy.optimize.linprog; floating output is never an acceptance criterion.
Witnesses bind the full source state and target assumption by SHA-256. The LP is
the 729-indicator native relaxation: all 324 exactly-one equalities, nonnegative
values, zeros outside effective candidate domains, and the advertised target.
This does not test stronger local-marginal LP or semidefinite hierarchies.
"""
import argparse
import hashlib
import json
from fractions import Fraction as F
from pathlib import Path


def source_state(fixture):
    state = fixture.get("assumed_state", fixture.get("state", fixture))
    # Real-state research fixtures retain both pre-assumption and fully closed
    # post-assumption domains. Test the stronger, recorded assumed state.
    domains = state.get("assumed_domains", state.get("assumed_closed_domains", state.get("domains")))
    if domains is None:
        raise ValueError("Missing domains")
    values = state.get("values")
    if values is None:
        puzzle = fixture.get("puzzle")
        if puzzle is not None:
            if not isinstance(puzzle, str) or len(puzzle) != 81 or any(c not in ".0123456789" for c in puzzle):
                raise ValueError("Invalid givens string")
            values = [0 if char == "." else int(char) for char in puzzle]
        else:
            values = [0] * 81
    if len(domains) != 81 or len(values) != 81:
        raise ValueError("Expected 81 values and domains")
    effective = []
    for value, domain in zip(values, domains):
        if type(value) is not int or type(domain) is not int or not 0 <= value <= 9 or not 0 < domain <= 511:
            raise ValueError("Invalid value or domain")
        if value and domain != 1 << (value - 1):
            raise ValueError("Placed value must have a singleton effective domain")
        effective.append(domain)
    target = fixture["target"]
    cell, digit = target["cell"], target["digit"]
    assumed = target.get("assumed", target.get("assumption", target.get("value", target.get("inference", "elimination") != "placement")))
    if type(assumed) is not bool or type(cell) is not int or type(digit) is not int:
        raise ValueError("Malformed target assumption")
    if not 0 <= cell < 81 or not 1 <= digit <= 9:
        raise ValueError("Target must name a cell and a digit in 1..9")
    # Under the recorded opposite assumption, propagation may already have
    # placed this cell to a different digit. That remains a legitimate LP
    # premise; the deduction certificate separately checks the unforced target.
    if values[cell] and ((values[cell] == digit) != assumed):
        raise ValueError("Post-assumption target placement contradicts the assumption")
    bit = 1 << (digit - 1)
    if assumed and not effective[cell] & bit:
        raise ValueError("Assumed-true target candidate is absent")
    if assumed:
        effective[cell] = bit
    else:
        effective[cell] &= ~bit
    return values, domains, effective, {"cell": cell, "digit": digit, "assumed": assumed}


def binding(fixture):
    values, domains, _, target = source_state(fixture)
    payload = {"schema": "native-lp-source-v1", "values": values, "domains": domains, "target": target}
    return hashlib.sha256(json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def requirement_ids(cell, digit):
    r, c = divmod(cell, 9)
    return [cell, 81 + 9 * r + digit - 1, 162 + 9 * c + digit - 1,
            243 + 9 * (3 * (r // 3) + c // 3) + digit - 1]


def verify(fixture, witness):
    _, _, effective, target = source_state(fixture)
    if witness.get("schema_version") != 1 or witness.get("state_hash") != binding(fixture):
        raise ValueError("LP witness version/state/target mismatch")
    if witness.get("target") != target:
        raise ValueError("LP witness target differs")
    dense = [F(0)] * 729
    seen = set()
    for item in witness["nonzero"]:
        if set(item) != {"cell", "digit", "numerator", "denominator"}:
            raise ValueError("Malformed nonzero entry")
        cell, digit, num, den = (item[k] for k in ["cell", "digit", "numerator", "denominator"])
        if any(type(v) is not int for v in [cell, digit, num, den]):
            raise ValueError("Rational entries must use integers")
        if not 0 <= cell < 81 or not 1 <= digit <= 9 or den <= 0 or num <= 0:
            raise ValueError("Invalid candidate or nonpositive rational entry")
        if not effective[cell] & (1 << (digit - 1)) or (cell, digit) in seen:
            raise ValueError("Unavailable or duplicate candidate entry")
        seen.add((cell, digit))
        dense[9 * cell + digit - 1] = F(num, den)
    totals = [F(0)] * 324
    for cell in range(81):
        for digit in range(1, 10):
            value = dense[9 * cell + digit - 1]
            if not 0 <= value <= 1:
                raise ValueError("Indicator outside [0,1]")
            for requirement in requirement_ids(cell, digit):
                totals[requirement] += value
    failed = [i for i, total in enumerate(totals) if total != 1]
    if failed:
        raise ValueError("Native equality failed: " + str(failed[:5]))
    target_value = dense[9 * target["cell"] + target["digit"] - 1]
    if target_value != int(target["assumed"]):
        raise ValueError("Target assumption is not satisfied")
    return {"equalities": 324, "nonzeros": len(seen),
            "max_denominator_bits": max(v.denominator.bit_length() for v in dense),
            "max_numerator_bits": max(v.numerator.bit_length() for v in dense)}


def verify_farkas(fixture, witness):
    """Exact negative-control check: a linear certificate rejecting a fixture."""
    _, _, effective, target = source_state(fixture)
    if witness.get("schema_version") != 1 or witness.get("state_hash") != binding(fixture) or witness.get("target") != target:
        raise ValueError("Farkas witness version/state/target mismatch")
    weights = {}
    for item in witness["terms"]:
        if not isinstance(item, list) or len(item) != 3 or any(type(v) is not int for v in item):
            raise ValueError("Malformed Farkas term")
        source, numerator, denominator = item
        if not 0 <= source < 324 or source in weights or denominator <= 0 or numerator == 0:
            raise ValueError("Invalid Farkas term")
        weights[source] = F(numerator, denominator)
    rhs = sum(weights.values(), F(0))
    if not rhs < 0 or rhs != F(witness["weighted_rhs"]):
        raise ValueError("Farkas right-hand side is not the recorded negative value")
    for cell in range(81):
        for digit in range(1, 10):
            if effective[cell] & (1 << (digit - 1)):
                if sum(weights.get(i, F(0)) for i in requirement_ids(cell, digit)) < 0:
                    raise ValueError("Farkas left-hand side has a negative available coefficient")
    return {"terms": len(weights), "rhs": str(rhs), "conclusion": "native LP infeasible"}


def rational_basic_solution(matrix, approximate, tolerance=1e-9):
    """Recover the floating basic solution's nonzero support over exact Q.

    This can fail safely for nearly-zero entries or unsuitable numerical bases.
    Any returned vector still needs the independent full-state verifier.
    """
    active = [i for i, value in enumerate(approximate) if value > tolerance]
    rows = [[F(int(row[j])) for j in active] + [F(1)] for row in matrix]
    rank = 0
    pivots = []
    for col in range(len(active)):
        pivot = next((r for r in range(rank, len(rows)) if rows[r][col]), None)
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        divisor = rows[rank][col]
        rows[rank] = [v / divisor for v in rows[rank]]
        for r in range(len(rows)):
            if r == rank or not rows[r][col]:
                continue
            multiplier = rows[r][col]
            rows[r] = [a - multiplier * b for a, b in zip(rows[r], rows[rank])]
        pivots.append(col)
        rank += 1
    if any(not any(row[:-1]) and row[-1] for row in rows):
        raise ValueError("Floating support has no exact rational solution")
    free = set(range(len(active))) - set(pivots)
    solution = [F(0)] * len(active)
    for col in free:
        solution[col] = F(float(approximate[active[col]])).limit_denominator(1_000_000)
    for row, col in zip(rows[:rank], pivots):
        solution[col] = row[-1] - sum(row[j] * solution[j] for j in free)
    result = [F(0)] * len(approximate)
    for index, value in zip(active, solution):
        result[index] = value
    if any(value < 0 for value in result):
        raise ValueError("Recovered exact basic solution is negative")
    return result


def search(fixture):
    import numpy as np
    from scipy.optimize import linprog
    _, _, effective, target = source_state(fixture)
    variables = [(cell, digit) for cell in range(81) for digit in range(1, 10)
                 if effective[cell] & (1 << (digit - 1))]
    matrix = np.zeros((324, len(variables)), dtype=np.int8)
    for j, (cell, digit) in enumerate(variables):
        for requirement in requirement_ids(cell, digit):
            matrix[requirement, j] = 1
    solved = linprog(np.zeros(len(variables)), A_eq=matrix, b_eq=np.ones(324),
                     bounds=(0, None), method="highs-ds")
    if not solved.success:
        raise ValueError("No floating primal feasibility witness: " + solved.message)
    # Cheap rational reconstruction first, exact basis elimination on failure.
    rational = [F(float(v)).limit_denominator(10000) for v in solved.x]
    for attempt in range(2):
        witness = {"schema_version": 1, "state_hash": binding(fixture), "target": target,
                   "nonzero": [{"cell": c, "digit": d, "numerator": q.numerator,
                                "denominator": q.denominator}
                               for (c, d), q in zip(variables, rational) if q]}
        try:
            witness["verification"] = verify(fixture, witness)
            return witness
        except ValueError:
            if attempt:
                raise
            rational = rational_basic_solution(matrix, solved.x)
    raise AssertionError("unreachable")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=["search", "verify", "verify-farkas"])
    parser.add_argument("fixture", type=Path)
    parser.add_argument("witness", type=Path)
    args = parser.parse_args()
    fixture = json.loads(args.fixture.read_text())
    if args.mode == "search":
        witness = search(fixture)
        args.witness.write_text(json.dumps(witness, indent=2) + "\n")
    else:
        witness = json.loads(args.witness.read_text())
    verifier = verify_farkas if args.mode == "verify-farkas" else verify
    print(json.dumps(verifier(fixture, witness), sort_keys=True))


if __name__ == "__main__":
    main()
