"""Discover small native list-coloring cores in unique old-engine stalls.

Optional PySAT is an untrusted discovery oracle. Accepted deductions must still
have an actual Rust exhaustive certificate and independent exact checks.
Per-cell activation assumptions select precisely an induced native peer graph.
"""
import argparse
from collections import Counter
import itertools
import json
from pathlib import Path
import random
import time

from pysat.solvers import Solver
from real_verify import close, sees

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent


class GraphOracle:
    def __init__(self, domains):
        self.domains = domains
        clauses = []
        self.cells = [c for c, m in enumerate(domains) if m & (m - 1)]
        # Singleton cells are harmless to omit once native propagation closed.
        for c in self.cells:
            digits = [d for d in range(9) if domains[c] & (1 << d)]
            variables = [9 * c + d + 1 for d in digits]
            clauses.append([-1000 - c] + variables)
            clauses.extend([-a, -b] for a, b in itertools.combinations(variables, 2))
        for a, b in itertools.combinations(self.cells, 2):
            if sees(a, b):
                clauses.extend([-9 * a - d - 1, -9 * b - d - 1]
                               for d in range(9) if domains[a] & domains[b] & (1 << d))
        self.solver = Solver(name="m22", bootstrap_with=clauses)

    def delete(self):
        self.solver.delete()

    def core(self, cells, literal, rng, tries=6):
        assumptions = [1000 + c for c in cells] + [literal]
        if self.solver.solve(assumptions=assumptions):
            return None
        available = self.solver.get_core()
        best = sorted(a - 1000 for a in available if a >= 1000)
        target = (abs(literal) - 1) // 9
        if target not in best:
            best.append(target)
            best.sort()
        initial = best.copy()
        for _ in range(tries):
            current = initial.copy()
            order = current.copy()
            rng.shuffle(order)
            for c in order:
                if c == target:
                    continue
                candidate = [v for v in current if v != c]
                if not self.solver.solve(assumptions=[1000 + v for v in candidate] + [literal]):
                    current = candidate
            if len(current) < len(best):
                best = current
        return best


def known_three_digit(domains):
    """Reject any bare three-palette coloring contradiction after local closure.

    This is a deliberate discovery filter, not a complete known-method census.
    """
    oracle = GraphOracle(domains)
    try:
        for digits in itertools.combinations(range(9), 3):
            mask = sum(1 << d for d in digits)
            cells = [c for c in oracle.cells if domains[c] & ~mask == 0]
            if len(cells) >= 5 and not oracle.solver.solve(assumptions=[1000 + c for c in cells]):
                core = sorted(a - 1000 for a in oracle.solver.get_core())
                return {"palette": [d + 1 for d in digits], "cells": core}
        return None
    finally:
        oracle.delete()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("fixtures", type=Path, nargs="+")
    p.add_argument("--seconds", type=float, default=180)
    p.add_argument("--seed", type=int, default=720201)
    p.add_argument("--tries", type=int, default=8)
    p.add_argument("--output-dir", type=Path, default=ROOT / "target/different-sat")
    args = p.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    rng = random.Random(args.seed)
    started = time.monotonic()
    counts, records, positives = Counter(), [], []
    for path in args.fixtures:
        if time.monotonic() - started >= args.seconds:
            break
        case = json.loads(path.read_text())
        if case.get("status") != "full-old-stall":
            continue
        domains = case["base_domains"]
        assert close(domains)[1] == []
        oracle = GraphOracle(domains)
        solution = list(map(int, case["solution"]))
        targets = [(c, d, d == solution[c]) for c in oracle.cells
                   for d in range(1, 10) if domains[c] & (1 << (d - 1))]
        rng.shuffle(targets)
        targets.sort(key=lambda t: (bin(domains[t[0]]).count("1"), not t[2]))
        try:
            for cell, digit, placement in targets:
                if time.monotonic() - started >= args.seconds:
                    break
                counts["targets"] += 1
                raw = domains.copy()
                bit = 1 << (digit - 1)
                raw[cell] = raw[cell] & ~bit if placement else bit
                try:
                    assumed, trace = close(raw)
                except ValueError:
                    counts["local-refuted"] += 1
                    continue
                counts["local-survivors"] += 1
                known = known_three_digit(assumed)
                if known:
                    counts["bare-three-digit-refutation"] += 1
                    records.append({"source": str(path), "cell": cell, "digit": digit,
                                    "placement": placement, "rejected_q3": known})
                    continue
                counts["no-bare-three-digit-refutation"] += 1
                literal = -(9 * cell + digit) if placement else 9 * cell + digit
                cells = oracle.core(oracle.cells, literal, rng, args.tries)
                if cells is None:
                    counts["no-selected-core"] += 1
                    continue
                counts["core-size-" + str(len(cells))] += 1
                record = {"source": str(path), "cell": cell, "digit": digit,
                          "placement": placement, "cells": cells,
                          "assumption_trace": trace, "assumed_domains": assumed}
                records.append(record)
                if len(cells) <= 32:
                    name = "candidate-" + str(len(positives)) + ".json"
                    result = {**case, "target": {"cell": cell, "digit": digit, "value": not placement},
                              "conclusion": {"cell": cell, "digit": digit, "value": placement,
                                             "kind": "placement" if placement else "elimination"},
                              "graph_cells": cells, "assumption_trace": trace,
                              "assumed_domains": assumed,
                              "assumed_state": {"version": 1, "values": [m.bit_length() if m & (m-1) == 0 else 0 for m in assumed], "domains": assumed},
                              "discovery_only": "Requires exact proof/template/novelty acceptance"}
                    (args.output_dir / name).write_text(json.dumps(result, indent=2) + "\n")
                    positives.append(name)
                    print("SMALL_CORE", name, len(cells), cell, digit, placement, flush=True)
        finally:
            oracle.delete()
        print("CASE", path.name, dict(counts), flush=True)
    summary = {"schema_version": 1, "seed": args.seed, "seconds": time.monotonic() - started,
               "counts": dict(counts), "positive_discovery_files": positives, "records": records,
               "scope": "SAT discovery only; known-three-digit filter is not historical novelty proof"}
    (args.output_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps({k: v for k, v in summary.items() if k != "records"}), flush=True)


if __name__ == "__main__":
    main()
