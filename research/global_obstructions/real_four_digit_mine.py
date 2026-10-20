"""Try givens-only realizations of the discovered thirteen-cell four-color core.

Candidate givens come from the respecting complete Sudoku, outside the core and
outside its four-digit palette. Minimal per-digit peer covers preserve the desired
core restrictions. Floating LP is only a discovery screen; a positive requires
the separate exact rational verifier before being reported as a certified result.
"""
import argparse
import itertools
import json
from pathlib import Path
import tempfile
import time

from real_verify import ALL, HOUSES, close, sees

HERE = Path(__file__).resolve().parent


def minimal_covers(completion, core, guardian, digit):
    possible = [c for c, value in enumerate(completion) if value == digit and c not in core]
    required = set(core) - ({guardian} if digit == completion[guardian] else set())
    covers = []
    for mask in range(1 << len(possible)):
        chosen = {possible[i] for i in range(len(possible)) if mask & (1 << i)}
        if any(earlier < chosen for earlier in covers):
            continue
        if all(any(sees(cell, given) for given in chosen) for cell in required):
            covers = [earlier for earlier in covers if not chosen < earlier]
            covers.append(chosen)
    return [sorted(cells) for cells in covers]


def initial_domains(givens):
    return [1 << (d - 1) if d else ALL & ~sum({1 << (givens[p] - 1)
            for house in HOUSES if cell in house for p in house if givens[p]})
            for cell, d in enumerate(givens)]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-trials", type=int, default=5040)
    parser.add_argument("--candidate-output", type=Path,
                        default=Path(tempfile.gettempdir()) / "sudoku-global-four-digit-discovery.json")
    args = parser.parse_args()
    import numpy as np
    from scipy.optimize import linprog
    source = json.loads((HERE / "structural_four_digit_13.json").read_text())
    completion = list(map(int, source["completion"]))
    core, palette = source["support"], source["palette"]
    guardian, extra = source["target"]["cell"], source["target"]["digit"]
    choices = {digit: minimal_covers(completion, core, guardian, digit)
               for digit in range(1, 10) if digit not in palette}
    matrix = np.zeros((324, 729), dtype=np.int8)
    for cell in range(81):
        row, col = divmod(cell, 9)
        box = row // 3 * 3 + col // 3
        for digit in range(9):
            matrix[[cell, 81 + row * 9 + digit, 162 + col * 9 + digit,
                    243 + box * 9 + digit], cell * 9 + digit] = 1
    stats = {"schema_version": 1, "selection": "All inclusion-minimal independent per-digit covers",
             "cover_counts": {str(k): len(v) for k, v in choices.items()},
             "trial_limit": args.max_trials, "trials": 0, "floating_lp_feasible": 0,
             "local_survivors": 0, "positives": [],
             "completion": source["completion"], "core": core, "palette": palette,
             "target": {"cell": guardian, "digit": extra, "value": False}}
    started = time.perf_counter()
    for groups in itertools.islice(itertools.product(*choices.values()), args.max_trials):
        stats["trials"] += 1
        selected = set().union(*map(set, groups))
        givens = [d if cell in selected else 0 for cell, d in enumerate(completion)]
        domains = initial_domains(givens)
        assumed = domains.copy()
        assumed[guardian] &= ~(1 << (extra - 1))
        variables = [cell * 9 + digit for cell, mask in enumerate(assumed)
                     for digit in range(9) if mask & (1 << digit)]
        result = linprog(np.zeros(len(variables)), A_eq=matrix[:, variables],
                         b_eq=np.ones(324), bounds=(0, None), method="highs-ds")
        if result.success:
            stats["floating_lp_feasible"] += 1
            try:
                base, base_trace = close(domains)
                branch = base.copy()
                branch[guardian] &= ~(1 << (extra - 1))
                branch, branch_trace = close(branch)
            except ValueError:
                continue
            stats["local_survivors"] += 1
            positive = {"schema_version": 1, "name": "generated-four-digit-realization",
                        "puzzle": "".join(str(d) if d else "." for d in givens),
                        "initial_domains": domains, "base_domains": base,
                        "base_trace": base_trace, "assumed_domains": branch,
                        "assumption_trace": branch_trace, "completion": source["completion"],
                        "graph_cells": core, "palette": palette, "target": stats["target"],
                        "trial": stats["trials"],
                        "status": "Uncertified discovery candidate; exact LP witness required"}
            stats["positives"].append(positive)
            args.candidate_output.write_text(json.dumps(positive, indent=2) + "\n")
            print("Candidate found at trial", stats["trials"], positive["puzzle"], flush=True)
            break
        if stats["trials"] % 500 == 0:
            print("trials", stats["trials"], "floating LP survivors", stats["floating_lp_feasible"], flush=True)
    stats["seconds"] = time.perf_counter() - started
    stats["negative_result_scope"] = "Discovery screen only; floating infeasibility is not an exact certificate. The tested clue-cover family is restricted."
    (HERE / "real_four_digit_mining.json").write_text(json.dumps(stats, indent=2) + "\n")
    print(json.dumps({k: v for k, v in stats.items() if k != "positives"}), flush=True)


if __name__ == "__main__":
    main()
