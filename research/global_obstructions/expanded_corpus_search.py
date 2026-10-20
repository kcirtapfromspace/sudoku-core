"""Search published corpus states for native four-digit coloring obstructions.

Discovery only. A retained candidate still needs an exact rational LP witness,
independent fixture replay, and a production coloring certificate. The frozen
baseline is 0532d63, with nonforcing and full-dispatcher outcomes kept separate.
"""
import argparse
from collections import Counter
import hashlib
import itertools
import json
import os
from pathlib import Path
import subprocess
import time

for variable in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(variable, "1")

from real_four_digit_mine import initial_domains
from real_verify import close, exact_completions, sees
from structural_discover import four_colorable

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def state(domains):
    return {"version": 1, "domains": domains,
            "values": [m.bit_length() if m & (m - 1) == 0 else 0 for m in domains]}


def old_probe(binary, item, domains, full):
    for cell, mask in enumerate(domains):
        if mask & (mask - 1):
            if any(other != cell and sees(cell, other) and
                   other_mask & (other_mask - 1) == 0 and mask & other_mask
                   for other, other_mask in enumerate(domains)):
                raise ValueError("Probe input has unpropagated placed singleton candidates")
    request = {"id": item["id"], "puzzle": item["puzzle"],
               "state": state(domains), "unique": True, "full": full}
    response = subprocess.run([str(binary)], input=json.dumps(request) + "\n",
                              text=True, capture_output=True, timeout=90)
    if response.returncode:
        raise ValueError("Frozen baseline execution failed: " + response.stderr[-400:])
    return json.loads(response.stdout)


def joint_close(binary, item, domains, completion):
    # Singleton masks become placed values in the portable state. Propagate
    # them before handing masks to older detectors that expect legal peers.
    domains, trace = close(domains)
    for _ in range(100):
        old = old_probe(binary, item, domains, False)
        for _, cell, digits, placement in old["trace"]:
            if (digits != [completion[cell]] if placement else completion[cell] in digits):
                raise ValueError("Frozen baseline step excludes the independently counted completion")
        if old["trace"]:
            trace.append({"kind": "existing-engine", "steps": old["trace"],
                          "state_after": old["state"]})
        reduced, local_trace = close(old["state"]["domains"])
        trace.extend(local_trace)
        if not all(mask & (1 << (digit - 1)) for mask, digit in zip(reduced, completion)):
            raise ValueError("Closure excluded the independently counted completion")
        if not old["trace"] and not local_trace:
            return reduced, trace
        domains = reduced
    raise ValueError("Joint closure limit")


def patterns(domains, counts):
    for digits in itertools.combinations(range(1, 10), 4):
        palette = sum(1 << (d - 1) for d in digits)
        pure = [c for c, mask in enumerate(domains)
                if mask & ~palette == 0 and mask & (mask - 1)]
        for guardian, mask in enumerate(domains):
            extra = mask & ~palette
            if bin(extra).count("1") != 1 or bin(mask & palette).count("1") < 2:
                continue
            pool, seen, stack = set(pure) | {guardian}, {guardian}, [guardian]
            while stack:
                cell = stack.pop()
                neighbors = {other for other in pool - seen if sees(cell, other)}
                seen.update(neighbors)
                stack.extend(neighbors)
            if not 5 <= len(seen) <= 24:
                counts["components_outside_5_to_24"] += 1
                continue
            cells = sorted(seen)
            counts["component_presentations"] += 1
            try:
                colorable, _ = four_colorable(cells, node_limit=100000)
            except RuntimeError:
                counts["uniform_coloring_cutoffs"] += 1
                continue
            if colorable:
                counts["uniform_four_colorable"] += 1
                continue
            counts["uniform_non_four_colorable"] += 1
            yield list(digits), cells, guardian, extra.bit_length()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, required=True,
                        help="Frozen real_existing_probe built by real_existing_check.build_engine")
    parser.add_argument("--limit", type=int, default=106)
    parser.add_argument("--corpus", type=Path,
                        default=HERE.parent / "nonmodular/real_game_corpus.json")
    parser.add_argument("--output", type=Path, default=ROOT / "target/expanded-corpus-search.json")
    args = parser.parse_args()
    if not __debug__:
        raise SystemExit("Run without -O")
    import numpy as np
    from scipy.optimize import linprog
    matrix = np.zeros((324, 729), dtype=np.int8)
    for cell in range(81):
        row, column = divmod(cell, 9)
        box = row // 3 * 3 + column // 3
        for d in range(9):
            matrix[[cell, 81 + row * 9 + d, 162 + column * 9 + d,
                    243 + box * 9 + d], cell * 9 + d] = 1
    corpus_path = args.corpus
    corpus = json.loads(corpus_path.read_text())[:args.limit]
    counts, records, candidates = Counter(), [], []
    started = time.monotonic()
    for item in corpus:
        counts["puzzles_tested"] += 1
        record = {"id": item["id"]}
        try:
            completions = exact_completions(item["puzzle"])
            if len(completions) != 1:
                counts["not_unique"] += 1
                continue
            counts["independently_unique"] += 1
            completion = list(map(int, completions[0]))
            initial = initial_domains([int(c) if c in "123456789" else 0 for c in item["puzzle"]])
            base, trace = joint_close(args.binary, item, initial, completion)
            record["joint_trace_events"] = len(trace)
            if all(mask & (mask - 1) == 0 for mask in base):
                counts["joint_baseline_solved"] += 1
                record["status"] = "solved by older nonforcing plus local closure"
            else:
                counts["joint_baseline_stalled"] += 1
                record["status"] = "nonforcing/local fixed point"
                for palette, cells, cell, digit in patterns(base, counts):
                    branch = base.copy()
                    branch[cell] &= ~(1 << (digit - 1))
                    try:
                        branch, branch_trace = close(branch)
                    except ValueError:
                        counts["opposite_local_contradiction"] += 1
                        continue
                    variables = [c * 9 + d for c, mask in enumerate(branch)
                                 for d in range(9) if mask & (1 << d)]
                    lp = linprog(np.zeros(len(variables)), A_eq=matrix[:, variables],
                                 b_eq=np.ones(324), bounds=(0, None), method="highs-ds")
                    if not lp.success:
                        counts["opposite_lp_screen_rejected"] += 1
                        continue
                    counts["discovery_candidates"] += 1
                    full = old_probe(args.binary, item, base, True)
                    if full["trace"] or full["state"]["domains"] != base:
                        counts["repeated_baseline_changed"] += 1
                        continue
                    candidate = {"name": item["id"], "puzzle": item["puzzle"],
                                 "source": item, "solution": completions[0], "completion": completions[0],
                                 "solutions_capped_at_2": 1, "initial_domains": initial,
                                 "base_domains": base, "base_trace": trace,
                                 "assumed_domains": branch, "assumption_trace": branch_trace,
                                 "state": state(base), "assumed_state": state(branch),
                                 "graph_cells": cells, "palette": palette,
                                 "target": {"cell": cell, "digit": digit, "value": False},
                                 "conclusion": {"cell": cell, "digit": digit,
                                                "value": True, "kind": "placement"},
                                 "full_old_dispatcher": full,
                                 "status": "Discovery only; exact LP and proof certificate required"}
                    candidate["graph_edges"] = [[a, b] for i, a in enumerate(cells)
                                                for b in cells[i + 1:] if sees(a, b)]
                    candidates.append(candidate)
                    counts["full_old_stalls" if full["full_current_finding"] is None
                           else "full_old_finds_step"] += 1
                    print("CANDIDATE", item["id"], cell, digit,
                          "old full", full["full_current_finding"], flush=True)
        except (ValueError, subprocess.TimeoutExpired) as error:
            counts["execution_or_validation_failures"] += 1
            record["status"], record["error"] = "quarantined failure", str(error)
        records.append(record)
        result = {"schema_version": 1, "baseline_commit": "0532d63d9241a69c1fb57e9c34a1533effe1981d",
                  "corpus_sha256": hashlib.sha256(corpus_path.read_bytes()).hexdigest(),
                  "counts": dict(counts), "records": records, "candidates": candidates,
                  "seconds": time.monotonic() - started,
                  "scope": "First N specified corpus entries; joint old nonforcing/local closure; all four-digit palettes, single guardians and entire connected components of 5..24 unplaced cells; uniform non-four-colorability required. Numerical LP is discovery only."}
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2) + "\n")
        print(item["id"], record["status"], dict(counts), flush=True)


if __name__ == "__main__":
    main()
