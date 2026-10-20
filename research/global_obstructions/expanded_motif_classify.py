"""Exactly classify retained native graphs; optionally certify synthetic LP.

The optional SciPy solve proposes rational role weights. Every proposed weight,
bound and native equality is then checked with Fraction, so solver status alone
never establishes the positive LP claim. This does not certify givens-derived
old-method stalling or uniqueness.
"""
import os
for _name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_name] = "1"
import argparse
from fractions import Fraction
from functools import lru_cache
import itertools
import json
from pathlib import Path

from structural_discover import HOUSES, PEERS, four_colorable

HERE = Path(__file__).resolve().parent


def classify(motif, exact_role=False):
    cells = motif["support"]
    completion = list(map(int, motif["completion"]))
    palette = motif["palette"]
    guardian = motif["guardian"]
    assert len(cells) == len(set(cells)) and cells == sorted(cells)
    assert len(palette) == 4 and len(set(palette)) == 4
    assert all({completion[c] for c in h} == set(range(1,10)) for h in HOUSES)
    assert guardian in cells and completion[guardian] not in palette
    assert all(completion[c] in palette for c in cells if c != guardian)
    adjacency = [sum(1 << j for j,c in enumerate(cells) if c in PEERS[v]) for v in cells]

    @lru_cache(None)
    def alpha(mask):
        if not mask:
            return 0
        vertex = max((i for i in range(len(cells)) if mask & (1 << i)),
                     key=lambda i: bin(adjacency[i] & mask).count("1"))
        without = mask & ~(1 << vertex)
        return max(alpha(without), 1+alpha(without & ~adjacency[vertex]))

    independence = alpha((1 << len(cells))-1)
    colorable, calls = four_colorable(cells, node_limit=1000000)
    assert colorable is False
    assert all(four_colorable(set(cells)-{c}, node_limit=1000000)[0] for c in cells)
    result = {"motif_id":motif["id"], "source_id":motif["source_id"],
              "support_size":len(cells), "native_edges":sum(bin(a).count("1") for a in adjacency)//2,
              "alpha":independence, "four_digit_rank_gap":len(cells)-4*independence,
              "chromatic_number":5, "vertex_critical":True,
              "exploratory_normalized_color_calls":calls,
              "maximum_house_support":max(len(set(cells)&h) for h in HOUSES),
              "support":cells, "palette":palette, "guardian":guardian,
              "completion":motif["completion"]}
    assert result["native_edges"] == motif["edges"]
    if exact_role:
        import numpy as np
        from scipy.optimize import linprog
        matrix = np.zeros((27,81), dtype=np.int8)
        for i,h in enumerate(HOUSES):
            for c in h:
                matrix[i,c] = 1
        answer = linprog(np.zeros(81), A_eq=matrix, b_eq=np.full(27,4),
                         bounds=[(1,1) if c in cells else (0,1) for c in range(81)], method="highs-ds")
        assert answer.success
        role = [Fraction(float(t)).limit_denominator(100000) for t in answer.x]
        assert all(0 <= t <= 1 for t in role)
        assert all(role[c] == 1 for c in cells)
        assert all(sum(role[c] for c in h) == 4 for h in HOUSES)
        point = [[role[c]/4 if d in palette else (1-role[c])/5
                  for d in range(1,10)] for c in range(81)]
        assert all(sum(row) == 1 for row in point)
        assert all(sum(point[c][d] for c in h) == 1 for h in HOUSES for d in range(9))
        assert all(point[c][d-1] == 0 for c in cells for d in range(1,10) if d not in palette)
        assert point[guardian][completion[guardian]-1] == 0
        result["exact_synthetic_role_lp"] = {"native_equalities":324,
            "role_weights":[[t.numerator,t.denominator] for t in role],
            "nonzero_variables":sum(t != 0 for row in point for t in row)}
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inputs", type=Path, nargs="+", help="expanded search metadata files")
    parser.add_argument("--exact-role-lp", action="store_true")
    parser.add_argument("--output", type=Path, default=HERE/"expanded_motif_classification.json")
    args = parser.parse_args()
    results = []
    for path in args.inputs:
        for motif in json.loads(path.read_text())["motifs"]:
            result = classify(motif, args.exact_role_lp)
            result["source_artifact"] = str(path)
            results.append(result)
    output = {"schema_version":1,"scope":"Exact synthetic graph and optional native LP claims; no unique-puzzle or old-method-stall acceptance", "motifs":results}
    args.output.write_text(json.dumps(output,separators=(",",":"))+"\n")
    print(json.dumps({"motifs":len(results),"capacity_rank_obstructions":sum(r["four_digit_rank_gap"]>0 for r in results),"exact_role_lp_verified":sum("exact_synthetic_role_lp" in r for r in results)}))


if __name__ == "__main__":
    main()
