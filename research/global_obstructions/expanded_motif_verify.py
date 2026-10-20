"""Standard-library exact checks of expanded synthetic obstruction witnesses.

No old-engine execution, givens-puzzle uniqueness, or historical originality
claim follows from this certificate. All native edges and maximal independent
sets are reconstructed, not trusted from a discovery program.
"""
from fractions import Fraction
import hashlib
import itertools
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
HOUSES = [set(range(9*r,9*r+9)) for r in range(9)]
HOUSES += [set(range(c,81,9)) for c in range(9)]
HOUSES += [set((3*br+r)*9+3*bc+c for r in range(3) for c in range(3))
           for br in range(3) for bc in range(3)]


def sees(a,b):
    return a//9 == b//9 or a%9 == b%9 or (a//27,a%9//3) == (b//27,b%9//3)


def maximal_independent_sets(cells):
    size = len(cells)
    all_vertices = (1 << size)-1
    # A clique in the complement is an independent set in the native graph.
    complement = [sum(1 << j for j,b in enumerate(cells) if a != b and not sees(a,b)) for a in cells]
    answer = []

    def walk(chosen, possible, excluded):
        if not possible and not excluded:
            answer.append(sorted(cells[j] for j in range(size) if chosen & (1 << j)))
            return
        union = possible|excluded
        pivot = max((j for j in range(size) if union & (1 << j)),
                    key=lambda j: bin(complement[j]&possible).count("1")) if union else 0
        branches = possible & ~complement[pivot]
        while branches:
            bit = branches & -branches
            branches ^= bit
            j = bit.bit_length()-1
            walk(chosen|bit, possible & complement[j], excluded & complement[j])
            possible ^= bit
            excluded |= bit

    walk(0,all_vertices,0)
    return sorted(answer)


def verify():
    if not __debug__:
        raise SystemExit("Run without -O: this exact verifier uses assertions")
    cases = json.loads((HERE/"expanded_motif_classification.json").read_text())["motifs"]
    role_count = 0
    for case in cases:
        cells, palette = case["support"], case["palette"]
        sets = maximal_independent_sets(cells)
        assert max(map(len,sets)) == case["alpha"]
        completion = list(map(int,case["completion"]))
        assert all({completion[c] for c in h} == set(range(1,10)) for h in HOUSES)
        assert all(completion[c] in palette for c in cells if c != case["guardian"])
        assert completion[case["guardian"]] not in palette
        witness = case.get("exact_synthetic_role_lp")
        if witness:
            role = [Fraction(*pair) for pair in witness["role_weights"]]
            assert len(role) == 81 and all(0 <= t <= 1 for t in role)
            assert all(role[c] == 1 for c in cells)
            point = [[role[c]/4 if d in palette else (1-role[c])/5 for d in range(1,10)] for c in range(81)]
            assert all(sum(row) == 1 for row in point)
            assert all(sum(point[c][d] for c in h) == 1 for h in HOUSES for d in range(9))
            assert all(point[c][d-1] == 0 for c in cells for d in range(1,10) if d not in palette)
            assert point[case["guardian"]][completion[case["guardian"]]-1] == 0
            role_count += 1
    fractional = json.loads((HERE/"expanded_motif_fractional_coloring.json").read_text())["motifs"]
    for case in fractional:
        cells = case["support"]
        sets = maximal_independent_sets(cells)
        assert len(sets) == case["maximal_independent_sets_count"]
        assert hashlib.sha256(json.dumps(sets,separators=(",",":")).encode()).hexdigest() == case["maximal_independent_sets_sha256"]
        weights = [Fraction(*p) for p in case["dual_vertex_weights"]]
        assert len(weights) == len(cells) and all(w >= 0 for w in weights)
        by_cell = dict(zip(cells,weights))
        assert all(sum(by_cell[c] for c in independent) <= 1 for independent in sets)
        cover = [(entry["cells"],Fraction(*entry["weight"])) for entry in case["fractional_cover"]]
        assert all(w >= 0 and set(st) <= set(cells) and
                   all(not sees(a,b) for a,b in itertools.combinations(st,2)) for st,w in cover)
        assert all(sum(w for st,w in cover if c in st) >= 1 for c in cells)
        assert sum(w for st,w in cover) == sum(weights) == Fraction(*case["fractional_chromatic_number"]) == Fraction(25,6)
        assert sum(weights) > 4
    weighted_supports = {tuple(case["support"]) for case in fractional}
    for case in cases:
        assert case["support_size"] > 4*case["alpha"] or tuple(case["support"]) in weighted_supports
        # The lower bound is the checked rank/weighted-rank argument; the
        # respecting completion uses at most the four palette colors plus the
        # guardian's extra color, providing the independent upper bound.
        assert len({int(case["completion"][c]) for c in case["support"]}) <= 5
        assert case["chromatic_number"] == 5
    print(json.dumps({"synthetic_native_lp_witnesses":role_count,
                      "exact_independence_numbers":len(cases),
                      "exact_fractional_chromatic_25_over_6":len(fractional),
                      "exact_chromatic_number_5":len(cases),
                      "unique_puzzle_acceptance":False}))


if __name__ == "__main__":
    verify()
