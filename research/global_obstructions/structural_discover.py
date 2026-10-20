"""Reproduce the bounded four-digit synthetic support search.

Discovery uses optional NumPy/SciPy floating-point screening. Its output is not
a certificate. structural_verify.py checks the retained example independently
with exact arithmetic and needs only Python's standard library.
"""

import argparse
import itertools
import json
from pathlib import Path
import random

HERE = Path(__file__).resolve().parent
HOUSES = [set(range(9 * r, 9 * r + 9)) for r in range(9)]
HOUSES += [set(range(c, 81, 9)) for c in range(9)]
HOUSES += [set(9 * (3 * br + r) + 3 * bc + c for r in range(3) for c in range(3))
           for br in range(3) for bc in range(3)]
PEERS = [set().union(*(h for h in HOUSES if c in h)) - {c} for c in range(81)]


def four_colorable(vertices, node_limit=100000):
    """Exploratory DSATUR, with sound unused-color symmetry normalization.

    All initial vertex domains are the same four colors. Therefore colors not
    yet assigned are interchangeable. This optimization is absent from the
    independent full-branch validation and production checker.
    """
    vertices = tuple(sorted(vertices))
    n = len(vertices)
    adjacency = [sum(1 << j for j, v in enumerate(vertices) if v in PEERS[c]) for c in vertices]
    calls = 0
    dead = set()

    def search(domains, assigned, used):
        nonlocal calls
        calls += 1
        if calls > node_limit:
            raise RuntimeError("Coloring budget exhausted; no conclusion may be drawn")
        if assigned == (1 << n) - 1:
            return True
        key = tuple(domains)
        if key in dead:
            return False
        vertex = min((i for i in range(n) if not assigned & (1 << i)),
                     key=lambda i: (bin(domains[i]).count("1"),
                                    -bin(adjacency[i] & ~assigned).count("1")))
        choices = domains[vertex]
        unused = choices & ~used
        if unused:
            choices = (choices & used) | (unused & -unused)
        while choices:
            bit = choices & -choices
            choices -= bit
            updated = domains[:]
            updated[vertex] = bit
            neighbors = adjacency[vertex] & ~assigned
            contradiction = False
            while neighbors:
                neighbor = neighbors & -neighbors
                neighbors -= neighbor
                index = neighbor.bit_length() - 1
                updated[index] &= ~bit
                if not updated[index]:
                    contradiction = True
                    break
            if not contradiction and search(updated, assigned | (1 << vertex), used | bit):
                return True
        dead.add(key)
        return False

    return search([15] * n, 0, 0), calls


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cyclic", action="store_true", help="Run the earlier cyclic-completion negative control, without LP screening")
    parser.add_argument("--budget", type=int, default=48384, help="Maximum number of candidate presentations")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    assert 1 <= args.budget <= 48384
    fixture = json.loads((HERE / "structural_four_digit_13.json").read_text())
    completion = ([(3 * r + r // 3 + c) % 9 + 1 for r in range(9) for c in range(9)]
                  if args.cyclic else list(map(int, fixture["completion"])))
    assert all({completion[c] for c in house} == set(range(1, 10)) for house in HOUSES)
    if not args.cyclic:
        import numpy as np
        from scipy.optimize import linprog
        coefficients = np.zeros((27, 81))
        for i, house in enumerate(HOUSES):
            for c in house:
                coefficients[i, c] = 1
    generator = random.Random(8341)
    tested = graph_obstructions = linear_rejections = 0
    result = None
    for palette in itertools.combinations(range(1, 10), 4):
        base = {c for c, d in enumerate(completion) if d in palette}
        guardians = [c for c, d in enumerate(completion) if d not in palette]
        generator.shuffle(guardians)
        for guardian in guardians[:6]:
            choices = [tuple(sorted(base & house)) for house in HOUSES if guardian in house]
            removals = list(itertools.product(*choices))
            generator.shuffle(removals)
            for removed in removals:
                if tested == args.budget:
                    break
                vertices = base - set(removed) | {guardian}
                assert all(len(vertices & house) <= 4 for house in HOUSES)
                colorable, calls = four_colorable(vertices)
                tested += 1
                if colorable:
                    continue
                graph_obstructions += 1
                if args.cyclic:
                    raise AssertionError("Cyclic negative control unexpectedly had a graph obstruction")
                linear = linprog(np.zeros(81), A_eq=coefficients, b_eq=np.full(27, 4.),
                                 bounds=[(1, 1) if c in vertices else (0, 1) for c in range(81)],
                                 method="highs")
                if not linear.success:
                    assert linear.status == 2, linear.message
                    linear_rejections += 1
                    continue
                # This is a discovery screen, not exact acceptance of the result.
                original = sorted(vertices)
                for cell in sorted(vertices):
                    if not four_colorable(vertices - {cell})[0]:
                        vertices.remove(cell)
                result = {"palette": list(palette), "guardian": guardian,
                          "support": sorted(vertices), "original_support": original,
                          "removed_completion_palette_cells": list(removed),
                          "approximate_role_witness": list(linear.x),
                          "completion": "".join(map(str, completion))}
                break
            if result is not None or tested == args.budget:
                break
        if result is not None or tested == args.budget:
            break
    output = {"schema_version": 1, "seed": 8341, "cyclic_control": args.cyclic,
              "candidate_presentations_tested": tested,
              "non_four_colorable_presentations": graph_obstructions,
              "floating_linear_rejections": linear_rejections,
              "result": result,
              "scope": "126 palettes, six seeded guardian positions per palette, 64 ordered removal choices; presentations may describe duplicate graphs"}
    if args.cyclic and args.budget == 48384:
        assert graph_obstructions == 0
    if not args.cyclic and args.budget >= 45091:
        assert tested == 45091 and graph_obstructions == 9 and linear_rejections == 8
        assert result["support"] == fixture["support"]
    encoded = json.dumps(output, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(encoded)
    print(encoded)


if __name__ == "__main__":
    main()
