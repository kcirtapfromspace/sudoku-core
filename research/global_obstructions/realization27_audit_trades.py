#!/usr/bin/env python3
"""Exhaustively verify the native-27 outside-clue uniqueness obstruction.

Enumerates every labeled four-coloring of the 26 non-guardian cells.  A signature
is the 27 native-house color masks. Equal signatures give an interchangeable
completion on these cells, regardless of all cells outside the support.
"""
import hashlib
import itertools
import json
from collections import Counter, defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent


def peer(a, b):
    ar, ac = divmod(a, 9)
    br, bc = divmod(b, 9)
    return a != b and (ar == br or ac == bc or
                       (ar // 3, ac // 3) == (br // 3, bc // 3))


def houses():
    return ([[9*r+c for c in range(9)] for r in range(9)] +
            [[9*r+c for r in range(9)] for c in range(9)] +
            [[9*r+c for r in range(3*br, 3*br+3)
                         for c in range(3*bc, 3*bc+3)]
             for br in range(3) for bc in range(3)])


def enumerate_colorings(cells, q):
    adjacency = {v: [w for w in cells if peer(v, w)] for v in cells}
    answer, assignment = [], {}

    def visit():
        if len(assignment) == len(cells):
            answer.append(tuple(assignment[c] for c in cells))
            return
        v = max((w for w in cells if w not in assignment),
                key=lambda w: (len({assignment[x] for x in adjacency[w] if x in assignment}),
                               len(adjacency[w]), -w))
        banned = {assignment[x] for x in adjacency[v] if x in assignment}
        for color in range(q):
            if color not in banned:
                assignment[v] = color
                visit()
                del assignment[v]

    visit()
    return answer


def main():
    if not __debug__:
        raise ValueError('Run without -O')
    path = HERE / 'nonrank_graph_native27_candidate.json'
    fixture = json.loads(path.read_text())
    support = fixture['graph_cells']
    guardian = fixture['target']['cell']
    palette = fixture['palette']
    assert guardian == 41 and len(palette) == 4 and len(support) == 27
    assert not enumerate_colorings(support, 4)
    cells = [c for c in support if c != guardian]
    colorings = enumerate_colorings(cells, 4)
    assert len(colorings) == len(set(colorings)) == 1320
    normalized = set()
    groups = defaultdict(list)
    units = houses()
    for colors in colorings:
        renaming = {}
        normalized.add(tuple(renaming.setdefault(d, len(renaming)) for d in colors))
        assignment = dict(zip(cells, colors))
        signature = tuple(sum(1 << assignment[c] for c in unit if c in assignment) for unit in units)
        # Proper native coloring means sum and bitwise OR coincide in every house.
        assert all(len({assignment[c] for c in unit if c in assignment}) ==
                   sum(c in assignment for c in unit) for unit in units)
        groups[signature].append(colors)
    histogram = Counter(map(len, groups.values()))
    assert len(normalized) == 55 and dict(histogram) == {2: 192, 3: 96, 5: 72, 6: 48}
    assert min(histogram) == 2
    completion = list(map(int, fixture['completion']))
    color_index = {d: i for i, d in enumerate(palette)}
    base_colors = tuple(color_index[completion[c]] for c in cells)
    assert base_colors in colorings
    bucket = next(group for group in groups.values() if base_colors in group)
    alternate_colors = next(colors for colors in bucket if colors != base_colors)
    alternate = completion.copy()
    for c, color in zip(cells, alternate_colors):
        alternate[c] = palette[color]
    assert alternate[guardian] == completion[guardian] == fixture['target']['digit']
    assert all(alternate[c] == completion[c] for c in range(81) if c not in cells)
    assert all(set(alternate[c] for c in unit) == set(range(1, 10)) for unit in units)
    difference = [c for c in cells if completion[c] != alternate[c]]
    assert difference
    result = {
        'schema_version': 1, 'verified': True,
        'fixture_sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
        'support': support, 'guardian': guardian, 'palette': palette,
        'four_colorings_of_all_27': 0,
        'labeled_four_colorings_without_guardian': len(colorings),
        'color_permutation_orbits': len(normalized),
        'native_house_signature_count': len(groups),
        'bucket_size_histogram': dict(sorted(histogram.items())),
        'smallest_bucket': min(histogram),
        'completion': ''.join(map(str, completion)),
        'alternate_completion': ''.join(map(str, alternate)),
        'example_trade_cells': difference,
        'scope': 'Every full completion with the 26 cells in the four-digit palette and guardian in palette+extra has a different completion agreeing outside these 26 cells. Thus clues confined outside the 27-cell support cannot give a unique puzzle admitting such a completion. This does not cover clues inside the support or arbitrary unsound extra pencilmark restrictions.'}
    (HERE / 'realization27_audit_trades.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
