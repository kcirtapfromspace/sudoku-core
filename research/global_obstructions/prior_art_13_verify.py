#!/usr/bin/env python3
"""Small exact CI check for the 13-cell motif's established-family reduction.

Standard library only. This verifies a mathematical mapping and current fixture
premises; it cannot verify a historical first-appearance claim.
"""
import itertools
import json
from pathlib import Path

if not __debug__:
    raise SystemExit('prior_art_13_verify.py must run without Python -O')

HERE = Path(__file__).resolve().parent
SUPPORT = [9, 15, 16, 17, 25, 27, 34, 36, 43, 46, 47, 51, 53]
BAGS = [[9], [15, 16, 17], [25], [34, 43], [51, 53], [46, 47], [27, 36]]
EXTRA = {(15, 51), (16, 34), (16, 43), (17, 53), (27, 34), (36, 43)}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def edge(a, b):
    return tuple(sorted((a, b)))


def peer(a, b):
    r, c = divmod(a, 9)
    s, d = divmod(b, 9)
    return a != b and (r == s or c == d or (r // 3, c // 3) == (s // 3, d // 3))


def houses():
    return ([[r * 9 + c for c in range(9)] for r in range(9)]
            + [[r * 9 + c for r in range(9)] for c in range(9)]
            + [[r * 9 + c for r in range(br * 3, br * 3 + 3)
                for c in range(bc * 3, bc * 3 + 3)]
               for br in range(3) for bc in range(3)])


def verify(fixture=None):
    if fixture is None:
        fixture = json.loads((HERE / 'structural_four_digit_13.json').read_text())
    require(fixture['support'] == SUPPORT, 'fixture support changed')
    require(len(set(SUPPORT)) == 13, 'support must be thirteen distinct cells')
    require(sorted(itertools.chain.from_iterable(BAGS)) == SUPPORT,
            'seven bags must partition support exactly')
    require(list(map(len, BAGS)) == [1, 3, 1, 2, 2, 2, 2], 'bag sizes changed')
    native = {edge(a, b) for a, b in itertools.combinations(SUPPORT, 2) if peer(a, b)}
    skeleton = set()
    for bag in BAGS:
        skeleton.update(edge(a, b) for a, b in itertools.combinations(bag, 2))
    for i, bag in enumerate(BAGS):
        skeleton.update(edge(a, b) for a in bag for b in BAGS[(i + 1) % 7])
    require(len(native) == 35 and len(skeleton) == 29, 'edge counts changed')
    require(skeleton <= native, 'a hyperhole edge is not a native peer edge')
    require(native - skeleton == EXTRA, 'six extra native edges do not match')

    # The 7-hyperhole's stable sets choose nonconsecutive bags. Check the actual
    # 13-vertex graph as well, without relying on the literature's formula.
    triple = [9, 25, 46]
    require(all(edge(a, b) not in native for a, b in itertools.combinations(triple, 2)),
            'independent triple lower-bound witness invalid')
    four_sets = list(itertools.combinations(SUPPORT, 4))
    require(len(four_sets) == 715, 'four-set enumeration incomplete')
    require(all(any(edge(a, b) in skeleton for a, b in itertools.combinations(s, 2))
                for s in four_sets), 'skeleton has an independent four-set')
    triples = [s for s in itertools.combinations(SUPPORT, 3)
               if all(edge(a, b) not in native for a, b in itertools.combinations(s, 2))]
    require(len(triples) == 32, 'native independent-triple count changed')

    completion = [int(d) for d in fixture['completion']]
    require(len(completion) == 81, 'completion must have 81 digits')
    require(all({completion[v] for v in h} == set(range(1, 10)) for h in houses()),
            'completion is not a standard Sudoku solution')
    require(len({completion[v] for v in SUPPORT}) == 5, 'completion is not a five-color witness')
    require(all(completion[a] != completion[b] for a, b in native),
            'completion fails a native edge')
    require(all(fixture['base_domains'][v] & (1 << (completion[v] - 1)) for v in range(81)),
            'completion violates a base domain')

    palette = fixture['palette']
    require(palette == [4, 6, 7, 8], 'palette changed')
    palette_mask = sum(1 << (d - 1) for d in palette)
    target = fixture['target']
    require(target['cell'] == 16 and target['digit'] == 1 and target['value'] is False,
            'guardian-false assumption changed')
    require(fixture['base_domains'][16] == palette_mask | 1, 'guardian domain changed')
    require(all(fixture['base_domains'][v] == palette_mask for v in SUPPORT if v != 16),
            'non-guardian support base domain changed')
    assumed = fixture['assumed_domains']
    require(len(assumed) == 81 and all(type(m) is int and 0 < m <= 511 for m in assumed),
            'assumed domain snapshot invalid')
    require(all(assumed[v] == palette_mask for v in SUPPORT),
            'selected assumed domains must be the common four-digit palette')

    # Exact strong-link availability in the *native houses*, rather than merely
    # in our selected support. These depend on the guardian-false snapshot.
    for box, left, right in [(2, BAGS[1], BAGS[2]),
                             (5, BAGS[3], BAGS[4]),
                             (3, BAGS[5], BAGS[6])]:
        native_house = houses()[18 + box]
        for digit in palette:
            positions = {v for v in native_house if assumed[v] & (1 << (digit - 1))}
            require(positions == set(left + right), f'box{box + 1} lacks exhaustive grouped link')
            require(set(left).isdisjoint(right), 'strong groups overlap')
    # Four weak links: row2, column8, row6, column1 respectively.
    for left, right, house_index in [(BAGS[0], BAGS[1], 1),
                                      (BAGS[2], BAGS[3], 9 + 7),
                                      (BAGS[4], BAGS[5], 5),
                                      (BAGS[6], BAGS[0], 9)]:
        require(set(left + right) <= set(houses()[house_index]), 'weak groups lack stated native house')
        require(all(peer(a, b) for a in left for b in right), 'weak grouped link not complete')

    # Independently enumerate the seven group propositions. A group is true
    # exactly when it contains the digit. Box pairs are exhaustive and disjoint,
    # so their group truth values differ. Native weak links forbid both true.
    link_models = []
    for values in itertools.product([False, True], repeat=7):
        strong = all(values[a] != values[b] for a, b in [(1, 2), (3, 4), (5, 6)])
        weak = all(not (values[a] and values[b]) for a, b in [(0, 1), (2, 3), (4, 5), (6, 0)])
        if strong and weak:
            link_models.append(values)
    require(bool(link_models), 'link implication would be vacuous')
    require(all(values[0] is False for values in link_models),
            'grouped loop failed to eliminate the start candidate')
    return {'vertices': 13, 'native_edges': 35, 'hyperhole_edges': 29,
            'alpha': 3, 'chi': 5, 'four_sets_checked': 715,
            'palette_digits_with_grouped_loop': palette,
            'seven_group_boolean_models_checked': 128,
            'historical_priority': 'not established by this verifier'}


if __name__ == '__main__':
    print(json.dumps(verify(), sort_keys=True))
