#!/usr/bin/env python3
"""Reproduce the exact 13-cell motif's graph invariants and named comparisons.

Uses NetworkX for isomorphism only; graph construction, independence counts,
coloring, and the canonical individualization/refinement key are explicit.
This is a characterization aid, not a historical-novelty proof.
"""
import hashlib
import itertools
import json
from pathlib import Path

import networkx as nx

HERE = Path(__file__).resolve().parent
CELLS = [9, 15, 16, 17, 25, 27, 34, 36, 43, 46, 47, 51, 53]


def peer(a, b):
    r, c = divmod(a, 9)
    s, d = divmod(b, 9)
    return a != b and (r == s or c == d or (r // 3, c // 3) == (s // 3, d // 3))


def graph(cells):
    g = nx.Graph()
    g.add_nodes_from(range(len(cells)))
    g.add_edges_from((i, j) for i in g for j in g if i < j and peer(cells[i], cells[j]))
    return g


def coloring(g, q):
    colors = {}

    def visit():
        if len(colors) == len(g):
            return dict(colors)
        v = max((v for v in g if v not in colors),
                key=lambda v: (len({colors[w] for w in g[v] if w in colors}), g.degree(v), -v))
        banned = {colors[w] for w in g[v] if w in colors}
        for color in range(q):
            if color not in banned:
                colors[v] = color
                answer = visit()
                if answer is not None:
                    return answer
                del colors[v]
        return None

    return visit()


def canonical(g):
    """Isomorphism-invariant minimum among complete IR leaves.

    Refinement orders degree/signature classes canonically; individualization
    branches over every member, never prunes nonisomorphic possibilities.
    This is a canonical convention, not nauty's particular convention.
    """
    def refine(partition):
        while True:
            refined = []
            for block in partition:
                groups = {}
                for v in block:
                    sig = tuple(sum(w in g[v] for w in other) for other in partition)
                    groups.setdefault(sig, []).append(v)
                refined.extend(tuple(groups[s]) for s in sorted(groups))
            if refined == partition:
                return refined
            partition = refined

    def visit(partition):
        partition = refine(partition)
        if all(len(p) == 1 for p in partition):
            order = tuple(p[0] for p in partition)
            code = ''.join('1' if g.has_edge(order[i], order[j]) else '0'
                           for j in range(1, len(g)) for i in range(j))
            return code, order
        k = next(i for i, p in enumerate(partition) if len(p) > 1)
        return min(visit(partition[:k] + [(v,), tuple(w for w in partition[k] if w != v)]
                         + partition[k + 1:]) for v in partition[k])

    degrees = sorted(set(dict(g.degree()).values()))
    return visit([tuple(v for v in g if g.degree(v) == d) for d in degrees])


def compare(g):
    candidates = {}
    # C13^3, sometimes called a web / cycle power (notation differs by author).
    power = nx.Graph()
    power.add_nodes_from(range(13))
    power.add_edges_from((i, j) for i in range(13) for j in range(i + 1, 13)
                         if min((i - j) % 13, (j - i) % 13) <= 3)
    candidates['cycle_13_power_3'] = power
    candidates['paley_13_undirected'] = nx.Graph(nx.paley_graph(13))
    out = {}
    for name, h in candidates.items():
        matcher = nx.algorithms.isomorphism.GraphMatcher(g, h)
        match = matcher.is_isomorphic()
        out[name] = {'vertices': len(h), 'edges': h.number_of_edges(),
                     'degree_sequence': sorted(dict(h.degree()).values()),
                     'isomorphic': match,
                     'mapping_if_isomorphic': dict(matcher.mapping) if match else None}
    return out


def geometry_canonical(cells):
    """Minimum sorted support over all 3,359,232 standard spatial symmetries.

    Digit renaming is immaterial to support. Row/column maps may permute bands
    or stacks and rows/columns within them. Transposition exchanges the axes.
    Map entries are OLD index -> NEW index (zero based).
    """
    perms = list(itertools.permutations(range(3)))
    maps = []
    for outer in perms:
        for inside in itertools.product(perms, repeat=3):
            maps.append(tuple(3 * outer[r // 3] + inside[r // 3][r % 3]
                              for r in range(9)))
    best = None
    witness = None
    for transpose in [False, True]:
        pairs = [divmod(v, 9)[::-1] if transpose else divmod(v, 9) for v in cells]
        for rows in maps:
            rparts = [(9 * rows[r], c) for r, c in pairs]
            for cols in maps:
                support = tuple(sorted(r + cols[c] for r, c in rparts))
                if best is None or support < best:
                    best = support
                    witness = {'transpose_first': transpose,
                               'row_old_to_new': rows, 'column_old_to_new': cols}
    return {'canonical_support': best, 'mapping': witness,
            'symmetries_considered': 2 * len(maps) ** 2}


def hyperhole_reduction(g):
    bags = [[9], [15, 16, 17], [25], [34, 43], [51, 53], [46, 47], [27, 36]]
    idx = {v: i for i, v in enumerate(CELLS)}
    h = nx.Graph()
    h.add_nodes_from(g)
    for bag in bags:
        h.add_edges_from((idx[a], idx[b]) for a, b in itertools.combinations(bag, 2))
    for i, bag in enumerate(bags):
        h.add_edges_from((idx[a], idx[b]) for a in bag for b in bags[(i + 1) % 7])
    assert all(g.has_edge(a, b) for a, b in h.edges())
    assert coloring(h, 4) is None and coloring(h, 5) is not None
    edge_critical = []
    for a, b in h.edges():
        k = h.copy()
        k.remove_edge(a, b)
        witness = coloring(k, 4)
        assert witness is not None
        edge_critical.append({'removed_cells': [CELLS[a], CELLS[b]],
                              'four_coloring': {CELLS[v]: c + 1 for v, c in witness.items()}})
    # Verify all strong/weak grouped links on the actual guardian-false snapshot.
    fixture = json.loads((HERE / 'structural_four_digit_13.json').read_text())
    domains = fixture['assumed_domains']
    palette = fixture['palette']
    strong_boxes = [(2, bags[1] + bags[2]), (5, bags[3] + bags[4]),
                    (3, bags[5] + bags[6])]
    for box, positions in strong_boxes:
        box_cells = [r * 9 + c for r in range(box // 3 * 3, box // 3 * 3 + 3)
                     for c in range(box % 3 * 3, box % 3 * 3 + 3)]
        for digit in palette:
            assert {v for v in box_cells if domains[v] & (1 << (digit - 1))} == set(positions)
    return {'bags_as_board_cells_in_cyclic_order': bags,
            'bag_sizes': list(map(len, bags)), 'edges': h.number_of_edges(),
            'graph6_in_support_order': nx.to_graph6_bytes(h, header=False).decode().strip(),
            'extra_native_edges': [[CELLS[a], CELLS[b]] for a, b in sorted(g.edges())
                                   if not h.has_edge(a, b)],
            'edge_critical': True, 'edge_deletion_four_colorings': edge_critical,
            'grouped_loop_actual_strong_box_premises_verified': True,
            'five_coloring_as_board_cells': {CELLS[v]: c + 1 for v, c in coloring(h, 5).items()}}


def main():
    g = graph(CELLS)
    try:
        nx.inverse_line_graph(g)
        is_line_graph = True
    except nx.NetworkXError:
        is_line_graph = False
    code, order = canonical(g)
    independent = [s for k in range(len(g) + 1) for s in itertools.combinations(g, k)
                   if all(not g.has_edge(a, b) for a, b in itertools.combinations(s, 2))]
    chi = next(q for q in range(1, 14) if coloring(g, q) is not None)
    edge_deletions = []
    for a, b in sorted(g.edges()):
        h = g.copy()
        h.remove_edge(a, b)
        edge_deletions.append({'cells': [CELLS[a], CELLS[b]],
                               'four_colorable': coloring(h, 4) is not None})
    # Vertex-criticality includes explicit 4-colorings of all 13 deletions.
    vertex_deletions = []
    for v in g:
        h = g.copy()
        h.remove_node(v)
        answer = coloring(h, 4)
        vertex_deletions.append({'removed_cell': CELLS[v],
                                'four_coloring': {CELLS[w]: c + 1 for w, c in answer.items()}
                                if answer is not None else None})
    autos = list(nx.algorithms.isomorphism.GraphMatcher(g, g).isomorphisms_iter())
    orbits = []
    unseen = set(g)
    while unseen:
        v = min(unseen)
        orbit = sorted({a[v] for a in autos})
        orbits.append([CELLS[w] for w in orbit])
        unseen.difference_update(orbit)
    report = {
        'schema_version': 1,
        'cells_zero_based': CELLS,
        'cells_human': [f'r{a // 9 + 1}c{a % 9 + 1}' for a in CELLS],
        'edges_as_board_cells': [[CELLS[a], CELLS[b]] for a, b in sorted(g.edges())],
        'vertices': len(g), 'edges': g.number_of_edges(),
        'degree_sequence': sorted(dict(g.degree()).values()),
        'independence_polynomial_coefficients': [sum(len(s) == k for s in independent)
                                                for k in range(4)],
        'alpha': max(map(len, independent)), 'chi': chi,
        'omega': max(map(len, nx.find_cliques(g))),
        'maximal_cliques': [[CELLS[v] for v in sorted(c)] for c in nx.find_cliques(g)],
        'triangle_count': sum(nx.triangles(g).values()) // 3,
        'automorphism_count': len(autos), 'automorphism_orbits_as_cells': orbits,
        'graph6_in_support_order': nx.to_graph6_bytes(g, header=False).decode().strip(),
        'canonical_adjacency_bits_IR_v1': code,
        'canonical_order_as_board_cells_IR_v1': [CELLS[v] for v in order],
        'canonical_bits_sha256': hashlib.sha256(code.encode()).hexdigest(),
        'vertex_deletions': vertex_deletions,
        'edge_deletions': edge_deletions,
        'edge_critical': all(row['four_colorable'] for row in edge_deletions),
        'is_line_graph': is_line_graph,
        'comparisons': compare(g),
        'hyperhole_spanning_reduction': hyperhole_reduction(g),
        'native_geometry': {
            'row_occupancies': [sum(v // 9 == r for v in CELLS) for r in range(9)],
            'column_occupancies': [sum(v % 9 == c for v in CELLS) for c in range(9)],
            'box_occupancies': [sum((v // 27) * 3 + v % 9 // 3 == b for v in CELLS)
                                for b in range(9)],
            **geometry_canonical(CELLS),
        },
    }
    path = HERE / 'prior_art_13_invariants.json'
    path.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({k: report[k] for k in ['vertices', 'edges', 'degree_sequence', 'alpha',
                     'chi', 'omega', 'triangle_count', 'automorphism_count', 'automorphism_orbits_as_cells',
                     'edge_critical', 'graph6_in_support_order', 'canonical_bits_sha256', 'comparisons']}, indent=2))


if __name__ == '__main__':
    main()
