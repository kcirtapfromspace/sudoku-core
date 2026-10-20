#!/usr/bin/env python3
"""Reproduce native-27 graph identity and a short research-only refutation.

NetworkX is used for graph invariants/isomorphism; the refutation checker and
fractional-cover arithmetic are explicit. Color normalization is justified by
uniform four-color domains and a native K4, not a production certificate feature.
"""
import hashlib
import itertools
import json
from collections import Counter
from fractions import Fraction
from pathlib import Path

import networkx as nx

HERE = Path(__file__).resolve().parent


def peer(a, b):
    r, c = divmod(a, 9)
    s, d = divmod(b, 9)
    return a != b and (r == s or c == d or (r//3, c//3) == (s//3, d//3))


def canonical(g):
    def refine(partition):
        while True:
            next_partition = []
            for block in partition:
                groups = {}
                for v in block:
                    signature = tuple(sum(w in g[v] for w in other) for other in partition)
                    groups.setdefault(signature, []).append(v)
                next_partition.extend(tuple(groups[s]) for s in sorted(groups))
            if next_partition == partition:
                return partition
            partition = next_partition

    def visit(partition):
        partition = refine(partition)
        if all(len(block) == 1 for block in partition):
            order = tuple(block[0] for block in partition)
            code = ''.join('1' if g.has_edge(order[i], order[j]) else '0'
                           for j in range(1, len(g)) for i in range(j))
            return code, order
        k = next(i for i, block in enumerate(partition) if len(block) > 1)
        return min(visit(partition[:k] + [(v,), tuple(w for w in partition[k] if w != v)]
                         + partition[k+1:]) for v in partition[k])

    return visit([tuple(v for v in g if g.degree(v) == degree)
                  for degree in sorted(set(dict(g.degree()).values()))])


def coloring(g, q):
    assigned = {}

    def visit():
        if len(assigned) == len(g):
            return dict(assigned)
        v = max((w for w in g if w not in assigned), key=lambda w:
                (len({assigned[x] for x in g[w] if x in assigned}), g.degree(w), -w))
        for d in range(q):
            if all(assigned.get(x) != d for x in g[v]):
                assigned[v] = d
                answer = visit()
                if answer is not None:
                    return answer
                del assigned[v]
        return None

    return visit()


def compact_tree(g, anchors):
    initial = dict.fromkeys(g, 15)
    assigned = set()

    def assign(domains, cell, color):
        domains[cell] = 1 << color
        for other in g[cell]:
            domains[other] &= ~(1 << color)

    for color, cell in enumerate(anchors):
        assign(initial, cell, color)
        assigned.add(cell)

    def search(domains, assigned):
        forced = []
        while not any(mask == 0 for mask in domains.values()):
            single = next((v for v in g if v not in assigned and
                           bin(domains[v]).count('1') == 1), None)
            if single is None:
                break
            color = domains[single].bit_length() - 1
            forced.append({'cell': single, 'color': color})
            assign(domains, single, color)
            assigned.add(single)
        empty = next((v for v in g if domains[v] == 0), None)
        if empty is not None:
            return {'forced': forced, 'conflict': empty}
        remaining = [v for v in g if v not in assigned]
        if not remaining:
            raise ValueError('Unexpected four-coloring')
        pivot = min(remaining, key=lambda v: (bin(domains[v]).count('1'), -g.degree(v), v))
        branches = []
        for color in range(4):
            if domains[pivot] & (1 << color):
                next_domains = domains.copy()
                assign(next_domains, pivot, color)
                branches.append({'color': color, 'child': search(next_domains, assigned | {pivot})})
        return {'forced': forced, 'pivot': pivot, 'branches': branches}

    proof = {'schema': 'research-normalized-four-color-tree-v1', 'anchors': anchors,
             'tree': search(initial, assigned)}
    # Independent replay: no MRV decisions or forced-step discovery is trusted.
    stats = Counter()

    def check(node, masks, assigned, depth):
        stats['nodes'] += 1
        stats['maximum_branch_depth'] = max(stats['maximum_branch_depth'], depth)
        for step in node['forced']:
            c, d = step['cell'], step['color']
            assert c in masks and c not in assigned and masks[c] == 1 << d
            assigned.add(c)
            for v in g[c]:
                masks[v] &= ~(1 << d)
            stats['forced_assignments'] += 1
        if 'conflict' in node:
            assert masks[node['conflict']] == 0
            stats['conflict_leaves'] += 1
            return
        c = node['pivot']
        assert c not in assigned and masks[c]
        wanted = {d for d in range(4) if masks[c] & (1 << d)}
        got = [b['color'] for b in node['branches']]
        assert set(got) == wanted and len(got) == len(wanted)
        stats['branching_nodes'] += 1
        for branch in node['branches']:
            d = branch['color']
            following = masks.copy()
            following[c] = 1 << d
            for v in g[c]:
                following[v] &= ~(1 << d)
            check(branch['child'], following, assigned | {c}, depth + 1)

    assert len(anchors) == 4 and all(g.has_edge(a, b) for a, b in itertools.combinations(anchors, 2))
    check(proof['tree'], initial, assigned.copy(), 0)
    proof['verified_stats'] = dict(stats)
    return proof


def main():
    if not __debug__:
        raise ValueError('Run without -O')
    fixture = json.loads((HERE / 'nonrank_graph_native27_candidate.json').read_text())
    cells = fixture['graph_cells']
    g = nx.Graph()
    g.add_nodes_from(cells)
    g.add_edges_from((a, b) for a, b in itertools.combinations(cells, 2) if peer(a, b))
    cover = dict.fromkeys(cells, Fraction(0))
    total = Fraction(0)
    for term in fixture['fractional_core_cover']:
        stable = term['cells']
        weight = Fraction(*term['weight'])
        assert weight > 0 and all(not peer(a, b) for a, b in itertools.combinations(stable, 2))
        total += weight
        for c in stable:
            cover[c] += weight
    assert total == 4 and all(x == 1 for x in cover.values())
    maximal_cliques = list(nx.find_cliques(g))
    omega = max(map(len, maximal_cliques))
    assert omega == 4
    proof = compact_tree(g, [31, 40, 41, 50])
    five = coloring(g, 5)
    assert five is not None
    deletion_witnesses = []
    for c in cells:
        witness = coloring(g.subgraph([v for v in cells if v != c]), 4)
        assert witness is not None
        deletion_witnesses.append({'removed_cell': c, 'four_coloring': witness})
    code, order = canonical(g)
    claws = [[c] + list(triple) for c in cells for triple in itertools.combinations(g[c], 3)
             if all(not g.has_edge(a, b) for a, b in itertools.combinations(triple, 2))]
    result = {'schema_version': 1, 'vertices': len(g), 'edges': g.number_of_edges(),
              'degree_histogram': dict(sorted(Counter(dict(g.degree()).values()).items())),
              'omega': omega, 'K4_count': sum(1 for x in nx.enumerate_all_cliques(g) if len(x) == 4),
              'alpha': max(map(len, nx.find_cliques(nx.complement(g)))),
              'chi': 5, 'chi_fractional_exact': [4, 1],
              'automorphisms': sum(1 for _ in nx.algorithms.isomorphism.GraphMatcher(g, g).isomorphisms_iter()),
              'native_graph_vertex_connectivity': nx.node_connectivity(g),
              'induced_claws': len(claws), 'example_induced_claw': claws[0],
              'true_twin_pairs': [[a, b] for a, b in itertools.combinations(cells, 2)
                                  if set(g[a]) | {a} == set(g[b]) | {b}],
              'graph6_in_support_order': nx.to_graph6_bytes(g, nodes=cells, header=False).decode().strip(),
              'canonical_convention': 'Full degree-initialized ordered individualization/refinement, binary upper triangle; not nauty canonical form.',
              'canonical_order_board_cells': order,
              'canonical_binary_sha256': hashlib.sha256(code.encode()).hexdigest(),
              'five_coloring': five, 'vertex_deletion_four_colorings': deletion_witnesses,
              'compact_research_proof_stats': proof['verified_stats'],
              'limits': 'Invariant mismatches exclude particular named identities, not historical novelty. Normalized compact proof is research-only; production replays every color explicitly.'}
    (HERE / 'realization27_prior_art_invariants.json').write_text(json.dumps(result, indent=2) + '\n')
    (HERE / 'realization27_prior_art_compact_proof.json').write_text(json.dumps(proof, indent=2) + '\n')
    print(json.dumps({k: v for k, v in result.items() if k not in ['vertex_deletion_four_colorings', 'five_coloring']}, indent=2))


if __name__ == '__main__':
    main()
