#!/usr/bin/env python3
"""Independent, stdlib-only audit of the synthetic native-27 Rust replay.

Rebuilds peer edges from coordinates; does not import production or search code.
The certificate proves a local conclusion, not a givens-only puzzle realization.
"""
import copy
import hashlib
import json
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent


def require(condition, message):
    if not condition:
        raise ValueError(message)


def native(a, b):
    ar, ac = divmod(a, 9)
    br, bc = divmod(b, 9)
    return a != b and (ar == br or ac == bc or
                       (ar // 3, ac // 3) == (br // 3, bc // 3))


def audit(replay, fixture=None):
    require(replay['version'] == 1, 'replay version')
    state, proof = replay['state'], replay['proof']
    require(state['version'] == proof['version'] == 1, 'version')
    require(len(state['values']) == len(state['domains']) == 81, 'state size')
    digest = hashlib.sha256(b'sudoku-global-coloring-v1-native')
    for value, mask in zip(state['values'], state['domains']):
        require(type(value) is int and 0 <= value <= 9, 'value')
        require(type(mask) is int and 0 < mask <= 511, 'domain')
        require(value == 0 or mask == 1 << (value - 1), 'placed singleton')
        digest.update(bytes([value]))
        digest.update(mask.to_bytes(2, 'little'))
    require(digest.hexdigest() == proof['state_hash'], 'state binding')
    require(set(proof) == {'version', 'state_hash', 'cells', 'cell', 'digit',
                           'value', 'witness', 'nodes'}, 'proof fields')
    cells, nodes = proof['cells'], proof['nodes']
    require(0 < len(cells) <= 32 and cells == sorted(set(cells)), 'support')
    require(all(type(c) is int and 0 <= c < 81 for c in cells), 'cell indices')
    require(0 < len(nodes) <= 32768, 'hard node cap')
    target, digit, result = proof['cell'], proof['digit'], proof['value']
    require(target in cells and type(digit) is int and 1 <= digit <= 9 and
            type(result) is bool, 'target')
    bit = 1 << (digit - 1)
    require(state['values'][target] == 0 and state['domains'][target] & bit, 'target availability')
    adjacency = {c: [d for d in cells if native(c, d)] for c in cells}
    witness = proof['witness']
    require(len(witness) == len(cells), 'witness length')
    coloring = dict(zip(cells, witness))
    for c in cells:
        d = coloring[c]
        require(type(d) is int and 1 <= d <= 9 and state['domains'][c] & (1 << (d - 1)), 'base witness domain')
        require(all(coloring[v] != d for v in adjacency[c]), 'base witness peers')
    require((coloring[target] == digit) == result, 'base endpoint')
    if fixture is not None:
        require(state == fixture['state'], 'fixture exact base snapshot')
        require(cells == fixture['graph_cells'], 'fixture support')
        require((target, digit, not result) == (fixture['target']['cell'],
                fixture['target']['digit'], fixture['target']['value']), 'fixture polarity')
        completion = list(map(int, fixture['completion']))
        require(len(completion) == 81, 'completion size')
        for c, d in enumerate(completion):
            require(state['domains'][c] & (1 << (d - 1)), 'completion respects base')
            require(all(completion[v] != d for v in range(c) if native(c, v)), 'native completion')
        require(completion[target] == digit, 'completion conclusion endpoint')
    domains = {c: state['domains'][c] for c in cells}
    domains[target] = domains[target] & ~bit if result else bit
    seen, stats = set(), Counter()
    max_depth = 0

    def visit(index, masks, assigned, depth):
        nonlocal max_depth
        require(type(index) is int and 0 <= index < len(nodes) and index not in seen, 'tree child/reuse')
        seen.add(index)
        max_depth = max(max_depth, depth)
        node = nodes[index]
        require(type(node) is dict and len(node) == 1, 'node')
        if 'Conflict' in node:
            data = node['Conflict']
            require(set(data) == {'cell'} and data['cell'] in masks and masks[data['cell']] == 0, 'empty leaf')
            stats['conflict_leaves'] += 1
            return
        require(set(node) == {'Split'}, 'node type')
        data = node['Split']
        require(set(data) == {'cell', 'branches'}, 'split fields')
        cell, branches = data['cell'], data['branches']
        require(cell in masks and cell not in assigned and masks[cell] != 0, 'pivot')
        wanted = {d for d in range(1, 10) if masks[cell] & (1 << (d - 1))}
        require(type(branches) is list and len(branches) == len(wanted), 'branch count')
        got = []
        stats['split_nodes'] += 1
        stats['singleton_split_nodes' if len(wanted) == 1 else 'branching_split_nodes'] += 1
        for branch in branches:
            require(set(branch) == {'digit', 'child'}, 'branch fields')
            d, child = branch['digit'], branch['child']
            require(type(d) is int and d in wanted and type(child) is int and 0 <= child < index, 'branch/postorder')
            got.append(d)
            next_masks = masks.copy()
            next_masks[cell] = 1 << (d - 1)
            for v in adjacency[cell]:
                next_masks[v] &= ~(1 << (d - 1))
            visit(child, next_masks, assigned | {cell}, depth + 1)
        require(set(got) == wanted and len(set(got)) == len(got), 'complete exact coverage')

    visit(len(nodes) - 1, domains, set(), 0)
    require(len(seen) == len(nodes), 'unreachable node')
    return {'support_cells': len(cells), 'native_edges': sum(map(len, adjacency.values())) // 2,
            'nodes': len(nodes), 'maximum_assignment_depth': max_depth,
            'target': {'cell': target, 'digit': digit, 'result': result},
            'state_hash': digest.hexdigest(), **stats}


def main():
    require(__debug__, 'Run without -O; all checks must execute')
    replay_path = HERE / 'replay_native27.json'
    fixture_path = HERE / 'nonrank_graph_native27_candidate.json'
    replay = json.loads(replay_path.read_text())
    fixture = json.loads(fixture_path.read_text())
    result = audit(replay, fixture)
    corruptions = []
    for key, value in [('state_hash', 'bad'), ('value', False), ('cell', 81), ('digit', 0)]:
        r = copy.deepcopy(replay)
        r['proof'][key] = value
        corruptions.append(r)
    r = copy.deepcopy(replay)
    r['state']['domains'][0] ^= 1
    corruptions.append(r)
    r = copy.deepcopy(replay)
    r['proof']['witness'][0] = 0
    corruptions.append(r)
    split = next(i for i, n in enumerate(replay['proof']['nodes'])
                 if 'Split' in n and len(n['Split']['branches']) > 1)
    leaf = next(i for i, n in enumerate(replay['proof']['nodes']) if 'Conflict' in n)
    r = copy.deepcopy(replay)
    r['proof']['nodes'][split]['Split']['branches'].pop()
    corruptions.append(r)
    r = copy.deepcopy(replay)
    branches = r['proof']['nodes'][split]['Split']['branches']
    branches[1]['digit'] = branches[0]['digit']
    corruptions.append(r)
    r = copy.deepcopy(replay)
    r['proof']['nodes'][split]['Split']['branches'][0]['child'] = split
    corruptions.append(r)
    r = copy.deepcopy(replay)
    r['proof']['nodes'][leaf]['Conflict']['cell'] = 81
    corruptions.append(r)
    r = copy.deepcopy(replay)
    r['proof']['cells'].append(80)
    corruptions.append(r)
    r = copy.deepcopy(replay)
    r['proof']['anchors'] = [31, 40, 41, 50]
    corruptions.append(r)
    for i, corrupt in enumerate(corruptions):
        try:
            audit(corrupt, fixture)
        except (ValueError, KeyError, TypeError, IndexError):
            continue
        raise ValueError('Accepted mutation {}'.format(i))
    result.update({'schema_version': 1, 'verified': True, 'rejected_mutations': len(corruptions),
        'fixture_sha256': hashlib.sha256(fixture_path.read_bytes()).hexdigest(),
        'replay_sha256': hashlib.sha256(replay_path.read_bytes()).hexdigest(),
        'limits': 'Synthetic domains; does not prove a givens-only or unique-puzzle realization.'})
    (HERE / 'realization27_audit_export.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
