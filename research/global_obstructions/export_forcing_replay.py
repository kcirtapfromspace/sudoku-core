"""Translate the independently verified seven-node research tree to Rust replay.

Run from the repository root. This is a format translation, not a proof oracle;
the Rust replay checker independently checks every propagation event and leaf.
"""
import hashlib
import itertools
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent


def steps(trace):
    return [{('House' if e['kind'] == 'house-gac' else 'ThreeHouse'): {
        ('sector' if e['kind'] == 'house-gac' else 'center'): e.get('sector', e.get('center')),
        'edits': e['edits']}} for e in trace]


def atoms(center):
    r, c = divmod(center, 9)
    box = {rr*9+cc for rr in range(r//3*3, r//3*3+3) for cc in range(c//3*3, c//3*3+3)}
    row, col = set(range(r*9, r*9+9)), set(range(c, 81, 9))
    return [[center], sorted(row & box - {center}), sorted(col & box - {center}),
            sorted(row-box), sorted(col-box), sorted(box-row-col)]


def export():
    fixture = json.loads((HERE / 'different_lead.json').read_text())
    certificate = json.loads((HERE / 'deeper_forcing_certificate.json').read_text())
    nodes = []

    def convert(n, root=False):
        branches = [{'digit': b['digit'], 'child': convert(b['node'])} for b in n.get('branches', [])]
        leaf = None
        if n['kind'] == 'hall-contradiction':
            leaf = {'Hall': {'sector': n['sector'], 'cells': n['witness']['cells']}}
        if n['kind'] == 'three-house-contradiction':
            partition, blocked = atoms(n['center']), {}
            for digit in range(9):
                dm = 1 << digit
                if not dm & n['domains'][n['center']]:
                    continue
                for aa in itertools.combinations([1 << k for k in range(9) if k != digit], 2):
                    a = aa[0] | aa[1]
                    for bb in itertools.combinations([1 << k for k in range(9) if not (dm | a) & (1 << k)], 2):
                        b = bb[0] | bb[1]
                        sig = [dm, a, b, 511 ^ (dm | a), 511 ^ (dm | b), 511 ^ (dm | a | b)]
                        blocker = next(x for x in n['cover']['blockers']
                                       if bin(x['neighbor_union'] & sig[x['atom']]).count('1') < len(x['cells']))
                        atom, digits = blocker['atom'], sig[blocker['atom']]
                        subset = sum(1 << partition[atom].index(cell) for cell in blocker['cells'])
                        blocked[atom, digits] = {'atom': atom, 'digits': digits, 'cells': subset}
            leaf = {'ThreeHouse': {'center': n['center'], 'blockers': list(blocked.values())}}
        nodes.append({'steps': steps((fixture['assumption_trace'] if root else []) + n['closure_trace']),
                      'leaf': leaf, 'cell': n.get('cell'), 'branches': branches})
        return len(nodes)-1

    convert(certificate['proof'], True)
    state = fixture['state']
    digest = hashlib.sha256(b'sudoku-forcing-v1-native')
    for value, domain in zip(state['values'], state['domains']):
        digest.update(bytes([value]))
        digest.update(domain.to_bytes(2, 'little'))
    proof = {'version': 1, 'state_hash': digest.hexdigest(), 'cell': 7, 'digit': 3, 'value': False, 'nodes': nodes}
    output = HERE.parents[1] / 'tests/fixtures/forcing_mith515.json'
    output.write_text(json.dumps({'state': state, 'proof': proof}, indent=2) + '\n')
    print('Exported', len(nodes), 'nodes and', sum(len(n['steps']) for n in nodes), 'propagation steps to', output)


if __name__ == '__main__':
    export()
