#!/usr/bin/env python3
"""Bounded shallow reductions to already-known contradictions.

Assume one additional cell value, close exact native house/crossing triples,
then recognize a two-color odd cycle, Loki tridagon or old13 hyperhole. A
composition is reported only when EVERY candidate of one pivot is covered.
No-hit/budget exhaustion are explicitly inconclusive.
"""
import argparse
import hashlib
import itertools
import json
import time
from collections import deque
from pathlib import Path

from different_audit_known_patterns import masks_from_fixture, peer, tridagon, hyperhole
from real_verify import close, replay_trace


def odd_cycle(domains):
    for digits in itertools.combinations(range(1, 10), 2):
        mask = sum(1 << (d-1) for d in digits)
        cells = [c for c in range(81) if domains[c] & ~mask == 0]
        adjacency = {c: [d for d in cells if peer(c, d)] for c in cells}
        side, parent = {}, {}
        for first in cells:
            if first in side:
                continue
            side[first], parent[first] = 0, None
            queue = deque([first])
            while queue:
                a = queue.popleft()
                for b in adjacency[a]:
                    if b not in side:
                        side[b], parent[b] = side[a] ^ 1, a
                        queue.append(b)
                    elif side[b] == side[a]:
                        left, right = [a], [b]
                        while parent[left[-1]] is not None:
                            left.append(parent[left[-1]])
                        while parent[right[-1]] is not None:
                            right.append(parent[right[-1]])
                        common = next(c for c in left if c in right)
                        cycle = left[:left.index(common)+1] + list(reversed(right[:right.index(common)]))
                        assert len(cycle) % 2 == 1 and len(cycle) == len(set(cycle))
                        assert all(peer(cycle[i], cycle[(i+1) % len(cycle)]) for i in range(len(cycle)))
                        return {'family': 'known-two-color-odd-cycle', 'palette': digits,
                                'cycle_cells': cycle, 'actual_domains': [domains[c] for c in cycle]}
    return None


def shallow(case, seconds=40, max_pivots=32, use_hyperhole=False, preferred_pivots=None):
    domains = masks_from_fixture(case)
    cells = case.get('graph_cells', list(range(81)))
    pivots = sorted((c for c in cells if bin(domains[c]).count('1') > 1),
                    key=lambda c: (bin(domains[c]).count('1'),
                                   -sum(peer(c, d) for d in cells), c))[:max_pivots]
    if preferred_pivots:
        pivots = [c for c in preferred_pivots if c in pivots] + [c for c in pivots if c not in preferred_pivots]
    started, tested = time.monotonic(), []
    for pivot in pivots:
        branches = []
        complete = True
        for digit in range(1, 10):
            bit = 1 << (digit-1)
            if not domains[pivot] & bit:
                continue
            if time.monotonic()-started >= seconds:
                return {'classification': 'inconclusive-budget', 'pivots_tested': tested,
                        'seconds': time.monotonic()-started}
            assumed = domains.copy()
            assumed[pivot] = bit
            try:
                final, trace = close(assumed)
            except ValueError as error:
                if not str(error).startswith(('Hall contradiction', 'Three-house contradiction')):
                    raise
                branches.append({'digit': digit, 'classification': 'local-contradiction',
                                 'contradiction': str(error)})
                continue
            assert replay_trace(assumed, trace) == final
            known = odd_cycle(final)
            if known is None:
                hits, _ = tridagon(final)
                known = hits[0] if hits else None
            if known is None and use_hyperhole:
                known, _, _ = hyperhole(final, budget=10000)
            branch = {'digit': digit, 'classification': 'known-pattern' if known else 'unclassified',
                      'closure_events': len(trace), 'known_pattern': known,
                      'closed_domains_sha256': hashlib.sha256(json.dumps(final, separators=(',',':')).encode()).hexdigest()}
            if known is not None:
                branch.update(closure_trace=trace, closed_domains=final)
            branches.append(branch)
            if known is None:
                complete = False
                break
        tested.append({'pivot': pivot, 'branches': branches})
        if complete:
            expected = [d for d in range(1, 10) if domains[pivot] & (1 << (d-1))]
            assert [b['digit'] for b in branches] == expected
            return {'classification': 'known-shallow-composition', 'pivot': pivot,
                    'pivots_tested': tested, 'seconds': time.monotonic()-started,
                    'limits': 'A standard case split plus exact local propagation/known-pattern lemmas, not an identity of the full starting graph with a named pattern.'}
    return {'classification': 'no-checked-shallow-composition', 'pivots_tested': tested,
            'seconds': time.monotonic()-started,
            'limits': 'Inconclusive: excludes only searched single additional pivot reductions, not all forcing or named-rule compositions.'}


def main():
    if not __debug__:
        raise ValueError('Run without -O')
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('fixture', type=Path)
    parser.add_argument('--seconds', type=float, default=40)
    parser.add_argument('--max-pivots', type=int, default=32)
    parser.add_argument('--hyperhole', action='store_true')
    parser.add_argument('--all-board-pivots', action='store_true', help='Test pivots outside the selected proof support too')
    parser.add_argument('--prefer-cells', default='', help='Comma-separated zero-based pivots tried first')
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    preferred = [int(c) for c in args.prefer_cells.split(',') if c]
    case = json.loads(args.fixture.read_text())
    if args.all_board_pivots:
        case['graph_cells'] = list(range(81))
    result = shallow(case, args.seconds, args.max_pivots, args.hyperhole, preferred)
    result['pivot_scope'] = 'all81 board cells' if args.all_board_pivots else 'selected support cells'
    result['fixture'] = str(args.fixture)
    result['fixture_sha256'] = hashlib.sha256(args.fixture.read_bytes()).hexdigest()
    if args.output:
        args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k:v for k,v in result.items() if k != 'pivots_tested'}, indent=2))


if __name__ == '__main__':
    main()
