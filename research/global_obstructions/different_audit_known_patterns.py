#!/usr/bin/env python3
"""Bounded, exact recognizer for two already-established Sudoku obstructions.

Recognizes the published 12-cell tridagon by a native graph isomorphism and
the earlier 13-cell odd hyperhole by an edge-preserving injection. Actual
domain lists must be subsets of the indicated three-/four-digit palette.
No-hit and budget exhaustion are never evidence of novelty.
"""
import argparse
import hashlib
import itertools
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent


def peer(a, b):
    r, c = divmod(a, 9)
    s, d = divmod(b, 9)
    return a != b and (r == s or c == d or (r//3, c//3) == (s//3, d//3))


def graph(cells):
    return {a: {b for b in cells if peer(a, b)} for a in cells}


def inject(pattern, host, budget, induced=False):
    """Return a checked vertex injection, or a bounded-search status."""
    mapping, used = {}, set()
    visits = [0]

    def visit():
        visits[0] += 1
        if visits[0] > budget:
            raise TimeoutError
        if len(mapping) == len(pattern):
            return dict(mapping)
        available = {}
        for v in pattern:
            if v in mapping:
                continue
            available[v] = [w for w in host if w not in used and len(host[w]) >= len(pattern[v])
                            and all(((w in host[y]) == (v in pattern[x]) if induced
                                     else (v not in pattern[x] or w in host[y]))
                                    for x, y in mapping.items())]
        pivot = min(available, key=lambda v: (len(available[v]),
                    -sum(x in mapping for x in pattern[v]), -len(pattern[v]), v))
        for w in available[pivot]:
            mapping[pivot] = w
            used.add(w)
            result = visit()
            if result is not None:
                return result
            used.remove(w)
            del mapping[pivot]
        return None

    try:
        answer = visit()
        status = 'match' if answer is not None else 'exhaustive-no-match'
    except TimeoutError:
        answer, status = None, 'budget-exhausted'
    if answer is not None:
        assert len(answer) == len(set(answer.values())) == len(pattern)
        assert all(answer[b] in host[answer[a]] for a in pattern for b in pattern[a])
        if induced:
            assert all((b in pattern[a]) == (answer[b] in host[answer[a]])
                       for a in pattern for b in pattern)
    return answer, status, visits[0]


def masks_from_fixture(case):
    state = case.get('assumed_state', case.get('state', case))
    domains = case.get('assumed_domains', case.get('assumed_closed_domains',
                       state.get('domains', case.get('base_domains'))))
    if domains is None or len(domains) != 81 or any(type(m) is not int or not 0 < m <= 511 for m in domains):
        raise ValueError('Expected 81 nonempty nine-bit domains')
    domains = domains.copy()
    if 'target' in case:
        target = case['target']
        assumed = target.get('assumed', target.get('value'))
        if type(assumed) is not bool:
            raise ValueError('Target must record the assumption polarity')
        bit = 1 << (target['digit'] - 1)
        domains[target['cell']] = domains[target['cell']] & bit if assumed else domains[target['cell']] & ~bit
        if domains[target['cell']] == 0:
            raise ValueError('Target gives an immediate empty domain')
    return domains


def tridagon(domains, max_matches=1):
    known = json.loads((HERE / 'real_Loki.json').read_text())
    pattern = graph(known['graph_cells'])
    matches, tested = [], 0
    perms = list(itertools.permutations(range(3)))
    for digits in itertools.combinations(range(1, 10), 3):
        mask = sum(1 << (d-1) for d in digits)
        for bands in itertools.combinations(range(3), 2):
            for stacks in itertools.combinations(range(3), 2):
                options = []
                for band, stack in itertools.product(bands, stacks):
                    diagonals = [[9*(3*band+r)+3*stack+perm[r] for r in range(3)] for perm in perms]
                    options.append([ds for ds in diagonals if all(domains[c] & ~mask == 0 for c in ds)])
                if not all(options):
                    continue
                for selected in itertools.product(*options):
                    cells = sorted(sum(selected, []))
                    tested += 1
                    mapping, status, visits = inject(pattern, graph(cells), 10000, induced=True)
                    if mapping is None:
                        continue
                    matches.append({'family': 'published-Loki-12-cell-tridagon', 'palette': digits,
                        'cells': cells, 'known_to_candidate_mapping': sorted(mapping.items()),
                        'actual_domains': [domains[c] for c in cells],
                        'primary_source': known['source'], 'mapping_visits': visits})
                    if len(matches) >= max_matches:
                        return matches, tested
    return matches, tested


def hyperhole(domains, budget=25000):
    bags = [[9], [15,16,17], [25], [34,43], [51,53], [46,47], [27,36]]
    pattern = {v: set() for v in sum(bags, [])}
    for i, bag in enumerate(bags):
        for a in bag:
            pattern[a].update(set(bag) - {a})
            pattern[a].update(bags[(i-1) % 7] + bags[(i+1) % 7])
    tested, exhausted = 0, 0
    for digits in itertools.combinations(range(1, 10), 4):
        mask = sum(1 << (d-1) for d in digits)
        eligible = [c for c in range(81) if domains[c] & ~mask == 0]
        host = graph(eligible)
        # The known hyperhole has minimum degree four; prune host's four-core.
        while True:
            low = [v for v in host if len(host[v]) < 4]
            if not low:
                break
            host = {v: ns - set(low) for v, ns in host.items() if v not in low}
        pending = set(host)
        while pending:
            todo, component = [min(pending)], set()
            while todo:
                v = todo.pop()
                if v not in component:
                    component.add(v)
                    todo.extend(host[v] - component)
            pending -= component
            if len(component) < 13:
                continue
            tested += 1
            subhost = {v: host[v] & component for v in component}
            mapping, status, visits = inject(pattern, subhost, budget, induced=False)
            exhausted += status == 'budget-exhausted'
            if mapping is not None:
                cells = sorted(mapping.values())
                return {'family': 'known-13-cell-7-hyperhole', 'palette': digits, 'cells': cells,
                        'known_to_candidate_mapping': sorted(mapping.items()),
                        'actual_domains': [domains[c] for c in cells],
                        'required_edges': sum(map(len, pattern.values()))//2,
                        'mapping_visits': visits,
                        'primary_principle': 'https://arxiv.org/pdf/1907.11905v1'}, tested, exhausted
    return None, tested, exhausted


def classify(case, close_branch=False, check_hyperhole=True):
    domains = masks_from_fixture(case)
    closure_events = 0
    if close_branch:
        # This optional discovery step uses the already separately audited closure.
        # A saved hit is verified by coordinate edges and actual final masks above.
        from real_verify import close
        try:
            domains, trace = close(domains)
            closure_events = len(trace)
        except ValueError as error:
            if not str(error).startswith(('Hall contradiction', 'Three-house contradiction')):
                raise
            return {'classification': 'local-contradiction', 'details': str(error)}
    known, tested = tridagon(domains)
    h, ht, exhausted = (None, 0, 0) if known or not check_hyperhole else hyperhole(domains)
    return {'schema_version': 1, 'classification': 'known-pattern' if known or h else 'no-checked-known-match',
            'tridagon_matches': known, 'hyperhole_match': h, 'tridagon_supports_tested': tested,
            'hyperhole_components_tested': ht, 'hyperhole_mapping_budget_exhaustions': exhausted,
            'local_closure_events': closure_events,
            'effective_domains_sha256': hashlib.sha256(json.dumps(domains, separators=(',',':')).encode()).hexdigest(),
            'limits': 'No-hit is not novelty. Tests two explicit known motifs; does not classify all MSLS, Exocet, forcing, impossible-pattern catalogues or arbitrary conditional reductions.'}


def main():
    if not __debug__:
        raise ValueError('Run without -O')
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('fixture', type=Path)
    parser.add_argument('--close-branch', action='store_true')
    parser.add_argument('--skip-hyperhole', action='store_true')
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = classify(json.loads(args.fixture.read_text()), args.close_branch, not args.skip_hyperhole)
    result['fixture'] = str(args.fixture)
    result['fixture_sha256'] = hashlib.sha256(args.fixture.read_bytes()).hexdigest()
    if args.output:
        args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
