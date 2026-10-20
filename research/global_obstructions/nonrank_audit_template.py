#!/usr/bin/env python3
"""Independent exact audit of full-board single-digit-template LP certificates.

No solver/numerical packages are imported. Positive checking needs only the
listed templates; negative checking exhaustively reconstructs all46656 native
templates before candidate/singleton filtering. Fractional feasibility is not
a Sudoku solution, and no uniqueness or closure trace is checked here.
"""
import argparse
import copy
import hashlib
import itertools
import json
from fractions import Fraction
from pathlib import Path

if not __debug__:
    raise SystemExit('nonrank_audit_template.py must run without Python -O')

RELAXATION = 'native-single-digit-template-v1'
HERE = Path(__file__).resolve().parent


def require(condition, message):
    if not condition:
        raise ValueError(message)


def source(fixture):
    """Reconstruct snapshot and opposite assumption without importing search."""
    state = fixture.get('assumed_state', fixture.get('state', fixture))
    require(isinstance(state, dict), 'candidate state must be an object')
    domains = state.get('assumed_domains', state.get('assumed_closed_domains', state.get('domains')))
    require(type(domains) is list and len(domains) == 81, 'expected81 source domains')
    values = state.get('values')
    if values is None:
        puzzle = fixture.get('puzzle')
        if puzzle is None:
            values = [0] * 81
        else:
            require(type(puzzle) is str and len(puzzle) == 81
                    and all(c in '.0123456789' for c in puzzle), 'invalid source givens')
            values = [0 if c == '.' else int(c) for c in puzzle]
    require(type(values) is list and len(values) == 81, 'expected81 source values')
    for value, mask in zip(values, domains):
        require(type(value) is int and 0 <= value <= 9 and type(mask) is int
                and 0 < mask <= 511, 'invalid value or domain')
        require(not value or mask == 1 << (value - 1), 'placed value/domain mismatch')
    target = fixture['target']
    cell, digit = target.get('cell'), target.get('digit')
    assumed = target.get('assumed', target.get('assumption', target.get('value')))
    require(type(cell) is int and 0 <= cell < 81 and type(digit) is int
            and 1 <= digit <= 9 and type(assumed) is bool, 'invalid target assumption')
    require(not values[cell] or (values[cell] == digit) == assumed,
            'placed target contradicts assumption')
    effective = list(domains)
    bit = 1 << (digit - 1)
    require(not assumed or bool(effective[cell] & bit), 'assumed candidate unavailable')
    effective[cell] = bit if assumed else effective[cell] & ~bit
    require(bool(effective[cell]), 'assumption empties target domain')
    normalized = {'cell': cell, 'digit': digit, 'assumed': assumed}
    payload = {'schema': 'native-template-lp-source-v1', 'values': values,
               'domains': effective, 'target': normalized}
    digest = hashlib.sha256(json.dumps(payload, sort_keys=True, separators=(',', ':')).encode()).hexdigest()
    return values, effective, normalized, digest


def common(fixture, witness):
    values, effective, target, digest = source(fixture)
    require(type(witness.get('schema_version')) is int and witness['schema_version'] == 1
            and witness.get('relaxation') == RELAXATION,
            'witness version or relaxation mismatch')
    require(witness.get('state_hash') == digest and witness.get('target') == target,
            'state/assumption binding mismatch')
    return values, effective, target


def legal(cells):
    return (type(cells) is list and len(cells) == 9
            and all(type(v) is int and 0 <= v < 81 for v in cells)
            and cells == sorted(set(cells))
            and {v // 9 for v in cells} == set(range(9))
            and {v % 9 for v in cells} == set(range(9))
            and {(v // 27) * 3 + v % 9 // 3 for v in cells} == set(range(9)))


def verify(fixture, witness):
    values, domains, target = common(fixture, witness)
    masses = [Fraction(0)] * 9
    coverage = [Fraction(0)] * 81
    projected = [[Fraction(0) for _ in range(9)] for _ in range(81)]
    seen = set()
    require(type(witness.get('terms')) is list and bool(witness['terms']), 'missing primal terms')
    for term in witness['terms']:
        require(type(term) is dict and set(term) == {'digit', 'cells', 'numerator', 'denominator'},
                'malformed primal term')
        d, cells, num, den = (term[k] for k in ['digit', 'cells', 'numerator', 'denominator'])
        require(type(d) is int and 1 <= d <= 9 and legal(cells), 'invalid digit/template')
        require(type(num) is int and type(den) is int and num > 0 and den > 0,
                'template coefficient must be a positive exact rational')
        require((d, tuple(cells)) not in seen, 'duplicate digit/template')
        require(all(domains[v] & (1 << (d - 1)) for v in cells), 'template uses unavailable cell')
        require(all(v in cells for v, value in enumerate(values) if value == d),
                'template omits a placed value of its digit')
        seen.add((d, tuple(cells)))
        weight = Fraction(num, den)
        masses[d - 1] += weight
        for v in cells:
            coverage[v] += weight
            projected[v][d - 1] += weight
    require(all(q == 1 for q in masses), 'digit distribution mass is not exactly1')
    require(all(q == 1 for q in coverage), 'cell coverage is not exactly1')
    # Redundant native324 equality audits; independent from template construction.
    houses = ([[r * 9 + c for c in range(9)] for r in range(9)]
              + [[r * 9 + c for r in range(9)] for c in range(9)]
              + [[r * 9 + c for r in range(br * 3, br * 3 + 3)
                  for c in range(bc * 3, bc * 3 + 3)] for br in range(3) for bc in range(3)])
    require(all(sum(projected[v][d] for v in house) == 1
                for house in houses for d in range(9)), 'projected native house equality fails')
    require(projected[target['cell']][target['digit'] - 1] == int(target['assumed']),
            'projected target assumption fails')
    require(all(not value or projected[v][value - 1] == 1 for v, value in enumerate(values)),
            'projected placed value fails')
    return {'result': 'exact feasible full single-digit-template LP', 'digit_masses': 9,
            'cell_coverages': 81, 'projected_native_equalities': 324, 'terms': len(seen)}


def templates():
    """Each column permutation picks one cell per row/column; test all boxes."""
    for columns in itertools.permutations(range(9)):
        if all(len({columns[r] // 3 for r in range(first, first + 3)}) == 3
               for first in [0, 3, 6]):
            yield tuple(9 * r + c for r, c in enumerate(columns))


def verify_negative(fixture, witness):
    """Exact90-row Farkas certificate: nine digit rows then81 cell rows."""
    _, domains, _ = common(fixture, witness)
    weights = [Fraction(0)] * 90
    seen = set()
    for term in witness.get('weights', []):
        require(type(term) is list and len(term) == 3
                and all(type(v) is int for v in term), 'malformed dual weight')
        row, num, den = term
        require(0 <= row < 90 and row not in seen and num != 0 and den > 0,
                'invalid or duplicate dual row')
        weights[row] = Fraction(num, den)
        seen.add(row)
    rhs = sum(weights, Fraction(0))
    require(rhs < 0 and rhs == Fraction(witness['weighted_rhs']), 'dual RHS must be exact and negative')
    forced = [{v for v, mask in enumerate(domains) if mask == 1 << d} for d in range(9)]
    counts = [0] * 9
    total = 0
    for cells in templates():
        total += 1
        positions = set(cells)
        for d in range(9):
            if forced[d] <= positions and all(domains[v] & (1 << d) for v in cells):
                require(weights[d] + sum(weights[9 + v] for v in cells) >= 0,
                        'negative coefficient on an admissible full template')
                counts[d] += 1
    require(total == 46656, 'native template enumeration count mismatch')
    return {'result': 'exact infeasible full single-digit-template LP', 'rhs': str(rhs),
            'native_templates_enumerated': total, 'admissible_templates_per_digit': counts}


def verify_native_negative(fixture, witness):
    """Separate exact324-indicator Farkas check; no native search imports."""
    values, domains, target, _ = source(fixture)
    state = fixture.get('assumed_state', fixture.get('state', fixture))
    raw = state.get('assumed_domains', state.get('assumed_closed_domains', state.get('domains')))
    payload = {'schema': 'native-lp-source-v1', 'values': values, 'domains': raw, 'target': target}
    digest = hashlib.sha256(json.dumps(payload, sort_keys=True, separators=(',', ':')).encode()).hexdigest()
    require(type(witness.get('schema_version')) is int and witness['schema_version'] == 1
            and witness.get('state_hash') == digest and witness.get('target') == target,
            'native Farkas version/state/target mismatch')
    weights = [Fraction(0)] * 324
    seen = set()
    for term in witness.get('terms', []):
        require(type(term) is list and len(term) == 3
                and all(type(v) is int for v in term), 'malformed native dual term')
        row, num, den = term
        require(0 <= row < 324 and row not in seen and num != 0 and den > 0,
                'invalid native dual row')
        weights[row] = Fraction(num, den)
        seen.add(row)
    rhs = sum(weights, Fraction(0))
    require(rhs < 0 and rhs == Fraction(witness['weighted_rhs']), 'native dual RHS is not negative')
    count = 0
    for cell in range(81):
        r, c = divmod(cell, 9)
        box = r // 3 * 3 + c // 3
        for d in range(9):
            if domains[cell] & (1 << d):
                coefficient = (weights[cell] + weights[81 + 9 * r + d]
                               + weights[162 + 9 * c + d] + weights[243 + 9 * box + d])
                require(coefficient >= 0, 'negative native candidate coefficient')
                count += 1
    return {'result': 'exact infeasible native324 LP', 'terms': len(seen),
            'available_columns_checked': count, 'rhs': str(rhs)}


def shell(fixture):
    _, _, target, digest = source(fixture)
    return {'schema_version': 1, 'relaxation': RELAXATION, 'state_hash': digest, 'target': target}


def mutations(fixture, witness):
    cases = []
    for key, value in [('state_hash', 'wrong'), ('schema_version', True), ('relaxation', 'wrong')]:
        changed = copy.deepcopy(witness)
        changed[key] = value
        cases.append(changed)
    changed = copy.deepcopy(witness)
    changed['target']['assumed'] = not changed['target']['assumed']
    cases.append(changed)
    for key, value in [('digit', 0), ('digit', True), ('numerator', -1), ('denominator', 0)]:
        changed = copy.deepcopy(witness)
        changed['terms'][0][key] = value
        cases.append(changed)
    changed = copy.deepcopy(witness)
    changed['terms'][0]['cells'][0] = changed['terms'][0]['cells'][1]
    cases.append(changed)
    changed = copy.deepcopy(witness)
    changed['terms'].append(copy.deepcopy(changed['terms'][0]))
    cases.append(changed)
    changed = copy.deepcopy(witness)
    changed['terms'].pop()
    cases.append(changed)
    changed = copy.deepcopy(witness)
    changed['terms'][0]['numerator'] += 1
    cases.append(changed)
    for changed in cases:
        try:
            verify(fixture, changed)
        except (ValueError, TypeError, KeyError):
            continue
        raise ValueError('accepted a tampered template witness')
    return len(cases)


def known_pattern(fixture):
    """Recheck the recorded exact tridagon graph maps using coordinates only."""
    result = json.loads((HERE / 'nonrank_audit_first_candidate_results.json').read_text())
    core = result['classic_tridagon_core']
    require(len(core) == len(set(core)) == 12
            and all(fixture['assumed_domains'][v] == 200 for v in core),
            'recorded four-seven-eight tridagon domains changed')

    def peer(a, b):
        return a != b and (a // 9 == b // 9 or a % 9 == b % 9
                           or (a // 27, a % 9 // 3) == (b // 27, b % 9 // 3))

    maps = result['exact_graph_isomorphisms_to_published_controls']
    require(len(maps) == 3, 'expected three explicit published-control maps')
    for item in maps:
        reference = json.loads((HERE / item['reference']).read_text())
        mapping = {int(a): b for a, b in item['mapping_newcell_to_referencecell'].items()}
        require(set(mapping) == set(core) and set(mapping.values()) == set(reference['graph_cells']),
                'tridagon map is not a bijection')
        require(all(peer(a, b) == peer(mapping[a], mapping[b])
                    for a, b in itertools.combinations(core, 2)),
                'tridagon map changes native adjacency')
        require(item['reference_source'] == reference['source'], 'published reference provenance changed')
    return {'known_family': 'established twelve-cell tridagon', 'exact_graph_maps_checked': 3}


def controls():
    """Fractional feasible infrastructure control and old hyperhole negative."""
    fixture = {'name': 'two-grid-mixture-control', 'values': [0] * 81,
               'domains': [511] * 81, 'target': {'cell': 0, 'digit': 9, 'value': False}}
    grid = [(3 * r + r // 3 + c) % 9 + 1 for r in range(9) for c in range(9)]
    second = [d % 9 + 1 for d in grid]
    witness = shell(fixture)
    witness['terms'] = [{'digit': d, 'cells': [v for v, value in enumerate(g) if value == d],
                         'numerator': 1, 'denominator': 2}
                        for g in [grid, second] for d in range(1, 10)]
    positive = verify(fixture, witness)
    positive['rejected_mutations'] = mutations(fixture, witness)
    old = json.loads((HERE / 'structural_four_digit_13.json').read_text())
    dual = shell(old)
    dual['weights'] = [[d - 1, 3, 1] for d in old['palette']]
    dual['weights'] += [[9 + v, -1, 1] for v in old['support']]
    dual['weighted_rhs'] = '-1'
    negative = verify_negative(old, dual)
    result = {'positive_control': positive, 'negative_old13_control': negative}
    for name in ['real_4401', 'real_4662', 'real_Loki']:
        path = HERE / ('nonrank_audit_template_' + name + '.json')
        if path.exists():
            original = json.loads((HERE / (name + '.json')).read_text())
            control = json.loads(path.read_text())
            checked = verify(original, control)
            checked['rejected_mutations'] = mutations(original, control)
            result[name] = checked
    current = HERE / 'nonrank_full_old_stall.json'
    if current.exists():
        original = json.loads(current.read_text())
        control = json.loads((HERE / 'nonrank_template_primal.json').read_text())
        checked = verify(original, control)
        checked['rejected_mutations'] = mutations(original, control)
        checked['known_pattern'] = known_pattern(original)
        result['modern_full_old_stall'] = checked
    rejected = HERE / 'nonrank_audit_27_raw_state.json'
    if rejected.exists():
        original = json.loads(rejected.read_text())
        native = json.loads((HERE / 'nonrank_audit_27_native_farkas.json').read_text())
        template = json.loads((HERE / 'nonrank_audit_27_template_farkas.json').read_text())
        result['rejected27_native'] = verify_native_negative(original, native)
        result['rejected27_template'] = verify_negative(original, template)
    native27 = HERE / 'nonrank_graph_native27_candidate.json'
    if native27.exists():
        original = json.loads(native27.read_text())
        control = json.loads((HERE / 'nonrank_audit_native27_template_primal.json').read_text())
        checked = verify(original, control)
        checked['rejected_mutations'] = mutations(original, control)
        result['synthetic_native27_positive'] = checked
    return fixture, witness, dual, result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('fixture', type=Path, nargs='?')
    parser.add_argument('witness', type=Path, nargs='?')
    parser.add_argument('--negative', action='store_true')
    parser.add_argument('--self-test', action='store_true')
    args = parser.parse_args()
    if args.self_test:
        print(json.dumps(controls()[-1], sort_keys=True))
    else:
        require(args.fixture is not None and args.witness is not None, 'fixture and witness required')
        f, w = json.loads(args.fixture.read_text()), json.loads(args.witness.read_text())
        if 'witness' in w:
            w = w['witness']  # Accept the numerical search's envelope; audit only its actual certificate.
        print(json.dumps((verify_negative if args.negative else verify)(f, w), sort_keys=True))


if __name__ == '__main__':
    main()
