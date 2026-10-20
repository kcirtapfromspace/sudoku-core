#!/usr/bin/env python3
"""Independent exact check of the different-search negative-novelty control.

Checks the recorded conditional composition, positive/negative recognizer
controls, and independent primal/Rust replay validators. No new-family claim.
"""
import copy
import itertools
import hashlib
import json
from pathlib import Path

from audit_global_certificates import verify as verify_rust
from different_audit_known_patterns import classify, masks_from_fixture, peer, graph, inject, tridagon
from different_audit_geometry import overlap
from different_audit_conditional import odd_cycle
from nonrank_audit_template import verify as verify_template
from real_four_digit_mine import initial_domains
from real_verify import close, exact_completions, graph_colorings, replay_trace, triple_feasible

HERE = Path(__file__).resolve().parent


def require(value, message):
    if not value:
        raise ValueError(message)


def check_pattern(domains, record):
    family = record['family']
    palette = record['palette']
    require(len(set(palette)) == len(palette) and all(type(d) is int and 1 <= d <= 9 for d in palette), 'palette')
    mask = sum(1 << (d-1) for d in palette)
    if family == 'known-two-color-odd-cycle':
        cells = record['cycle_cells']
        require(len(palette) == 2 and len(cells) >= 3 and len(cells) % 2 == 1, 'odd-cycle size')
        require(len(set(cells)) == len(cells), 'odd-cycle repeats')
        require(all(peer(cells[i], cells[(i+1) % len(cells)]) for i in range(len(cells))), 'odd-cycle native edge')
    elif family == 'published-Loki-12-cell-tridagon':
        require(len(palette) == 3, 'tridagon palette')
        known = json.loads((HERE/'real_Loki.json').read_text())
        pairs = record['known_to_candidate_mapping']
        require(len(pairs) == 12, 'mapping length')
        mapping = dict(pairs)
        require(set(mapping) == set(known['graph_cells']) and len(set(mapping.values())) == 12, 'tridagon bijection')
        for a in mapping:
            for b in mapping:
                require(peer(a, b) == peer(mapping[a], mapping[b]), 'full native isomorphism')
        cells = sorted(mapping.values())
        require(cells == record['cells'], 'recorded support')
        # Establish non-three-colorability independently of the family name.
        require(graph_colorings([7]*81, known['graph_cells'])[0] == 0, 'known core must be non-three-colorable')
    elif family == 'known-13-cell-7-hyperhole':
        require(len(palette) == 4, 'hyperhole palette')
        bags = [[9], [15,16,17], [25], [34,43], [51,53], [46,47], [27,36]]
        mapping = dict(record['known_to_candidate_mapping'])
        require(set(mapping) == set(sum(bags, [])) and len(set(mapping.values())) == 13, 'hyperhole injection')
        for i, bag in enumerate(bags):
            for a in bag:
                for b in set(bag + bags[(i+1) % 7]) - {a}:
                    require(peer(mapping[a], mapping[b]), 'required hyperhole edge')
        cells = sorted(mapping.values())
        require(cells == record['cells'], 'hyperhole recorded support')
    else:
        raise ValueError('Unrecognized known lemma')
    require(all(type(c) is int and 0 <= c < 81 and domains[c] > 0 and domains[c] & ~mask == 0 for c in cells), 'actual list containment')
    require(record['actual_domains'] == [domains[c] for c in cells], 'recorded domain binding')


def check_composition(case, record):
    require(record['classification'] == 'known-shallow-composition', 'composition status')
    domains = masks_from_fixture(case)
    pivot = record['pivot']
    require(type(pivot) is int and 0 <= pivot < 81 and bin(domains[pivot]).count('1') >= 2, 'pivot')
    pivots = [p for p in record['pivots_tested'] if p['pivot'] == pivot]
    require(len(pivots) == 1, 'one chosen pivot')
    branches = pivots[0]['branches']
    expected = [d for d in range(1, 10) if domains[pivot] & (1 << (d-1))]
    require([b['digit'] for b in branches] == expected, 'every pivot candidate exactly once')
    families = []
    for branch in branches:
        assumed = domains.copy()
        assumed[pivot] = 1 << (branch['digit']-1)
        if branch['classification'] == 'local-contradiction':
            try:
                close(assumed)
            except ValueError as error:
                require(str(error) == branch['contradiction'] and str(error).startswith(
                    ('Hall contradiction', 'Three-house contradiction')), 'exact local contradiction')
            else:
                raise ValueError('Purported contradictory branch is locally feasible')
            families.append('local-contradiction')
        else:
            require(branch['classification'] == 'known-pattern', 'uncovered branch')
            replayed = replay_trace(assumed, branch['closure_trace'])
            recomputed, _ = close(assumed)
            require(replayed == recomputed == branch['closed_domains'], 'exact conditional closure')
            require(hashlib.sha256(json.dumps(replayed, separators=(',',':')).encode()).hexdigest() ==
                    branch['closed_domains_sha256'], 'conditional state hash')
            check_pattern(replayed, branch['known_pattern'])
            families.append(branch['known_pattern']['family'])
    return {'pivot': pivot, 'digits': expected, 'branch_families': families}


def recognizer_controls():
    known = json.loads((HERE/'real_Loki.json').read_text())
    domains = [511]*81
    for c in known['graph_cells']:
        domains[c] = 7
    result = classify({'domains': domains})
    require(result['classification'] == 'known-pattern' and result['tridagon_matches'], 'positive tridagon recognition')
    check_pattern(domains, result['tridagon_matches'][0])
    # One extra candidate is a real guardian, so no direct contradiction exists.
    domains[known['graph_cells'][0]] |= 8
    require(classify({'domains': domains})['classification'] == 'no-checked-known-match', 'guardian false positive')
    require(classify({'domains': [511]*81})['classification'] == 'no-checked-known-match', 'blank false positive')
    old = json.loads((HERE/'structural_four_digit_13.json').read_text())
    result = classify(old)
    require(result['hyperhole_match'] is not None, 'positive hyperhole recognition')
    check_pattern(masks_from_fixture(old), result['hyperhole_match'])
    return 4


def check_manifests():
    """Validate committed compact metadata/maps, not absent discovery traces.

    Whole conditional proofs stay under target; only the stable control above
    is independently replayed from givens in the default CI run.
    """
    results = []
    corpus = {c['id'] for c in json.loads((HERE/'expanded_modern_corpus.json').read_text())}
    known_cells = json.loads((HERE/'real_Loki.json').read_text())['graph_cells']
    for name, expected in [('different_audit_manifest.json', 327),
                           ('different_audit_geometryfree_manifest.json', 20),
                           ('different_audit_mith_manifest.json', 135)]:
        path = HERE/name
        require(path.exists(), 'required campaign manifest missing')
        manifest = json.loads(path.read_text())
        rows = manifest['records']
        require(manifest['completed'] == manifest['requested'] == len(rows) == expected, 'complete manifest count')
        observed = {}
        sources = set()
        hits = 0
        for row in rows:
            require(row['source_fixture'] not in sources, 'duplicate discovery record')
            sources.add(row['source_fixture'])
            require(row['host_id'] in corpus, 'unrecognized published host')
            cells = row['core_cells']
            require(cells == sorted(set(cells)) and len(cells) == row['core_size'] and
                    all(type(c) is int and 0 <= c < 81 for c in cells), 'native selected support')
            target = row['target_assumption']
            require(type(target['cell']) is int and 0 <= target['cell'] < 81 and
                    type(target['digit']) is int and 1 <= target['digit'] <= 9 and
                    type(target['value']) is bool, 'manifest target')
            pairs = row['known12_geometry_mapping']
            if pairs is not None:
                mapping = dict(pairs)
                require(set(mapping) == set(known_cells) and len(mapping) == 12 and
                        len(set(mapping.values())) == 12 and set(mapping.values()) <= set(cells), 'geometry bijection')
                require(all(peer(a,b) == peer(mapping[a],mapping[b]) for a in mapping for b in mapping), 'geometry native edges')
                hits += 1
            else:
                require(not overlap(cells)['known_tridagon_geometry_present'], 'missed standard four-box geometry')
            status = row['conditional_status']
            require(status in ['known-shallow-composition','no-checked-shallow-composition','inconclusive-budget'], 'manifest status')
            observed[status] = observed.get(status, 0) + 1
            digest = row['conditional_proof_sha256']
            if status == 'known-shallow-composition':
                require(type(digest) is str and len(digest) == 64 and all(c in '0123456789abcdef' for c in digest), 'semantic proof hash')
            else:
                require(digest is None, 'unclassified entry cannot claim a checked proof')
        require(observed == manifest['conditional_status_counts'] and hits == manifest['known12_geometry_overlaps'], 'manifest summary counts')
        results.append({'manifest': name, 'records': len(rows), 'validated_native_geometry_maps': hits,
                        'conditional_status_counts': observed,
                        'trust_boundary': 'Maps and aggregate metadata checked; bulk target proofs are not embedded in this manifest.'})
    direct = json.loads((HERE/'different_audit_geometryfree_direct.json').read_text())
    require(direct['total_screened'] == 20 and direct['direct_full_state_Loki12_count'] == 18, 'direct screen metadata')
    for entry in direct['key_cases']:
        for record in entry['tridagon_matches']:
            dummy = [511]*81
            for c, mask in zip(record['cells'], record['actual_domains']):
                dummy[c] = mask
            check_pattern(dummy, record)
    return results


def check_lead_negative_menu(case, record):
    domains = masks_from_fixture(case)
    require(record['classification'] == 'no-checked-shallow-composition', 'lead menu must finish without timeout')
    require(record['effective_domains_sha256'] == hashlib.sha256(json.dumps(domains,
            separators=(',',':')).encode()).hexdigest(), 'lead opposite snapshot')
    expected = {c for c in range(81) if bin(domains[c]).count('1') > 1}
    tested = [p['pivot'] for p in record['pivots_tested']]
    require(len(tested) == len(set(tested)) and set(tested) == expected, 'every full-board eligible pivot')
    for pivot in record['pivots_tested']:
        require(pivot['branches'], 'empty failed-pivot evidence')
        failed = [b for b in pivot['branches'] if b['classification'] == 'unclassified']
        require(len(failed) == 1, 'one menu-unclassified branch per failed pivot')
        branch = failed[0]
        bit = 1 << (branch['digit']-1)
        require(domains[pivot['pivot']] & bit, 'counterbranch candidate')
        assumed = domains.copy()
        assumed[pivot['pivot']] = bit
        final, _ = close(assumed)
        require(branch['closed_domains_sha256'] == hashlib.sha256(json.dumps(final,
                separators=(',',':')).encode()).hexdigest(), 'counterbranch exact closure')
        require(odd_cycle(final) is None and not tridagon(final)[0], 'known menu would cover counterbranch')
    return len(tested)


def check_lead():
    fixture = HERE/'different_lead.json'
    case = json.loads(fixture.read_text())
    record = json.loads((HERE/'different_audit_lead_classification.json').read_text())
    require(record['fixture_sha256'] == hashlib.sha256(fixture.read_bytes()).hexdigest(), 'lead fixture hash')
    pivots = check_lead_negative_menu(case, record)
    cells = case['graph_cells']
    require(len(cells) == 20 and sum(peer(a,b) for i,a in enumerate(cells) for b in cells[:i]) == 51, 'lead graph invariants')
    require(graph_colorings([7]*81, cells)[0] == 0, 'lead ordinary three-color refutation')
    # Explicit proper four-color witness, checked only against coordinate edges.
    witness = {3:0,13:1,9:0,15:2,16:3,7:1,2:2,19:1,43:0,33:1,
               30:2,57:1,66:3,64:2,54:3,68:0,76:2,28:0,31:3,74:0}
    require(set(witness) == set(cells) and set(witness.values()) == set(range(4)), 'four-color witness support')
    require(all(witness[a] != witness[b] for a in cells for b in cells if peer(a,b)), 'four-color native conflicts')
    known = json.loads((HERE/'real_Loki.json').read_text())['graph_cells']
    mapping, status, visits = inject(graph(known), graph(cells), 1000000, induced=False)
    require(mapping is None and status == 'exhaustive-no-match', 'complete abstract Loki12 edge-injection search')
    require(not overlap(cells)['known_tridagon_geometry_present'], 'complete standard four-box geometry search')
    require(exact_completions(case['puzzle'], 2) == [case['solution']], 'lead independent uniqueness')
    initial = initial_domains([int(c) if c in '123456789' else 0 for c in case['puzzle']])
    require(initial == case['initial_domains'], 'lead exact givens candidates')
    base = replay_trace(initial, case['base_trace'], list(map(int,case['solution'])))
    require(base == case['base_domains'] == case['state']['domains'], 'lead base premise trace')
    raw = base.copy()
    target = case['target']
    bit = 1 << (target['digit']-1)
    raw[target['cell']] = raw[target['cell']] & bit if target['value'] else raw[target['cell']] & ~bit
    branch = replay_trace(raw,case['assumption_trace'])
    require(branch == case['assumed_domains'] == case['assumed_state']['domains'], 'lead opposite premise trace')
    require(close(branch)[1] == [], 'lead Hall/all81 fixedpoint')
    for selected in itertools.combinations(range(27),3):
        require(triple_feasible(branch,selected), 'lead native three-house conjunction')
    template = verify_template(case, json.loads((HERE/'different_lead_template_primal.json').read_text()))
    replay = json.loads((HERE/'replay_different_lead.json').read_text())
    require(replay['state'] == case['state'], 'lead production exact base state')
    require((replay['proof']['cell'],replay['proof']['digit'],replay['proof']['value']) ==
            (case['conclusion']['cell'],case['conclusion']['digit'],case['conclusion']['value']), 'lead result endpoint')
    result = {'schema_version':1, 'verified':True,
              'fixture_sha256':hashlib.sha256(fixture.read_bytes()).hexdigest(),
              'full_board_non_singleton_pivots_rechecked':pivots,
              'menu_no_hit_scope':'One additional cell split followed by Hall/all81 closure, then local contradiction, Loki12 containment or a uniform two-color odd cycle. Conditional hyperholes, all 630 patterns and arbitrary forcing are not excluded.',
              'ordinary_graph_chromatic_number':4, 'native_edges':51,
              'proper_four_coloring':witness,
              'abstract_Loki12_edge_injection_status':status, 'abstract_mapping_search_visits':visits,
              'standard_four_box_geometry_search_complete':True,
              'independent_solution_count_capped_at_2':1,
              'givens_and_opposite_traces_replayed':True,
              'native_three_house_conjunctions_checked':2925,
              'exact_template_relaxation':template, 'independent_rust_replay':verify_rust(replay),
              'historical_classification':'Unclassified research lead; no historical novelty claim.'}
    (HERE/'different_audit_lead_result.json').write_text(json.dumps(result,indent=2)+'\n')
    return result


def main():
    require(__debug__, 'Run without -O')
    fixture = HERE/'different_control.json'
    case = json.loads(fixture.read_text())
    record = json.loads((HERE/'different_audit_control_classification.json').read_text())
    require(record['fixture_sha256'] == hashlib.sha256(fixture.read_bytes()).hexdigest(), 'fixture hash')
    composition = check_composition(case, record)
    corruptions = []
    bad = copy.deepcopy(record)
    chosen = next(p for p in bad['pivots_tested'] if p['pivot'] == bad['pivot'])
    chosen['branches'].pop()
    corruptions.append(bad)
    bad = copy.deepcopy(record)
    bad['pivot'] = 81
    corruptions.append(bad)
    bad = copy.deepcopy(record)
    branch = next(b for p in bad['pivots_tested'] for b in p['branches'] if b['classification'] == 'known-pattern')
    branch['known_pattern']['known_to_candidate_mapping'][0][1] = 81
    corruptions.append(bad)
    bad = copy.deepcopy(record)
    branch = next(b for p in bad['pivots_tested'] for b in p['branches'] if b['classification'] == 'known-pattern')
    branch['closed_domains'][0] ^= 1
    corruptions.append(bad)
    for bad in corruptions:
        try:
            check_composition(case, bad)
        except (ValueError, AssertionError, IndexError, KeyError, TypeError):
            continue
        raise ValueError('Tampered composition accepted')
    require(exact_completions(case['puzzle'], 2) == [case['solution']], 'independent uniqueness')
    initial = initial_domains([int(c) if c in '123456789' else 0 for c in case['puzzle']])
    require(initial == case['initial_domains'], 'givens-only initial candidates')
    base = replay_trace(initial, case['base_trace'], list(map(int, case['solution'])))
    require(base == case['base_domains'] == case['state']['domains'], 'exact base premise trace')
    raw = base.copy()
    target = case['target']
    bit = 1 << (target['digit']-1)
    raw[target['cell']] = raw[target['cell']] & bit if target['value'] else raw[target['cell']] & ~bit
    branch = replay_trace(raw, case['assumption_trace'])
    require(branch == case['assumed_domains'] == case['assumed_state']['domains'], 'exact opposite-target trace')
    require(close(branch)[1] == [], 'Hall/all81 fixedpoint')
    triples = 0
    for selected in itertools.combinations(range(27), 3):
        require(triple_feasible(branch, selected), 'infeasible three-house conjunction')
        triples += 1
    template = verify_template(case, json.loads((HERE/'different_template_primal.json').read_text()))
    replay = json.loads((HERE/'replay_different_control.json').read_text())
    require(replay['state'] == case['state'], 'production base state')
    require((replay['proof']['cell'], replay['proof']['digit'], replay['proof']['value']) ==
            (case['conclusion']['cell'], case['conclusion']['digit'], case['conclusion']['value']), 'production conclusion')
    rust = verify_rust(replay)
    output = {'schema_version': 1, 'verified': True,
              'fixture_sha256': hashlib.sha256(fixture.read_bytes()).hexdigest(),
              'independent_solution_count_capped_at_2': 1,
              'givens_and_opposite_traces_replayed': True,
              'native_three_house_conjunctions_checked': triples,
              'checked_composition': composition, 'rejected_composition_mutations': len(corruptions),
              'positive_negative_recognizer_controls': recognizer_controls(),
              'compact_manifests': check_manifests(),
              'exact_template_relaxation': template, 'independent_rust_replay': rust,
              'lead_audit': check_lead(),
              'classification': 'Valid exact-control deduction, standard conditional tridagon composition; not evidence of a new family.'}
    (HERE/'different_audit_control_result.json').write_text(json.dumps(output, indent=2)+'\n')
    print(json.dumps(output, indent=2))


if __name__ == '__main__':
    main()
