#!/usr/bin/env python3
"""Bounded classification of SAT discoveries; bulky proof traces stay in target.

For each host puzzle, reuse pivots already giving checked known compositions.
Compact manifests distinguish geometry overlap from an actual checked reduction.
"""
import argparse
from collections import Counter, defaultdict
import hashlib
import json
from pathlib import Path
import time

from different_audit_conditional import shallow
from different_audit_geometry import overlap


def proof_payload(record):
    if record['classification'] != 'known-shallow-composition':
        return None
    return {'pivot': record['pivot'], 'branches': next(p['branches'] for p in record['pivots_tested']
                                                       if p['pivot'] == record['pivot'])}


def main():
    if not __debug__:
        raise ValueError('Run without -O')
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directories', nargs='+', type=Path)
    parser.add_argument('--seconds-per-case', type=float, default=3)
    parser.add_argument('--reuse-first', type=Path)
    parser.add_argument('--output-dir', type=Path, default=Path('target/different-audit-batch'))
    parser.add_argument('--manifest', type=Path, default=Path('research/global_obstructions/different_audit_manifest.json'))
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    previous = {}
    preferred = defaultdict(Counter)
    if args.reuse_first:
        prior = json.loads((args.reuse_first/'summary.json').read_text())
        for row in prior['records']:
            old = args.reuse_first / Path(row['file']).name
            previous[row['file']] = json.loads(old.read_text())
            if row['classification'] == 'known-shallow-composition':
                preferred[row['puzzle']][row['pivot']] += 1
    items = [(len(c['graph_cells']), p, c) for d in args.directories
             for p in d.glob('candidate-*.json') for c in [json.loads(p.read_text())]]
    items.sort(key=lambda item: (item[0], str(item[1])))
    rows, started = [], time.monotonic()
    for size, path, case in items:
        geometry = overlap(case['graph_cells'])
        record = previous.get(str(path))
        if record is None:
            pivots = [p for p, _ in preferred[case['name']].most_common()]
            record = shallow(case, seconds=args.seconds_per_case, preferred_pivots=pivots)
            if record['classification'] == 'known-shallow-composition':
                preferred[case['name']][record['pivot']] += 1
        destination = args.output_dir / (path.parent.name + '-' + path.name)
        destination.write_text(json.dumps(record, indent=2) + '\n')
        payload = proof_payload(record)
        rows.append({'source_fixture': str(path),
                     'fixture_sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
                     'host_id': case['name'], 'target_assumption': case['target'],
                     'core_cells': case['graph_cells'], 'core_size': size,
                     'known12_geometry_mapping': geometry['known_to_candidate_mapping'],
                     'conditional_status': record['classification'],
                     'conditional_pivot': record.get('pivot'),
                     'conditional_proof_sha256': hashlib.sha256(json.dumps(payload, sort_keys=True,
                         separators=(',', ':')).encode()).hexdigest() if payload else None,
                     'proof_file': str(destination)})
        manifest = {'schema_version': 1, 'records': rows,
                    'completed': len(rows), 'requested': len(items),
                    'host_counts': dict(Counter(r['host_id'] for r in rows)),
                    'conditional_status_counts': dict(Counter(r['conditional_status'] for r in rows)),
                    'known12_geometry_overlaps': sum(r['known12_geometry_mapping'] is not None for r in rows),
                    'seconds': time.monotonic()-started,
                    'limits': 'Geometry overlap alone does not prove mixed-list equivalence. No checked reduction and exhausted budgets are inconclusive, not novelty. Sources and full traces are discovery bulk under target.'}
        args.manifest.write_text(json.dumps(manifest, indent=2) + '\n')
        if len(rows) % 25 == 0:
            print(len(rows), manifest['conditional_status_counts'], flush=True)
    print('DONE', manifest['completed'], manifest['conditional_status_counts'], flush=True)


if __name__ == '__main__':
    main()
