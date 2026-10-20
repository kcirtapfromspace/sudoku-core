#!/usr/bin/env python3
"""Flag known Loki12 geometry inside discovered supports, ignoring digit lists.

This is a conservative overlap flag, not a proof that a mixed-list deduction
is equivalent to a tridagon. A real reduction must also check domain premises.
"""
import argparse
import json
from pathlib import Path
from different_audit_known_patterns import tridagon


def overlap(cells):
    dummy = [511] * 81
    for c in cells:
        dummy[c] = 7
    hits, tested = tridagon(dummy)
    return {'known_tridagon_geometry_present': bool(hits),
            'known_to_candidate_mapping': hits[0]['known_to_candidate_mapping'] if hits else None,
            'embedded_cells': hits[0]['cells'] if hits else None,
            'supports_tested': tested,
            'limits': 'Ignores actual digit lists: a geometric overlap is not a valid tridagon proof or historical identity by itself.'}


def main():
    if not __debug__:
        raise ValueError('Run without -O')
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directory', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    rows = []
    for path in sorted(args.directory.glob('candidate-*.json')):
        case = json.loads(path.read_text())
        rows.append({'file': str(path), 'core_cells': len(case['graph_cells']),
                     **overlap(case['graph_cells'])})
    result = {'schema_version': 1, 'records': rows,
              'total': len(rows),
              'known_tridagon_geometry_present': sum(r['known_tridagon_geometry_present'] for r in rows)}
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k:v for k,v in result.items() if k != 'records'}))


if __name__ == '__main__':
    main()
