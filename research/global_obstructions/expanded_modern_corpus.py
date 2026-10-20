#!/usr/bin/env python3
"""Validate the pinned modern200 givens corpus; optionally refetch/re-extract it.

Offline: python3 research/global_obstructions/expanded_modern_corpus.py
Source replay: add --fetch-sources (three immutable raw GitHub files).

Extraction uses only the FIRST whitespace-separated token on each selected
source line. It must be exactly81 cells, containing dots/digits and a blank;
ratings, dates, solutions, expanded grids and sukaku are never puzzle columns.
This does not count solutions or establish uniqueness; the engine corpus run
performs that independent check.
"""
import argparse
import hashlib
import json
import re
import urllib.request
from pathlib import Path

if not __debug__:
    raise SystemExit('expanded_modern_corpus.py must run without Python -O')

HERE = Path(__file__).resolve().parent
PIN = '15c9d217b08d2d3ebf0f147d4f7c0beb5c2fcdb2'
REPOSITORY = 'https://github.com/denis-berthier/Sudoku-classif'
RAW = f'https://raw.githubusercontent.com/denis-berthier/Sudoku-classif/{PIN}/'
SOURCES = {
    'Monhard-until-2023-08-15-TE23/TE2/puzzles.txt': (
        'c0c0640b7764a28793dd7383a6171359c38d52dd6edf10717ab53ba841b0cb19',
        100, 'monhard-te2', 11.7,
        'Hendrik Monhard, Aug2023 minimal puzzle collection; maintained/classified by Denis Berthier in Sudoku-classif; selected source SER>=11.7'),
    'Monhard-until-2023-08-15-TE23/TE3/puzzles.txt': (
        '0fc9e6f84f6845c249452c4d88502ea037d27647cba791c3e0bd27485358131b',
        50, 'monhard-te3', 11.7,
        'Hendrik Monhard, Aug2023 minimal puzzle collection; maintained/classified by Denis Berthier in Sudoku-classif; selected source SER>=11.7'),
    'mith-158276-TE3/Sample500/minimals.txt': (
        '70ab1ad899500e440e15eb70bb7e91a5b57805f6ccee7eafa658bfd173cd532a',
        50, 'mith-minimals', None,
        'Philip Newman (mith), post2022 T&E3 collection; minimal reconstructions and classification by Denis Berthier in Sudoku-classif Sample500/minimals.txt'),
}
FIELDS = {'id', 'puzzle', 'source_url', 'source_sha256', 'source_index_1based', 'attribution'}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def native_houses():
    return ([[r * 9 + c for c in range(9)] for r in range(9)]
            + [[r * 9 + c for r in range(9)] for c in range(9)]
            + [[r * 9 + c for r in range(br * 3, br * 3 + 3)
                for c in range(bc * 3, bc * 3 + 3)] for br in range(3) for bc in range(3)])


def validate(corpus_bytes, metadata, old_records):
    require(metadata['schema_version'] == 1 and metadata['commit'] == PIN
            and metadata['repository'] == REPOSITORY, 'source pin or schema changed')
    require(hashlib.sha256(corpus_bytes).hexdigest() == metadata['corpus_sha256'],
            'corpus SHA256 mismatch')
    records = json.loads(corpus_bytes)
    require(type(records) is list and len(records) == metadata['count'] == 200,
            'expected exactly200 records')
    require([source['path'] for source in metadata['sources']] == list(SOURCES),
            'expected exactly the three pinned source files, in order')
    expected_provenance = []
    for source in metadata['sources']:
        path = source['path']
        digest, count, label, _, attribution = SOURCES[path]
        rows = source['selected_rows']
        require(source['sha256'] == digest and source['count'] == len(rows) == count,
                'source hash or selection count changed')
        require(len(set(rows)) == count and all(type(row) is int and 1 <= row <= source['source_lines']
                                               for row in rows), 'source line indices invalid')
        expected_provenance.extend((path, row, digest, label, attribution) for row in rows)
    identifiers, puzzles = set(), set()
    for record, (path, row, digest, label, attribution) in zip(records, expected_provenance):
        require(set(record) == FIELDS, 'record schema changed')
        require(record['id'] == f'modern-{label}-{row:06d}' and record['id'] not in identifiers,
                'record ID mismatch or duplicate')
        require(record['source_url'] == RAW + path and record['source_sha256'] == digest
                and record['source_index_1based'] == row and record['attribution'] == attribution,
                'record provenance differs from manifest')
        puzzle = record['puzzle']
        require(type(puzzle) is str and re.fullmatch(r'[.1-9]{81}', puzzle) and '.' in puzzle,
                'puzzle must be81 plain givens with blanks')
        require(puzzle not in puzzles, 'duplicate literal puzzle')
        for house in native_houses():
            givens = [puzzle[v] for v in house if puzzle[v] != '.']
            require(len(givens) == len(set(givens)), 'native house contains conflicting givens')
        identifiers.add(record['id'])
        puzzles.add(puzzle)
    old = {record['puzzle'].replace('0', '.') for record in old_records}
    require(len(old_records) == 106, 'old corpus count changed; review exclusion baseline')
    require(not (puzzles & old), 'literal overlap with old106 corpus')
    return records


def refetch_and_extract(metadata):
    """Rebuild all six fields from pinned bytes, selections and attribution specs."""
    rebuilt = []
    for source in metadata['sources']:
        path = source['path']
        digest, _, label, minimum_ser, attribution = SOURCES[path]
        with urllib.request.urlopen(RAW + path, timeout=30) as response:
            data = response.read()
        require(hashlib.sha256(data).hexdigest() == digest, 'refetched source SHA256 mismatch')
        lines = data.decode('utf-8').splitlines()
        require(len(lines) == source['source_lines'], 'refetched source line count mismatch')
        for row in source['selected_rows']:
            columns = lines[row - 1].split()
            puzzle = columns[0].replace('0', '.') if columns else ''
            require(re.fullmatch(r'[.1-9]{81}', puzzle) and '.' in puzzle,
                    'first source column is not a plain81-cell givens puzzle')
            if minimum_ser is not None:
                rating = re.match(r'[0-9]+(?:\.[0-9]+)?', columns[1]) if len(columns) > 1 else None
                require(rating is not None and float(rating.group()) >= minimum_ser,
                        'selected source record does not meet SER threshold')
            rebuilt.append(dict(id=f'modern-{label}-{row:06d}', puzzle=puzzle,
                                source_url=RAW + path, source_sha256=digest,
                                source_index_1based=row, attribution=attribution))
    return rebuilt


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fetch-sources', action='store_true',
                        help='refetch the three pinned files and re-extract exact selected first columns')
    args = parser.parse_args()
    corpus_bytes = (HERE / 'expanded_modern_corpus.json').read_bytes()
    metadata = json.loads((HERE / 'expanded_modern_corpus_meta.json').read_text())
    old = json.loads((HERE.parent / 'nonmodular/real_game_corpus.json').read_text())
    records = validate(corpus_bytes, metadata, old)
    if args.fetch_sources:
        rebuilt = refetch_and_extract(metadata)
        require(rebuilt == records, 'exact source extraction differs from corpus')
        require((json.dumps(rebuilt, indent=2) + '\n').encode() == corpus_bytes,
                'source reconstruction does not reproduce exact corpus bytes')
    print(json.dumps({'records': len(records), 'distinct_literal_puzzles': len(records),
                      'old106_literal_overlap': 0, 'native_givens': 'pass',
                      'pinned_source_reextraction': 'pass' if args.fetch_sources else 'not requested',
                      'uniqueness': 'not checked here; checked separately by the engine corpus run',
                      'corpus_sha256': metadata['corpus_sha256']}, sort_keys=True))


if __name__ == '__main__':
    main()
