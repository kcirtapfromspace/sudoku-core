#!/usr/bin/env python3
"""Fetch and parse the pinned primary author's 592 public Imp630 rules.

The author's manual explicitly excludes38 of the630 IDs as no-free-cell
patterns. We account for their IDs but do not fabricate their missing geometry.
Original GPL source files remain in target; output is derived mathematical data
and per-file provenance, not a vendored copy of the CLIPS implementation.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
from collections import Counter
import hashlib
import json
from pathlib import Path
import re
import urllib.request

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
COMMIT = '7eefd71566b37e55c15f330d474140bfc4f4371d'
REPO = 'denis-berthier/CSP-Rules-V2.1'
COUNTS = {10:31,12:38,13:290,14:159,15:102,16:10}
EXCLUDED = {13:[207,211,217,218,219],14:[23,33,36,73,92,112,113],
            15:[11,16,18,30,32,36,38,40,42,45,47,49,51,53,56,57,66,68,69,80,83,85,100,101],
            16:[9,10]}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def expected_paths():
    return [('EL{}c{}'.format(n,i), 'SudoRules-V20.1/EXOTIC/Imp630/{}c/EL{}c{}.clp'.format(n,n,i))
            for n, count in COUNTS.items() for i in range(1,count+1) if i not in EXCLUDED.get(n,[])]


def parse_rule(identifier, raw):
    source = raw.decode('utf-8-sig')
    strings = re.findall(r'\(assert\s+\(3-digit-pattern\s+"([.X]{81})"\)\)',source)
    require(len(strings) == 2, identifier+': expected normal and transposed activation')
    cells = [i for i,c in enumerate(strings[0]) if c == 'X']
    transpose = sorted((i % 9)*9+i//9 for i in cells)
    require([i for i,c in enumerate(strings[1]) if c == 'X'] == transpose, identifier+': transposition')
    require(len(cells) == int(re.match(r'EL(\d+)c',identifier).group(1)), identifier+': advertised cell count')
    for suffix, expected in [('',cells),('s',transpose)]:
        marker = '(defrule 3-digit-pattern-rule-'+identifier+suffix+'\n'
        require(source.count(marker) == 1, identifier+': rule marker')
        lhs = source.split(marker,1)[1].split('=>',1)[0]
        candidates = [line for line in lhs.splitlines() if '(candidate ' in line]
        observed = Counter()
        for line in candidates:
            require('(status cand)' in line, identifier+': candidate status semantics')
            nb = re.search(r'\(number \?nb([123])\)',line)
            row = re.search(r'\(row \?row(\d+)',line)
            col = re.search(r'\(column \?col(\d+)',line)
            require(nb and row and col, identifier+': unexpected candidate syntax')
            cell = 9*(int(row.group(1))-1)+int(col.group(1))-1
            observed[(cell,int(nb.group(1)))] += 1
        require(observed == Counter((c,d) for c in expected for d in [1,2,3]), identifier+': full three-candidate premises')
        require('find the additional candidates' in source and '?real-guardians' in source,
                identifier+': guardian construction absent')
    return cells


def fetch_and_parse(cache, workers=8):
    paths = expected_paths()
    require(len(paths) == 592 and sum(COUNTS.values()) == 630 and
            sum(map(len,EXCLUDED.values())) == 38, '630=592+38 coverage')
    cache.mkdir(parents=True,exist_ok=True)

    def one(item):
        identifier,path = item
        destination = cache/(identifier+'.clp')
        if not destination.exists():
            url='https://raw.githubusercontent.com/'+REPO+'/'+COMMIT+'/'+path
            last = None
            for _ in range(3):
                try:
                    destination.write_bytes(urllib.request.urlopen(url,timeout=15).read())
                    break
                except Exception as error:
                    last=error
            else:
                raise last
        raw=destination.read_bytes()
        return {'id':identifier,'path':path,'sha256':hashlib.sha256(raw).hexdigest(),
                'cells':parse_rule(identifier,raw),'bytes':len(raw)}

    with ThreadPoolExecutor(max_workers=workers) as executor:
        rows=list(executor.map(one,paths))
    require(len({r['id'] for r in rows}) == 592 and len({tuple(r['cells']) for r in rows}) == 592,
            'duplicate public pattern ID/geometry')
    return rows


def main():
    require(__debug__, 'Run without -O')
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cache',type=Path,default=ROOT/'target/catalog630-rules')
    parser.add_argument('--workers',type=int,default=8)
    args=parser.parse_args()
    rows=fetch_and_parse(args.cache,args.workers)
    manual=ROOT/'target/catalog630-UMRN.pdf'
    result={'schema_version':1,'repository':REPO,'commit':COMMIT,
            'primary_url':'https://github.com/'+REPO+'/tree/'+COMMIT+'/SudoRules-V20.1/EXOTIC/Imp630',
            'manual_url':'https://github.com/'+REPO+'/blob/'+COMMIT+'/Docs/2023-UMRN-r1.pdf',
            'manual_sha256':hashlib.sha256(manual.read_bytes()).hexdigest() if manual.exists() else None,
            'total_catalogue_ids':630,'public_geometries_parsed':592,
            'excluded_no_free_cell_ids':{'{}c'.format(n):ids for n,ids in EXCLUDED.items()},
            'original_counts':COUNTS,'patterns':rows,
            'semantics':'Each public rule requires all three distinct palette candidates at every mapped pattern cell, and forms an OR of additional candidates. Missing palette candidates are not allowed by this implementation.',
            'coverage_limit':'All630 IDs accounted; only592 public geometries parsed.38 author-excluded no-free-cell geometries unavailable here. Strict-rule inapplicability to unique givens puzzles is a separate theorem, not a fabricated parse of those38.'}
    (HERE/'catalog630_catalogue.json').write_text(json.dumps(result,separators=(',',':'))+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='patterns'},indent=2))


if __name__=='__main__':
    main()
