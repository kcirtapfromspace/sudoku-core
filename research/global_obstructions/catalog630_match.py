#!/usr/bin/env python3
"""Exact strict-public-rule embeddings and independent guardian checks.

Pair relations preserve rows, columns, bands and stacks (or their transpose),
so a complete injection extends to an ordinary Sudoku geometry symmetry.
This tests source rule premises, not arbitrary OR-chain resolution.
"""
import argparse
from collections import Counter
import hashlib
import itertools
import json
from pathlib import Path
import time

HERE=Path(__file__).resolve().parent


def relation(a,b,transpose=False):
    ar,ac=divmod(a,9);br,bc=divmod(b,9)
    if transpose:ar,ac,br,bc=ac,ar,bc,br
    return int(ar==br)+2*int(ac==bc)+4*int(ar//3==br//3)+8*int(ac//3==bc//3)


def peers(a,b):
    r,c=divmod(a,9);s,d=divmod(b,9)
    return a!=b and (r==s or c==d or (r//3,c//3)==(s//3,d//3))


def proper_coloring(cells,q):
    adjacency={c:[v for v in cells if peers(c,v)]for c in cells}
    assigned={}
    def visit():
        if len(assigned)==len(cells):return dict(assigned)
        choices={c:[d for d in range(q)if all(assigned.get(v)!=d for v in adjacency[c])]
                 for c in cells if c not in assigned}
        c=min(choices,key=lambda v:(len(choices[v]),-len(adjacency[v]),v))
        for d in choices[c]:
            assigned[c]=d
            answer=visit()
            if answer is not None:return answer
            del assigned[c]
        return None
    return visit()


def embeddings(pattern,host,transpose=False):
    if len(pattern)>len(host):return
    pr={(a,b):relation(a,b,transpose)for a in pattern for b in pattern}
    hr={(a,b):relation(a,b)for a in host for b in host}
    ps={a:Counter(pr[a,b]for b in pattern if b!=a)for a in pattern}
    hs={a:Counter(hr[a,b]for b in host if b!=a)for a in host}
    available={a:[b for b in host if all(hs[b][r]>=n for r,n in ps[a].items())]for a in pattern}
    mapping,used={},set()
    def visit():
        if len(mapping)==len(pattern):
            yield dict(mapping);return
        opts={a:[b for b in available[a]if b not in used and all(pr[a,x]==hr[b,y]for x,y in mapping.items())]
              for a in pattern if a not in mapping}
        a=min(opts,key=lambda x:(len(opts[x]),-sum(pr[x,z]!=0 for z in mapping),x))
        for b in opts[a]:
            mapping[a]=b;used.add(b)
            yield from visit()
            used.remove(b);del mapping[a]
    yield from visit()


def check_occurrence(pattern,domains,entry):
    pairs=entry['mapping'];mapping=dict(pairs)
    if len(pairs)!=len(pattern)or set(mapping)!=set(pattern)or len(set(mapping.values()))!=len(pattern):raise ValueError('pattern injection')
    if not all(relation(a,b,entry['transpose'])==relation(mapping[a],mapping[b])for a in pattern for b in pattern):raise ValueError('geometry')
    digits=entry['palette']
    if len(digits)!=3 or any(type(d)is not int or d<1 or d>9 for d in digits)or len(set(digits))!=3:raise ValueError('palette')
    if any(type(c)is not int or c<0 or c>=81 for c in mapping.values()):raise ValueError('board cell')
    mask=sum(1<<(d-1)for d in digits)
    if not all(domains[c]&mask==mask for c in mapping.values()):raise ValueError('strict actual lists')
    expected=sorted([c,d]for c in mapping.values()for d in range(1,10)if domains[c]&(1<<(d-1))and not mask&(1<<(d-1)))
    if entry['guardians']!=expected:raise ValueError('complete guardian OR')
    return True


def scan(patterns,domains,allowed_cells=None):
    allowed=set(range(81))if allowed_cells is None else set(allowed_cells)
    out=[];dedup=set();triple_hosts=[];mapping_calls=0
    minimum=min(len(p['cells'])for p in patterns)
    for digits in itertools.combinations(range(1,10),3):
        mask=sum(1<<(d-1)for d in digits)
        host=[c for c in sorted(allowed)if domains[c]&mask==mask]
        if len(host)<minimum:continue
        triple_hosts.append({'palette':digits,'cells':host})
        for p in patterns:
            for transpose in [False,True]:
                mapping_calls+=1
                for mapping in embeddings(p['cells'],host,transpose):
                    signature=(p['id'],digits,tuple(sorted(mapping.values())))
                    if signature in dedup:continue
                    dedup.add(signature)
                    guardians=sorted([c,d]for c in mapping.values()for d in range(1,10)
                                     if domains[c]&(1<<(d-1))and not mask&(1<<(d-1)))
                    entry={'id':p['id'],'palette':digits,'transpose':transpose,
                           'mapping':sorted(mapping.items()),'guardians':guardians}
                    check_occurrence(p['cells'],domains,entry)
                    out.append(entry)
    return {'occurrences':out,'eligible_triple_hosts':triple_hosts,'orientation_search_calls':mapping_calls,
            'search_complete':True,'scope':'All592 public strict all-three-candidate rules under standard Sudoku geometry symmetries and transpose. This is not an OR-chain closure.'}


def main():
    if not __debug__:raise ValueError('Run without -O')
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=HERE/'catalog630_match_result.json')
    args=parser.parse_args()
    catalogue=json.loads((HERE/'catalog630_catalogue.json').read_text())
    case_path=HERE/'different_lead.json';case=json.loads(case_path.read_text())
    base=case['base_domains'];raw=base.copy();target=case['target'];bit=1<<(target['digit']-1)
    raw[target['cell']]=raw[target['cell']]&bit if target['value']else raw[target['cell']]&~bit
    sources=[('selected20_raw_opposite',raw,case['graph_cells']),
             ('selected20_closed_opposite',case['assumed_domains'],case['graph_cells']),
             ('full_base',base,None),('full_raw_opposite',raw,None),
             ('full_closed_opposite',case['assumed_domains'],None)]
    started=time.monotonic();runs={}
    for name,domains,allowed in sources:
        runs[name]=scan(catalogue['patterns'],domains,allowed)
        runs[name]['domains_sha256']=hashlib.sha256(json.dumps(domains,separators=(',',':')).encode()).hexdigest()
        print(name,len(runs[name]['occurrences']),'occurrences',flush=True)
    result={'schema_version':1,'fixture_sha256':hashlib.sha256(case_path.read_bytes()).hexdigest(),
            'catalogue_sha256':hashlib.sha256((HERE/'catalog630_catalogue.json').read_bytes()).hexdigest(),
            'public_patterns':len(catalogue['patterns']),'runs':runs,'seconds':time.monotonic()-started}
    args.output.write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':main()
