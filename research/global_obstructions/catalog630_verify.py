#!/usr/bin/env python3
"""Check derived catalogue coverage, exact native maps and optional full rescans.

The38 omitted geometries are author-reported no-free-cell exclusions, not
parsed source. Their strict-native-rule inapplicability is a theorem with a
stated provenance premise; this checker does not turn it into geometric data.
"""
import argparse
import copy
import hashlib
import json
from pathlib import Path
import catalog630_match as matching
import catalog630_sources as source

HERE=Path(__file__).resolve().parent


def require(ok,message):
    if not ok:raise ValueError(message)


def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def native_map(pattern,domains,record):
    pairs=record['mapping'];mapping=dict(pairs)
    require(len(pairs)==len(pattern)and set(mapping)==set(pattern),'mapped pattern cells')
    require(len(set(mapping.values()))==len(pattern),'mapping injection')
    require(type(record['transpose'])is bool,'orientation')
    require(all(type(x)is int and 0<=x<81 for x in mapping.values()),'native cell range')
    # Independent coordinates; preservation of equality/non-equality of row,
    # column, band and stack extends to permutations of each unused group.
    for a in pattern:
        ar,ac=divmod(a,9)
        if record['transpose']:ar,ac=ac,ar
        for b in pattern:
            br,bc=divmod(b,9)
            if record['transpose']:br,bc=bc,br
            x,y=divmod(mapping[a],9);u,v=divmod(mapping[b],9)
            require((ar==br,ac==bc,ar//3==br//3,ac//3==bc//3)==
                    (x==u,y==v,x//3==u//3,y//3==v//3),'native geometry')
    digits=record['palette']
    require(len(digits)==3 and len(set(digits))==3 and all(type(d)is int and 1<=d<=9 for d in digits),'palette')
    bits=sum(1<<(d-1)for d in digits)
    require(all(domains[c]&bits==bits for c in mapping.values()),'three actual candidate premises')
    guardians=sorted([c,d]for c in mapping.values()for d in range(1,10)
                     if domains[c]&(1<<(d-1))and d not in digits)
    require(record['guardians']==guardians,'exact complete guardian OR')
    return True


def verify(rescan=False,source_cache=None):
    cat_path=HERE/'catalog630_catalogue.json';fixture=HERE/'different_lead.json'
    catalogue=json.loads(cat_path.read_text());case=json.loads(fixture.read_text())
    result=json.loads((HERE/'catalog630_match_result.json').read_text())
    require(result['fixture_sha256']==digest(fixture)and result['catalogue_sha256']==digest(cat_path),'source binding')
    require(catalogue['commit']==source.COMMIT and catalogue['repository']==source.REPO,'pinned source')
    rows=catalogue['patterns'];expected=dict(source.expected_paths())
    require(len(rows)==592 and len({r['id']for r in rows})==592,'592 distinct parsed IDs')
    require({r['id']:r['path']for r in rows}==expected,'complete public ID/path coverage')
    require(catalogue['original_counts']=={str(k):v for k,v in source.COUNTS.items()},'630 original family sizes')
    require(catalogue['excluded_no_free_cell_ids']=={str(k)+'c':v for k,v in source.EXCLUDED.items()},'38 author exclusions')
    require(sum(source.COUNTS.values())==630 and sum(map(len,source.EXCLUDED.values()))==38,'630 accounting')
    for row in rows:
        cells=row['cells']
        require(cells==sorted(set(cells))and all(type(c)is int and 0<=c<81 for c in cells),'pattern cells')
        require(len(cells)==int(row['id'][2:].split('c')[0]),'advertised cell count')
        # This independently proves the OR clause, even without trusting CLIPS.
        require(matching.proper_coloring(cells,3)is None,'catalogue graph must be non-three-colorable')
        if source_cache:
            raw=(source_cache/(row['id']+'.clp')).read_bytes()
            require(hashlib.sha256(raw).hexdigest()==row['sha256'],'pinned source bytes')
            require(source.parse_rule(row['id'],raw)==cells,'source parsing')
    raw=case['base_domains'].copy();target=case['target'];bit=1<<(target['digit']-1)
    raw[target['cell']]=raw[target['cell']]&bit if target['value']else raw[target['cell']]&~bit
    snapshots={'selected20_raw_opposite':(raw,case['graph_cells']),
               'selected20_closed_opposite':(case['assumed_domains'],case['graph_cells']),
               'full_base':(case['base_domains'],None),'full_raw_opposite':(raw,None),
               'full_closed_opposite':(case['assumed_domains'],None)}
    byid={row['id']:row['cells']for row in rows};checks=0
    require(set(result['runs'])==set(snapshots),'snapshot coverage')
    for name,(domains,allowed)in snapshots.items():
        run=result['runs'][name]
        require(run['domains_sha256']==hashlib.sha256(json.dumps(domains,separators=(',',':')).encode()).hexdigest(),'snapshot binding')
        require(run['search_complete']is True,'complete search marker')
        for occurrence in run['occurrences']:
            native_map(byid[occurrence['id']],domains,occurrence)
            require(allowed is None or set(dict(occurrence['mapping']).values())<=set(allowed),'selected support scope')
            checks+=1
        if rescan:require(json.loads(json.dumps(matching.scan(rows,domains,allowed)))=={k:v for k,v in run.items()if k!='domains_sha256'},'exact complete rescan')
    # A unique Sudoku needs >=8 given labels: swapping two absent labels would
    # produce a distinct solution. This fixture has eight given labels. The
    # source-reported38 no-free-cell shapes would forbid givens of all3 palette
    # labels under strict premises, under every native symmetry. Hence absent.
    given_labels={int(c)for c in case['puzzle']if c in '123456789'}
    require(len(given_labels)>=8,'strict no-free-cell exclusion premise')
    positive=result['runs']['full_base']['occurrences']
    require(len(positive)==1 and positive[0]['id']=='EL12c38'and positive[0]['guardians']==[[68,7]],'known OR1 control')
    mutations=0
    for name,change in [('missing_guardian',lambda r:r['guardians'].clear()),
                        ('invalid_palette',lambda r:r['palette'].__setitem__(0,0)),
                        ('duplicate_map',lambda r:r['mapping'].__setitem__(0,r['mapping'][1])),
                        ('wrong_orientation',lambda r:r.__setitem__('transpose',not r['transpose']))]:
        bad=copy.deepcopy(positive[0]);change(bad)
        try:native_map(byid[bad['id']],case['base_domains'],bad)
        except (ValueError,IndexError):mutations+=1
        else:raise ValueError('Mutation accepted: '+name)
    return {'parsed_public_patterns':len(rows),'all_catalogue_ids_accounted':630,
            'author_excluded_geometries_not_parsed':38,'non_three_colorable_graphs_checked':len(rows),
            'occurrence_maps_checked':checks,'negative_mutations_rejected':mutations,
            'full_rescan':rescan,'given_digit_labels':sorted(given_labels),
            'coverage':'Strict all-three-candidate native-geometry rules only; not OR-chain closure or guarded sparse variants of omitted38.'}


def main():
    require(__debug__,'Run without -O')
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--rescan',action='store_true')
    p.add_argument('--source-cache',type=Path);a=p.parse_args()
    print(json.dumps(verify(a.rescan,a.source_cache),sort_keys=True))


if __name__=='__main__':main()
