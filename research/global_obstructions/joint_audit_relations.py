#!/usr/bin/env python3
"""Independent guardian provenance and signed binary implication acceptance.

No strong_baseline module is imported. Complete coverage, when requested, is
recomputed rather than inferred from a supplied `complete` flag.
"""
import argparse
import copy
from functools import lru_cache
import hashlib
import json
from pathlib import Path
import catalog630_match
from catalog630_verify import native_map

HERE=Path(__file__).resolve().parent
CAT_PATH=HERE/'catalog630_catalogue.json'
CAT=json.loads(CAT_PATH.read_text());CAT_HASH=hashlib.sha256(CAT_PATH.read_bytes()).hexdigest()
PATTERNS={p['id']:p for p in CAT['patterns']}


@lru_cache(None)
def impossible(cells):return catalog630_match.proper_coloring(cells,3)is None


@lru_cache(maxsize=100000)
def occurrence_valid(identifier,origin,palette,transpose,mapping,guardians):
    pattern=PATTERNS[identifier]['cells']
    native_map(pattern,list(origin),{'palette':list(palette),'transpose':transpose,
               'mapping':[list(p)for p in mapping],'guardians':[list(g)for g in guardians]})
    need(impossible(tuple(pattern)),'native impossible-pattern premise')
    return True


def need(ok,message):
    if not ok:raise ValueError(message)


def digest(value):return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(',',':')).encode()).hexdigest()
def peer(a,b):return a!=b and (a//9==b//9 or a%9==b%9 or (a//27,a%9//3)==(b//27,b%9//3))


def house(h):
    need(type(h)is int and 0<=h<27,'native house')
    if h<9:return [9*h+c for c in range(9)]
    if h<18:return [9*r+h-9 for r in range(9)]
    box=h-18;return [9*(3*(box//3)+r)+3*(box%3)+c for r in range(3)for c in range(3)]


def masks(domains):need(type(domains)is list and len(domains)==81 and all(type(m)is int and 0<m<=511 for m in domains),'native domains')


def left(clause,domains):
    result=[]
    for c,d in clause['guardians']:
        if domains[c]==1<<(d-1):return None
        if domains[c]&(1<<(d-1)):result.append((c,d))
    return result


def signature(clause):
    return clause.get('catalogue_id',clause.get('id')),tuple(clause['palette']),tuple(sorted(dict(clause['mapping']).values()))


def relation_state(domains,relations,require_complete=False):
    need(__debug__,'Run without -O');masks(domains)
    need(relations['version']==1 and relations['catalogue_sha256']==CAT_HASH,'relation version/catalogue')
    for key,origin in relations['origins'].items():
        masks(origin);need(digest(origin)==key,'origin hash')
        need(all(not m&~old for m,old in zip(domains,origin)),'monotone clause ancestry')
    ids=set()
    for clause in relations['clauses']:
        need(clause['id']not in ids,'unique clause identity');ids.add(clause['id'])
        origin=relations['origins'][clause['origin']]
        occurrence_valid(clause['catalogue_id'],tuple(origin),tuple(clause['palette']),clause['transpose'],
                         tuple(tuple(p)for p in clause['mapping']),tuple(tuple(g)for g in clause['guardians']))
        need(clause['id']==digest({k:v for k,v in clause.items()if k!='id'}),'immutable clause identity')
    if require_complete:
        valid=[]
        for row in relations['coverage']:
            if row['complete']is not True or row['catalogue_sha256']!=CAT_HASH:continue
            origin=relations['origins'][row['origin']]
            expected=catalog630_match.scan(CAT['patterns'],origin)['occurrences']
            actual=[c for c in relations['clauses']if c['origin']==row['origin']]
            need({signature(c)for c in actual}=={signature(c)for c in expected},'coverage lost a persistent clause')
            valid.append(row)
        need(valid,'no independently complete ancestor scan')
    return {c['id']:c for c in relations['clauses']}


def decode(literal,domains):
    need(type(literal)is int and literal!=0 and abs(literal)<=729,'signed nonzero candidate')
    c,d=divmod(abs(literal)-1,9);d+=1
    need(domains[c]&(1<<(d-1)),'literal is currently available')
    return c,d


def edge(step,domains,clauses):
    u,v=step['from'],step['to'];a,d=decode(u,domains);b,e=decode(v,domains);reason=step['reason'];kind=reason['kind']
    if kind=='native-conflict':
        need(u>0 and v<0 and ((a==b and d!=e)or(d==e and peer(a,b))),'native weak implication')
        need({tuple(reason['a']),tuple(reason['b'])}=={(a,d),(b,e)},'native reason endpoints');return
    need(u<0 and v>0,'positive OR implication signs')
    if kind=='cell-or':
        c=reason['cell'];need(type(c)is int and 0<=c<81,'cell OR')
        literals={(c,n)for n in range(1,10)if domains[c]&(1<<(n-1))}
    elif kind=='house-or':
        n=reason['digit'];need(type(n)is int and 1<=n<=9,'house OR digit')
        literals={(c,n)for c in house(reason['house'])if domains[c]&(1<<(n-1))}
    elif kind=='guardian-or':
        remaining=left(clauses[reason['clause']],domains)
        need(remaining is not None,'guardian OR already satisfied');literals=set(remaining)
    else:raise ValueError('Unknown edge reason')
    need(len(literals)==2 and {(a,d),(b,e)}==literals,'exact binary OR endpoints')


def step(domains,relations,proof):
    clauses=relation_state(domains,relations);kind=proof['kind']
    if kind!='guardian-contradiction':
        need(type(proof['cell'])is int and 0<=proof['cell']<81 and type(proof['digit'])is int and 1<=proof['digit']<=9,'native inference endpoint')
    if kind in ('guardian-unit','guardian-contradiction'):
        remaining=left(clauses[proof['clause']],domains);need(remaining is not None,'satisfied clause cannot force')
        if kind=='guardian-contradiction':need(not remaining,'guardian contradiction');return True
        need(remaining==[(proof['cell'],proof['digit'])]and proof['value']is True,'exact guardian unit');return True
    need(kind=='binary-forcing','supported proof kind')
    c,d,value=proof['cell'],proof['digit'],proof['value'];need(type(value)is bool,'conclusion polarity')
    literal=9*c+d;decode(literal,domains);start,goal=(-literal,literal)if value else(literal,-literal)
    path=proof['path'];need(path and path[0]['from']==start and path[-1]['to']==goal,'failed-literal path endpoints')
    for i,entry in enumerate(path):
        if i:need(path[i-1]['to']==entry['from'],'connected implication path')
        edge(entry,domains,clauses)
    return True


def self_test():
    case=json.loads((HERE/'different_lead.json').read_text());domains=case['base_domains']
    occurrence=json.loads((HERE/'catalog630_match_result.json').read_text())['runs']['full_base']['occurrences'][0]
    origin=digest(domains)
    clause={'origin':origin,'catalogue_id':occurrence['id'],**{k:occurrence[k]for k in ('palette','transpose','mapping','guardians')}}
    clause['id']=digest(clause)
    relations={'version':1,'context_id':'known-lead-control','catalogue_sha256':CAT_HASH,'origins':{origin:domains},
               'clauses':[clause],'coverage':[{'origin':origin,'catalogue_sha256':CAT_HASH,'complete':True}]}
    relation_state(domains,relations,True)
    unit={'kind':'guardian-unit','clause':clause['id'],'cell':68,'digit':7,'value':True}
    step(domains,relations,unit);step(case['assumed_domains'],relations,unit)
    false=case['assumed_domains'].copy();false[68]&=511^(1<<6)
    step(false,relations,{'kind':'guardian-contradiction','clause':clause['id']})
    rejected=0
    for edit in [lambda x:x['clauses'].clear(),lambda x:x['clauses'][0]['guardians'].clear(),
                 lambda x:x['clauses'][0].__setitem__('id','0'*64),lambda x:x['coverage'][0].__setitem__('complete',False)]:
        bad=copy.deepcopy(relations);edit(bad)
        try:relation_state(domains,bad,True)
        except (ValueError,KeyError):rejected+=1
        else:raise ValueError('Corrupt relation accepted')
    bad=copy.deepcopy(unit);bad['digit']=8
    try:step(domains,relations,bad)
    except ValueError:rejected+=1
    else:raise ValueError('Wrong guardian unit accepted')
    # A nonvacuous native pair control: cells0,1 are {1,2}, cell2 {1,3}.
    # The five-edge implication chain proves cell2!=1; many completions exist.
    small=[511]*81;small[0]=small[1]=3;small[2]=5
    empty={'version':1,'catalogue_sha256':CAT_HASH,'origins':{},'clauses':[],'coverage':[]}
    path=[{'from':19,'to':-1,'reason':{'kind':'native-conflict','a':[2,1],'b':[0,1]}},
          {'from':-1,'to':2,'reason':{'kind':'cell-or','cell':0}},
          {'from':2,'to':-11,'reason':{'kind':'native-conflict','a':[0,2],'b':[1,2]}},
          {'from':-11,'to':10,'reason':{'kind':'cell-or','cell':1}},
          {'from':10,'to':-19,'reason':{'kind':'native-conflict','a':[1,1],'b':[2,1]}}]
    binary={'kind':'binary-forcing','cell':2,'digit':1,'value':False,'path':path}
    step(small,empty,binary)
    for edit in [lambda x:x['path'][0].__setitem__('to',1),
                 lambda x:x['path'][1]['reason'].__setitem__('cell',2),
                 lambda x:x.__setitem__('value',True)]:
        bad=copy.deepcopy(binary);edit(bad)
        try:step(small,empty,bad)
        except ValueError:rejected+=1
        else:raise ValueError('Wrong binary implication accepted')
    # Independently exercise the exact native-house OR endpoint convention.
    for h in range(27):
        found=False
        for d in range(1,10):
            candidates=[c for c in house(h)if domains[c]&(1<<(d-1))]
            if len(candidates)==2:
                x,y=[9*c+d for c in candidates]
                edge({'from':-x,'to':y,'reason':{'kind':'house-or','house':h,'digit':d}},domains,{})
                found=True;break
        if found:break
    need(found,'house OR sanity control')
    return {'complete_ancestor_scans_recomputed':1,'persistent_unit_snapshots_checked':2,
            'guardian_false_refutations_checked':1,'corruptions_rejected':rejected,
            'binary_forcing_paths_checked':1,'binary_path_edges_checked':5,'native_house_or_edges_checked':1,
            'scope':'Exact clause provenance and implication-step soundness; not a proof that a full forcing menu reached a fixed point.'}


def main():
    need(__debug__,'Run without -O');p=argparse.ArgumentParser(description=__doc__);p.add_argument('--self-test',action='store_true');a=p.parse_args()
    need(a.self_test,'Use --self-test or import exact acceptance functions')
    print(json.dumps(self_test(),sort_keys=True))


if __name__=='__main__':main()
