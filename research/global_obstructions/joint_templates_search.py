"""Exact bounded compatibility propagation between full native digit templates.

This is a known template-CSP/overlay reasoning scheme, not a new Sudoku rule.
A numerical solver is never used. Refutation events identify unsupported full
one-digit templates; an independent checker must reconstruct every claim.
"""
import argparse
from functools import lru_cache
import hashlib
import itertools
import json
from pathlib import Path
import time

HERE=Path(__file__).resolve().parent


def canonical(value):return json.dumps(value,sort_keys=True,separators=(',',':')).encode()
def digest(value):return hashlib.sha256(canonical(value)).hexdigest()


def source(case):
    state=case.get('assumed_state',case.get('state',case))
    domains=state.get('domains',case.get('assumed_domains',case.get('base_domains')))
    if domains is None:raise ValueError('Missing domains')
    domains=domains.copy()
    if len(domains)!=81 or any(type(m) is not int or not 0<m<=511 for m in domains):raise ValueError('Invalid81 masks')
    t=case['target'];target={'cell':t['cell'],'digit':t['digit'],'assumed':t.get('assumed',t.get('value'))}
    c,d,v=target.values()
    if type(c)is not int or not 0<=c<81 or type(d)is not int or not 1<=d<=9 or type(v)is not bool:raise ValueError('Malformed target')
    bit=1<<(d-1)
    domains[c]=domains[c]&bit if v else domains[c]&~bit
    if not domains[c]:raise ValueError('Immediate empty target domain')
    payload={'schema':'joint-template-source-v1','domains':domains,'target':target}
    return domains,target,digest(payload)


@lru_cache(None)
def universe():
    answer=[]
    for columns in itertools.permutations(range(9)):
        if all(len({columns[r]//3 for r in range(start,start+3)})==3 for start in (0,3,6)):
            cells=tuple(9*r+c for r,c in enumerate(columns));answer.append((sum(1<<c for c in cells),cells))
    assert len(answer)==46656
    return tuple(answer)


def initial_templates(domains):
    tables=[]
    for d in range(9):
        allowed=sum(1<<c for c,m in enumerate(domains) if m>>d&1)
        forced=sum(1<<c for c,m in enumerate(domains) if m==1<<d)
        tables.append([(mask,cells) for mask,cells in universe() if not mask&~allowed and mask&forced==forced])
    return tables


def snapshot(tables,active):
    lists=[[list(cells) for i,(_,cells) in enumerate(table) if active[d]>>i&1] for d,table in enumerate(tables)]
    return {'counts':[len(x) for x in lists],'hashes':[digest(x) for x in lists]}


class Cutoff(Exception):pass


def search(case,seconds=30,max_per_digit=2000,max_pair_checks=20000000,max_support_checks=20000000,order=3):
    if order not in (2,3):raise ValueError('Support order must be2 or3')
    started=time.monotonic();domains,target,binding=source(case);tables=initial_templates(domains)
    active=[(1<<len(t))-1 for t in tables];initial=snapshot(tables,active);events=[];stats={'pair_mask_comparisons':0,'support_checks':0,'pair_rounds':0,'triple_rounds':0}
    proof={'schema_version':1,'kind':'joint-template-refutation-v1','source_hash':binding,
      'source_domains':domains,'target':target,'initial':initial,'pairwise_fixed_point':None,'events':events,'contradiction_digit':None}
    caches={}
    def budget():
        if time.monotonic()-started>=seconds:raise Cutoff('wall-time limit')
        if stats['pair_mask_comparisons']>max_pair_checks:raise Cutoff('pair-comparison limit')
        if stats['support_checks']>max_support_checks:raise Cutoff('support-check limit')
    def compatibility(a,b):
        if (a,b) not in caches:
            budget();matrix=[]
            for amask,_ in tables[a]:
                allowed=0
                for j,(bmask,_) in enumerate(tables[b]):
                    if not amask&bmask:allowed|=1<<j
                matrix.append(allowed);stats['pair_mask_comparisons']+=len(tables[b]);budget()
            caches[a,b]=matrix
        return caches[a,b]
    def indices(mask):
        while mask:
            bit=mask&-mask;mask-=bit;yield bit.bit_length()-1
    def remove(a,partners,bad):
        if not bad:return False
        events.append({'digit':a+1,'partners':[d+1 for d in partners],
          'templates':[list(tables[a][i][1]) for i in indices(bad)]})
        active[a]&=~bad
        if not active[a]:proof['contradiction_digit']=a+1
        return True
    def pair_close():
        while True:
            stats['pair_rounds']+=1;changed=False
            for a in range(9):
                for b in range(9):
                    if a==b:continue
                    matrix=compatibility(a,b);bad=0
                    for i in indices(active[a]):
                        stats['support_checks']+=1
                        if not matrix[i]&active[b]:bad|=1<<i
                    budget();changed|=remove(a,[b],bad)
                    if proof['contradiction_digit'] is not None:return
            if not changed:return
    try:
        budget()
        if any(len(t)>max_per_digit for t in tables):raise Cutoff('templates-per-digit limit')
        if any(not t for t in tables):proof['contradiction_digit']=next(d+1 for d,t in enumerate(tables) if not t)
        else:pair_close()
        if proof['contradiction_digit'] is None:
            proof['pairwise_fixed_point']={'event_index':len(events),**snapshot(tables,active)}
            if order>=3:
                while True:
                    stats['triple_rounds']+=1;changed=False
                    for a in range(9):
                        for b,c in itertools.combinations([d for d in range(9) if d!=a],2):
                            ab=compatibility(a,b);ac=compatibility(a,c);bc=compatibility(b,c);bad=0
                            for i in indices(active[a]):
                                bs=ab[i]&active[b];cs=ac[i]&active[c];supported=False
                                for j in indices(bs):
                                    stats['support_checks']+=1
                                    if bc[j]&cs:supported=True;break
                                if not supported:bad|=1<<i
                            budget();changed|=remove(a,[b,c],bad)
                            if proof['contradiction_digit'] is not None:break
                        if proof['contradiction_digit'] is not None:break
                    if proof['contradiction_digit'] is not None or not changed:break
                    pair_close()
                    if proof['contradiction_digit'] is not None:break
        status='refuted' if proof['contradiction_digit'] is not None else 'support-fixed-point'
        reason=None
    except Cutoff as e:status='cutoff';reason=str(e)
    proof['final']=snapshot(tables,active)
    return {'schema_version':1,'status':status,'cutoff_reason':reason,'order_limit':order,
      'seconds':time.monotonic()-started,'limits':{'seconds':seconds,'max_per_digit':max_per_digit,'max_pair_checks':max_pair_checks,'max_support_checks':max_support_checks},
      'stats':stats,'proof':proof,'scope':'Known template-CSP support propagation; a fixedpoint or cutoff is not a Sudoku solution or novelty claim.'}


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('fixture',type=Path);p.add_argument('--seconds',type=float,default=30);p.add_argument('--order',type=int,choices=(2,3),default=3)
    p.add_argument('--max-per-digit',type=int,default=2000);p.add_argument('--max-pair-checks',type=int,default=20000000);p.add_argument('--max-support-checks',type=int,default=20000000);p.add_argument('--output',type=Path)
    a=p.parse_args();result=search(json.loads(a.fixture.read_text()),a.seconds,a.max_per_digit,a.max_pair_checks,a.max_support_checks,a.order)
    if a.output:a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='proof'}|{'initial_counts':result['proof']['initial']['counts'],'final_counts':result['proof']['final']['counts'],'events':len(result['proof']['events']),'pairwise_fixed_point':result['proof']['pairwise_fixed_point'] is not None}))

if __name__=='__main__':main()
