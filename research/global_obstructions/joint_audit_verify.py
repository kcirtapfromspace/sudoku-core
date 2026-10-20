#!/usr/bin/env python3
"""Independent exact acceptance of ordered joint-template refutations.

No searcher import, numerical solver, or supplied template universe is trusted.
Enumeration recursively selects a cell in each row using column and box masks.
"""
import argparse
import copy
import hashlib
import itertools
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent


def need(condition,message):
    if not condition:raise ValueError(message)


def canonical(obj):return json.dumps(obj,sort_keys=True,separators=(',',':')).encode()
def sha(obj):return hashlib.sha256(canonical(obj)).hexdigest()


def source_state(case):
    state=case.get('assumed_state',case.get('state',case))
    domains=state.get('domains')
    if domains is None:domains=case.get('assumed_domains',case.get('base_domains'))
    need(type(domains)is list and len(domains)==81,'81 source masks')
    need(all(type(m)is int and 0<m<=511 for m in domains),'nonempty nine-bit source masks')
    if 'assumed_state'in case and 'assumed_domains'in case:
        need(domains==case['assumed_domains'],'conflicting fixture snapshots')
    values=state.get('values')
    if values is not None:
        need(len(values)==81 and all(type(v)is int and 0<=v<=9 for v in values),'source values')
        need(all(not v or domains[c]==1<<(v-1)for c,v in enumerate(values)),'placed values and domains')
    t=case['target'];target={'cell':t['cell'],'digit':t['digit'],'assumed':t.get('assumed',t.get('value'))}
    c,d,v=target['cell'],target['digit'],target['assumed']
    need(type(c)is int and 0<=c<81 and type(d)is int and 1<=d<=9 and type(v)is bool,'target type/range')
    effective=domains.copy();effective[c]&=(1<<(d-1))if v else 511^(1<<(d-1))
    need(effective[c]>0,'nonempty opposite literal')
    return effective,target


def templates(domains,digit):
    bit=1<<(digit-1);rows=[]
    for r in range(9):
        forced=[9*r+c for c in range(9)if domains[9*r+c]==bit]
        if len(forced)>1:return []
        rows.append(forced or [9*r+c for c in range(9)if domains[9*r+c]&bit])
    answer=[];path=[]
    def visit(r,columns,boxes):
        if r==9:answer.append(tuple(path));return
        for cell in rows[r]:
            col=cell%9;box=3*(r//3)+col//3
            if columns&(1<<col)or boxes&(1<<box):continue
            path.append(cell);visit(r+1,columns|(1<<col),boxes|(1<<box));path.pop()
    visit(0,0,0)
    return answer


def summary(active):
    return {'counts':[len(x)for x in active],
            'hashes':[sha([list(t)for t in sorted(x)])for x in active]}


def masks_for(active):return [{t:sum(1<<c for c in t)for t in domain}for domain in active]


def pair_fixed(active):
    need(all(active),'pairwise checkpoint has an empty domain')
    masks=masks_for(active)
    for a in range(9):
        for b in range(9):
            if a==b:continue
            for am in masks[a].values():
                need(any(not am&bm for bm in masks[b].values()),'checkpoint lacks pairwise support')


def verify(case,wrapper):
    need(__debug__,'Run without -O')
    need(wrapper.get('status')=='refuted'and wrapper.get('cutoff_reason')is None,'only completed refutations accepted')
    proof=wrapper['proof'];need(proof['schema_version']==1 and proof['kind']=='joint-template-refutation-v1','certificate version/kind')
    domains,target=source_state(case)
    need(proof['source_domains']==domains and proof['target']==target,'effective snapshot/target binding')
    need(proof['source_hash']==sha({'schema':'joint-template-source-v1','domains':domains,'target':target}),'source hash')
    active=[set(templates(domains,d))for d in range(1,10)]
    need(proof['initial']==summary(active),'complete initial template universe')
    initial_counts=[len(x)for x in active];events=proof['events'];checkpoint=proof['pairwise_fixed_point']
    index=None
    if checkpoint is not None:
        index=checkpoint['event_index'];need(type(index)is int and 0<=index<=len(events),'checkpoint event index')
    triple_events=[i for i,e in enumerate(events)if len(e['partners'])==2]
    need(type(wrapper['order_limit'])is int and wrapper['order_limit']in (2,3),'declared support order')
    if triple_events:need(wrapper['order_limit']==3,'triple event outside declared order')
    if triple_events:need(index is not None and index<=min(triple_events),'checkpoint before ternary phase')
    removals=0
    for i in range(len(events)+1):
        if i==index:
            need({k:checkpoint[k]for k in ('counts','hashes')}==summary(active),'checkpoint active universe')
            need(all(len(e['partners'])==1 for e in events[:i]),'pair-only checkpoint history')
            pair_fixed(active)
        if i==len(events):break
        need(all(active),'events after first contradiction')
        event=events[i];a=event['digit'];partners=event['partners']
        need(type(a)is int and 1<=a<=9,'event digit')
        need(type(partners)is list and len(partners)in (1,2)and len(set(partners))==len(partners),'partner arity/distinctness')
        need(all(type(d)is int and 1<=d<=9 and d!=a for d in partners),'partner digits')
        batch=event['templates'];need(type(batch)is list and len(batch)>0,'nonempty event batch')
        tuples=[]
        for cells in batch:
            need(type(cells)is list and len(cells)==9 and all(type(c)is int for c in cells),'template syntax')
            t=tuple(cells);need(t in active[a-1],'template currently active');tuples.append(t)
        need(len(tuples)==len(set(tuples)),'no duplicate batch deletion')
        partner_masks=[[sum(1<<c for c in t)for t in active[d-1]]for d in partners]
        for t in tuples:
            tm=sum(1<<c for c in t)
            if len(partners)==1:
                need(not any(not tm&u for u in partner_masks[0]),'claimed pair deletion has support')
            else:
                # Exhaustive exact join; no supplied incompatibility table trusted.
                need(not any(not tm&u and not (tm|u)&v for u in partner_masks[0]for v in partner_masks[1]),'claimed triple deletion has support')
        active[a-1].difference_update(tuples);removals+=len(tuples)
    need(proof['final']==summary(active),'final active universe')
    contradiction=proof['contradiction_digit']
    need(type(contradiction)is int and 1<=contradiction<=9 and not active[contradiction-1],'empty contradiction domain')
    return {'accepted':True,'initial_counts':initial_counts,'events':len(events),'templates_removed':removals,
            'pairwise_fixed_point_checked':checkpoint is not None,'triple_events':len(triple_events),
            'separation':bool(checkpoint is not None and triple_events),
            'source_hash':proof['source_hash'],'scope':'Known template-CSP support propagation; exact two-versus-three support separation only when flagged.'}


def mutation_checks(case,wrapper):
    proof=wrapper['proof'];mutations=[]
    def add(name,edit):
        bad=copy.deepcopy(wrapper);edit(bad);mutations.append((name,bad))
    add('cutoff',lambda x:x.__setitem__('status','cutoff'))
    add('binding',lambda x:x['proof'].__setitem__('source_hash','0'*64))
    add('target',lambda x:x['proof']['target'].__setitem__('assumed',not x['proof']['target']['assumed']))
    add('source_domain',lambda x:x['proof']['source_domains'].__setitem__(0,0))
    add('initial_count',lambda x:x['proof']['initial']['counts'].__setitem__(0,x['proof']['initial']['counts'][0]+1))
    add('initial_hash',lambda x:x['proof']['initial']['hashes'].__setitem__(0,'0'*64))
    add('event_digit',lambda x:x['proof']['events'][0].__setitem__('digit',0))
    add('partner_self',lambda x:x['proof']['events'][0].__setitem__('partners',[x['proof']['events'][0]['digit']]))
    add('duplicate_deletion',lambda x:x['proof']['events'][0]['templates'].append(x['proof']['events'][0]['templates'][0]))
    add('empty_batch',lambda x:x['proof']['events'][0].__setitem__('templates',[]))
    add('wrong_final',lambda x:x['proof']['final']['counts'].__setitem__(0,999))
    add('wrong_contradiction',lambda x:x['proof'].__setitem__('contradiction_digit',0))
    domains,_=source_state(case);initial=[templates(domains,d)for d in range(1,10)]
    supported=None
    for a in range(9):
        for t in initial[a]:
            tm=set(t)
            for b in range(9):
                if a!=b and any(not tm.intersection(u)for u in initial[b]):
                    supported={'digit':a+1,'partners':[b+1],'templates':[list(t)]};break
            if supported:break
        if supported:break
    need(supported is not None,'supported deletion mutation control')
    add('delete_supported_template',lambda x:x['proof']['events'].__setitem__(0,supported))
    if proof['pairwise_fixed_point']is not None:
        add('checkpoint_count',lambda x:x['proof']['pairwise_fixed_point']['counts'].__setitem__(0,0))
        add('checkpoint_missing',lambda x:x['proof'].__setitem__('pairwise_fixed_point',None))
    for name,bad in mutations:
        try:verify(case,bad)
        except (ValueError,KeyError,TypeError,IndexError):pass
        else:raise ValueError('Corruption accepted: '+name)
    return len(mutations)


def main():
    need(__debug__,'Run without -O')
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('fixture',type=Path,nargs='?')
    parser.add_argument('certificate',type=Path,nargs='?');parser.add_argument('--mutations',action='store_true')
    parser.add_argument('--self-test',action='store_true');args=parser.parse_args()
    if args.self_test:
        need(len(templates([511]*81,1))==46656,'complete blank-board enumeration')
        records=json.loads((HERE/'joint_templates_controls.json').read_text())['records'];answers=[]
        need(len(records)==2,'both stable controls required')
        for record in records:
            fixture=HERE/record['fixture'];certificate=HERE/record['proof']
            need(hashlib.sha256(fixture.read_bytes()).hexdigest()==record['fixture_sha256'],'control fixture bytes')
            need(hashlib.sha256(certificate.read_bytes()).hexdigest()==record['proof_sha256'],'control proof bytes')
            case=json.loads(fixture.read_text());proof=json.loads(certificate.read_text());answer=verify(case,proof)
            answer['negative_mutations_rejected']=mutation_checks(case,proof)
            answer['fixture']=record['fixture'];answer['certificate']=record['proof'];answers.append(answer)
        need([r['separation']for r in answers]==[False,True],'pair/triple control classification')
        print(json.dumps({'blank_templates':46656,'records':answers,'mutations':sum(r['negative_mutations_rejected']for r in answers)},sort_keys=True));return
    need(args.fixture is not None and args.certificate is not None,'fixture and certificate required')
    case=json.loads(args.fixture.read_text());proof=json.loads(args.certificate.read_text())
    result=verify(case,proof)
    if args.mutations:result['negative_mutations_rejected']=mutation_checks(case,proof)
    print(json.dumps(result,sort_keys=True))


if __name__=='__main__':main()
