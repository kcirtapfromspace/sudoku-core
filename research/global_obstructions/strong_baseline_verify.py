"""Stdlib chronological acceptance of strong-baseline positive reasoning.

Complete catalogue coverage is independently rescanned. Fixedpoint labels also
require a fresh configured menu pass; focused/weaker profiles stay explicit.
"""
import argparse
from copy import deepcopy
import hashlib
import json
from pathlib import Path

from strong_baseline_relations import digest,CATALOGUE_HASH
from joint_audit_relations import relation_state,step as relation_step,left,masks
from joint_audit_patterns import pattern_step
from deeper_forcing_verify import valid_trace,verify_node
from real_verify import close,replay_trace,HOUSES


def empty(context):
    return {'version':1,'context_id':context,'catalogue_sha256':CATALOGUE_HASH,'origins':{},'clauses':[],'coverage':[]}


def apply(domains,event):
    c,d=event['cell'],event['digit']
    assert type(c)is int and 0<=c<81 and type(d)is int and 1<=d<=9
    assert type(event['value'])is bool and domains[c]&(1<<(d-1))
    domains[c]=(1<<(d-1))if event['value']else domains[c]&~(1<<(d-1))
    assert domains[c]>0


def forcing(initial,relations,node,options,stats):
    domains,local,refuted=replay(initial,deepcopy(relations),node['trace'],options,stats)
    assert domains==node['domains']
    if node['kind']=='refuted':
        assert refuted
        return 0
    assert node['kind']=='split'and not refuted
    literals=[]
    if node['split_kind']=='cell':
        c=node['pivot'];assert type(c)is int and 0<=c<81
        assert domains[c]&(domains[c]-1)
        literals=[[c,d]for d in range(1,10)if domains[c]&(1<<(d-1))]
    else:
        assert node['split_kind']=='guardian'
        clauses=relation_state(domains,local)
        remaining=left(clauses[node['pivot']],domains)
        assert remaining is not None and 2<=len(remaining)<=options['forcing_guardian_width']
        literals=[list(x)for x in remaining]
    assert [b['literal']for b in node['branches']]==literals
    depths=[]
    for branch in node['branches']:
        c,d=branch['literal'];assumed=domains.copy();assumed[c]=1<<(d-1)
        depths.append(forcing(assumed,local,branch['proof'],options,stats))
    stats['forcing_splits']+=1
    return 1+max(depths)


def replay(initial,relations,trace,options,stats):
    domains=initial.copy();masks(domains);relation_state(domains,relations)
    refuted=False
    for event in trace:
        assert not refuted,'A refuted branch must terminate'
        kind=event['kind'];stats['events']+=1
        if kind in('catalogue-scan','history-catalogue-scan'):
            origin=event['origin_domains'];masks(origin)
            if kind=='catalogue-scan':assert origin==domains
            else:assert all(not m&~old for m,old in zip(domains,origin))
            key=digest(origin)
            row=event['scan']
            if 'origin'in row:relations['origins'][key]=origin.copy()
            known={c['id']for c in relations['clauses']}
            for clause in event['new_clauses']:
                assert clause['origin']==key and clause['id']not in known
                relations['clauses'].append(deepcopy(clause));known.add(clause['id'])
            assert row['status']in('complete','cutoff')
            if 'origin'in row:
                assert row['origin']==key and row['complete']==(row['status']=='complete')
                assert row['catalogue_sha256']==CATALOGUE_HASH
                relations['coverage'].append({k:v for k,v in row.items()if k!='status'})
            else:
                assert not event['new_clauses']and row.get('inherited_complete_coverage')is True
            relation_state(domains,relations)
        elif kind=='local':
            cert=event['certificate'];valid_trace(cert['closure_trace'])
            if cert['kind']=='fixedpoint':
                domains=replay_trace(domains,cert['closure_trace'])
                assert domains==cert['domains']and close(domains)==(domains,[])
            else:
                counters={k:0 for k in['nodes','splits','hall_leaves','three_house_leaves','closure_events','hall_blockers','boundary_signatures_covered']}
                assert verify_node(domains,cert,counters)==0
                domains=cert['domains'].copy();refuted=True
            stats['local_calls']+=1
        elif kind in('guardian-unit','guardian-contradiction','binary-forcing'):
            relation_step(domains,relations,event)
            if kind=='guardian-contradiction':refuted=True
            else:apply(domains,event)
            stats['relation_inferences']+=1
        elif kind=='pattern-elimination':
            pattern_step(domains,event);apply(domains,event);stats['pattern_inferences']+=1
        else:
            assert kind=='forcing-elimination'and event['value']is False
            c,d=event['cell'],event['digit'];assert type(c)is int and 0<=c<81 and type(d)is int and 1<=d<=9
            assert domains[c]&(1<<(d-1))and domains[c]&(domains[c]-1)
            assumed=domains.copy();assumed[c]=1<<(d-1)
            depth=forcing(assumed,relations,event['proof'],options,stats)
            assert depth<=options['forcing_depth']
            apply(domains,event);stats['forcing_inferences']+=1
    return domains,relations,refuted


def verify_result(result,check_fixedpoint=True):
    assert result['schema_version']==1
    initial=result['input_domains'];masks(initial)
    relations=deepcopy(result.get('input_relations',empty(result['relations']['context_id'])))
    relation_state(initial,relations)
    stats={k:0 for k in['events','local_calls','relation_inferences','pattern_inferences','forcing_inferences','forcing_splits']}
    domains,derived,refuted=replay(initial,relations,result['trace'],result['options'],stats)
    assert domains==result['domains']
    # Trace must establish every returned fact; no hidden future/branch fact.
    assert {c['id']for c in derived['clauses']}=={c['id']for c in result['relations']['clauses']}
    assert derived['origins']==result['relations']['origins']
    # Every recorded scan result must correspond to its chronological event.
    expected=[{k:v for k,v in row.items()if k!='status'}for row in derived['coverage']]
    assert expected==result['relations']['coverage']
    relation_state(domains,result['relations'])
    status=result['status'];assert status in('complete','refuted','fixedpoint','cutoff')
    assert (status=='refuted')==refuted
    if status in('complete','fixedpoint'):
        assert not any(p['status']=='cutoff'for p in result['phase_results'])
    complete_catalogue=False
    if result['options']['catalogue']and status in('complete','fixedpoint'):
        relation_state(domains,result['relations'],True);complete_catalogue=True
    if status=='complete':
        assert all(m&(m-1)==0 for m in domains)
        assert all(len({domains[c]for c in h})==9 for h in HOUSES)
    fresh=False
    if status=='fixedpoint'and check_fixedpoint:
        from strong_baseline_core import close_state,Options
        options=Options(**result['options']);options.history_snapshots=()
        repeat=close_state(domains,result['relations'],options)
        assert repeat['status']=='fixedpoint'and repeat['domains']==domains
        assert not any(e['kind']in('guardian-unit','binary-forcing','pattern-elimination','forcing-elimination')for e in repeat['trace'])
        fresh=True
    if status=='cutoff':
        assert result['error']or any(p['status']=='cutoff'for p in result['phase_results'])or any(e.get('scan',{}).get('status')=='cutoff'for e in result['trace'])
    full_scope=not result['options']['forcing_targets']
    profile=result['options']
    strong_profile=(full_scope and profile['forcing_depth']>=2 and all(profile[k]for k in
                    ['catalogue','unit_binary','junior','relaxed_junior','msls','forcing_fresh_catalogue']))
    return {**stats,'status':status,'positive_steps_exact':True,'catalogue_completeness_rescanned':complete_catalogue,
            'configured_fixedpoint_repeated':fresh,'full_candidate_scope':full_scope,
            'configured_phase_fixedpoint_candidate':status=='fixedpoint'and fresh and full_scope,
            'mandatory_strong_profile_enabled':strong_profile,
            'configured_strong_profile_fixedpoint':status=='fixedpoint'and fresh and strong_profile,
            'strong_separation_candidate':False,
            'required_external_bindings':['verified original givens history coverage','same-state corrected retained-old union fixedpoint'],
            'limits':'Configured profile only. Historicalgivens coverage and union with retainedoldengine require callerbinding; ORkwhips are not fullyimplemented.'}


def self_test():
    here=Path(__file__).resolve().parent
    control=json.loads((here/'strong_baseline_controls.json').read_text())
    case=json.loads((here/'different_lead.json').read_text())
    assert control['fixture_sha256']==hashlib.sha256((here/'different_lead.json').read_bytes()).hexdigest()
    for name in ['known20_refutation','focused_branch','persistent_branch']:
        assert control[name]['input_domains']==case['assumed_domains']
    assert control['persistent_origin']['input_domains']==case['base_domains']
    checks={name:verify_result(control[name])for name in ['known20_refutation','focused_branch','persistent_origin','persistent_branch']}
    assert checks['known20_refutation']['status']=='refuted'
    assert checks['persistent_branch']['status']=='refuted'
    assert checks['focused_branch']['mandatory_strong_profile_enabled']is False
    assert checks['persistent_origin']['mandatory_strong_profile_enabled']is False
    assert all(not c['strong_separation_candidate']for c in checks.values())
    trials=[]
    def alter(label,name,fn):
        copy=deepcopy(control[name]);fn(copy);trials.append((label,copy))
    def forced(r):return next(e for e in r['trace']if e['kind']=='forcing-elimination')['proof']
    def unit(r):return next(e for e in r['trace']if e['kind']=='guardian-unit')
    alter('omitted forcing branch','focused_branch',lambda r:forced(r)['branches'].pop())
    alter('duplicate forcing branch','focused_branch',lambda r:forced(r)['branches'].append(deepcopy(forced(r)['branches'][0])))
    alter('wrong split literal','focused_branch',lambda r:forced(r)['branches'][0].update(literal=[9,4]))
    alter('unsupported depth claim','focused_branch',lambda r:r['options'].update(forcing_depth=0))
    alter('wrong guardian unit','persistent_branch',lambda r:unit(r).update(digit=8))
    alter('dropped inherited clause','persistent_branch',lambda r:r['input_relations']['clauses'].clear())
    alter('mutated inherited guardians','persistent_branch',lambda r:r['input_relations']['clauses'][0]['guardians'].clear())
    alter('wrong clause identity','persistent_origin',lambda r:r['relations']['clauses'][0].update(id='0'*64))
    def forge(r):
        r['relations']['clauses'].clear()
        for e in r['trace']:
            if e['kind']=='catalogue-scan':e['new_clauses'].clear()
    alter('forged complete coverage with missing occurrence','persistent_origin',forge)
    alter('cutoff phase hidden as fixedpoint','persistent_origin',lambda r:r['phase_results'].append({'phase':'catalogue','status':'cutoff'}))
    alter('refutation relabeled fixedpoint','known20_refutation',lambda r:r.update(status='fixedpoint'))
    alter('wrong output mask','known20_refutation',lambda r:r['domains'].__setitem__(0,511))
    alter('nonmonotone inherited origin','persistent_branch',lambda r:r['input_relations']['origins'].__setitem__(next(iter(r['input_relations']['origins'])),[1]*81))
    alter('guardians lost only in returned state','persistent_origin',lambda r:r['relations']['clauses'].clear())
    for label,result in trials:
        try:verify_result(result)
        except(AssertionError,ValueError,KeyError,IndexError):continue
        raise AssertionError('Tamper accepted: '+label)
    return {'positive_controls':4,'known20_refuted':True,'persistent_clause_survives_palette_deletion':True,
            'focused_and_disabled_profiles_not_strong':True,'mutations_rejected':len(trials),'checks':checks}


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('result',type=Path,nargs='?');p.add_argument('--no-repeat-fixedpoint',action='store_true');p.add_argument('--self-test',action='store_true');a=p.parse_args()
    if a.self_test:print(json.dumps(self_test(),sort_keys=True))
    else:
        assert a.result is not None
        print(json.dumps(verify_result(json.loads(a.result.read_text()),not a.no_repeat_fixedpoint),sort_keys=True))


if __name__=='__main__':
    if not __debug__:raise SystemExit('Run without -O')
    main()
