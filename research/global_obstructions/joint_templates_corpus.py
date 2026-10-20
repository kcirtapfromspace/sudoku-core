"""Bounded joint-template discovery on frozen unique full-old-engine stalls.

These source states predate the strengthened baseline; discoveries here are
controls, never promoted as surviving that newer profile or as original rules.
"""
import argparse
from collections import Counter
import copy
import hashlib
import json
from pathlib import Path
import time

from joint_templates_search import canonical,digest,search
from real_verify import close,exact_completions,replay_trace
from real_four_digit_mine import initial_domains

ROOT=Path(__file__).resolve().parent.parent.parent


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('folders',type=Path,nargs='*');p.add_argument('--inputs',type=Path,help='Public JSON cases list or object with cases list')
    p.add_argument('--seconds',type=float,default=180);p.add_argument('--query-seconds',type=float,default=3)
    p.add_argument('--accepted-status',action='append',default=[])
    p.add_argument('--profile-label',default='frozen-legacy-old-snapshot')
    p.add_argument('--max-targets-per-host',type=int,default=1000);p.add_argument('--output-dir',type=Path,default=ROOT/'target/joint-templates-corpus')
    a=p.parse_args();a.output_dir.mkdir(parents=True,exist_ok=True);started=time.monotonic();counts=Counter();records=[];hosts=[];files=[]
    def checkpoint():
        result={'schema_version':1,'counts':dict(counts),'seconds':time.monotonic()-started,
         'limits':{'seconds':a.seconds,'query_seconds':a.query_seconds,'max_targets_per_host':a.max_targets_per_host},
         'profile_label':a.profile_label,'accepted_statuses':a.accepted_status or ['full-old-stall'],
         'hosts':hosts,'records':records,'saved_controls':files,
         'scope':'Frozen full-old stalls only; stronger baseline not yet applied. Target opposites exact local closed before template search. Budget cutoffs are inconclusive; proofs require independent acceptance.'}
        (a.output_dir/'summary.json').write_text(json.dumps(result,indent=2)+'\n');return result
    items=[(str(p),json.loads(p.read_text()),hashlib.sha256(p.read_bytes()).hexdigest())for p in sorted({p for f in a.folders for p in f.glob('modern*.json')})]
    if a.inputs:
        payload=json.loads(a.inputs.read_text());cases=payload.get('cases')if isinstance(payload,dict)else payload
        if not isinstance(cases,list):raise ValueError('--inputs expects a cases list')
        items.extend((str(a.inputs)+'#'+c['name'],c,digest(c))for c in cases)
    if not items:p.error('Provide source folders or --inputs')
    for source_name,case,source_hash in items:
        if time.monotonic()-started>=a.seconds:counts['campaign_cutoff_before_host']+=1;break
        if case.get('status')not in (a.accepted_status or ['full-old-stall']):continue
        assert exact_completions(case['puzzle'],2)==[case['solution']]
        initial=initial_domains([int(v)if v in '123456789'else 0 for v in case['puzzle']])
        assert initial==case['initial_domains']
        domains=case['base_domains'];assert replay_trace(initial,case['base_trace'],list(map(int,case['solution'])))==domains
        assert domains==case['state']['domains'] and close(domains)[1]==[]
        assert all(m>>(int(v)-1)&1 for m,v in zip(domains,case['solution']))
        hostcounts=Counter();endpoints=[(c,d,d==int(case['solution'][c])) for c,m in enumerate(domains) if m&(m-1) for d in range(1,10) if m>>(d-1)&1]
        endpoints.sort(key=lambda v:(bin(domains[v[0]]).count('1'),not v[2],v[0],v[1]));saved=set()
        host={'id':case['name'],'source':source_name,'source_sha256':source_hash,'targets_available':len(endpoints),'completed':False,'counts':hostcounts};hosts.append(host)
        for c,d,placement in endpoints[:a.max_targets_per_host]:
            if time.monotonic()-started>=a.seconds:counts['campaign_cutoff_during_host']+=1;break
            counts['target_opposites']+=1;hostcounts['target_opposites']+=1
            raw=domains.copy();bit=1<<(d-1);raw[c]=raw[c]&~bit if placement else bit
            record={'host':case['name'],'cell':c,'digit':d,'conclusion_value':placement,
              'raw_opposite_domains_hash':digest(raw)}
            try:closed,trace=close(raw)
            except ValueError as error:
                if not str(error).startswith(('Hall contradiction','Three-house contradiction')):raise
                record['status']='local-refuted';counts['local_refuted']+=1;hostcounts['local_refuted']+=1;records.append(record);continue
            query=copy.deepcopy(case);query['target']={'cell':c,'digit':d,'value':not placement}
            query['conclusion']={'cell':c,'digit':d,'value':placement,'kind':'placement' if placement else 'elimination'}
            query['assumed_domains']=closed;query['assumed_state']={'version':1,'domains':closed,'values':[m.bit_length() if not m&(m-1) else 0 for m in closed]};query['assumption_trace']=trace
            result=search(query,seconds=min(a.query_seconds,max(0.001,a.seconds-(time.monotonic()-started))))
            record.update(status=result['status'],closed_opposite_domains_hash=digest(closed),closure_events=len(trace),
              initial_template_counts=result['proof']['initial']['counts'],pairwise_checkpoint=result['proof']['pairwise_fixed_point'] is not None,
              event_count=len(result['proof']['events']),cutoff_reason=result['cutoff_reason'])
            kind='triple-refuted' if result['status']=='refuted' and record['pairwise_checkpoint'] else 'pair-refuted' if result['status']=='refuted' else result['status']
            record['classification']=kind;counts[kind]+=1;hostcounts[kind]+=1
            if kind in ('pair-refuted','triple-refuted') and kind not in saved:
                stem=case['name']+'-'+kind;fixture=a.output_dir/(stem+'-fixture.json');proof=a.output_dir/(stem+'-proof.json')
                query['joint_templates_acceptance_role']='Known-method discovery control from frozen older baseline; strengthened profile not yet applied.'
                fixture.write_text(json.dumps(query,indent=2)+'\n');proof.write_text(json.dumps(result,indent=2)+'\n');saved.add(kind)
                files.append({'fixture':str(fixture),'proof':str(proof),'kind':kind,'host':case['name']});record['saved_proof']=str(proof)
            records.append(record)
        else:host['completed']=len(endpoints)<=a.max_targets_per_host
        counts['hosts_visited']+=1;counts['hosts_completed']+=host['completed'];checkpoint();print('HOST',case['name'],dict(hostcounts),flush=True)
    result=checkpoint();print(json.dumps({k:v for k,v in result.items() if k not in ('records','hosts','saved_controls')}),flush=True)

if __name__=='__main__':main()
