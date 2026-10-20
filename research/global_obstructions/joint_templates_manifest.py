"""Export/check every bounded joint-template endpoint outcome, without novelty.

Counts are bookkeeping over explicit endpoints, not acceptance of every
refutation. The independent joint_audit checker accepts retained proof controls.
"""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
from joint_templates_search import canonical,digest
from real_verify import close,exact_completions,replay_trace
from real_four_digit_mine import initial_domains

HERE=Path(__file__).resolve().parent
RESULTS=HERE/'joint_templates_corpus_results.json'
TRIALS=HERE/'joint_templates_corpus_trials.jsonl'
RETAINED_RESULTS=HERE/'joint_templates_retained_results.json'
RETAINED_TRIALS=HERE/'joint_templates_retained_trials.jsonl'

def public_seeds(retained=False):
    public=HERE/'discovery_pipeline_inputs.json';cases={c['name']:c for c in json.loads(public.read_text())['cases']}
    if not retained:return cases
    from copy import deepcopy
    rows=json.loads((HERE/'discovery_pipeline_retained_results.json').read_text())['records'];answer={}
    for row in rows:
        if row['status']!='retained-configured-fixedpoint':continue
        c=deepcopy(cases[row['id']]);c['base_trace']+=row['additional_trace'];c['base_domains']=row['domains'];c['state']=row['last_probe']['state'];c['status']=row['status'];answer[c['name']]=c
    return answer


def export(raw,retained=False):
    results_path=RETAINED_RESULTS if retained else RESULTS;trials_path=RETAINED_TRIALS if retained else TRIALS
    result=json.loads(results_path.read_text());bulk=json.loads(raw.read_text())
    assert result['raw_summary_sha256']==hashlib.sha256(raw.read_bytes()).hexdigest()
    rows=bulk['records'];trials_path.write_bytes(b''.join(canonical(q)+b'\n'for q in rows))
    result['endpoint_manifest_sha256']=hashlib.sha256(trials_path.read_bytes()).hexdigest()
    result['endpoint_records']=len(rows)
    result['query_limits']={'support_order':3,'max_per_digit':2000,'max_pair_checks':20000000,'max_support_checks':20000000,'seconds':bulk['limits']['query_seconds']}
    public=HERE/'discovery_pipeline_inputs.json';seeds=public_seeds(retained)
    result['public_input_sha256']=hashlib.sha256(public.read_bytes()).hexdigest()
    for host in result['hosts']:
        host['public_seed_sha256']=digest(seeds[host['id']]);host.pop('base_domains',None)
    if retained:result['public_retained_sha256']=hashlib.sha256((HERE/'discovery_pipeline_retained_results.json').read_bytes()).hexdigest()
    results_path.write_text(json.dumps(result,indent=2)+'\n');return verify(retained=retained)


def verify(recount=False,replay=0,retained=False):
    results_path=RETAINED_RESULTS if retained else RESULTS;trials_path=RETAINED_TRIALS if retained else TRIALS
    data=trials_path.read_bytes();result=json.loads(results_path.read_text());rows=[json.loads(q)for q in data.splitlines()]
    assert hashlib.sha256(data).hexdigest()==result['endpoint_manifest_sha256']
    assert len(rows)==result['endpoint_records']==result['counts']['target_opposites']
    public=HERE/'discovery_pipeline_inputs.json';seeds=json.loads(public.read_text());assert seeds['schema_version']==1
    assert hashlib.sha256(public.read_bytes()).hexdigest()==result['public_input_sha256']
    seeds=public_seeds(retained)
    if retained:assert result['public_retained_sha256']==hashlib.sha256((HERE/'discovery_pipeline_retained_results.json').read_bytes()).hexdigest()
    hosts={h['id']:h for h in result['hosts']};seen=set();counts=Counter();hostcounts={h:Counter()for h in hosts};replayed=0
    assert set(hosts)==set(seeds)
    for h in hosts.values():
        c=seeds[h['id']];assert digest(c)==h['public_seed_sha256']
        assert c['puzzle']==h['puzzle'] and c['solution']==h['solution']
        assert digest(c['state'])==h['state_sha256']
        d=c['base_domains'];assert d==c['state']['domains']
        initial=initial_domains([int(v)if v in '123456789'else 0 for v in c['puzzle']])
        assert initial==c['initial_domains'] and replay_trace(initial,c['base_trace'],list(map(int,c['solution'])))==d
        assert len(d)==81 and all(type(m)is int and 0<m<=511 for m in d)
        assert all(m>>(int(v)-1)&1 for m,v in zip(d,h['solution']))
        if recount:assert exact_completions(h['puzzle'],2)==[h['solution']]
    for row in rows:
        host=hosts[row['host']];c,d,v=row['cell'],row['digit'],row['conclusion_value'];raw=seeds[row['host']]['base_domains'].copy()
        assert type(c)is int and 0<=c<81 and type(d)is int and 1<=d<=9 and type(v)is bool
        assert raw[c]&(raw[c]-1) and raw[c]>>(d-1)&1
        assert v==(int(host['solution'][c])==d)
        key=row['host'],c,d;assert key not in seen;seen.add(key)
        raw[c]=raw[c]&~(1<<(d-1)) if v else 1<<(d-1)
        assert digest(raw)==row['raw_opposite_domains_hash']
        counts['target_opposites']+=1;hostcounts[row['host']]['target_opposites']+=1
        kind='local_refuted'if row['status']=='local-refuted'else row['classification']
        assert kind in ('local_refuted','pair-refuted','triple-refuted','support-fixed-point','cutoff')
        if kind!='local_refuted':
            assert len(row['closed_opposite_domains_hash'])==64 and len(row['initial_template_counts'])==9
            assert row['pairwise_checkpoint']==(kind in ('triple-refuted','support-fixed-point')) or kind=='cutoff'
            assert (row['cutoff_reason']is not None)==(kind=='cutoff')
        counts[kind]+=1;hostcounts[row['host']][kind]+=1
        if replayed<replay:
            try:closed,trace=close(raw)
            except ValueError as error:
                assert kind=='local_refuted'and str(error).startswith(('Hall contradiction','Three-house contradiction'))
            else:
                assert kind!='local_refuted'and digest(closed)==row['closed_opposite_domains_hash']and len(trace)==row['closure_events']
            replayed+=1
    assert {k:result['counts'][k]for k in counts}==dict(counts)
    for name,h in hosts.items():
        assert dict(hostcounts[name])==h['counts']
        expected={(name,c,d)for c,m in enumerate(seeds[name]['base_domains'])if m&(m-1)for d in range(1,10)if m>>(d-1)&1}
        assert len(expected)==h['targets_available']
        if h['completed']:assert {k for k in seen if k[0]==name}==expected
    return {'records':len(rows),'hosts':len(hosts),'counts':dict(counts),'unique_recounts':len(hosts)if recount else 0,
     'opposite_closures_replayed':replayed,'scope':'Endpoint/count/binding audit only. Independent proof acceptance separate; fixedpoints/cutoffs not Sudoku solutions or novelty.'}


def main():
    if not __debug__:raise SystemExit('Run without -O')
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('mode',choices=('export','verify'));p.add_argument('--raw',type=Path,default=Path('target/joint-templates-corpus/summary.json'))
    p.add_argument('--retained',action='store_true');p.add_argument('--count',action='store_true');p.add_argument('--replay-closures',type=int,default=0);a=p.parse_args()
    print(json.dumps(export(a.raw,a.retained)if a.mode=='export'else verify(a.count,a.replay_closures,a.retained)))

if __name__=='__main__':main()
