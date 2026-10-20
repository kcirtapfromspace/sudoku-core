#!/usr/bin/env python3
"""Independently accept retained bulk controls; full corpus replay is not claimed.

Bulk fixtures/proofs deliberately remain under target. Stable CI controls are
covered by joint_audit_verify.py --self-test without those local bulk files.
"""
import argparse
import hashlib
import json
from pathlib import Path
import time
import joint_audit_verify as exact

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent.parent


def run(traces=False):
    manifest_path=HERE/'joint_templates_corpus_results.json';manifest=json.loads(manifest_path.read_text())
    records=[];hosts={};started=time.monotonic()
    if traces:
        from real_verify import exact_completions,replay_trace,HOUSES
    for row in manifest['retained_bulk_controls']:
        fixture=ROOT/row['fixture'];proof_path=ROOT/row['proof']
        exact.need(hashlib.sha256(fixture.read_bytes()).hexdigest()==row['fixture_sha256'],'bulk fixture hash')
        exact.need(hashlib.sha256(proof_path.read_bytes()).hexdigest()==row['proof_sha256'],'bulk proof hash')
        case=json.loads(fixture.read_text());proof=json.loads(proof_path.read_text())
        answer=exact.verify(case,proof);exact.need(answer['separation']==(row['kind']=='triple-refuted'),'bulk classification')
        if traces:
            if row['host']not in hosts:
                exact.need(exact_completions(case['puzzle'],2)==[case['solution']],'unique givens completion')
                values=[int(x)if x in '123456789'else 0 for x in case['puzzle']]
                initial=[]
                for c,d in enumerate(values):
                    excluded={values[q]for h in HOUSES if c in h for q in h if values[q]}
                    initial.append(1<<(d-1)if d else sum(1<<(k-1)for k in range(1,10)if k not in excluded))
                exact.need(initial==case['initial_domains'],'native initial givens masks')
                base=replay_trace(initial,case['base_trace'],list(map(int,case['solution'])))
                exact.need(base==case['base_domains'],'exact base trace')
                hosts[row['host']]={'puzzle':case['puzzle'],'solution':case['solution'],'base_domains':base}
            else:
                old=hosts[row['host']]
                exact.need(all(case[k]==old[k]for k in old),'shared host binding')
            branch=case['base_domains'].copy();target=case['target'];bit=1<<(target['digit']-1)
            branch[target['cell']]&=bit if target['value']else 511^bit
            exact.need(replay_trace(branch,case['assumption_trace'])==case['assumed_domains'],'exact assumption closure trace')
            answer['givens_and_closure_replayed']=True
        records.append({**row,'verification':answer})
    exact.need(len(records)==16 and len({r['host']for r in records})==8,'all retained hosts/controls')
    return {'schema_version':1,'manifest_sha256':hashlib.sha256(manifest_path.read_bytes()).hexdigest(),
            'accepted_controls':len(records),'pair_controls':sum(not r['verification']['separation']for r in records),
            'triple_separation_controls':sum(r['verification']['separation']for r in records),
            'unique_hosts_checked':len(hosts),'givens_and_closure_replayed':traces,
            'records':records,'seconds':time.monotonic()-started,
            'scope':'Independent acceptance of retained16 exact controls only. Does not independently replay all1200 corpus queries or stronger baseline completion.'}


def main():
    exact.need(__debug__,'Run without -O');p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--traces',action='store_true');p.add_argument('--output',type=Path);args=p.parse_args();result=run(args.traces)
    if args.output:args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items()if k!='records'},sort_keys=True))


if __name__=='__main__':main()
