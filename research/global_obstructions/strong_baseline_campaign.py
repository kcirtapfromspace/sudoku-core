"""One-worker bounded stronger-baseline controls on the eight unique hosts."""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import time
from strong_baseline_core import close_state,Options

HERE=Path(__file__).resolve().parent


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--seconds',type=float,default=60)
    p.add_argument('--catalogue-seconds',type=float,default=10)
    p.add_argument('--catalogue-nodes',type=int,default=100000)
    p.add_argument('--only',nargs='*',default=[])
    p.add_argument('--summary',type=Path,default=HERE/'strong_baseline_campaign_results.json')
    p.add_argument('--output',type=Path,default=Path('target/strong-baseline-corpus'))
    p.add_argument('--retained',type=Path,default=Path('target/discovery-retained-comparison.json'));a=p.parse_args()
    a.output.mkdir(parents=True,exist_ok=True)
    inputs=HERE/'discovery_pipeline_inputs.json';original=json.loads(inputs.read_text())['cases']
    comparison=json.loads(a.retained.read_text())if a.retained.exists()else None
    rows=comparison['records']if comparison else[{'id':c['name'],'fixture':c,'status':'original-old-stall'}for c in original]
    if a.only:rows=[r for r in rows if r['id']in a.only]
    rows.sort(key=lambda r:(r['id']!='modern-mith-minimals-000515',r['id']))
    records=[];began=time.monotonic()
    for row in rows:
        case=row['fixture'];domains=case['state']['domains']
        print('START',row['id'],row['status'],flush=True)
        options=Options(seconds=a.seconds,context_id=row['id'],history_snapshots=(case['initial_domains'],),
                        catalogue_seconds=a.catalogue_seconds,catalogue_nodes=a.catalogue_nodes)
        result=close_state(domains,options=options)
        solution=list(map(int,case['solution']))
        assert all(m&(1<<(d-1))for m,d in zip(result['domains'],solution))
        assert result['status']!='refuted','A genuinely unique base state cannot be refuted'
        result['source_id']=row['id'];result['puzzle']=case['puzzle'];result['solution']=case['solution']
        dest=a.output/(row['id']+'.json');dest.write_text(json.dumps(result,separators=(',',':'))+'\n')
        records.append({'id':row['id'],'seed_status':row['status'],'status':result['status'],'error':result['error'],
                        'unplaced':sum(bool(m&(m-1))for m in result['domains']),
                        'input_state_hash':__import__('strong_baseline_relations').digest(domains),
                        'output_state_hash':__import__('strong_baseline_relations').digest(result['domains']),
                        'trace_events':len(result['trace']),'relations':len(result['relations']['clauses']),
                        'seconds':result['seconds'],'result_path':str(dest),
                        'result_sha256':hashlib.sha256(dest.read_bytes()).hexdigest()})
        summary={'schema_version':1,'inputs_sha256':hashlib.sha256(inputs.read_bytes()).hexdigest(),
                 'profile':{'seconds':a.seconds,'catalogue_seconds':a.catalogue_seconds,'catalogue_nodes':a.catalogue_nodes},
                 'retained_comparison_sha256':hashlib.sha256(a.retained.read_bytes()).hexdigest()if comparison else None,
                 'counts':dict(Counter(r['status']for r in records)),'records':records,'seconds':time.monotonic()-began,
                 'scope':'Configured offline profile with supplied original givens-domain history. Cutoffs are inconclusive, never strengthened fixedpoints. Original unique hosts and corrected old-engine seeds are separate input provenance.'}
        (a.output/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
        print('DONE',row['id'],result['status'],round(result['seconds'],3),flush=True)
    a.summary.write_text(json.dumps(summary,indent=2)+'\n')


if __name__=='__main__':main()
