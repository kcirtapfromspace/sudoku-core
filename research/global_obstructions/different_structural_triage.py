"""Normalize a bounded diverse sample of geometry-hitting SAT discoveries."""
import os
for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):os.environ[name]='1'
import argparse
from collections import Counter
import json
from pathlib import Path
from real_verify import close


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--summary',type=Path,default=Path('target/different-structural-alternatives/summary.json'))
    p.add_argument('--limit',type=int,default=20);p.add_argument('--output-dir',type=Path,default=Path('target/different-structural-screened'))
    a=p.parse_args();a.output_dir.mkdir(parents=True,exist_ok=True);data=json.loads(a.summary.read_text());counts=Counter();records=[];cache={};seen=set()
    options=sorted(data['records'],key=lambda x:(len(x['cells']),len(x['union_digits']),x['source'],x['cell'],x['digit']))
    # Prefer different target endpoints; redundant proof cores can be screened later.
    for raw in options:
        if len(records)>=a.limit:break
        if len(raw['cells'])>32:continue
        key=(raw['source'],raw['cell'],raw['digit'],raw['placement'])
        if key in seen:continue
        seen.add(key);counts['unique_target_opposites_checked']+=1
        case=cache.setdefault(raw['source'],json.loads(Path(raw['source']).read_text()))
        domains=case['base_domains'].copy();bit=1<<(raw['digit']-1)
        domains[raw['cell']]=domains[raw['cell']]&~bit if raw['placement'] else bit
        try:assumed,trace=close(domains)
        except ValueError as error:
            counts['full_opposite_local_contradiction']+=1;continue
        counts['full_opposite_local_survivor']+=1
        state={'version':1,'domains':assumed,'values':[m.bit_length() if m&(m-1)==0 else 0 for m in assumed]}
        result={**case,'target':{'cell':raw['cell'],'digit':raw['digit'],'value':not raw['placement']},
          'conclusion':{'cell':raw['cell'],'digit':raw['digit'],'value':raw['placement'],'kind':'placement' if raw['placement'] else 'elimination'},
          'graph_cells':raw['cells'],'assumed_domains':assumed,'assumed_state':state,'assumption_trace':trace,
          'discovery_only':'No literal Loki geometry in selected cells; conditional known patterns/template feasibility/Rust certificate still unaccepted',
          'geometry_hitting_record':raw}
        name='candidate-%02d.json'%len(records);(a.output_dir/name).write_text(json.dumps(result,indent=2)+'\n');records.append({'file':str(a.output_dir/name),'source':raw['source'],'cells':len(raw['cells']),'target':result['target']})
        print('SCREENED',str(a.output_dir/name),len(raw['cells']),flush=True)
    result={'schema_version':1,'counts':dict(counts),'records':records,'scope':'Bounded smallest-core diverse-target sample; local closure passed only. No discovery acceptance.'}
    (a.output_dir/'summary.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result),flush=True)

if __name__=='__main__':main()
