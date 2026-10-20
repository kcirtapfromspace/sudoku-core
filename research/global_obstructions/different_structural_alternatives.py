"""Bounded SAT discovery after deliberately breaking a known tridagon geometry.

The old-stalled board and its candidate lists are frozen. At most two successive
known 12-cell geometries are hit by omitting one selected cell. A returned core
is discovery only: its native lists, proof tree, local/template relaxations and
conditional known-pattern classification still require independent acceptance.
"""
import os
for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):os.environ[name]='1'
import argparse
from collections import Counter
from functools import lru_cache
import json
from pathlib import Path
import random
import time
from different_sat_search import GraphOracle
from different_audit_known_patterns import tridagon


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('fixtures',type=Path,nargs='+')
    p.add_argument('--seconds',type=float,default=240);p.add_argument('--seed',type=int,default=221641)
    p.add_argument('--depth',type=int,default=2);p.add_argument('--output-dir',type=Path,default=Path('target/different-structural-alternatives'))
    a=p.parse_args();a.output_dir.mkdir(parents=True,exist_ok=True);rng=random.Random(a.seed);start=time.monotonic()
    counts=Counter();records=[];candidates=[]
    @lru_cache(None)
    def geometry(cells):
        fake=[511]*81
        for c in cells:fake[c]=7
        hits,_=tridagon(fake,max_matches=1)
        return tuple(hits[0]['cells']) if hits else ()
    for path in a.fixtures:
        if time.monotonic()-start>=a.seconds:break
        case=json.loads(path.read_text())
        if case.get('status')!='full-old-stall':continue
        domains=case['base_domains'];solution=list(map(int,case['solution']));oracle=GraphOracle(domains)
        targets=[(c,d,d==solution[c]) for c in oracle.cells for d in range(1,10) if domains[c]&(1<<(d-1))]
        rng.shuffle(targets);targets.sort(key=lambda t:(bin(domains[t[0]]).count('1'),not t[2]))
        try:
            for c,d,placement in targets:
                if time.monotonic()-start>=a.seconds:break
                literal=-(9*c+d) if placement else 9*c+d
                initial=oracle.core(oracle.cells,literal,rng,tries=1);counts['targets']+=1
                if initial is None:counts['initial-satisfiable']+=1;continue
                first=geometry(tuple(initial))
                if not first:
                    counts['initial-no-tridagon-geometry']+=1
                    continue # root's complementary direct core search owns these
                counts['initial-known-tridagon']+=1
                queue=[(frozenset(),initial,0)];seen=set()
                while queue and time.monotonic()-start<a.seconds:
                    omitted,core,depth=queue.pop()
                    hit=geometry(tuple(core))
                    if not hit:
                        counts['alternative-no-tridagon-geometry']+=1
                        record={'source':str(path),'cell':c,'digit':d,'placement':placement,
                          'omitted':sorted(omitted),'cells':core,'depth':depth,
                          'union_digits':sorted({v for q in core for v in range(1,10) if domains[q]&(1<<(v-1))}),
                          'discovery_only':'May still be a conditional known composition; no certificate accepted'}
                        records.append(record)
                        if len(core)<=32:
                            record['base_domains']=domains;record['solution']=case['solution'];record['puzzle']=case.get('puzzle')
                            name='candidate-%04d.json'%len(candidates);(a.output_dir/name).write_text(json.dumps(record,indent=2)+'\n');candidates.append(name)
                            if len(candidates)<=5 or len(candidates)%1000==0:
                                print('ALTERNATIVE',name,len(core),record['union_digits'],flush=True)
                        continue
                    if depth>=a.depth:counts['known-geometry-depth-cutoff']+=1;continue
                    choices=[q for q in hit if q!=c];rng.shuffle(choices)
                    for q in choices:
                        new=omitted|{q}
                        if new in seen:continue
                        seen.add(new);allowed=[q for q in oracle.cells if q not in new]
                        alternative=oracle.core(allowed,literal,rng,tries=1);counts['omission-queries']+=1
                        if alternative is None:counts['omission-satisfiable']+=1;continue
                        counts['omission-unsatisfiable']+=1
                        queue.append((new,alternative,depth+1))
                counts['targets-queue-budget-unfinished']+=bool(queue)
        finally:oracle.delete()
        counts['fixtures']+=1;print('CASE',path.name,dict(counts),flush=True)
    result={'schema_version':1,'seed':a.seed,'elapsed_seconds':time.monotonic()-start,'depth_limit':a.depth,
      'counts':dict(counts),'discovery_files':candidates,'records':records,
      'scope':'Geometry-hitting SAT discovery only; no historical novelty, full-relaxation or proof-certificate acceptance'}
    (a.output_dir/'summary.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k not in ('records','discovery_files')}),flush=True)

if __name__=='__main__':main()
