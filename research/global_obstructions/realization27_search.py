"""Bounded compatible-completion/givens realization search for the exact27 graph.

All27 support cells are excluded from givens. Exact full initial palette masks
and relaxed support lists are separate tiers. A maximal outside clue set is a
necessary uniqueness test: if it has a second completion, every subset fails.
No node/time cutoff is ever classified as uniqueness or an engine stall.
"""
import os
for _key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
             "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_key] = "1"
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import random
import subprocess
import sys
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
from real_verify import HOUSES, close, exact_completions
from real_four_digit_mine import initial_domains
from nonrank_puzzles_search import cached_engine, candidate_state, full_close


class Enumerator:
    def __init__(self):
        source = HERE / "realization27_search_solutions.cpp"
        binary = ROOT / "target/realization27-solutions"
        if not binary.exists() or binary.stat().st_mtime < source.stat().st_mtime:
            subprocess.run(["c++","-O3","-std=c++17",source,"-o",binary],check=True)
        self.process = subprocess.Popen([binary],stdin=subprocess.PIPE,
                                        stdout=subprocess.PIPE,text=True,bufsize=1)

    def run(self, domains=None, puzzle=None, seed=0, limit=2, nodes=1000000):
        text = ("0 " if puzzle is not None else "1 ")+str(seed)+" "+str(limit)+" "+str(nodes)+" "
        text += puzzle if puzzle is not None else " ".join(map(str,domains))
        self.process.stdin.write(text+"\n");self.process.stdin.flush()
        return json.loads(self.process.stdout.readline())

    def close(self):
        self.process.stdin.close();self.process.wait(timeout=10)


def valid_solution(solution,puzzle=None,domains=None):
    assert len(solution)==81 and set(solution)<=set("123456789")
    assert all({solution[c] for c in h}==set("123456789") for h in HOUSES)
    if puzzle is not None:
        assert all(d=="." or d==solution[c] for c,d in enumerate(puzzle))
    if domains is not None:
        assert all(m&(1<<(int(d)-1)) for m,d in zip(domains,solution))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed",type=int,default=61481)
    parser.add_argument("--samples",type=int,default=2000)
    parser.add_argument("--seconds",type=float,default=360)
    parser.add_argument("--max-positive-tests",type=int,default=4)
    parser.add_argument("--output-dir",type=Path,default=ROOT/"target/realization27-completions")
    args=parser.parse_args();args.output_dir.mkdir(parents=True,exist_ok=True)
    source_path=HERE/"nonrank_graph_native27_candidate.json"
    source=json.loads(source_path.read_text());core=source["graph_cells"]
    guardian=source["target"]["cell"];extra=source["target"]["digit"]
    palette=sum(1<<(d-1) for d in source["palette"])
    domains=[511]*81
    for c in core:domains[c]=palette
    domains[guardian]=1<<(extra-1)
    enum=Enumerator();counts=Counter();records=[];seen=set();rng=random.Random(args.seed)
    started=time.monotonic();trial_file=(args.output_dir/"trials.jsonl").open("w")
    for trial in range(args.samples+1):
        if time.monotonic()-started>args.seconds:break
        seed=0 if trial==0 else args.seed+trial-1
        if trial==0:
            solution=source["completion"];generation_nodes=0
        else:
            generated=enum.run(domains=domains,seed=seed,limit=1)
            counts["completion-generation-attempts"]+=1
            if generated["cutoff"]:
                counts["completion-generation-cutoffs"]+=1;continue
            assert len(generated["solutions"])==1
            solution=generated["solutions"][0];generation_nodes=generated["nodes"]
        valid_solution(solution,domains=domains)
        if solution in seen:
            counts["duplicate-completions"]+=1;continue
        seen.add(solution);counts["distinct-compatible-completions"]+=1
        puzzle="".join("." if c in core else solution[c] for c in range(81))
        tested=enum.run(puzzle=puzzle,limit=2)
        for answer in tested["solutions"]:valid_solution(answer,puzzle=puzzle)
        status="maximal-outside-nonunique"
        if tested["cutoff"] and len(tested["solutions"])<2:
            status="uniqueness-cutoff"
        elif len(tested["solutions"])==1:
            assert tested["solutions"]==[solution]
            status="maximal-outside-unique"
        else:assert len(tested["solutions"])==2
        counts[status]+=1
        record={"seed":seed,"completion":solution,"status":status,
                "generation_nodes":generation_nodes,"count_nodes":tested["nodes"],
                "count_cutoff":tested["cutoff"],"solutions_capped_at_2":len(tested["solutions"])}
        if trial==0 or status=="maximal-outside-unique":
            control={**record,"puzzle":puzzle,"independent_completions":tested["solutions"],
                     "graph_cells":core,"palette":source["palette"],"guardian":guardian,
                     "all_support_unclued":True,"source":"generated compatible completion"}
            (args.output_dir/("stored-control.json" if trial==0 else "unique-"+str(seed)+".json")).write_text(json.dumps(control,indent=2)+"\n")
        if status=="maximal-outside-unique" and counts["unique-endpoints-compared"]<args.max_positive_tests:
            print("UNIQUE_MAXIMAL",seed,puzzle,flush=True)
            # Random clue deletion can change earlier decisions; preserve the
            # same27 graph's list restrictions after exact local closure.
            givens=[c for c in range(81) if c not in core];order=givens.copy();rng.shuffle(order)
            for c in order:
                proposed="".join(solution[p] if p in givens and p!=c else "." for p in range(81))
                count=enum.run(puzzle=proposed,limit=2)
                counts["clue-deletion-trials"]+=1
                if count["cutoff"] or count["solutions"]!=[solution]:continue
                initial=initial_domains([int(d) if d!="." else 0 for d in proposed]);base,_=close(initial)
                if any(base[p]&~(palette|(1<<(extra-1) if p==guardian else 0)) for p in core):continue
                givens.remove(c);counts["accepted-clue-deletions"]+=1
            minimized="".join(solution[c] if c in givens else "." for c in range(81))
            assert exact_completions(minimized)==[solution]
            case={"id":"realization27-"+str(seed),"puzzle":minimized}
            result,initial,base,trace,findings,calls,comparison=full_close(
                cached_engine(),case,list(map(int,solution)),45,1000)
            out={"schema_version":1,"name":case["id"],"puzzle":minimized,"solution":solution,
                 "completion":solution,"solutions_capped_at_2":1,"initial_domains":initial,
                 "base_domains":base,"base_trace":trace,"state":candidate_state(base),
                 "graph_cells":core,"palette":source["palette"],"target":source["target"],
                 "status":result,"full_old_fixedpoint":comparison,"old_findings":findings,
                 "all_support_unclued":True,"old_calls":calls,"seed":seed}
            counts["unique-endpoints-compared"]+=1;counts["old-"+result]+=1
            (args.output_dir/("endpoint-"+str(seed)+".json")).write_text(json.dumps(out,indent=2)+"\n")
            print("ENDPOINT",seed,result,flush=True)
        trial_file.write(json.dumps(record,separators=(",",":"))+"\n");trial_file.flush();records.append(record)
        if trial%100==0:print("PROGRESS",trial,dict(counts),"seconds",round(time.monotonic()-started,2),flush=True)
    enum.close();trial_file.close()
    summary={"schema_version":1,"seed_start":args.seed,"requested_samples":args.samples,
             "seconds":time.monotonic()-started,"counts":dict(counts),"command":" ".join(sys.argv),
             "source_fixture_sha256":hashlib.sha256(source_path.read_bytes()).hexdigest(),
             "core":core,"palette":source["palette"],"guardian":guardian,"extra":extra,
             "scope":"All27 support cells excluded from givens; randomized domain-compatible completions, maximal54outside uniqueness test, then uniqueness-preserving clue deletion and full-old/local closure only for unique endpoints.",
             "limitation":"Finite randomized completion sample; nonunique maximum implies impossibility only for that completion, not for every compatible completion."}
    (args.output_dir/"summary.json").write_text(json.dumps(summary,indent=2)+"\n")
    print(json.dumps(summary),flush=True)


if __name__=="__main__":
    if not __debug__:raise SystemExit("Run without -O")
    main()
