"""Three-minute changed-premise search with one selected non-guardian given.

This intentionally abandons the all27-unclued premise. Clue removal preserves
independent uniqueness and the same selected graph's palette-list restriction.
Local, nonforcing and full-older outcomes remain separate from live stalls.
"""
import argparse
import json
from pathlib import Path
import random
import subprocess
import time

from realization27_search import Enumerator
from real_four_digit_mine import initial_domains
from real_verify import close,exact_completions,replay_trace
from nonrank_puzzles_search import cached_engine,candidate_state,full_close

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent.parent


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed",type=int,default=92137)
    parser.add_argument("--seconds",type=float,default=150)
    parser.add_argument("--restarts",type=int,default=40)
    parser.add_argument("--output-dir",type=Path,default=ROOT/"target/realization27-anchored")
    args=parser.parse_args();args.output_dir.mkdir(parents=True,exist_ok=True)
    original=json.loads((HERE/"nonrank_graph_native27_candidate.json").read_text())
    controls=json.loads((HERE/"realization27_results_internal.json").read_text())["controls"]
    core=original["graph_cells"];g=original["target"]["cell"];extra=original["target"]["digit"]
    palette=sum(1<<(d-1) for d in original["palette"]);solution=original["completion"]
    enum=Enumerator();started=time.monotonic();records=[];live=[];seen=set()
    totals={"clue_removals_tested":0,"accepted_clue_removals":0,"unique_count_cutoffs":0}
    for restart in range(args.restarts):
        if time.monotonic()-started>=args.seconds:break
        control=controls[restart%len(controls)];anchor=control["internal_given_cell"]
        givens=[c for c in range(81) if c not in core or c==anchor]
        order=[c for c in givens if c!=anchor];random.Random(args.seed+restart).shuffle(order)
        for cell in order:
            if time.monotonic()-started>=args.seconds:break
            proposal="".join(solution[c] if c in givens and c!=cell else "." for c in range(81))
            totals["clue_removals_tested"]+=1
            counted=enum.run(puzzle=proposal,limit=2)
            if counted["cutoff"]:
                totals["unique_count_cutoffs"]+=1;continue
            if counted["solutions"]!=[solution]:continue
            base,_=close(initial_domains([int(d) if d!="." else 0 for d in proposal]))
            if any(base[c]&~(palette|(1<<(extra-1) if c==g else 0)) for c in core):continue
            givens.remove(cell);totals["accepted_clue_removals"]+=1
        puzzle="".join(solution[c] if c in givens else "." for c in range(81))
        initial=initial_domains([int(d) if d!="." else 0 for d in puzzle]);base,trace=close(initial)
        target_live=base[g]&(base[g]-1)!=0
        record={"restart":restart,"seed":args.seed+restart,"anchor":anchor,"puzzle":puzzle,
                "givens":len(givens),"target_live_after_local_closure":target_live,
                "unplaced_support_after_local_closure":sum(base[c]&(base[c]-1)!=0 for c in core)}
        records.append(record)
        if target_live and puzzle not in seen:
            assert exact_completions(puzzle)==[solution]
            seen.add(puzzle)
            output={"schema_version":1,"name":"anchored-native27-"+str(restart),
                    "puzzle":puzzle,"solution":solution,"completion":solution,
                    "solutions_capped_at_2":1,"initial_domains":initial,"base_domains":base,
                    "base_trace":trace,"state":candidate_state(base),"graph_cells":core,
                    "palette":original["palette"],"target":original["target"],
                    "anchor":anchor,"changed_premise":"One support vertex is a given; original all27-live uniform pattern is not claimed."}
            branch=base.copy();branch[g]&=~(1<<(extra-1))
            try:
                assumed,at=close(branch)
                output.update(assumed_domains=assumed,assumption_trace=at,
                              assumed_state=candidate_state(assumed),opposite_local_status="survives")
            except ValueError:
                output["opposite_local_status"]="contradiction"
            path=args.output_dir/("live-"+str(restart)+".json")
            path.write_text(json.dumps(output,indent=2)+"\n")
            live.append(path)
            print("LIVE_UNIQUE",path,"givens",len(givens),"branch",output["opposite_local_status"],flush=True)
        print("RESTART",restart,"anchor",anchor,"givens",len(givens),"targetlive",target_live,
              "seconds",round(time.monotonic()-started,2),flush=True)
    enum.close()
    comparisons=[]
    for path in live[:4]:
        case=json.loads(path.read_text())
        try:
            status,initial,base,trace,findings,calls,comparison=full_close(
                cached_engine(),case,list(map(int,solution)),20,1000)
            assert replay_trace(initial,trace,list(map(int,solution)))==base
            case.update(old_status=status,old_trace=trace,old_domains=base,old_state=candidate_state(base),
                        old_findings=findings,old_calls=calls,full_old_comparison=comparison)
            comparisons.append({"file":str(path),"status":status,"target_value":case["old_state"]["values"][g],
                                "old_findings":len(findings),"old_calls":calls})
            path.write_text(json.dumps(case,indent=2)+"\n")
        except subprocess.TimeoutExpired:
            comparisons.append({"file":str(path),"status":"probe-timeout-not-stall"})
    summary={"schema_version":1,"scope":"Changed-premise tier: one selected support vertex is an actual given.",
             "seed_start":args.seed,"seconds":time.monotonic()-started,"records":records,"counts":totals,
             "distinct_unique_live_states":len(live),"full_older_comparisons":comparisons,
             "all27_unclued_request_fulfilled":False}
    (HERE/"realization27_results_anchored.json").write_text(json.dumps(summary,indent=2)+"\n")
    print(json.dumps({"restarts":len(records),"live_unique":len(live),"counts":totals,"comparisons":comparisons}),flush=True)


if __name__=="__main__":
    if not __debug__:raise SystemExit("Run without -O")
    main()
