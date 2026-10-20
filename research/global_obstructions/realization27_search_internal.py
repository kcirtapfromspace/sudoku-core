"""Changed-premise controls: one non-guardian core clue plus all54outside clues.

This does NOT realize the requested all27-unclued premise. It measures the
boundary of the universal outside-only trade obstruction and keeps solved
states separate from genuine full-older-engine stalls.
"""
import json
from pathlib import Path

from realization27_search import Enumerator,valid_solution
from nonrank_puzzles_search import cached_engine,candidate_state,full_close,old_probe
from real_verify import exact_completions

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent.parent


def main():
    source=json.loads((HERE/"nonrank_graph_native27_candidate.json").read_text())
    support=source["graph_cells"];guardian=source["target"]["cell"];solution=source["completion"]
    enum=Enumerator();records=[];controls=[]
    binary=ROOT/"target/nonrank-stdlib-baseline/instrumented/target/release/examples/real_existing_probe"
    if not binary.exists():binary=cached_engine()
    for cell in support:
        if cell==guardian:continue
        puzzle="".join("." if c in support and c!=cell else solution[c] for c in range(81))
        counted=enum.run(puzzle=puzzle,limit=2)
        assert not counted["cutoff"]
        for completed in counted["solutions"]:valid_solution(completed,puzzle=puzzle)
        record={"internal_given_cell":cell,"digit":int(solution[cell]),
                "solutions_capped_at_2":len(counted["solutions"]),"cutoff":False}
        if len(counted["solutions"])==1:
            assert exact_completions(puzzle)==[solution]
            case={"id":"changed-premise-core-given-"+str(cell),"puzzle":puzzle}
            status,initial,domains,trace,findings,calls,comparison=full_close(
                binary,case,list(map(int,solution)),45,1000)
            # If local closure has already solved it, a None dispatcher result
            # is recorded as solved, never as a live stall.
            if status=="solved":comparison=old_probe(binary,puzzle,domains,45)
            record.update(old_status=status,old_findings=len(findings),old_calls=calls)
            controls.append({"puzzle":puzzle,"solution":solution,"solutions_capped_at_2":1,
                "internal_given_cell":cell,"digit":int(solution[cell]),"givens":55,
                "initial_domains":initial,"base_domains":domains,"base_trace":trace,
                "state":candidate_state(domains),"old_findings":findings,"old_calls":calls,
                "full_old_comparison":comparison,"status":status,
                "changed_premise":"One selected support cell is a given; only26support cells remain unclued.",
                "strict_opposite_filter":"Guardian is placed by local closure; no LP-surviving opposite branch or global-method gap is claimed."})
        records.append(record)
    enum.close()
    result={"schema_version":1,"scope":"Separate changed-premise diagnostic; not the requested all27-unclued realization.",
            "support":support,"guardian":guardian,"records":records,"unique_cases":len(controls),
            "controls":controls,"node_cutoffs":0,
            "full_old_stalls":sum(c["status"]=="full-old-stall" for c in controls)}
    (HERE/"realization27_results_internal.json").write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps({"tests":len(records),"unique_cases":len(controls),
                      "statuses":[c["status"] for c in controls],"full_old_stalls":result["full_old_stalls"]}))


if __name__=="__main__":
    if not __debug__:raise SystemExit("Run without -O")
    main()
