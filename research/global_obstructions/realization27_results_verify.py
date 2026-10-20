"""Stdlib verification of realization27 records and explicit nonunique controls.

The full manifest's outcome labels remain discovery records. The two rejected
controls have explicit respecting completions; the stronger all-completion
trade theorem is checked by the separate exhaustive structural certificate.
"""
import argparse
import hashlib
import itertools
import json
from pathlib import Path

from real_four_digit_mine import initial_domains
from real_verify import HOUSES, exact_completions, replay_trace, sees

HERE=Path(__file__).resolve().parent


def completion(solution,puzzle=None):
    assert len(solution)==81 and set(solution)<=set("123456789")
    assert all({solution[c] for c in house}==set("123456789") for house in HOUSES)
    if puzzle is not None:
        assert all(d=="." or d==solution[c] for c,d in enumerate(puzzle))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--recount-sample",type=int,default=0)
    args=parser.parse_args()
    summary=json.loads((HERE/"realization27_results.json").read_text())
    metadata=summary["public_manifest"];data=(HERE/metadata["file"]).read_bytes()
    assert len(data)==metadata["bytes"] and hashlib.sha256(data).hexdigest()==metadata["sha256"]
    rows=[json.loads(line) for line in data.splitlines()]
    assert len(rows)==metadata["records"]==summary["counts"]["distinct-compatible-completions"]==4999
    assert len({r[0] for r in rows})==len({r[1] for r in rows})==4999
    core=summary["core"];palette=summary["palette"];guardian=summary["guardian"];extra=summary["extra"]
    for seed,solution,count in rows:
        assert count==2 and (seed==0 or summary["seed_start"]<=seed<summary["seed_start"]+5000)
        completion(solution)
        assert int(solution[guardian])==extra
        assert all(int(solution[c]) in palette for c in core if c!=guardian)
    controls_data=(HERE/summary["exact_controls"]["file"]).read_bytes()
    assert hashlib.sha256(controls_data).hexdigest()==summary["exact_controls"]["sha256"]
    controls=json.loads(controls_data)
    assert core==controls["core"] and palette==controls["palette"]
    for cell in range(81):
        if cell not in core:
            assert any(sees(cell,peer) for peer in core)
    exact=controls["exact_initial_mask_realization"]
    assert sum(d!="." for d in exact["puzzle"])==exact["givens"]==44
    assert all(exact["puzzle"][c]=="." for c in core)
    initial=initial_domains([int(d) if d!="." else 0 for d in exact["puzzle"]])
    assert initial==exact["initial_domains"]
    mask=sum(1<<(d-1) for d in palette)
    assert all(initial[c]==(mask|(1<<(extra-1)) if c==guardian else mask) for c in core)
    assert all(d=="." or int(d) not in palette for d in exact["puzzle"])
    reference=exact["respecting_completions"][0];orbit=set()
    for permutation in itertools.permutations(map(str,palette)):
        changed=reference.translate(str.maketrans(dict(zip(map(str,palette),permutation))))
        completion(changed,exact["puzzle"]);orbit.add(changed)
    assert len(orbit)==exact["solution_orbit_lower_bound"]==24
    for key in ("exact_initial_mask_realization","maximal_outside_stored_completion"):
        control=controls[key];solutions=control["respecting_completions"]
        assert len(set(solutions))==2
        for solution in solutions:completion(solution,control["puzzle"])
        assert len(exact_completions(control["puzzle"],2))==2
    for _,solution,_ in rows[:args.recount_sample]:
        puzzle="".join("." if c in core else solution[c] for c in range(81))
        assert len(exact_completions(puzzle,2))==2
    internal_path=HERE/"realization27_results_internal.json"
    if internal_path.exists():
        internal=json.loads(internal_path.read_text())
        assert len(internal["records"])==26 and internal["unique_cases"]==4
        assert internal["node_cutoffs"]==0 and internal["full_old_stalls"]==0
        assert sum(r["solutions_capped_at_2"]==1 for r in internal["records"])==4
        for control in internal["controls"]:
            verify_anchored(control,core,guardian)
            assert control["givens"]==55
    for name in ("realization27_results_anchored_control.json",
                 "realization27_results_anchored_minimal_control.json"):
        minimized_path=HERE/name
        if minimized_path.exists():
            verify_anchored(json.loads(minimized_path.read_text()),core,guardian)
    raw_path=HERE/"realization27_results_raw_old.json"
    if raw_path.exists():
        raw=json.loads(raw_path.read_text())
        assert raw["baseline_commit"]=="0532d63d9241a69c1fb57e9c34a1533effe1981d"
        assert {r["givens"] for r in raw["records"]}=={21,24} and len(raw["records"])==2
        for record in raw["records"]:verify_raw_old(record)
    print("4999 compatible records, exact controls, anchored uniqueness/traces and direct raw full-older-engine replay: PASS. Discovery labels are not all independently recounted.")


def verify_anchored(control,core,guardian):
    puzzle,solution=control["puzzle"],control["solution"]
    assert puzzle[guardian]=="."
    anchors=[c for c in core if puzzle[c]!="."]
    assert anchors==[control["internal_given_cell"]]
    assert exact_completions(puzzle)==[solution]
    initial=initial_domains([int(d) if d!="." else 0 for d in puzzle])
    assert initial==control["initial_domains"]
    if "initial_state" in control:
        assert control["initial_state"]=={"version":1,
            "values":[int(d) if d!="." else 0 for d in puzzle],"domains":initial}
        assert control["graph_cells"]==core
        assert control["conclusion"]=={"cell":guardian,"digit":1,"value":True,"kind":"placement"}
    final=replay_trace(initial,control["base_trace"],list(map(int,solution)))
    assert final==control["base_domains"]==control["state"]["domains"]
    assert all(mask&(mask-1)==0 for mask in final)
    assert control["status"]=="solved"
    comparison=control["full_old_comparison"]
    assert comparison["state"]==control["state"]
    assert comparison["full_dispatcher_executed"] is True
    assert comparison["full_current_finding"] is None and comparison["trace"]==[]
    if "local_domains" in control:
        assert replay_trace(initial,control["local_trace"],list(map(int,solution)))==control["local_domains"]
        assert control["local_domains"][guardian]==1


def verify_raw_old(record):
    source=HERE/record["fixture"]
    assert hashlib.sha256(source.read_bytes()).hexdigest()==record["fixture_sha256"]
    case=json.loads(source.read_text());solution=record["solution"]
    assert record["puzzle"]==case["puzzle"] and solution==case["solution"]
    assert exact_completions(record["puzzle"])==[solution]
    assert record["initial_state"]==case["initial_state"]
    assert record["initial_state"]["values"]==[int(d) if d!="." else 0 for d in record["puzzle"]]
    assert record["external_hall_or_three_house_closure"] is False
    assert record["promoted_unplaced_singletons"] is False
    assert len(record["full_older_calls"])==1
    call=record["full_older_calls"][0];returned=call["output"]
    assert call["input_state"]==record["initial_state"]
    assert returned["full_dispatcher_executed"] is True
    assert returned["arithmetic"]["error"] is None and returned["three_house"]["error"] is None
    assert returned["state"]==record["end_state"] and returned["trace"]==record["all_findings"]
    assert all(event["kind"]=="existing-engine" for event in record["trace"])
    assert [s for e in record["trace"] for s in e["steps"]]==record["all_findings"]
    final=replay_trace(record["initial_state"]["domains"],record["trace"],list(map(int,solution)))
    assert final==record["end_state"]["domains"]
    values=record["initial_state"]["values"].copy()
    for _,cell,digits,placement in record["all_findings"]:
        if placement:
            assert values[cell]==0 and digits==[int(solution[cell])]
            values[cell]=digits[0]
    assert values==record["end_state"]["values"]==list(map(int,solution))
    assert record["status"]=="solved-by-full-older-engine" and returned["full_current_finding"] is None


if __name__=="__main__":
    if not __debug__:raise SystemExit("Run without -O")
    main()
