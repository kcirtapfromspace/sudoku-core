"""Run the complete frozen older engine directly on retained raw givens states.

No external Hall/81 closure or singleton promotion is performed. All subsequent
placements/eliminations are replayed explicitly, preserving values and masks.
The saved unique completion checks every endpoint; solved and stalled differ.
"""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import subprocess
import time

from real_verify import exact_completions,replay_trace,sees

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent.parent


def apply(state,step,solution):
    result={"version":1,"values":state["values"].copy(),"domains":state["domains"].copy()}
    _,cell,digits,placement=step
    assert result["values"][cell]==0
    if placement:
        assert digits==[int(solution[cell])] and result["domains"][cell]&(1<<(digits[0]-1))
        result["values"][cell]=digits[0];result["domains"][cell]=1<<(digits[0]-1)
        for peer in range(81):
            if peer!=cell and sees(cell,peer) and result["values"][peer]==0:
                result["domains"][peer]&=~(1<<(digits[0]-1))
    else:
        assert int(solution[cell]) not in digits
        for digit in digits:
            assert result["domains"][cell]&(1<<(digit-1))
            result["domains"][cell]&=~(1<<(digit-1))
    assert all(mask&(1<<(int(digit)-1)) for mask,digit in zip(result["domains"],solution))
    assert all(not value or value==int(solution[c]) for c,value in enumerate(result["values"]))
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary",type=Path,default=ROOT/"target/nonrank-stdlib-baseline/instrumented/target/release/examples/real_existing_probe")
    parser.add_argument("--seconds",type=float,default=100)
    parser.add_argument("--output",type=Path,default=HERE/"realization27_results_raw_old.json")
    args=parser.parse_args();started=time.monotonic();records=[]
    for name in ("realization27_results_anchored_control.json","realization27_results_anchored_minimal_control.json"):
        path=HERE/name;case=json.loads(path.read_text());puzzle=case["puzzle"];solution=case["solution"]
        assert exact_completions(puzzle)==[solution]
        state=case["initial_state"];assert state["values"]==[int(d) if d!="." else 0 for d in puzzle]
        trace=[];steps=[];calls=[];status=None
        for iteration in range(1000):
            if time.monotonic()-started>=args.seconds:
                status="time-budget-not-stall";break
            request={"id":"raw-initial-old","puzzle":puzzle,"state":state,"unique":True,"full":True}
            returned=subprocess.run([args.binary],input=json.dumps(request)+"\n",text=True,
                                    capture_output=True,timeout=min(20,args.seconds-(time.monotonic()-started)))
            assert returned.returncode==0,returned.stderr
            older=json.loads(returned.stdout);assert older["full_dispatcher_executed"] is True
            assert older["arithmetic"]["error"] is None and older["three_house"]["error"] is None
            request_state=state
            for step in older["trace"]:
                state=apply(state,step,solution);steps.append(step)
            assert state==older["state"]
            if older["trace"]:
                trace.append({"kind":"existing-engine","steps":older["trace"],"state_after":state})
            calls.append({"input_state":request_state,"output":older})
            if all(state["values"]):
                status="solved-by-full-older-engine";break
            finding=older["full_current_finding"]
            if finding is None:
                assert any(value==0 for value in state["values"])
                status="live-full-older-stall";break
            state=apply(state,finding,solution);steps.append(finding)
            trace.append({"kind":"existing-engine","steps":[finding],"state_after":state})
        else:status="step-budget-not-stall"
        assert replay_trace(case["initial_domains"],trace,list(map(int,solution)))==state["domains"]
        assert status=="solved-by-full-older-engine", "Retained anchored control no longer reproduces its solved outcome"
        record={"fixture":name,"fixture_sha256":hashlib.sha256(path.read_bytes()).hexdigest(),
                "puzzle":puzzle,"solution":solution,"givens":case["givens"],
                "initial_state":case["initial_state"],"end_state":state,"status":status,
                "external_hall_or_three_house_closure":False,"promoted_unplaced_singletons":False,
                "full_older_calls":calls,"trace":trace,"all_findings":steps,
                "technique_counts":dict(Counter(step[0] for step in steps))}
        records.append(record);print(name,status,record["technique_counts"],flush=True)
    output={"schema_version":1,"baseline_commit":"0532d63d9241a69c1fb57e9c34a1533effe1981d",
            "binary_sha256":hashlib.sha256(args.binary.read_bytes()).hexdigest(),"seconds":time.monotonic()-started,
            "scope":"Direct full frozen older engine from actual givens-only values and ordinary initial candidates. Its own configured methods are allowed; no external closure or singleton promotion.","records":records}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(output,indent=2)+"\n")


if __name__=="__main__":
    if not __debug__:raise SystemExit("Run without -O")
    main()
