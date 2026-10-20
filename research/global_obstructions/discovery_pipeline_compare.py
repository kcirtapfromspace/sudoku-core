"""Compare unique research states with the candidate-preserving frozen ladder.

The retained dispatcher is a bounded discovery baseline. Its endpoints are
checked against independently counted unique completions; its internal named
method proofs are not independently reconstructed by this script.
"""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import subprocess
import time

from nonrank_puzzles_search import candidate_state, old_probe, valid_step
from real_four_digit_mine import initial_domains
from real_verify import close, exact_completions, replay_trace, sees

HERE = Path(__file__).resolve().parent


def sha(value):
    return hashlib.sha256(value).hexdigest()


def verify_input(case):
    if exact_completions(case["puzzle"], 2) != [case["solution"]]:
        raise ValueError("A research host must be independently unique")
    initial = initial_domains([int(c) if c in "123456789" else 0 for c in case["puzzle"]])
    if initial != case["initial_domains"]:
        raise ValueError("Givens do not match initial candidates")
    solution = list(map(int, case["solution"]))
    if replay_trace(initial, case["base_trace"], solution) != case["base_domains"]:
        raise ValueError("Premise trace does not reach base state")
    if case["state"] != candidate_state(case["base_domains"]):
        raise ValueError("Base state mismatch")
    if close(case["base_domains"])[1]:
        raise ValueError("Expected Hall/81 closure")
    return solution


def finish_old(binary, case, deadline, call_seconds=30, steps=100, checked_solution=None):
    # A staged campaign authenticates strong-rule states with its own grammar;
    # they cannot be replayed as an old-engine-only premise trace here.
    solution = verify_input(case) if checked_solution is None else checked_solution
    if solution != list(map(int, case["solution"])):
        raise ValueError("Checked completion identity mismatch")
    if not all(m & (1 << (d - 1)) for m, d in zip(case["base_domains"], solution)):
        raise ValueError("Seed excludes the checked completion")
    domains = case["base_domains"].copy()
    trace, calls = [], 0
    for _ in range(steps):
        if all(not m & (m - 1) for m in domains):
            return "solved", domains, trace, calls, None
        remaining = deadline - time.monotonic()
        if remaining <= 0:
            return "campaign-cutoff", domains, trace, calls, None
        try:
            old = old_probe(binary, case["puzzle"], domains, min(call_seconds, remaining))
            calls += 1
        except subprocess.TimeoutExpired:
            return "probe-timeout", domains, trace, calls, None
        if old["trace"]:
            for step in old["trace"]:
                valid_step(step, solution)
            event = {"kind": "existing-engine", "steps": old["trace"], "state_after": old["state"]}
            domains = replay_trace(domains, [event], solution)
            trace.append(event)
        else:
            if old["state"]["domains"] != domains:
                raise ValueError("Probe changed domains without recording endpoints")
            finding = old["full_current_finding"]
            if finding is None:
                if not old["full_dispatcher_executed"]:
                    raise ValueError("Older full dispatcher was skipped")
                complete = (not old["arithmetic"]["budget_exhausted"]
                            and old["arithmetic"]["hint"] is None
                            and old["three_house"]["exhausted"]
                            and old["three_house"]["centers_checked"] == 81
                            and old["three_house"]["hint"] is None)
                return ("retained-configured-fixedpoint" if complete else "detector-cutoff"), domains, trace, calls, old
            valid_step(finding, solution)
            _, cell, digits, placement = finding
            if placement:
                bit = 1 << (digits[0] - 1)
                domains[cell] = bit
                for peer in range(81):
                    if peer != cell and sees(cell, peer):
                        domains[peer] &= ~bit
            else:
                for digit in digits:
                    domains[cell] &= ~(1 << (digit - 1))
            trace.append({"kind": "existing-engine", "steps": [finding],
                          "state_after": candidate_state(domains)})
        domains, local = close(domains)
        trace.extend(local)
        if not all(m & (1 << (d - 1)) for m, d in zip(domains, solution)):
            raise ValueError("An endpoint removed the unique completion")
    return "step-cutoff", domains, trace, calls, None


def main():
    if not __debug__:
        raise SystemExit("Run without -O")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs", type=Path, default=HERE / "discovery_pipeline_inputs.json")
    parser.add_argument("--build-folder", type=Path, required=True)
    parser.add_argument("--seconds", type=float, default=240)
    parser.add_argument("--call-seconds", type=float, default=30)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    binary = (args.build_folder / "instrumented/target/release/examples/real_existing_probe").resolve()
    manifest = json.loads((args.build_folder / "retained_manifest.json").read_text())
    if sha(binary.read_bytes()) != manifest["binary_sha256"]:
        raise ValueError("Retained binary identity mismatch")
    inputs = json.loads(args.inputs.read_text())
    started = time.monotonic()
    records = []
    for case in inputs["cases"]:
        status, domains, trace, calls, old = finish_old(binary, case, started + args.seconds, args.call_seconds)
        updated = dict(case, base_domains=domains, base_trace=case["base_trace"] + trace,
                       state=candidate_state(domains), status=status)
        verify_input(updated)
        record = {"id": case["name"], "status": status, "calls": calls,
                  "additional_trace_events": len(trace), "unplaced": sum(bool(m & (m - 1)) for m in domains),
                  "fixture": updated, "last_probe": old}
        records.append(record)
        result = {"schema_version": 1, "input_sha256": sha(args.inputs.read_bytes()),
                  "instrument": manifest, "counts": dict(Counter(r["status"] for r in records)),
                  "seconds": time.monotonic() - started, "records": records,
                  "limits": {"campaign_seconds": args.seconds, "call_seconds": args.call_seconds, "steps_per_host": 100},
                  "scope": "Retained-candidate frozen configured ladder only; NOT the union fixed point with catalogue/chains/Exocet/MSLS/deeper forcing. All input and resulting premise traces replayed and hosts independently counted unique. Named inner methods remain comparison instrumentation."}
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2) + "\n")
        print(case["name"], status, record["unplaced"], flush=True)
    print(json.dumps({"counts": result["counts"], "seconds": result["seconds"]}))


if __name__ == "__main__":
    main()
