"""Continue the pinned modern corpus with complete frozen-engine comparisons."""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import time

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from nonrank_puzzles_search import BASELINE, ROOT, candidate_state, full_close
import nonrank_puzzles_search as frozen
from real_verify import exact_completions, replay_trace


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--start", type=int, default=106)
    p.add_argument("--stop", type=int, default=150)
    p.add_argument("--seconds", type=float, default=720)
    p.add_argument("--output", type=Path, default=ROOT / "target/different-corpus")
    p.add_argument("--binary", type=Path, default=ROOT / "target/nonrank-stdlib-baseline/instrumented/target/release/examples/real_existing_probe")
    a = p.parse_args()
    a.output = a.output.resolve()
    a.output.mkdir(parents=True, exist_ok=True)
    data = (HERE / "expanded_modern_corpus.json").read_bytes()
    cases = json.loads(data)
    counts, records, started = Counter(), [], time.monotonic()
    for index in range(a.start, min(a.stop, len(cases))):
        if time.monotonic() - started >= a.seconds:
            break
        case, beginning = cases[index], time.monotonic()
        print("START", index, case["id"], flush=True)
        record = {"corpus_index": index, "id": case["id"], "puzzle": case["puzzle"]}
        fixedpoint_calls = []
        original_probe = frozen.old_probe
        def recorded_probe(binary, puzzle, domains, timeout):
            answer = original_probe(binary, puzzle, domains, timeout)
            assert answer["full_dispatcher_executed"] is True
            if not answer["trace"] and answer["full_current_finding"] is None:
                fixedpoint_calls.append({"input_state": candidate_state(domains), "response": answer})
            return answer
        frozen.old_probe = recorded_probe
        try:
            solutions = exact_completions(case["puzzle"], 2)
            if len(solutions) != 1:
                raise ValueError("not unique: " + str(len(solutions)))
            sol = list(map(int, solutions[0]))
            status, initial, domains, trace, findings, calls, comparison = full_close(
                a.binary, case, sol, 45, 1000)
            assert replay_trace(initial, trace, sol) == domains
            if status == "full-old-stall":
                assert len(fixedpoint_calls) >= 2
                assert fixedpoint_calls[-1]["input_state"] == fixedpoint_calls[-2]["input_state"]
                assert all(not c["response"]["arithmetic"]["budget_exhausted"]
                           and c["response"]["three_house"]["exhausted"]
                           and c["response"]["three_house"]["centers_checked"] == 81
                           for c in fixedpoint_calls[-2:])
                assert comparison["full_dispatcher_executed"] is True
                assert comparison["arithmetic"]["error"] is None
                assert comparison["three_house"]["error"] is None
            output = {"schema_version": 1, "name": case["id"], "puzzle": case["puzzle"],
                      "source": case, "solution": solutions[0], "completion": solutions[0],
                      "solutions_capped_at_2": 1, "initial_domains": initial,
                      "base_domains": domains, "base_trace": trace,
                      "state": candidate_state(domains), "baseline_commit": BASELINE,
                      "status": status, "old_findings": findings, "old_calls": calls,
                      "full_old_fixedpoint": comparison,
                      "full_old_fixedpoint_calls": fixedpoint_calls[-2:] if status == "full-old-stall" else [],
                      "all_old_endpoints_completion_validated": True}
            destination = a.output / (case["id"] + ".json")
            destination.write_text(json.dumps(output, separators=(",", ":")) + "\n")
            record.update(status=status, old_findings=len(findings), calls=calls,
                          unplaced=sum(bool(m & (m-1)) for m in domains),
                          path=str(destination.relative_to(ROOT)),
                          state_sha256=hashlib.sha256(json.dumps(output["state"], sort_keys=True).encode()).hexdigest())
        except subprocess.TimeoutExpired:
            record["status"] = "probe-timeout"
        except (ValueError, AssertionError) as error:
            record.update(status="invalid-or-error", error=repr(error))
        finally:
            frozen.old_probe = original_probe
        record["seconds"] = round(time.monotonic() - beginning, 5)
        counts[record["status"]] += 1
        records.append(record)
        summary = {"schema_version": 1, "command": " ".join(sys.argv),
                   "baseline_commit": BASELINE,
                   "binary_sha256": hashlib.sha256(a.binary.read_bytes()).hexdigest(),
                   "corpus_sha256": hashlib.sha256(data).hexdigest(),
                   "counts": dict(counts), "records": records,
                   "seconds": time.monotonic() - started}
        (a.output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
        print("DONE", index, record["status"], record["seconds"], flush=True)


if __name__ == "__main__":
    main()
