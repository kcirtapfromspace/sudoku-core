"""Bounded unique-givens discovery, alternating both offline baseline parts.

Only an unchanged, completed strict catalogue/chain/JE/MSLS/forcing closure and
an unchanged retained full older dispatcher qualify as a configured union fixed
point. A cutoff never qualifies. Detailed positive certificates stay in output
artifacts and require independent strong-baseline acceptance.
"""
import argparse
from collections import Counter
import json
from pathlib import Path
import time

from discovery_pipeline_compare import finish_old, sha
from nonrank_puzzles_search import candidate_state
from real_four_digit_mine import initial_domains
from real_verify import close, exact_completions
from strong_baseline_core import close_state, Options

HERE = Path(__file__).resolve().parent


def main():
    if not __debug__:
        raise SystemExit("Run without -O")
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--indices", type=int, nargs="+", required=True)
    p.add_argument("--build-folder", type=Path, required=True)
    p.add_argument("--seconds", type=float, default=240)
    p.add_argument("--host-seconds", type=float, default=40)
    p.add_argument("--call-seconds", type=float, default=10)
    p.add_argument("--output-dir", type=Path, required=True)
    a = p.parse_args()
    binary = (a.build_folder / "instrumented/target/release/examples/real_existing_probe").resolve()
    instrument = json.loads((a.build_folder / "retained_manifest.json").read_text())
    assert sha(binary.read_bytes()) == instrument["binary_sha256"]
    corpus_path = HERE / "expanded_modern_corpus.json"
    corpus = json.loads(corpus_path.read_text())
    a.output_dir.mkdir(parents=True, exist_ok=True)
    started = time.monotonic()
    records = []
    for index in a.indices:
        if time.monotonic() - started >= a.seconds:
            break
        item = corpus[index]
        began = time.monotonic()
        deadline = min(started + a.seconds, began + a.host_seconds)
        print("START", index, item["id"], flush=True)
        solutions = exact_completions(item["puzzle"], 2)
        row = {"corpus_index": index, "id": item["id"], "puzzle": item["puzzle"]}
        if len(solutions) != 1:
            row["status"] = "not-unique"
        else:
            solution = solutions[0]
            initial = initial_domains([int(c) if c in "123456789" else 0 for c in item["puzzle"]])
            domains, trace = close(initial)
            case = {"schema_version": 1, "name": item["id"], "puzzle": item["puzzle"], "source": item,
                    "solution": solution, "completion": solution, "solutions_capped_at_2": 1,
                    "initial_domains": initial, "base_domains": domains, "base_trace": trace,
                    "state": candidate_state(domains), "status": "local-seed"}
            stages, relations, status = [], None, "round-cutoff"
            for _ in range(20):
                before = case["base_domains"].copy()
                old_status, domains, trace, calls, last = finish_old(
                    binary, case, deadline, a.call_seconds, checked_solution=list(map(int, solution)))
                stages.append({"phase": "retained-old", "status": old_status,
                               "input_domains": before, "domains": domains, "trace": trace,
                               "calls": calls, "last_probe": last})
                case.update(base_domains=domains, base_trace=case["base_trace"] + trace,
                            state=candidate_state(domains), status=old_status)
                if old_status == "solved":
                    status = "solved-retained-old"
                    break
                if old_status != "retained-configured-fixedpoint":
                    status = old_status
                    break
                remaining = deadline - time.monotonic()
                if remaining <= 0:
                    status = "campaign-cutoff"
                    break
                options = Options(context_id=item["id"], seconds=remaining,
                                  catalogue_seconds=min(30, remaining), catalogue_nodes=1000000,
                                  history_snapshots=(initial,))
                stronger = close_state(domains, inherited_relations=relations, options=options)
                relations = stronger["relations"]
                stages.append({"phase": "strong", "result": stronger})
                domains = stronger["domains"]
                assert all(m & (1 << (int(d) - 1)) for m, d in zip(domains, solution))
                assert stronger["status"] != "refuted", "Unique base state refuted"
                if stronger["status"] == "cutoff":
                    status = "strong-cutoff"
                    break
                if all(not m & (m - 1) for m in domains):
                    status = "solved-strong"
                    break
                if stronger["status"] == "fixedpoint" and domains == before:
                    status = "configured-union-fixedpoint"
                    break
                # Strong steps have their own exact grammar. Keep the old trace
                # at its old prefix; the stage log, not an existing-engine event,
                # authenticates the next seed for continued comparison.
                case.update(base_domains=domains, base_trace=[], state=candidate_state(domains))
            artifact = {"schema_version": 1, "id": item["id"], "puzzle": item["puzzle"],
                        "solution": solution, "source": item, "initial_domains": initial,
                        "local_trace": close(initial)[1], "stages": stages, "status": status,
                        "domains": domains,
                        "scope": "Discovery artifact. Strong-phase proof replay and independent audit are required before acceptance; old named inner proofs remain comparison instrumentation."}
            destination = a.output_dir / (item["id"] + ".json")
            destination.write_text(json.dumps(artifact, separators=(",", ":")) + "\n")
            row.update(status=status, solution=solution,
                       unplaced=sum(bool(m & (m - 1)) for m in domains),
                       stages=len(stages), artifact_path=str(destination),
                       artifact_sha256=sha(destination.read_bytes()))
        row["seconds"] = time.monotonic() - began
        records.append(row)
        summary = {"schema_version": 1, "corpus_sha256": sha(corpus_path.read_bytes()),
                   "instrument": instrument, "requested_indices": a.indices,
                   "limits": {"campaign_seconds": a.seconds, "host_seconds": a.host_seconds,
                              "old_call_seconds": a.call_seconds, "catalogue_nodes": 1000000},
                   "counts": dict(Counter(r["status"] for r in records)), "records": records,
                   "seconds": time.monotonic() - started,
                   "scope": "New unique author-source hosts. Completed solved, configured union fixedpoint and every cutoff are separate. No historical novelty claim and no unverified promotion."}
        (a.output_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
        print("DONE", item["id"], row["status"], round(row["seconds"], 2), flush=True)
    print(json.dumps({"counts": summary["counts"], "seconds": summary["seconds"]}))


if __name__ == "__main__":
    main()
