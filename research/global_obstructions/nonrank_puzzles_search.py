"""Find genuine full frozen-engine fixedpoints on the pinned modern200 corpus.

One worker; no timeout, cutoff, invalid step, or solved board is called a stall.
Every older endpoint is checked against an independently counted completion.
All stored domains survive placements; Hall and all81 crossing triples are
closed between older findings. Complete premise traces stay under target.
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
import shutil
import subprocess
import sys
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
sys.path.insert(0, str(HERE))
from real_four_digit_mine import initial_domains
from real_verify import close, exact_completions, replay_trace, sees

BASELINE = "0532d63d9241a69c1fb57e9c34a1533effe1981d"


def candidate_state(domains):
    return {"version": 1, "domains": domains,
            "values": [mask.bit_length() if mask & (mask - 1) == 0 else 0
                       for mask in domains]}


def cached_engine():
    """Reuse a frozen probe or build it with stdlib Python and Rust only."""
    previous = ROOT / "target/expanded-cover-old-engine/real_existing_probe"
    if previous.exists():
        return previous
    cache = ROOT / "target/nonrank-old-engine"
    binary = cache / "real_existing_probe"
    if not binary.exists():
        from real_existing_check import build_engine
        cache.mkdir(parents=True, exist_ok=True)
        folder = cache / "build"
        folder.mkdir()
        shutil.copy2(build_engine(folder), binary)
    return binary


def normalized(domains):
    for cell, mask in enumerate(domains):
        assert 0 < mask <= 511
        if mask & (mask - 1) == 0:
            assert all(not domains[peer] & mask for peer in range(81)
                       if peer != cell and sees(cell, peer))


def old_probe(binary, puzzle, domains, timeout):
    normalized(domains)
    request = {"id": "nonrank-puzzles", "puzzle": puzzle,
               "state": candidate_state(domains), "unique": True, "full": True}
    returned = subprocess.run([binary], input=json.dumps(request) + "\n", text=True,
                              capture_output=True, timeout=timeout)
    if returned.returncode:
        raise ValueError("Historical probe failed: " + returned.stderr[-300:])
    result = json.loads(returned.stdout)
    assert result["arithmetic"]["error"] is None and result["three_house"]["error"] is None
    return result


def valid_step(step, solution):
    _, cell, digits, placement = step
    assert digits == [solution[cell]] if placement else solution[cell] not in digits


def full_close(binary, case, solution, timeout, step_limit):
    initial = initial_domains([int(c) if c in "123456789" else 0 for c in case["puzzle"]])
    domains, trace = close(initial)
    findings, calls = [], 0
    for iteration in range(step_limit):
        assert all(mask & (1 << (digit - 1)) for mask, digit in zip(domains, solution))
        if all(mask & (mask - 1) == 0 for mask in domains):
            return "solved", initial, domains, trace, findings, calls, None
        older = old_probe(binary, case["puzzle"], domains, timeout)
        calls += 1
        if older["trace"]:
            for step in older["trace"]:
                valid_step(step, solution)
            findings.extend(older["trace"])
            trace.append({"kind": "existing-engine", "steps": older["trace"],
                          "state_after": older["state"]})
            reduced, local = close(older["state"]["domains"])
            assert reduced != domains
            domains = reduced
            trace.extend(local)
            continue
        assert older["state"]["domains"] == domains
        finding = older["full_current_finding"]
        if finding is None:
            # A second complete successful call confirms the live fixedpoint.
            repeated = old_probe(binary, case["puzzle"], domains, timeout)
            calls += 1
            assert not repeated["trace"] and repeated["full_current_finding"] is None
            assert repeated["state"]["domains"] == domains
            assert close(domains)[1] == []
            assert any(mask & (mask - 1) for mask in domains)
            return "full-old-stall", initial, domains, trace, findings, calls, repeated
        valid_step(finding, solution)
        findings.append(finding)
        _, cell, digits, placement = finding
        before = domains.copy()
        if placement:
            domains[cell] = 1 << (digits[0] - 1)
            for peer in range(81):
                if peer != cell and sees(cell, peer):
                    domains[peer] &= ~(1 << (digits[0] - 1))
        else:
            for digit in digits:
                assert domains[cell] & (1 << (digit - 1))
                domains[cell] &= ~(1 << (digit - 1))
        assert domains != before
        trace.append({"kind": "existing-engine", "steps": [finding],
                      "state_after": candidate_state(domains)})
        domains, local = close(domains)
        trace.extend(local)
    return "step-limit", initial, domains, trace, findings, calls, None


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--limit", type=int, default=200)
    parser.add_argument("--seconds", type=float, default=1080)
    parser.add_argument("--call-timeout", type=float, default=45)
    parser.add_argument("--step-limit", type=int, default=1000)
    parser.add_argument("--order", choices=("hard-first", "source"), default="hard-first")
    parser.add_argument("--output-dir", type=Path, default=ROOT / "target/nonrank-puzzles-modern200")
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    corpus_path = HERE / "expanded_modern_corpus.json"
    corpus_bytes = corpus_path.read_bytes()
    corpus = json.loads(corpus_bytes)
    if args.order == "hard-first":
        corpus = corpus[100:] + corpus[:100]
    corpus = corpus[:args.limit]
    binary = cached_engine()
    started, counts, records = time.monotonic(), Counter(), []
    for item in corpus:
        if time.monotonic() - started >= args.seconds:
            break
        print("START", item["id"], "elapsed", round(time.monotonic()-started, 2), flush=True)
        record = {"id": item["id"], "puzzle": item["puzzle"]}
        counts["started"] += 1
        beginning = time.monotonic()
        try:
            completions = exact_completions(item["puzzle"], 2)
            if len(completions) != 1:
                counts["not-unique"] += 1
                record["status"] = "not-unique"
                records.append(record)
                continue
            counts["independently-unique"] += 1
            solution = list(map(int, completions[0]))
            status, initial, domains, trace, findings, calls, comparison = full_close(
                binary, item, solution, args.call_timeout, args.step_limit)
            replayed = replay_trace(initial, trace, solution)
            assert replayed == domains
            output = {"schema_version": 1, "name": item["id"], "puzzle": item["puzzle"],
                      "source": item, "solution": completions[0], "completion": completions[0],
                      "solutions_capped_at_2": 1, "initial_domains": initial,
                      "base_domains": domains, "base_trace": trace,
                      "state": candidate_state(domains), "baseline_commit": BASELINE,
                      "status": status, "old_findings": findings, "old_calls": calls,
                      "full_old_fixedpoint": comparison,
                      "all_old_endpoints_completion_validated": True}
            destination = args.output_dir / (item["id"] + ".json")
            destination.write_text(json.dumps(output, indent=2) + "\n")
            counts[status] += 1
            record.update(status=status, old_findings=len(findings), calls=calls,
                          unplaced=sum(mask & (mask-1) != 0 for mask in domains),
                          domains_sha256=hashlib.sha256(json.dumps(domains).encode()).hexdigest(),
                          trace_sha256=hashlib.sha256(json.dumps(trace).encode()).hexdigest())
            if status == "full-old-stall":
                print("FULL_STALL", destination, "unplaced", record["unplaced"], flush=True)
        except subprocess.TimeoutExpired:
            counts["probe-timeout"] += 1
            record["status"] = "probe-timeout"
        except (ValueError, AssertionError) as error:
            counts["invalid-or-error"] += 1
            record.update(status="invalid-or-error", error=repr(error))
        record["seconds"] = time.monotonic() - beginning
        records.append(record)
        print("DONE", record["id"], record["status"], "seconds", round(record["seconds"], 2), flush=True)
        summary = {"schema_version": 1, "baseline_commit": BASELINE,
                   "binary_sha256": hashlib.sha256(binary.read_bytes()).hexdigest(),
                   "corpus_sha256": hashlib.sha256(corpus_bytes).hexdigest(),
                   "command": " ".join(sys.argv), "counts": dict(counts),
                   "records": records, "seconds": time.monotonic()-started,
                   "scope": "Complete old dispatcher including forcing alternating Hall/all81 crossing closure. Recorded statuses are discovery/engine-comparison outcomes; no new deduction or LP acceptance claimed."}
        (args.output_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")


if __name__ == "__main__":
    if not __debug__:
        raise SystemExit("Run without -O")
    main()
