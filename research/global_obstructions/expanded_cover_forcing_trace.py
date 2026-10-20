"""Replay a complete frozen-old forcing closure on a certified unique fixture.

The independent exact-cover count establishes that the recorded completion is
the sole one. Every older placement/elimination is then checked against that
completion, every domain update is retained, and all local-GAC steps are replayed
by the independent matching checker. This validates the deduction endpoints;
it does not independently reconstruct each older named technique's inner proof.
"""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path

from expanded_cover_search import BASELINE, Engine, candidate_state, cached_engine
from real_four_digit_mine import initial_domains
from real_verify import close, exact_completions, replay_trace, sees

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fixture", type=Path,
                        default=HERE / "expanded_unique_nonforcing.json")
    parser.add_argument("--output", type=Path,
                        default=ROOT / "target/expanded-cover-old-forcing-29.json")
    args = parser.parse_args()
    args.fixture = args.fixture.resolve()
    case = json.loads(args.fixture.read_text())
    puzzle, solution = case["puzzle"], case["solution"]
    assert exact_completions(puzzle) == [solution]
    initial = initial_domains([int(c) if c in "123456789" else 0 for c in puzzle])
    assert initial == case["initial_domains"]
    trace = case["base_trace"].copy()
    domains = replay_trace(initial, trace, list(map(int, solution)))
    assert domains == case["base_domains"]
    target = case["target"]
    guardian, extra = target["cell"], target["digit"]
    engine = Engine(cached_engine())
    first_full, findings = None, []
    result = None
    for _ in range(100):
        if candidate_state(domains)["values"][guardian] != 0:
            assert domains[guardian] == 1 << (extra - 1)
            result = "guardian-placed-by-older-methods"
            break
        before = domains.copy()
        older = engine.run(puzzle, domains, True, True)
        if older["trace"]:
            findings.extend(older["trace"])
            trace.append({"kind": "existing-engine", "steps": older["trace"],
                          "state_after": older["state"]})
            domains, local = close(older["state"]["domains"])
            trace.extend(local)
        else:
            finding = older["full_current_finding"]
            if finding is None:
                result = "full-old-stall-with-guardian-live"
                break
            if first_full is None:
                first_full = finding
            findings.append(finding)
            _, cell, digits, placement = finding
            if placement:
                assert digits == [int(solution[cell])]
                domains[cell] = 1 << (digits[0] - 1)
                for peer in range(81):
                    if peer != cell and sees(cell, peer):
                        domains[peer] &= ~(1 << (digits[0] - 1))
            else:
                assert int(solution[cell]) not in digits
                for digit in digits:
                    domains[cell] &= ~(1 << (digit - 1))
            trace.append({"kind": "existing-engine", "steps": [finding],
                          "state_after": candidate_state(domains)})
            domains, local = close(domains)
            trace.extend(local)
        assert domains != before
        assert all(mask & (1 << (int(digit) - 1))
                   for mask, digit in zip(domains, solution))
    else:
        raise ValueError("Full closure exceeded monotone bound")
    engine.close()
    replayed = replay_trace(initial, trace, list(map(int, solution)))
    assert replayed == domains
    output = {"schema_version": 1, "source_fixture": str(args.fixture.relative_to(ROOT)),
              "source_fixture_sha256": hashlib.sha256(args.fixture.read_bytes()).hexdigest(),
              "baseline_commit": BASELINE,
              "binary_sha256": hashlib.sha256(cached_engine().read_bytes()).hexdigest(),
              "puzzle": puzzle, "solution": solution, "solutions_capped_at_2": 1,
              "initial_domains": initial, "start_state": case["state"],
              "target": target, "first_full_finding": first_full,
              "findings_after_start": findings,
              "technique_counts_after_start": dict(Counter(step[0] for step in findings)),
              "trace_from_givens": trace, "end_domains": domains,
              "end_state": candidate_state(domains), "result": result,
              "independent_completion_validation": True,
              "independent_local_trace_replay": True,
              "limitation": "Older technique endpoints are completion-validated; their individual named inner proofs are not independently reconstructed here."}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps({"result": result, "first_full_finding": first_full,
                      "techniques": output["technique_counts_after_start"],
                      "output": str(args.output)}))


if __name__ == "__main__":
    main()
