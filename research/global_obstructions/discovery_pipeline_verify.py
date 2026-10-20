"""Replay the public retained-candidate comparison and unique source premises.

This checks givens, exact solution counts, monotone domains and outer endpoint
traces. It deliberately does not certify the frozen detectors' internal proofs.
"""
import argparse
from collections import Counter
import json
from pathlib import Path

from discovery_pipeline_compare import sha, verify_input
from nonrank_puzzles_search import candidate_state
from real_verify import replay_trace

HERE = Path(__file__).resolve().parent


def structural_invariant(puzzle):
    """Unequal signatures prove different native symmetry orbits; converse false."""
    given = [c in "123456789" for c in puzzle]
    rows = sorted(sum(given[9*r+c] for c in range(9)) for r in range(9))
    cols = sorted(sum(given[9*r+c] for r in range(9)) for c in range(9))
    boxes = sorted(sum(given[9*(3*b+r)+3*s+c] for r in range(3) for c in range(3))
                   for b in range(3) for s in range(3))
    digits = sorted(puzzle.count(str(d)) for d in range(1, 10))
    return tuple(digits), tuple(boxes), tuple(sorted((tuple(rows), tuple(cols))))


def export(source):
    original = json.loads(source.read_text())
    inputs_path = HERE / "discovery_pipeline_inputs.json"
    inputs = {c["name"]: c for c in json.loads(inputs_path.read_text())["cases"]}
    records = []
    for row in original["records"]:
        source_case = inputs[row["id"]]
        fixture = row["fixture"]
        prefix = len(source_case["base_trace"])
        if fixture["base_trace"][:prefix] != source_case["base_trace"]:
            raise ValueError("Comparison changed its input history")
        records.append({"id": row["id"], "status": row["status"], "calls": row["calls"],
                        "additional_trace": fixture["base_trace"][prefix:],
                        "domains": fixture["base_domains"], "unplaced": row["unplaced"],
                        "last_probe": row["last_probe"]})
    result = {k: original[k] for k in ("schema_version", "input_sha256", "instrument", "counts", "seconds", "limits", "scope")}
    result["records"] = records
    destination = HERE / "discovery_pipeline_retained_results.json"
    destination.write_text(json.dumps(result, separators=(",", ":")) + "\n")


def verify():
    inputs_path = HERE / "discovery_pipeline_inputs.json"
    inputs = json.loads(inputs_path.read_text())
    result = json.loads((HERE / "discovery_pipeline_retained_results.json").read_text())
    assert result["input_sha256"] == sha(inputs_path.read_bytes())
    cases = {c["name"]: c for c in inputs["cases"]}
    assert len(cases) == len(inputs["cases"]) == 8
    assert {r["id"] for r in result["records"]} == set(cases)
    assert len(result["records"]) == len(cases)
    assert result["counts"] == dict(Counter(r["status"] for r in result["records"]))
    steps = 0
    for row in result["records"]:
        case = cases[row["id"]]
        solution = verify_input(case)
        domains = replay_trace(case["base_domains"], row["additional_trace"], solution)
        assert domains == row["domains"]
        assert all(not new & ~old for new, old in zip(domains, case["base_domains"]))
        assert row["unplaced"] == sum(bool(m & (m - 1)) for m in domains)
        if row["status"] == "solved":
            assert row["unplaced"] == 0
            assert domains == [1 << (d - 1) for d in solution]
        elif row["status"] == "retained-configured-fixedpoint":
            old = row["last_probe"]
            assert row["unplaced"] > 0
            assert old["state"] == candidate_state(domains)
            assert not old["trace"] and old["full_dispatcher_executed"] is True
            assert old["full_current_finding"] is None
            assert old["arithmetic"]["hint"] is None and old["arithmetic"]["error"] is None
            assert old["arithmetic"]["budget_exhausted"] is False
            assert old["three_house"]["hint"] is None and old["three_house"]["error"] is None
            assert old["three_house"]["centers_checked"] == 81 and old["three_house"]["exhausted"] is True
        else:
            assert row["status"] in ("campaign-cutoff", "probe-timeout", "detector-cutoff", "step-cutoff")
        steps += sum(len(t.get("steps", [])) for t in row["additional_trace"])
    return {"accepted": True, "independently_unique_inputs": len(cases),
            "provably_distinct_native_puzzle_orbits_at_least": len({structural_invariant(c["puzzle"]) for c in cases.values()}),
            "counts": result["counts"], "completion_validated_additional_endpoints": steps,
            "scope": "Retained frozen dispatcher comparison only; not a joint stronger-baseline fixed point or a historical novelty certificate."}


def main():
    if not __debug__:
        raise SystemExit("Run without -O")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--export", type=Path)
    args = parser.parse_args()
    if args.export:
        export(args.export)
    print(json.dumps(verify(), sort_keys=True))


if __name__ == "__main__":
    main()
