"""Independent premise and stage replay for the fresh unique-host campaign."""
import argparse
from collections import Counter
import json
from pathlib import Path

from discovery_pipeline_compare import sha
from real_four_digit_mine import initial_domains
from real_verify import exact_completions, replay_trace

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "discovery_pipeline_expanded_results.json"
TRIALS = HERE / "discovery_pipeline_expanded_trials.jsonl"


def export(summary):
    result = json.loads(summary.read_text())
    lines = []
    for row in result["records"]:
        raw = Path(row["artifact_path"]).read_bytes()
        assert sha(raw) == row["artifact_sha256"]
        lines.append(raw)
    TRIALS.write_bytes(b"".join(lines))
    result["trial_manifest_sha256"] = sha(TRIALS.read_bytes())
    result["acceptance"] = "Negative discovery controls only; no qualified family or new engine rule."
    RESULTS.write_text(json.dumps(result, indent=2) + "\n")


def verify():
    data = TRIALS.read_bytes()
    result = json.loads(RESULTS.read_text())
    assert sha(data) == result["trial_manifest_sha256"]
    corpus_path = HERE / "expanded_modern_corpus.json"
    assert sha(corpus_path.read_bytes()) == result["corpus_sha256"]
    corpus = json.loads(corpus_path.read_text())
    lines = data.splitlines(keepends=True)
    assert len(lines) == len(result["records"])
    assert result["counts"] == dict(Counter(r["status"] for r in result["records"]))
    assert len({r["id"] for r in result["records"]}) == len(lines)
    endpoints, strong_steps = 0, 0
    for raw, row in zip(lines, result["records"]):
        assert sha(raw) == row["artifact_sha256"]
        trial = json.loads(raw)
        original = corpus[row["corpus_index"]]
        assert trial["source"] == original
        assert (trial["id"], trial["puzzle"]) == (original["id"], original["puzzle"])
        assert row["puzzle"] == trial["puzzle"]
        assert exact_completions(trial["puzzle"], 2) == [trial["solution"]]
        assert row["solution"] == trial["solution"]
        solution = list(map(int, trial["solution"]))
        initial = initial_domains([int(c) if c in "123456789" else 0 for c in trial["puzzle"]])
        assert initial == trial["initial_domains"]
        domains = replay_trace(initial, trial["local_trace"], solution)
        for stage in trial["stages"]:
            before = domains.copy()
            if stage["phase"] == "retained-old":
                assert stage["input_domains"] == domains
                domains = replay_trace(domains, stage["trace"], solution)
                assert domains == stage["domains"]
                endpoints += sum(len(e.get("steps", [])) for e in stage["trace"])
            else:
                assert stage["phase"] == "strong"
                from strong_baseline_verify import verify_result
                stronger = stage["result"]
                assert stronger["input_domains"] == domains
                assert list(initial) in list(stronger["options"]["history_snapshots"])
                assert stronger["options"]["context_id"] == trial["id"]
                accepted = verify_result(stronger)
                assert accepted["positive_steps_exact"]
                domains = stronger["domains"]
                strong_steps += accepted["events"]
            assert all(not new & ~old for new, old in zip(domains, before))
            assert all(m & (1 << (d - 1)) for m, d in zip(domains, solution))
        assert trial["domains"] == domains
        assert trial["status"] == row["status"]
        assert row["unplaced"] == sum(bool(m & (m - 1)) for m in domains)
        if row["status"] in ("solved-retained-old", "solved-strong"):
            assert domains == [1 << (d - 1) for d in solution]
        elif row["status"] == "configured-union-fixedpoint":
            # Both components must be completed and unchanged on this same state.
            old, strong = trial["stages"][-2:]
            assert old["phase"] == "retained-old" and old["status"] == "retained-configured-fixedpoint"
            assert old["input_domains"] == old["domains"] == domains
            assert strong["phase"] == "strong" and strong["result"]["status"] == "fixedpoint"
            assert strong["result"]["input_domains"] == strong["result"]["domains"] == domains
            probe = old["last_probe"]
            assert probe["full_dispatcher_executed"] is True and probe["full_current_finding"] is None
            assert not probe["trace"] and probe["state"]["domains"] == domains
            assert probe["arithmetic"]["error"] is None and probe["arithmetic"]["hint"] is None
            assert probe["arithmetic"]["budget_exhausted"] is False
            assert probe["three_house"]["error"] is None and probe["three_house"]["hint"] is None
            assert probe["three_house"]["exhausted"] is True and probe["three_house"]["centers_checked"] == 81
            options = strong["result"]["options"]
            assert all(options[k] is True for k in ("catalogue", "unit_binary", "junior", "relaxed_junior", "msls", "forcing_fresh_catalogue"))
            assert options["forcing_depth"] >= 2 and not options["forcing_targets"]
            assert accepted["catalogue_completeness_rescanned"] and accepted["configured_fixedpoint_repeated"]
        else:
            assert row["status"] in ("campaign-cutoff", "strong-cutoff", "probe-timeout", "detector-cutoff", "step-cutoff", "round-cutoff")
    return {"accepted": True, "independently_unique_hosts": len(lines),
            "counts": result["counts"], "completion_validated_old_endpoints": endpoints,
            "exact_strong_events": strong_steps,
            "scope": "Stage and premise audit. A completed older detector call is a configured comparison, not an exhaustive literature test."}


def main():
    if not __debug__:
        raise SystemExit("Run without -O")
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--export", type=Path)
    a = p.parse_args()
    if a.export:
        export(a.export)
    print(json.dumps(verify(), sort_keys=True))


if __name__ == "__main__":
    main()
