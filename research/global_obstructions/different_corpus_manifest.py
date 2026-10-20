"""Compact discovery manifest; no trial status is a mathematical novelty claim.

Export retained bulk logs, or check the public counts/source/completions using
only Python's standard library. --count independently recounts uniqueness.
--bulk also replays complete retained premise traces when available locally.
"""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from real_verify import exact_completions, replay_trace

TRIALS = HERE / "different_corpus_trials.jsonl"
RESULTS = HERE / "different_corpus_results.json"


def digest(data):
    return hashlib.sha256(data).hexdigest()


def export(folders):
    records, batches, interrupted = [], [], []
    for folder in folders:
        summary = json.loads((folder / "summary.json").read_text())
        batches.append({key: summary[key] for key in
                        ("command", "baseline_commit", "binary_sha256", "corpus_sha256", "counts", "seconds")})
        if (folder / "interruption.json").exists():
            interrupted.append(json.loads((folder / "interruption.json").read_text()))
        for row in summary["records"]:
            entry = dict(row)
            if row["status"] in ("solved", "full-old-stall", "step-limit"):
                fixture = json.loads((folder / (row["id"] + ".json")).read_text())
                entry.update(solution=fixture["solution"], independently_counted_completions_capped_at_2=1,
                             trace_sha256=digest(json.dumps(fixture["base_trace"], sort_keys=True).encode()))
                comparison = fixture["full_old_fixedpoint"]
                entry["retained_last_full_dispatcher_executed"] = comparison["full_dispatcher_executed"] if comparison else None
                if comparison:
                    entry.update(arithmetic_budget_exhausted=comparison["arithmetic"]["budget_exhausted"],
                                 arithmetic_beam_pruned=comparison["arithmetic"]["beam_pruned"],
                                 crossing_centers_checked=comparison["three_house"]["centers_checked"],
                                 crossing_scan_exhausted=comparison["three_house"]["exhausted"])
            entry.pop("path", None)
            records.append(entry)
    records.sort(key=lambda r: r["corpus_index"])
    encoded = "".join(json.dumps(r, sort_keys=True, separators=(",", ":")) + "\n" for r in records).encode()
    TRIALS.write_bytes(encoded)
    result = {"schema_version": 1, "trial_manifest_sha256": digest(encoded),
              "counts": dict(Counter(r["status"] for r in records)), "completed_trials": len(records),
              "batches": batches, "interrupted_excluded": interrupted,
              "meaning": "Discovery engine comparisons only, using the configured frozen dispatcher including forcing, alternating Hall/all81 crossing closure. Each completed puzzle was independently counted unique; all older endpoints were checked against its sole completion and premise traces replayed during the campaign. Arithmetic beam pruning is recorded separately: a configured engine stall is not an exhaustive failure of all imaginable arithmetic or forcing proofs. Public compact records bind counts, source puzzles and completions. They do not reproduce every old method's internal proof, prove novelty, or turn a timeout/cutoff/solved board into a stall.",
              "reproduction": "Run each batch.command from the repository root, then this script export against those output folders. Full trace logs remain under target; --bulk replays them. Fresh old-engine audit of any promoted deduction is separate."}
    RESULTS.write_text(json.dumps(result, indent=2) + "\n")
    return verify(False, folders)


def verify(count=False, folders=()):
    data = TRIALS.read_bytes()
    result = json.loads(RESULTS.read_text())
    assert digest(data) == result["trial_manifest_sha256"]
    records = [json.loads(line) for line in data.splitlines()]
    corpus_bytes = (HERE / "expanded_modern_corpus.json").read_bytes()
    corpus = json.loads(corpus_bytes)
    assert all(batch["corpus_sha256"] == digest(corpus_bytes) for batch in result["batches"])
    assert result["completed_trials"] == len(records)
    assert result["counts"] == dict(Counter(r["status"] for r in records))
    assert sum(sum(b["counts"].values()) for b in result["batches"]) == len(records)
    assert len({r["corpus_index"] for r in records}) == len(records)
    assert len({r["puzzle"] for r in records}) == len(records)
    bulk = {p.stem: p for folder in folders for p in folder.glob("modern*.json")}
    replayed = recounted = 0
    houses = [[9*r+c for c in range(9)] for r in range(9)]
    houses += [[9*r+c for r in range(9)] for c in range(9)]
    houses += [[9*(3*b+r)+3*s+c for r in range(3) for c in range(3)] for b in range(3) for s in range(3)]
    for row in records:
        source = corpus[row["corpus_index"]]
        assert (row["id"], row["puzzle"]) == (source["id"], source["puzzle"])
        if "solution" not in row:
            assert row["status"] not in ("solved", "full-old-stall")
            continue
        sol = list(map(int, row["solution"]))
        assert len(sol) == 81 and all({sol[c] for c in h} == set(range(1, 10)) for h in houses)
        assert all(p not in "123456789" or int(p) == sol[c] for c, p in enumerate(row["puzzle"]))
        if row["status"] == "solved":
            assert row["unplaced"] == 0
        if row["status"] == "full-old-stall":
            assert row["unplaced"] > 0 and row["retained_last_full_dispatcher_executed"] is True
            assert row["arithmetic_budget_exhausted"] is False
            assert row["crossing_centers_checked"] == 81 and row["crossing_scan_exhausted"] is True
        if count:
            assert exact_completions(row["puzzle"], 2) == [row["solution"]]
            recounted += 1
        if row["id"] in bulk:
            fixture = json.loads(bulk[row["id"]].read_text())
            assert fixture["puzzle"] == row["puzzle"] and fixture["solution"] == row["solution"]
            assert digest(json.dumps(fixture["state"], sort_keys=True).encode()) == row["state_sha256"]
            assert digest(json.dumps(fixture["base_trace"], sort_keys=True).encode()) == row["trace_sha256"]
            assert replay_trace(fixture["initial_domains"], fixture["base_trace"], sol) == fixture["base_domains"]
            if row["status"] == "full-old-stall":
                old = fixture["full_old_fixedpoint"]
                assert old["full_dispatcher_executed"] is True and old["full_current_finding"] is None
                assert not old["trace"] and old["state"]["domains"] == fixture["base_domains"]
                assert old["arithmetic"]["error"] is None and old["three_house"]["error"] is None
            replayed += 1
    return {"trials": len(records), "counts": result["counts"], "unique_recounts": recounted,
            "exact_trace_replays": replayed, "manifest_verified": True}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("mode", choices=("export", "verify"))
    p.add_argument("--bulk", type=Path, nargs="*", default=[])
    p.add_argument("--count", action="store_true")
    a = p.parse_args()
    print(json.dumps(export(a.bulk) if a.mode == "export" else verify(a.count, a.bulk), sort_keys=True))


if __name__ == "__main__":
    main()
