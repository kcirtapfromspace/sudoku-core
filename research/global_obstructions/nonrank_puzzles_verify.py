"""Stdlib acceptance of retained old-engine comparisons and pattern identity.

This never launches an engine, compiler or network request. It replays the
premises and validates saved successful full-dispatcher comparisons. Rebuilding
and calling the historical engine is the separate manual audit command.
"""
import hashlib
import json
from pathlib import Path
from collections import Counter

from real_four_digit_mine import initial_domains
from real_verify import close, exact_completions, graph_colorings, replay_trace, sees

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent


def main():
    results = json.loads((HERE / "nonrank_puzzles_results.json").read_text())
    metadata = results["manifest"]
    data = (HERE / metadata["file"]).read_bytes()
    assert len(data) == metadata["bytes"]
    assert hashlib.sha256(data).hexdigest() == metadata["sha256"]
    trials = [json.loads(line) for line in data.splitlines()]
    assert len(trials) == metadata["records"] == results["completed_cases"] == 6
    assert len({row[0] for row in trials}) == 6
    corpus_path = ROOT / results["source_corpus"]
    assert hashlib.sha256(corpus_path.read_bytes()).hexdigest() == results["source_corpus_sha256"]
    corpus = {entry["id"]: entry["puzzle"] for entry in json.loads(corpus_path.read_text())}
    statuses = metadata["status_codes"]
    assert all(len(row) == 6 and row[0] in corpus and row[1] == corpus[row[0]]
               and str(row[2]) in statuses for row in trials)
    assert dict(Counter(statuses[str(row[2])] for row in trials)) == results["status_counts"]
    assert results["status_counts"] == {"full-old-stall":5,"solved":1}
    assert all(row[4] > 0 if row[2] == 0 else row[4] == 0 for row in trials)

    case_path = HERE / "nonrank_full_old_stall.json"
    case = json.loads(case_path.read_text())
    audit = json.loads((HERE / "nonrank_puzzles_first_audit.json").read_text())
    assert audit["fixture"] == str(case_path.relative_to(ROOT))
    assert audit["fixture_sha256"] == hashlib.sha256(case_path.read_bytes()).hexdigest()
    assert audit["state_sha256"] == hashlib.sha256(json.dumps(case["state"], sort_keys=True,
                                separators=(",", ":")).encode()).hexdigest()
    assert audit["baseline_commit"] == "0532d63d9241a69c1fb57e9c34a1533effe1981d"
    assert exact_completions(case["puzzle"]) == [case["solution"]]
    solution = list(map(int, case["solution"]))
    initial = initial_domains([int(c) if c in "123456789" else 0 for c in case["puzzle"]])
    assert initial == case["initial_domains"]
    domains = replay_trace(initial, case["base_trace"], solution)
    assert domains == case["base_domains"] == case["state"]["domains"]
    first = next(row for row in trials if row[0] == case["source"]["id"])
    assert first[2] == 0 and first[3] == len(case["old_findings"])
    assert first[5] == hashlib.sha256(json.dumps(domains).encode()).hexdigest()
    assert close(domains)[1] == []
    assert sum(mask & (mask-1) != 0 for mask in domains) == audit["unplaced_cells"] == 50
    assert len(audit["full_older_engine_confirmations"]) == 2
    for comparison in audit["full_older_engine_confirmations"]:
        assert comparison["state"] == case["state"] and comparison["trace"] == []
        assert comparison["full_dispatcher_executed"] is True
        assert comparison["full_current_finding"] is None
        assert comparison["arithmetic"]["hint"] is None and comparison["arithmetic"]["error"] is None
        assert comparison["three_house"]["hint"] is None and comparison["three_house"]["error"] is None
        assert comparison["three_house"]["centers_checked"] == 81
    target = case["target"]
    raw = domains.copy()
    assert target == {"cell": 35, "digit": 9, "value": False}
    assert case["state"]["values"][35] == 0 and raw[35] & 256
    raw[35] &= ~256
    branch = replay_trace(raw, case["assumption_trace"])
    assert branch == case["assumed_domains"] and close(branch)[1] == []

    classification = json.loads((HERE / "nonrank_puzzles_first_classification.json").read_text())
    assert classification["fixture_sha256"] == audit["fixture_sha256"]
    raw_core, closed_core = classification["raw_list_core"], classification["closed_tridagon_core"]
    assert graph_colorings(raw, raw_core)[0] == 0
    assert graph_colorings(branch, closed_core)[0] == 0
    assert all(branch[cell] == sum(1 << (digit-1) for digit in (4,7,8)) for cell in closed_core)
    uniform = [sum(1 << (digit-1) for digit in (4,6,7,8))] * 81
    assert graph_colorings(uniform, raw_core)[0] == 1
    pairs = classification["isomorphism_to_published_Loki_control"]
    assert sorted(pair[0] for pair in pairs) == sorted(closed_core)
    known = json.loads((HERE / "real_Loki.json").read_text())
    assert sorted(pair[1] for pair in pairs) == sorted(known["graph_cells"])
    assert all(sees(a,c) == sees(b,d) for a,b in pairs for c,d in pairs)
    print("Public six-case manifest, independent uniqueness/premise replay, two retained fresh old-engine stall calls, live target, and Loki-tridagon graph isomorphism: PASS.")


if __name__ == "__main__":
    if not __debug__:
        raise SystemExit("Run without -O")
    main()
