"""Summarize retained broad-search logs without rerunning discovery."""
from collections import Counter
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
ENDPOINT_CODES = {"old-target-placed": 0, "joint-target-placed": 1,
                  "local-filter": 2, "post-old-local-filter": 3,
                  "nonforcing-arithmetic-gap": 4}


def manifest(filename, records, fields, code_maps):
    data = "".join(json.dumps(record, separators=(",", ":")) + "\n"
                   for record in records).encode()
    (HERE / filename).write_bytes(data)
    return {"file": filename, "records": len(records), "bytes": len(data),
            "sha256": hashlib.sha256(data).hexdigest(), "array_fields": fields,
            "code_maps": code_maps,
            "meaning": "Recorded discovery classifications only; this manifest does not independently certify every uniqueness claim or endpoint."}


def main():
    runs, all_puzzles, unique_raw, endpoints, gaps = [], set(), set(), {}, {}
    raw_status = {}
    for name in "abcdef":
        folder = ROOT / ("target/expanded-cover-milp-" + name)
        if not folder.exists():
            continue
        events_path = folder / "events.jsonl"
        events = [json.loads(line) for line in events_path.read_text().splitlines()]
        summary_path = folder / "summary.json"
        summary = json.loads(summary_path.read_text()) if summary_path.exists() else {}
        categories = Counter(event["result"] for event in events)
        logged = [event for event in events if "puzzle" in event]
        run = {"run": name, "summary": summary,
               "seed": summary.get("seed", 62841 if name == "a" else None),
               "events_sha256": hashlib.sha256(events_path.read_bytes()).hexdigest(),
               "fully_logged_candidate_sets": len(logged),
               "logged_categories": dict(categories),
               "role": "discarded preliminary pilot" if name == "a" else "bounded discovery run"}
        runs.append(run)
        if name == "a":
            continue
        for event in logged:
            all_puzzles.add(event["puzzle"])
            status = (-1 if event.get("enumeration_cutoff", False) and
                      event.get("solutions_capped", 0) < 2 else
                      1 if event.get("solutions_capped") == 1 else 0)
            assert raw_status.get(event["puzzle"], status) == status
            raw_status[event["puzzle"]] = status
            if event.get("solutions_capped") == 1 and not event.get("enumeration_cutoff", False):
                unique_raw.add(event["puzzle"])
        for path in folder.glob("*.json"):
            if path.name == "summary.json":
                continue
            case = json.loads(path.read_text())
            if "rejection_category" not in case:
                continue
            puzzle, category = case["puzzle"], case["rejection_category"]
            endpoints.setdefault(puzzle, set()).add(category)
            if "gap" in category or category == "full-stall":
                gaps.setdefault(puzzle, []).append({"file": str(path.relative_to(ROOT)),
                    "category": category, "givens": sum(c != "." for c in puzzle),
                    "forcing_result": case.get("strict_full_closure", {}).get("result"),
                    "first_full_finding": case.get("old_full_comparison", {}).get("full_current_finding")})
    endpoint_categories = Counter()
    conflicts = []
    for puzzle, categories in endpoints.items():
        for category in categories:
            endpoint_categories[category] += 1
        if len(categories) > 1:
            conflicts.append({"puzzle": puzzle, "categories": sorted(categories)})
    assert not conflicts
    trials_manifest = manifest("expanded_cover_trials.jsonl",
        [[puzzle, raw_status[puzzle]] for puzzle in sorted(all_puzzles)],
        ["puzzle", "discovery_uniqueness_status"],
        {"discovery_uniqueness_status": {-1: "enumeration-cutoff", 0: "discovery-nonunique",
                                           1: "discovery-unique"}})
    endpoint_rows = []
    for puzzle in sorted(endpoints):
        category = next(iter(endpoints[puzzle]))
        forcing = 0
        if puzzle in gaps:
            outcomes = {entry["forcing_result"] for entry in gaps[puzzle]
                        if entry["forcing_result"] is not None}
            assert outcomes == {"target-placed"}
            forcing = 1
        endpoint_rows.append([puzzle, ENDPOINT_CODES[category], forcing])
    endpoints_manifest = manifest("expanded_cover_endpoints.jsonl", endpoint_rows,
        ["puzzle", "discovery_category", "forcing_outcome"],
        {"discovery_category": {value: key for key, value in ENDPOINT_CODES.items()},
         "forcing_outcome": {0: "not-attempted", 1: "target-placed"}})
    result = {"experiment": "Expanded counterexample-guided clue-cover search",
              "baseline_commit": "0532d63d9241a69c1fb57e9c34a1533effe1981d",
              "runs": runs, "distinct_logged_raw_clue_sets_excluding_pilot": len(all_puzzles),
              "distinct_logged_unique_raw_clue_sets": len(unique_raw),
              "distinct_evaluated_unique_endpoints": len(endpoints),
              "distinct_endpoint_category_counts": dict(endpoint_categories),
              "public_manifests": {"trials": trials_manifest, "endpoints": endpoints_manifest},
              "category_conflicts": conflicts,
              "distinct_nonforcing_gap_puzzles": [{"puzzle": puzzle, "occurrences": occurrences}
                                                   for puzzle, occurrences in gaps.items()],
              "distinct_gap_forcing_closure_results": dict(Counter(
                  next((entry["forcing_result"] for entry in occurrences
                        if entry["forcing_result"] is not None), "not-requalified")
                  for occurrences in gaps.values())),
              "full_stall_status": "No full older-dispatcher stall found in these bounded runs.",
              "counting_notes": [
                  "Raw MILP clue sets and uniqueness-minimized endpoints are different populations.",
                  "Resumed alternate-solution cuts and requalified positive fixtures do not constitute newly generated puzzles.",
                  "Distinct counts deduplicate exact givens strings, not Sudoku symmetries.",
                  "Interrupted summaries may contain one started candidate not yet present in the fully logged counts.",
                  "The preliminary a pilot used a more restrictive debugging no-good rule and is excluded from aggregate evidence.",
                  "Per-MILP wall-time bounds make discovery trajectories hardware dependent; fixed candidate correctness checks are separate.",
                  "Floating LP/MILP feasibility is a discovery screen. Only separately exact-certified fixtures support a native-LP separation claim."
              ]}
    destination = HERE / "expanded_cover_summary.json"
    destination.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({key: value for key, value in result.items()
                      if key not in ("runs", "distinct_nonforcing_gap_puzzles", "counting_notes")}, indent=2))


if __name__ == "__main__":
    main()
