"""Independent exact verification of retained expanded-cover negative controls."""
import hashlib
import json
from pathlib import Path
from collections import Counter

from real_four_digit_mine import initial_domains
from real_verify import HOUSES, close, exact_completions, graph_colorings, replay_trace

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
if not __debug__:
    raise SystemExit("Run without -O: this checker uses assertions.")


def givens_domains(puzzle):
    return initial_domains([int(c) if c in "123456789" else 0 for c in puzzle])


def main():
    summary = json.loads((HERE / "expanded_cover_summary.json").read_text())
    manifests = {}
    for name, metadata in summary["public_manifests"].items():
        data = (HERE / metadata["file"]).read_bytes()
        assert len(data) == metadata["bytes"]
        assert hashlib.sha256(data).hexdigest() == metadata["sha256"]
        rows = [json.loads(line) for line in data.splitlines()]
        assert len(rows) == metadata["records"]
        assert len({row[0] for row in rows}) == len(rows)
        assert all(len(row) == len(metadata["array_fields"]) and
                   isinstance(row[0], str) and len(row[0]) == 81 and
                   set(row[0]) <= set(".123456789") for row in rows)
        manifests[name] = rows
    trials, endpoints = manifests["trials"], manifests["endpoints"]
    assert len(trials) == summary["distinct_logged_raw_clue_sets_excluding_pilot"]
    assert all(type(row[1]) is int and row[1] in (-1, 0, 1) for row in trials)
    assert sum(row[1] == 1 for row in trials) == summary["distinct_logged_unique_raw_clue_sets"]
    assert len(endpoints) == summary["distinct_evaluated_unique_endpoints"]
    categories = {0: "old-target-placed", 1: "joint-target-placed", 2: "local-filter",
                  3: "post-old-local-filter", 4: "nonforcing-arithmetic-gap"}
    assert all(type(row[1]) is int and row[1] in categories and
               type(row[2]) is int and row[2] in (0, 1) for row in endpoints)
    assert dict(Counter(categories[row[1]] for row in endpoints)) == summary["distinct_endpoint_category_counts"]
    assert all(row[2] == (1 if row[1] == 4 else 0) for row in endpoints)
    assert {"target-placed": sum(row[2] == 1 for row in endpoints)} == summary["distinct_gap_forcing_closure_results"]
    print("Public discovery manifests: deduplication, hashes, 767 raw/438 unique and 611 endpoint counts verified; classifications not recertified.")

    forced = json.loads((HERE / "expanded_cover_old_forcing_trace_29.json").read_text())
    source = ROOT / forced["source_fixture"]
    assert hashlib.sha256(source.read_bytes()).hexdigest() == forced["source_fixture_sha256"]
    assert exact_completions(forced["puzzle"]) == [forced["solution"]]
    assert givens_domains(forced["puzzle"]) == forced["initial_domains"]
    assert replay_trace(forced["initial_domains"], forced["trace_from_givens"],
                        list(map(int, forced["solution"]))) == forced["end_domains"]
    guardian, digit = forced["target"]["cell"], forced["target"]["digit"]
    assert forced["end_domains"][guardian] == 1 << (digit - 1)
    assert forced["result"] == "guardian-placed-by-older-methods"
    print("Unique 29-given control: five old Nishio endpoints completion-validated; full givens trace replayed.")

    rejected = json.loads((HERE / "expanded_cover_rejected_nonunique.json").read_text())
    solutions = rejected["independent_completions"]
    assert len(set(solutions)) == 2 and len(exact_completions(rejected["puzzle"], 2)) == 2
    assert givens_domains(rejected["puzzle"]) == rejected["initial_domains"]
    for solution in solutions:
        assert all({solution[c] for c in house} == set("123456789") for house in HOUSES)
        assert all(c == "." or c == solution[i] for i, c in enumerate(rejected["puzzle"]))
        assert all(mask & (1 << (int(digit) - 1))
                   for mask, digit in zip(rejected["base_domains"], solution))
        assert replay_trace(rejected["initial_domains"], rejected["base_trace"],
                            list(map(int, solution))) == rejected["base_domains"]
    assert close(rejected["base_domains"])[1] == []
    assert close(rejected["assumed_domains"])[1] == []
    assert graph_colorings(rejected["assumed_domains"], rejected["graph_cells"])[0] == 0
    print("Rejected 32-given control: two independent respecting completions verified; local fixedpoints checked.")


if __name__ == "__main__":
    main()
