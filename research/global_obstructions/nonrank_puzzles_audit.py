"""Fresh-build audit of a full older-engine fixedpoint and its givens trace.

This validates the engine comparison and mathematical premises. It does not
claim originality or independently reconstruct each older named inner proof.
"""
import os
for _key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
             "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_key] = "1"

import argparse
import hashlib
import json
from pathlib import Path
import subprocess

from nonrank_puzzles_search import normalized, old_probe, valid_step
from real_existing_check import build_engine
from real_four_digit_mine import initial_domains
from real_verify import close, exact_completions, replay_trace

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fixture", type=Path, default=HERE / "nonrank_full_old_stall.json")
    parser.add_argument("--fresh-folder", type=Path, default=ROOT / "target/nonrank-fresh-baseline")
    parser.add_argument("--binary", type=Path,
                        help="Reuse a manually built frozen probe instead of rebuilding")
    parser.add_argument("--archive-source", type=Path,
                        help="Its uninstrumented archived solver/mod.rs; required with --binary")
    parser.add_argument("--output", type=Path, default=HERE / "nonrank_puzzles_first_audit.json")
    args = parser.parse_args()
    case = json.loads(args.fixture.read_text())
    solutions = exact_completions(case["puzzle"], 2)
    assert solutions == [case["solution"]]
    solution = list(map(int, case["solution"]))
    initial = initial_domains([int(c) if c in "123456789" else 0 for c in case["puzzle"]])
    assert initial == case["initial_domains"]
    domains = replay_trace(initial, case["base_trace"], solution)
    assert domains == case["base_domains"] == case["state"]["domains"]
    assert close(domains)[1] == []
    normalized(domains)
    assert any(mask & (mask - 1) for mask in domains)
    for event in case["base_trace"]:
        if event["kind"] == "existing-engine":
            for finding in event["steps"]:
                valid_step(finding, solution)
    if "assumed_domains" in case:
        assumed = domains.copy()
        target = case["target"]
        bit = 1 << (target["digit"] - 1)
        assumed[target["cell"]] = (assumed[target["cell"]] & bit if target["value"]
                                    else assumed[target["cell"]] & ~bit)
        branch = replay_trace(assumed, case["assumption_trace"])
        assert branch == case["assumed_domains"]
        assert close(branch)[1] == []
    # Rebuild the archived pre-global source in a newly created directory.
    if args.binary is None:
        args.fresh_folder.mkdir(parents=True, exist_ok=False)
        binary = build_engine(args.fresh_folder)
        archive_source = args.fresh_folder / "baseline/src/solver/mod.rs"
    else:
        binary = args.binary.resolve()
        if args.archive_source is None:
            parser.error("--binary requires its --archive-source provenance")
        archive_source = args.archive_source.resolve()
    assert "global_coloring" not in archive_source.read_text()
    assert archive_source.read_bytes() == subprocess.check_output(
        ["git", "show", "0532d63:src/solver/mod.rs"], cwd=ROOT)
    checks = []
    for _ in range(2):
        old = old_probe(binary, case["puzzle"], domains, 120)
        assert old["trace"] == [] and old["state"] == case["state"]
        assert old["full_dispatcher_executed"] is True
        assert old["full_current_finding"] is None
        assert old["arithmetic"]["error"] is None
        assert old["three_house"]["error"] is None
        assert old["arithmetic"]["hint"] is None
        assert old["arithmetic"]["budget_exhausted"] is False
        assert old["three_house"]["hint"] is None
        assert old["three_house"]["centers_checked"] == 81
        assert old["three_house"]["exhausted"] is True
        checks.append(old)
    result = {"schema_version": 1, "fixture": str(args.fixture.resolve().relative_to(ROOT)),
              "fixture_sha256": hashlib.sha256(args.fixture.read_bytes()).hexdigest(),
              "baseline_commit": subprocess.check_output(
                  ["git", "rev-parse", "0532d63"], cwd=ROOT, text=True).strip(),
              "archived_solver_module_sha256": hashlib.sha256(archive_source.read_bytes()).hexdigest(),
              "fresh_binary_sha256": hashlib.sha256(binary.read_bytes()).hexdigest(),
              "fresh_build_this_run": args.binary is None,
              "state_sha256": hashlib.sha256(json.dumps(case["state"], sort_keys=True,
                  separators=(",", ":")).encode()).hexdigest(),
              "independent_solution_count_capped_at_2": 1,
              "exact_givens_trace_replay": True,
              "all_older_endpoints_completion_validated": True,
              "unplaced_cells": sum(mask & (mask-1) != 0 for mask in domains),
              "full_older_engine_confirmations": checks,
              "result": "live full older-engine fixedpoint confirmed twice by fresh build"}
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({key: value for key, value in result.items()
                      if key != "full_older_engine_confirmations"}))


if __name__ == "__main__":
    if not __debug__:
        raise SystemExit("Run without -O")
    main()
