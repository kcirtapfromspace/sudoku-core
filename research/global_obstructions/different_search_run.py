"""Screen explicit public list certificates on previously verified live stalls."""
import argparse
from collections import Counter
import json
from pathlib import Path
import subprocess
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from different_audit_known_patterns import tridagon
from nonrank_puzzles_search import candidate_state
from real_verify import close


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("fixtures", type=Path, nargs="+")
    p.add_argument("--binary", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    a.output.mkdir(parents=True, exist_ok=True)
    cases = [json.loads(f.read_text()) for f in a.fixtures]
    cases = [c for c in cases if c["status"] == "full-old-stall"]
    response = subprocess.run([a.binary], input="".join(json.dumps(c)+"\n" for c in cases),
                              text=True, capture_output=True, timeout=60)
    assert response.returncode == 0, response.stderr
    assert len(response.stdout.splitlines()) == len(cases)
    rows, counts = [], Counter()
    for case, line in zip(cases, response.stdout.splitlines()):
        result = json.loads(line)
        counts["states"] += 1
        counts["explicit_patterns"] += result["patterns_tested"]
        counts["explicit_pattern_cutoffs"] += result["budget_exhausted"]
        for hit in result["successes"]:
            counts["rust_certified_conclusions"] += 1
            inference = hit["conclusion"]
            domains = case["base_domains"].copy()
            bit = 1 << (inference["digit"]-1)
            c = inference["cell"]
            domains[c] = domains[c] & ~bit if inference["value"] else bit
            entry = {"source": case["name"], "conclusion": inference,
                     "cells": len(hit["graph_cells"])}
            try:
                assumed, trace = close(domains)
            except ValueError:
                counts["local-refuted"] += 1
                entry["classification"] = "local-refuted"
                rows.append(entry)
                continue
            known, tested = tridagon(assumed)
            entry["classification"] = "known-direct-tridagon" if known else "no-direct-tridagon-match"
            counts[entry["classification"]] += 1
            if not known:
                counts["retained_for_independent_audit"] += 1
                fixture = {**case, "conclusion": inference,
                           "target": {**inference, "value": not inference["value"]},
                           "graph_cells": hit["graph_cells"], "palette": hit["palette"],
                           "assumed_domains": assumed, "assumed_state": candidate_state(assumed),
                           "assumption_trace": trace,
                           "discovery_only": "Direct-known-pattern no-hit is not historical novelty or full-template acceptance."}
                filename = "candidate-"+str(counts["retained_for_independent_audit"])+".json"
                (a.output / filename).write_text(json.dumps(fixture, indent=2)+"\n")
                (a.output / (filename[:-5]+"-replay.json")).write_text(json.dumps(hit["replay"])+"\n")
                entry["path"] = str(a.output / filename)
            rows.append(entry)
    summary = {"counts": dict(counts), "records": rows,
               "command": " ".join(sys.argv), "scope": "Public Rust certificate discovery only; independent full-template and conditional known-pattern audit required."}
    (a.output / "summary.json").write_text(json.dumps(summary, indent=2)+"\n")
    print(json.dumps(dict(counts), sort_keys=True))


if __name__ == "__main__":
    main()
