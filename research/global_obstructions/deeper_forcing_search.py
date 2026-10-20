"""Exact finite-depth forcing under Hall/all81 crossing-house propagation.

This is established case-split reasoning, never an originality criterion.
Only all-values-covered trees are accepted; deadlines and node bounds are cutoffs.
"""
import argparse
from collections import Counter
from functools import lru_cache
import hashlib
import itertools
import json
from pathlib import Path
import sys
import time

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from real_verify import ALL, HOUSES, matching_support, scan_center, replay_trace


def count(mask):
    return bin(mask).count("1")


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def hall_subset(domains, cells, allowed=511):
    for n in range(1, len(cells)+1):
        for subset in itertools.combinations(cells, n):
            union = 0
            for c in subset:
                union |= domains[c] & allowed
            if count(union) < n:
                return {"cells": list(subset), "neighbor_union": union}
    raise AssertionError("A failed matching lacked a Hall subset")


@lru_cache(maxsize=20000)
def propagate(original):
    domains, trace = list(original), []
    for _ in range(100):
        changed = False
        for sector, cells in enumerate(HOUSES):
            support = matching_support([domains[c] for c in cells], ALL)
            if support is None:
                return {"domains": domains, "closure_trace": trace, "kind": "hall-contradiction",
                        "sector": sector, "witness": hall_subset(domains, cells)}
            edits = [{"cell": c, "before": domains[c], "after": m}
                     for c, m in zip(cells, support) if domains[c] != m]
            if edits:
                trace.append({"kind": "house-gac", "sector": sector, "edits": edits})
                for e in edits:
                    domains[e["cell"]] = e["after"]
                changed = True
        for center in range(81):
            result = scan_center(domains, center)
            if result["signatures"] == 0:
                return {"domains": domains, "closure_trace": trace,
                        "kind": "three-house-contradiction", "center": center}
            edits = [{"cell": e["cell"], "before": domains[e["cell"]],
                      "after": domains[e["cell"]] & ~e["mask"]} for e in result["eliminations"]]
            if edits:
                trace.append({"kind": "three-house-gac", "center": center, "edits": edits})
                for e in edits:
                    domains[e["cell"]] = e["after"]
                changed = True
        if not changed:
            return {"domains": domains, "closure_trace": trace, "kind": "fixedpoint"}
    raise AssertionError("Monotone closure exceeded its finite bound")


def pivots(domains):
    return sorted((c for c,m in enumerate(domains) if m & (m-1)),
                  key=lambda c: (count(domains[c]), -sum(count(domains[p])>1 for h in HOUSES if c in h for p in h), c))


class Search:
    def __init__(self, seconds, nodes, lemmas=False):
        self.deadline, self.max_nodes = time.monotonic()+seconds, nodes
        self.stats, self.cache, self.lemmas = Counter(), {}, lemmas

    def refute(self, initial, depth):
        key = (tuple(initial), depth)
        if key in self.cache:
            self.stats["cache_hits"] += 1
            return self.cache[key]
        self.stats["calls"] += 1
        if time.monotonic()>self.deadline or self.stats["calls"]>self.max_nodes:
            raise TimeoutError("Explicit time/node cutoff; no refutation claim")
        closed = propagate(tuple(initial))
        self.stats["propagation_calls"] += 1
        if closed["kind"] != "fixedpoint":
            self.stats[closed["kind"]] += 1
            self.cache[key] = closed
            return closed
        domains = closed["domains"]
        if self.lemmas:
            from different_audit_conditional import odd_cycle
            from different_audit_known_patterns import tridagon
            known = odd_cycle(domains)
            if known is None:
                hits, _ = tridagon(domains)
                known = hits[0] if hits else None
            if known:
                answer = {**closed, "kind": "known-lemma", "lemma": known}
                self.stats["known-lemma"] += 1
                self.cache[key] = answer
                return answer
        if depth == 0:
            self.cache[key] = None
            return None
        for pivot in pivots(domains):
            branches = []
            for digit in range(1, 10):
                if not domains[pivot] & (1 << (digit-1)):
                    continue
                assumed = domains.copy()
                assumed[pivot] = 1 << (digit-1)
                child = self.refute(assumed, depth-1)
                if child is None:
                    break
                branches.append({"digit": digit, "node": child})
            else:
                answer = {**closed, "kind": "split", "cell": pivot, "branches": branches}
                self.stats["complete_splits"] += 1
                self.cache[key] = answer
                return answer
        self.cache[key] = None
        return None


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("fixture", type=Path)
    p.add_argument("--depth", type=int, default=2)
    p.add_argument("--seconds", type=float, default=240)
    p.add_argument("--nodes", type=int, default=20000)
    p.add_argument("--lemmas", action="store_true")
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    fixture = json.loads(a.fixture.read_text())
    domains = fixture["assumed_domains"]
    assert propagate(tuple(domains))["kind"] == "fixedpoint"
    started = time.monotonic()
    search = Search(a.seconds,a.nodes,a.lemmas)
    try:
        proof = search.refute(domains, a.depth)
        status = "refuted" if proof else "exhaustive-no-refutation-at-depth"
    except TimeoutError as error:
        proof, status = None, "cutoff"
    result = {"schema_version": 1, "fixture": str(a.fixture),
              "fixture_sha256": hashlib.sha256(a.fixture.read_bytes()).hexdigest(),
              "closed_opposite_state_hash": digest(domains), "target": fixture["target"],
              "depth": a.depth, "propagator": "Hall GAC plus all81 crossing-three-house GAC; monotone fixedpoint",
              "known_lemmas": a.lemmas, "status": status, "stats": dict(search.stats),
              "seconds": time.monotonic()-started, "proof": proof,
              "scope": "Specified finite case-split calculus only. A generic DPLL/forcing composition is established reasoning; no originality claim."}
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps({k:v for k,v in result.items()if k!="proof"},sort_keys=True),flush=True)


if __name__ == "__main__":
    main()
