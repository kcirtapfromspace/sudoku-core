"""Bounded discovery across different completed grids and clue covers.

Optional NumPy/SciPy/NetworkX are discovery dependencies. No floating solver
status or old-engine trace is a proof certificate. Positive candidates require
the independent exact acceptance workflow before production use.
"""
import os
for _variable in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_variable] = "1"
import argparse
from collections import Counter
import hashlib
import itertools
import json
from pathlib import Path
import random
import subprocess
import time

import numpy as np
from scipy.optimize import linprog
import networkx as nx

from structural_discover import HOUSES, PEERS, four_colorable
from real_four_digit_mine import initial_domains, minimal_covers
from real_verify import close, exact_completions

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent


def matrices():
    native = np.zeros((324, 729), dtype=np.int8)
    role = np.zeros((27, 81), dtype=np.int8)
    for h, cells in enumerate(HOUSES):
        for c in cells:
            role[h, c] = 1
    for c in range(81):
        r, k = divmod(c, 9)
        b = r // 3 * 3 + k // 3
        for d in range(9):
            native[[c, 81 + r * 9 + d, 162 + k * 9 + d, 243 + b * 9 + d], c * 9 + d] = 1
    return native, role


def native_lp(domains, native):
    if any(not mask for mask in domains):
        return False
    variables = [9*c+d for c, mask in enumerate(domains) for d in range(9) if mask & (1 << d)]
    result = linprog(np.zeros(len(variables)), A_eq=native[:, variables],
                     b_eq=np.ones(324), bounds=(0, None), method="highs-ds")
    if not result.success and result.status != 2:
        raise RuntimeError(result.message)
    return bool(result.success)


def state(domains):
    return {"version": 1, "domains": domains,
            "values": [m.bit_length() if m & (m-1) == 0 else 0 for m in domains]}


class Discovery:
    def __init__(self, args):
        self.args = args
        self.started = time.monotonic()
        self.native, self.role = matrices()
        self.stats = Counter()
        self.motifs = []
        self.families = []
        self.candidates = []
        self.current_completion = None
        self.output = {"schema_version": 1, "seed": args.seed,
                       "baseline_commit": "0532d63", "scope": "Bounded discovery only; floating LP is a screen, not exact acceptance",
                       "arguments": json.loads(json.dumps(vars(args), default=str))}
        args.raw_dir.mkdir(parents=True, exist_ok=True)
        args.output.parent.mkdir(parents=True, exist_ok=True)

    def checkpoint(self):
        self.output.update(seconds=time.monotonic()-self.started, stats=dict(self.stats),
                           motifs=self.motifs, candidates=self.candidates)
        self.args.output.write_text(json.dumps(self.output, indent=2)+"\n")

    def artifact(self, suffix):
        return self.args.raw_dir / (self.args.output.stem+"_"+suffix+".json")

    def expired(self):
        return time.monotonic()-self.started >= self.args.seconds

    def old_engine(self, puzzle, domains, unique=False, full=False):
        request = {"id": "expanded-motif", "puzzle": puzzle, "state": state(domains),
                   "unique": unique, "full": full}
        self.stats["old_engine_calls"] += 1
        try:
            process = subprocess.run([self.args.binary], input=json.dumps(request)+"\n",
                                     text=True, capture_output=True, timeout=30, check=True)
            answer = json.loads(process.stdout)
        except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as error:
            self.stats["old_engine_execution_failures"] += 1
            if self.stats["old_engine_execution_failures"] <= 4:
                self.artifact("baseline_error_"+str(self.stats["old_engine_execution_failures"])).write_text(json.dumps({"request":request,"completion":self.current_completion,"stderr":str(getattr(error,"stderr",error))},indent=2)+"\n")
            return None
        invalid = [step for step in answer["trace"] if self.current_completion and
                   ((step[3] and step[2] != [int(self.current_completion[step[1]])]) or
                    (not step[3] and int(self.current_completion[step[1]]) in step[2]))]
        if invalid:
            self.stats["old_engine_invalid_steps"] += 1
            self.artifact("invalid_steps_"+str(self.stats["old_engine_invalid_steps"])).write_text(json.dumps({"request":request,"completion":self.current_completion,"answer":answer,"invalid_steps":invalid},indent=2)+"\n")
            return None
        if self.current_completion and not all(m & (1 << (int(d)-1)) for m,d in zip(answer["state"]["domains"],self.current_completion)):
            self.stats["old_engine_removed_completion"] += 1
            return None
        return answer

    def fixed_point(self, puzzle, domains, guardian, extra, unique=False):
        # state() represents singleton domains as placed values. Their peer
        # exclusions must already be propagated before the native old handlers,
        # which assume normalized placed-cell masks, are invoked.
        try:
            domains, trace = close(domains)
        except ValueError:
            self.stats["base_local_contradiction"] += 1
            return None
        if domains[guardian] == 1 << (extra-1):
            self.stats["local_forced_guardian"] += 1
            return None
        for _ in range(15):
            answer = self.old_engine(puzzle, domains, unique)
            if answer is None:
                return None
            domains = answer["state"]["domains"]
            trace.append({"kind": "existing-engine", "steps": answer["trace"], "state_after": answer["state"]})
            if domains[guardian] == 1 << (extra-1):
                self.stats["old_engine_forced_guardian"] += 1
                return None
            try:
                domains, local = close(domains)
            except ValueError:
                self.stats["base_local_contradiction"] += 1
                return None
            trace += local
            if not local:
                return domains, trace, answer
        self.stats["fixedpoint_iteration_budget"] += 1
        return None

    def evaluate(self, motif, givens, augment_seed=None, test_full=False):
        self.stats["clue_presentations_tested"] += 1
        guardian, extra = motif["guardian"], motif["extra"]
        domains = initial_domains(givens)
        assumed = domains.copy()
        assumed[guardian] &= ~(1 << (extra-1))
        if not native_lp(assumed, self.native):
            self.stats["raw_clue_lp_rejected"] += 1
            return None
        self.stats["raw_clue_lp_survivors"] += 1
        puzzle = "".join(str(d) if d else "." for d in givens)
        fixed = self.fixed_point(puzzle, domains, guardian, extra)
        if fixed is None:
            return None
        base, trace, answer = fixed
        assumed = base.copy()
        assumed[guardian] &= ~(1 << (extra-1))
        try:
            assumed, assumption_trace = close(assumed)
        except ValueError:
            self.stats["assumed_local_rejected"] += 1
            return None
        if not native_lp(assumed, self.native):
            self.stats["closed_lp_rejected"] += 1
            return None
        self.stats["closed_lp_survivors"] += 1
        completions = exact_completions(puzzle)
        unique = len(completions) == 1
        if unique:
            self.stats["unique_candidates"] += 1
            fixed = self.fixed_point(puzzle, base, guardian, extra, True)
            if fixed is None:
                self.stats["unique_old_structural_rejected"] += 1
                return None
            base, more_trace, answer = fixed
            trace += more_trace
            assumed = base.copy()
            assumed[guardian] &= ~(1 << (extra-1))
            try:
                assumed, assumption_trace = close(assumed)
            except ValueError:
                self.stats["unique_assumed_local_rejected"] += 1
                return None
            if not native_lp(assumed, self.native):
                self.stats["unique_closed_lp_rejected"] += 1
                return None
            # Full dispatcher includes forcing; this is deliberately distinct
            # from the frozen nonforcing closure above.
            answer = self.old_engine(puzzle, base, True, True)
            if answer is None:
                return None
            if answer["full_current_finding"] is not None:
                self.stats["unique_full_dispatcher_finding"] += 1
            else:
                self.stats["unique_full_dispatcher_stall"] += 1
        result = {"schema_version": 1, "name": "expanded-motif-discovery-only",
                  "puzzle": puzzle, "completion": motif["completion"], "solution": motif["completion"],
                  "initial_domains": domains, "base_domains": base, "base_trace": trace,
                  "assumed_domains": assumed, "assumption_trace": assumption_trace,
                  "state": state(base), "assumed_state": state(assumed),
                  "graph_cells": motif["support"], "palette": motif["palette"],
                  "target": {"cell": guardian, "digit": extra, "value": False},
                  "solutions_capped_at_2": len(completions), "motif_id": motif["id"],
                  "augment_seed": augment_seed, "givens_count": sum(bool(d) for d in givens),
                  "baseline": answer, "status": "Discovery only; independent exact acceptance is required"}
        self.stats["accepted_discovery_states"] += 1
        if unique:
            path = self.artifact("unique_"+str(self.stats["unique_candidates"]))
            path.write_text(json.dumps(result, indent=2)+"\n")
            print("UNIQUE DISCOVERY", path, "full finding", answer["full_current_finding"], flush=True)
        return result

    def clue_search(self, motif):
        completion = list(map(int, motif["completion"]))
        self.current_completion = motif["completion"]
        choices = [minimal_covers(completion, motif["support"], motif["guardian"], d)
                   for d in range(1, 10) if d not in motif["palette"]]
        motif["minimal_cover_counts"] = list(map(len, choices))
        total = 1
        for group in choices:
            total *= len(group)
        motif["clue_cover_presentations"] = total
        if not total:
            self.stats["motif_no_clue_cover"] += 1
            return
        if not self.args.covers:
            return
        rng = random.Random(self.args.seed+10000+motif["id"])
        if total <= 100000:
            presentations = list(itertools.product(*choices))
            rng.shuffle(presentations)
        else:
            # Bound memory on new completions with unusually many covers.
            presentations = []
            for index in rng.sample(range(total), min(self.args.covers,total)):
                groups = []
                for group in choices:
                    index, remainder = divmod(index,len(group))
                    groups.append(group[remainder])
                presentations.append(tuple(groups))
        best = None
        for groups in presentations[:self.args.covers]:
            if self.expired():
                break
            selected = set().union(*map(set, groups))
            givens = [d if c in selected else 0 for c, d in enumerate(completion)]
            result = self.evaluate(motif, givens)
            if result is None:
                continue
            if best is None or result["givens_count"] > best["givens_count"]:
                best = result
                self.artifact("candidate_"+str(motif["id"])).write_text(json.dumps(best, indent=2)+"\n")
            self.candidates.append({k: result[k] for k in ("puzzle", "givens_count", "solutions_capped_at_2", "motif_id")})
            print("SURVIVING COVER", motif["id"], "givens", result["givens_count"], "solutions", result["solutions_capped_at_2"], flush=True)
            # Diversify early covers rather than spending the entire run on
            # one embedding. Augmentation is repeated for bounded shuffled orders.
            for seed in range(self.args.augment_seeds):
                if self.expired():
                    break
                current = givens.copy()
                available = [c for c in range(81) if c != motif["guardian"] and not current[c]
                             and (self.args.allow_core_clues or c not in motif["support"])]
                rng = random.Random(self.args.seed+100000*motif["id"]+seed)
                rng.shuffle(available)
                if seed % 3 == 0:
                    available.sort(key=lambda c: completion[c] in motif["palette"])
                elif seed % 3 == 1:
                    available.sort(key=lambda c: completion[c] not in motif["palette"])
                for c in available:
                    if self.expired():
                        break
                    augmented = current.copy()
                    augmented[c] = completion[c]
                    trial = self.evaluate(motif, augmented, seed)
                    if trial is None:
                        continue
                    current = augmented
                    if best is None or trial["givens_count"] > best["givens_count"]:
                        best = trial
                        self.artifact("best").write_text(json.dumps(best, indent=2)+"\n")
                    if trial["solutions_capped_at_2"] == 1:
                        self.candidates.append({k: trial[k] for k in ("puzzle", "givens_count", "solutions_capped_at_2", "motif_id")})
                        return
            break
        if best is not None:
            motif["best_givens_count"] = best["givens_count"]
            motif["best_solutions_capped_at_2"] = best["solutions_capped_at_2"]
            self.artifact("candidate_"+str(motif["id"])).write_text(json.dumps(best, indent=2)+"\n")

    def run(self):
        if self.args.motif_input:
            for path in self.args.motif_input:
                source = json.loads(path.read_text())
                for item in source["motifs"]:
                    if self.expired():
                        break
                    motif = dict(item)
                    motif["parent_motif_source"] = str(path)
                    motif["parent_motif_id"] = motif["id"]
                    motif["id"] = len(self.motifs)
                    self.motifs.append(motif)
                    print("RECHECK MOTIF", motif["id"], motif["source_id"], len(motif["support"]), flush=True)
                    self.clue_search(motif)
                    self.checkpoint()
            self.checkpoint()
            print(json.dumps({"stats": dict(self.stats), "motifs":len(self.motifs), "seconds":time.monotonic()-self.started}), flush=True)
            return
        corpus = json.loads((ROOT / "research/nonmodular/real_game_corpus.json").read_text())
        selected = [corpus[i] for i in self.args.indices]
        original = json.loads((HERE / "structural_four_digit_13.json").read_text())
        source_graph = nx.Graph()
        source_graph.add_nodes_from(original["support"])
        source_graph.add_edges_from((a,b) for a,b in itertools.combinations(original["support"],2) if b in PEERS[a])
        for source in selected:
            if self.expired():
                break
            solution = exact_completions(source["puzzle"], 1)[0]
            completion = list(map(int, solution))
            self.stats["distinct_completed_grids"] += 1
            print("GRID", source["id"], solution, flush=True)
            rng = random.Random(self.args.seed+self.stats["distinct_completed_grids"])
            palettes = list(itertools.combinations(range(1,10),4))
            rng.shuffle(palettes)
            retained = 0
            local_families = []
            for palette in palettes:
                base = {c for c,d in enumerate(completion) if d in palette}
                guardians = [c for c,d in enumerate(completion) if d not in palette]
                rng.shuffle(guardians)
                for guardian in guardians[:self.args.guardians]:
                    choices = [tuple(sorted(base & h)) for h in HOUSES if guardian in h]
                    removals = list(itertools.product(*choices)); rng.shuffle(removals)
                    for removed in removals:
                        if self.expired():
                            break
                        self.stats["native_support_presentations"] += 1
                        vertices = base-set(removed)|{guardian}
                        try:
                            colorable, _ = four_colorable(vertices)
                        except RuntimeError:
                            self.stats["coloring_budget_rejected"] += 1
                            continue
                        if colorable:
                            continue
                        self.stats["non_four_colorable_supports"] += 1
                        role = linprog(np.zeros(81), A_eq=self.role, b_eq=np.full(27,4),
                                       bounds=[(1,1) if c in vertices else (0,1) for c in range(81)], method="highs-ds")
                        if not role.success:
                            self.stats["synthetic_role_lp_rejected"] += 1
                            continue
                        for c in sorted(vertices):
                            if not four_colorable(vertices-{c})[0]:
                                vertices.remove(c)
                        graph = nx.Graph()
                        graph.add_nodes_from(vertices)
                        graph.add_edges_from((a,b) for a,b in itertools.combinations(vertices,2) if b in PEERS[a])
                        signature = nx.weisfeiler_lehman_graph_hash(graph)
                        if any(nx.is_isomorphic(graph, g) for g in local_families):
                            self.stats["duplicate_graph_motifs_skipped"] += 1
                            continue
                        local_families.append(graph)
                        family = next((i for i,g in enumerate(self.families) if nx.is_isomorphic(graph,g)),len(self.families))
                        if family == len(self.families):
                            self.families.append(graph)
                        motif = {"id": len(self.motifs), "source_id": source["id"], "source_puzzle": source["puzzle"],
                                 "completion": solution, "palette": list(palette), "guardian": guardian,
                                 "extra": completion[guardian], "support": sorted(vertices),
                                 "edges": graph.number_of_edges(), "degrees": sorted(dict(graph.degree()).values()),
                                 "weisfeiler_lehman_hash": signature,
                                 "isomorphic_to_original13": nx.is_isomorphic(graph, source_graph),
                                 "graph_family": family,
                                 "removed_completion_cells": list(removed)}
                        self.motifs.append(motif)
                        retained += 1
                        print("NEW MOTIF", motif["id"], "n", len(vertices), "m", graph.number_of_edges(), "original13", motif["isomorphic_to_original13"], flush=True)
                        self.checkpoint()
                        self.clue_search(motif)
                        self.checkpoint()
                        if retained >= self.args.motifs_per_grid:
                            break
                    if retained >= self.args.motifs_per_grid or self.expired():
                        break
                if retained >= self.args.motifs_per_grid or self.expired():
                    break
            self.checkpoint()
        self.checkpoint()
        print(json.dumps({"stats": dict(self.stats), "motifs": len(self.motifs), "seconds":time.monotonic()-self.started}), flush=True)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--indices", type=lambda s: [int(x)-1 for x in s.split(",")], default=[0,1,2,3,4,5,6,7])
    p.add_argument("--guardians", type=int, default=6)
    p.add_argument("--motifs-per-grid", type=int, default=2)
    p.add_argument("--covers", type=int, default=80)
    p.add_argument("--augment-seeds", type=int, default=2)
    p.add_argument("--allow-core-clues", action="store_true")
    p.add_argument("--motif-input", type=lambda s: [Path(x) for x in s.split(",")])
    p.add_argument("--seconds", type=int, default=900)
    p.add_argument("--seed", type=int, default=2198341)
    p.add_argument("--binary", type=Path, default=Path("/tmp/global-real-engine-probe/target/release/examples/real_existing_probe"))
    p.add_argument("--output", type=Path, default=HERE / "expanded_motif_results.json")
    p.add_argument("--raw-dir", type=Path, default=ROOT / "target/expanded-motif-raw")
    Discovery(p.parse_args()).run()


if __name__ == "__main__":
    main()
