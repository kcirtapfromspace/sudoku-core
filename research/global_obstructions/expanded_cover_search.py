"""Counterexample-guided clue-set search with continuous-LP side constraints.

This is a bounded discovery experiment, not a proof of historical originality or
nonexistence. Every alternate complete Sudoku contributes a necessary hitting
set constraint for uniqueness. A mixed-integer optimizer changes whole clue
sets, including earlier removals/swaps, rather than only adding clues greedily.
One continuous native Sudoku LP point must coexist with the selected givens and
the guardian-false four-palette core. Positive discovery outputs still require
an independent exact rational witness before publication.
"""
import os
for _name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
              "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[_name] = "1"

import argparse
import hashlib
import json
from pathlib import Path
import random
import shutil
import subprocess
import sys
import time
import warnings

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, linprog, milp
from scipy.sparse import csc_matrix, lil_matrix

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
sys.path.insert(0, str(HERE))
from real_existing_check import build_engine
from real_four_digit_mine import initial_domains
from real_verify import close, exact_completions, graph_colorings, sees

BASELINE = "0532d63d9241a69c1fb57e9c34a1533effe1981d"


def candidate_state(domains):
    return {"version": 1, "domains": domains,
            "values": [mask.bit_length() if mask & (mask - 1) == 0 else 0
                       for mask in domains]}


def cached_engine():
    cache = ROOT / "target/expanded-cover-old-engine"
    binary = cache / "real_existing_probe"
    if not binary.exists():
        source = Path("/tmp/global-real-engine-probe/target/release/examples/real_existing_probe")
        cache.mkdir(parents=True, exist_ok=True)
        if source.exists():
            shutil.copy2(source, binary)
        else:
            build_folder = cache / "build"
            build_folder.mkdir(exist_ok=True)
            shutil.copy2(build_engine(build_folder), binary)
    return binary


class Engine:
    def __init__(self, binary):
        self.process = subprocess.Popen([binary], stdin=subprocess.PIPE,
                                        stdout=subprocess.PIPE, text=True, bufsize=1)
        self.calls = 0

    def run(self, puzzle, domains, unique=False, full=False):
        self.calls += 1
        # Conversion of a singleton mask to a placement is safe only after all
        # native peers have lost that digit. Callers close Hall/three-house GAC
        # first; audit this before exposing a state to CandidateFabric.
        for cell, mask in enumerate(domains):
            assert 0 < mask <= 511
            if mask & (mask - 1) == 0:
                assert all(not domains[peer] & mask for peer in range(81)
                           if peer != cell and sees(cell, peer)), (cell, mask)
        request = {"id": "expanded-cover", "puzzle": puzzle,
                   "state": candidate_state(domains), "unique": unique, "full": full}
        self.process.stdin.write(json.dumps(request) + "\n")
        self.process.stdin.flush()
        line = self.process.stdout.readline()
        if not line:
            raise RuntimeError("Historical baseline probe terminated")
        return json.loads(line)

    def close(self):
        self.process.stdin.close()
        self.process.wait(timeout=10)


class Enumerator:
    def __init__(self):
        binary = ROOT / "target/expanded-cover-solutions"
        source = HERE / "expanded_cover_solutions.cpp"
        if not binary.exists() or binary.stat().st_mtime < source.stat().st_mtime:
            subprocess.run(["c++", "-O3", "-std=c++17", source, "-o", binary], check=True)
        self.process = subprocess.Popen([binary], stdin=subprocess.PIPE,
                                        stdout=subprocess.PIPE, text=True, bufsize=1)

    def run(self, puzzle, limit=8):
        self.process.stdin.write(puzzle + " " + str(limit) + " 1000000\n")
        self.process.stdin.flush()
        return json.loads(self.process.stdout.readline())

    def close(self):
        self.process.stdin.close()
        self.process.wait(timeout=10)


class Search:
    def __init__(self, args):
        self.args, self.started = args, time.monotonic()
        self.script_sha256 = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
        shutil.copy2(__file__, args.output_dir / "discovery-script-snapshot.py")
        shutil.copy2(HERE / "expanded_cover_solutions.cpp", args.output_dir / "enumerator-snapshot.cpp")
        self.rng = random.Random(args.seed)
        self.source = json.loads((HERE / "real_four_digit_candidate.json").read_text())
        self.solution = self.source["completion"]
        self.core = self.source["graph_cells"]
        self.guardian = self.source["target"]["cell"]
        self.extra = self.source["target"]["digit"]
        self.palette = self.source["palette"]
        self.cells = [c for c in range(81) if c != self.guardian and
                      (not args.outside_core or c not in self.core)]
        self.index = {c: 729 + i for i, c in enumerate(self.cells)}
        self.n = 729 + len(self.cells)
        self.native = lil_matrix((324, 729), dtype=float)
        for cell in range(81):
            row, column = divmod(cell, 9)
            box = row // 3 * 3 + column // 3
            for digit in range(9):
                self.native[[cell, 81 + row * 9 + digit,
                             162 + column * 9 + digit, 243 + box * 9 + digit],
                            cell * 9 + digit] = 1
        self.native = self.native.tocsc()
        self.rows, self.lower, self.upper = [], [], []
        for equation in range(324):
            self.add(dict(zip(self.native.getrow(equation).indices,
                              self.native.getrow(equation).data)), 1, 1)
        for cell in self.cells:
            self.add({cell * 9 + int(self.solution[cell]) - 1: -1,
                      self.index[cell]: 1}, -np.inf, 0)
        # Every forbidden core digit must be blocked by an actual selected clue.
        # A selected core cell itself also restricts its domain to its clue value.
        for cell in self.core:
            for digit in range(1, 10):
                if digit in self.palette or (cell == self.guardian and digit == self.extra):
                    continue
                blockers = {self.index[p]: 1 for p in self.cells
                            if (p == cell and int(self.solution[p]) != digit) or
                            (p != cell and sees(cell, p) and int(self.solution[p]) == digit)}
                self.add(blockers, 1, np.inf)
        upper = np.ones(self.n)
        for cell in self.core:
            for digit in range(1, 10):
                if digit not in self.palette:
                    upper[cell * 9 + digit - 1] = 0
        self.bounds = Bounds(np.zeros(self.n), upper)
        self.integrality = np.array([0] * 729 + [1] * len(self.cells))
        self.engine = Engine(cached_engine())
        self.enumerator = Enumerator()
        self.stats = {"milp_calls": 0, "milp_timeouts": 0, "milp_infeasible": 0,
                      "candidate_clue_sets": 0, "unique_clue_sets": 0,
                      "alternate_solutions_seen": 0, "distinct_alternate_cuts": 0,
                      "local_filter_rejections": 0, "old_target_placement_rejections": 0,
                      "postclosure_lp_rejections": 0, "nonforcing_gaps": 0,
                      "arithmetic_gaps": 0, "full_dispatcher_stalls": 0}
        self.stats.update(enumeration_nodes=0, enumeration_cutoffs=0, resumed_puzzles=0)
        self.stats.update(minimization_trials=0, minimized_unique_endpoints=0,
                          duplicate_minimized_endpoints=0)
        self.stats.update(full_closure_attempts=0, full_closure_target_placed=0,
                          full_closure_local_failures=0, full_closure_lp_failures=0,
                          full_closure_step_limits=0)
        self.cuts, self.seen = set(), set()
        self.minimized_seen = set()
        self.events = (args.output_dir / "events.jsonl").open("w")
        self.verify_cache()
        self.minimum_clues_row = len(self.rows)
        self.add({self.index[cell]: 1 for cell in self.cells}, args.min_clues, np.inf)
        for directory in args.resume:
            for line in (directory / "events.jsonl").read_text().splitlines():
                event = json.loads(line)
                if "puzzle" not in event:
                    continue
                completions = self.enumerator.run(event["puzzle"], args.completions)
                self.learn(completions["solutions"])
                self.stats["resumed_puzzles"] += 1
        for index, path in enumerate(args.qualify):
            case = json.loads(path.read_text())
            assert case["completion"] == self.solution
            counted = self.enumerator.run(case["puzzle"], 2)
            assert not counted["cutoff"] and counted["solutions"] == [self.solution]
            category, record, error = self.evaluate_unique(
                case["puzzle"], case["initial_domains"], "qualified-" + str(index))
            record.update(rejection_category=category, rejection_detail=error,
                          qualified_from=str(path))
            destination = args.output_dir / (category + "-qualified-" + str(index) + ".json")
            destination.write_text(json.dumps(record, indent=2) + "\n")
            print("QUALIFIED", category, destination, flush=True)

    def add(self, row, lower, upper):
        self.rows.append(row)
        self.lower.append(lower)
        self.upper.append(upper)

    def verify_cache(self):
        for name, expected in (("4401", None), ("4662", "FinnedSwordfish")):
            case = json.loads((HERE / ("real_" + name + ".json")).read_text())
            result = self.engine.run(case["puzzle"], case["base_domains"], True, True)
            first = result["trace"][0][0] if result["trace"] else None
            assert first == expected and result["arithmetic"]["hint"] is None
            assert result["full_current_finding"][0] == "DynamicForcingChain"
        for name in ("4401", "4662", "Loki", "four_digit_candidate", "four_digit_unique"):
            case = json.loads((HERE / ("real_" + name + ".json")).read_text())
            counted = self.enumerator.run(case["puzzle"], 2)
            assert not counted["cutoff"]
            assert len(counted["solutions"]) == case["solutions_capped_at_2"]
            for completion in counted["solutions"]:
                assert all(case["puzzle"][c] == "." or case["puzzle"][c] == completion[c]
                           for c in range(81))

    def learn(self, completions):
        for completion in completions:
            if completion == self.solution:
                continue
            difference = tuple(cell for cell in self.cells
                               if completion[cell] != self.solution[cell])
            self.stats["alternate_solutions_seen"] += 1
            assert difference, "An alternate differs only at the unclued guardian"
            if difference not in self.cuts:
                self.cuts.add(difference)
                self.stats["distinct_alternate_cuts"] += 1
                self.add({self.index[cell]: 1 for cell in difference}, 1, np.inf)

    def log(self, event):
        event["seconds"] = time.monotonic() - self.started
        self.events.write(json.dumps(event) + "\n")
        self.events.flush()

    def lp(self, domains):
        variables = [cell * 9 + digit for cell, mask in enumerate(domains)
                     for digit in range(9) if mask & (1 << digit)]
        return linprog(np.zeros(len(variables)), A_eq=self.native[:, variables],
                       b_eq=np.ones(324), bounds=(0, 1), method="highs-ds").success

    def evaluate_unique(self, puzzle, initial, trial):
        record = {"schema_version": 1, "name": "expanded-cover-unique-" + str(trial),
                  "puzzle": puzzle, "solution": self.solution, "completion": self.solution,
                  "initial_domains": initial, "graph_cells": self.core,
                  "palette": self.palette, "target": self.source["target"],
                  "conclusion": self.source["conclusion"], "solutions_capped_at_2": 1,
                  "baseline_commit": BASELINE, "source_kind": "generated-clue-set-search",
                  "graph_edges": self.source["graph_edges"]}
        try:
            base, trace = close(initial)
        except ValueError:
            raise AssertionError("A subset of a respecting completion cannot be inconsistent")
        record.update(base_domains=base, base_trace=trace, state=candidate_state(base))
        try:
            branch = base.copy()
            branch[self.guardian] &= ~(1 << (self.extra - 1))
            if not branch[self.guardian]:
                raise ValueError("guardian already placed")
            assumed, at = close(branch)
        except ValueError as error:
            self.stats["local_filter_rejections"] += 1
            return "local-filter", record, str(error)
        record.update(assumed_domains=assumed, assumption_trace=at,
                      assumed_state=candidate_state(assumed))
        for _ in range(50):
            old = self.engine.run(puzzle, base, True, False)
            if old["state"]["values"][self.guardian] != 0:
                self.stats["old_target_placement_rejections"] += 1
                record["old_rejection"] = old
                return "old-target-placed", record, None
            if not old["trace"]:
                break
            trace.append({"kind": "existing-engine", "steps": old["trace"],
                          "state_after": old["state"]})
            base, more = close(old["state"]["domains"])
            trace.extend(more)
            if base[self.guardian] == 1 << (self.extra - 1):
                self.stats["old_target_placement_rejections"] += 1
                return "joint-target-placed", record, None
        else:
            raise RuntimeError("Joint closure monotone bound exceeded")
        record.update(base_domains=base, base_trace=trace, state=candidate_state(base))
        try:
            branch = base.copy()
            branch[self.guardian] &= ~(1 << (self.extra - 1))
            assumed, at = close(branch)
        except ValueError as error:
            self.stats["local_filter_rejections"] += 1
            return "post-old-local-filter", record, str(error)
        record.update(assumed_domains=assumed, assumption_trace=at,
                      assumed_state=candidate_state(assumed), old_fixedpoint=old)
        if not self.lp(assumed):
            self.stats["postclosure_lp_rejections"] += 1
            return "postclosure-lp", record, None
        assert graph_colorings(assumed, self.core)[0] == 0
        self.stats["nonforcing_gaps"] += 1
        full = self.engine.run(puzzle, base, True, True)
        record["old_full_comparison"] = full
        arithmetic_gap = full["arithmetic"]["hint"] is None
        if arithmetic_gap:
            self.stats["arithmetic_gaps"] += 1
        full_stall = full["full_current_finding"] is None
        if full_stall:
            self.stats["full_dispatcher_stalls"] += 1
        elif arithmetic_gap and self.args.full_steps:
            advanced = self.full_closure(record)
            if advanced is not None:
                record = advanced
                full_stall = True
                self.stats["full_dispatcher_stalls"] += 1
        tier = "full-stall" if full_stall else "nonforcing-arithmetic-gap" if arithmetic_gap else "nonforcing-gap"
        return tier, record, None

    def full_closure(self, record):
        self.stats["full_closure_attempts"] += 1
        base, trace = record["base_domains"].copy(), record["base_trace"].copy()
        steps = []
        for iteration in range(self.args.full_steps):
            if candidate_state(base)["values"][self.guardian] != 0:
                self.stats["full_closure_target_placed"] += 1
                record["strict_full_closure"] = {"result": "target-placed", "steps": steps}
                return None
            old = self.engine.run(record["puzzle"], base, True, True)
            if old["trace"]:
                trace.append({"kind": "existing-engine", "steps": old["trace"],
                              "state_after": old["state"]})
                base, local = close(old["state"]["domains"])
                trace.extend(local)
                steps.append({"kind": "nonforcing", "count": len(old["trace"])})
                continue
            finding = old["full_current_finding"]
            if finding is None:
                assert candidate_state(base)["values"][self.guardian] == 0
                try:
                    assumed = base.copy()
                    assumed[self.guardian] &= ~(1 << (self.extra - 1))
                    assumed, assumption_trace = close(assumed)
                except ValueError:
                    self.stats["full_closure_local_failures"] += 1
                    record["strict_full_closure"] = {"result": "local-filter", "steps": steps}
                    return None
                if not self.lp(assumed):
                    self.stats["full_closure_lp_failures"] += 1
                    record["strict_full_closure"] = {"result": "lp-filter", "steps": steps}
                    return None
                assert graph_colorings(assumed, self.core)[0] == 0
                result = dict(record)
                result.update(base_domains=base, base_trace=trace, state=candidate_state(base),
                              assumed_domains=assumed, assumption_trace=assumption_trace,
                              assumed_state=candidate_state(assumed), old_full_fixedpoint=old,
                              strict_full_closure={"result": "full-stall", "steps": steps})
                return result
            technique, cell, digits, placement = finding
            if placement:
                assert digits == [int(self.solution[cell])]
                base[cell] = 1 << (digits[0] - 1)
                for peer in range(81):
                    if peer != cell and sees(cell, peer):
                        base[peer] &= ~(1 << (digits[0] - 1))
            else:
                assert int(self.solution[cell]) not in digits
                for digit in digits:
                    base[cell] &= ~(1 << (digit - 1))
            trace.append({"kind": "existing-engine", "steps": [finding],
                          "state_after": candidate_state(base)})
            base, local = close(base)
            trace.extend(local)
            steps.append({"kind": "full-dispatcher", "finding": finding})
        self.stats["full_closure_step_limits"] += 1
        record["strict_full_closure"] = {"result": "step-limit", "steps": steps}
        return None

    def minimize(self, puzzle):
        selected = [cell for cell, digit in enumerate(puzzle) if digit != "."]
        order = selected.copy()
        self.rng.shuffle(order)
        palette_mask = sum(1 << (digit - 1) for digit in self.palette)
        for cell in order:
            proposed = "".join(self.solution[c] if c in selected and c != cell else "."
                               for c in range(81))
            domains = initial_domains([int(d) if d != "." else 0 for d in proposed])
            # Keep an actual givens-only realization of the same global core.
            if any(domains[c] & ~(palette_mask | (1 << (self.extra - 1)
                    if c == self.guardian else 0)) for c in self.core):
                continue
            self.stats["minimization_trials"] += 1
            counted = self.enumerator.run(proposed, 2)
            self.stats["enumeration_nodes"] += counted["nodes"]
            self.stats["enumeration_cutoffs"] += counted["cutoff"]
            if not counted["cutoff"] and counted["solutions"] == [self.solution]:
                selected.remove(cell)
        return "".join(self.solution[c] if c in selected else "." for c in range(81))

    def run(self):
        for trial in range(self.args.trials):
            if time.monotonic() - self.started > self.args.seconds:
                break
            matrix = lil_matrix((len(self.rows), self.n), dtype=float)
            self.lower[self.minimum_clues_row] = self.args.min_clues + trial % self.args.clue_span
            for row, coefficients in enumerate(self.rows):
                for column, value in coefficients.items():
                    matrix[row, column] = value
            objective = np.zeros(self.n)
            # Positive, independently jittered clue weights favor sparse unique puzzles.
            for cell in self.cells:
                objective[self.index[cell]] = self.rng.uniform(0.65, 1.35)
            self.stats["milp_calls"] += 1
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)
                found = milp(objective, integrality=self.integrality, bounds=self.bounds,
                             constraints=LinearConstraint(matrix.tocsc(), self.lower, self.upper),
                             options={"time_limit": self.args.milp_seconds, "threads": 1,
                                      "mip_rel_gap": 0.03})
            if found.x is None:
                key = "milp_infeasible" if found.status == 2 else "milp_timeouts"
                self.stats[key] += 1
                self.log({"trial": trial, "result": key, "status": found.status})
                if found.status == 2:
                    break
                continue
            selected = tuple(cell for cell in self.cells if found.x[self.index[cell]] > 0.5)
            # Exclude exactly this clue set, while allowing all its supersets.
            # The alternate-completion cuts, not this no-good, enforce uniqueness.
            self.add({self.index[cell]: (1 if cell in selected else -1)
                      for cell in self.cells}, -np.inf, len(selected) - 1)
            if selected in self.seen:
                raise AssertionError("No-good cut failed")
            self.seen.add(selected)
            self.stats["candidate_clue_sets"] += 1
            puzzle = "".join(self.solution[c] if c in selected else "." for c in range(81))
            counted = self.enumerator.run(puzzle, self.args.completions)
            completions = counted["solutions"]
            self.stats["enumeration_nodes"] += counted["nodes"]
            self.stats["enumeration_cutoffs"] += counted["cutoff"]
            assert all(len(s) == 81 for s in completions)
            self.learn(completions)
            category = "nonunique"
            if counted["cutoff"] and len(completions) < 2:
                category = "enumeration-cutoff"
            elif len(completions) == 1:
                assert completions[0] == self.solution
                self.stats["unique_clue_sets"] += 1
                for pass_index in range(max(1, self.args.minimize_passes)):
                    minimized = self.minimize(puzzle) if self.args.minimize_passes else puzzle
                    if minimized in self.minimized_seen:
                        self.stats["duplicate_minimized_endpoints"] += 1
                        continue
                    self.minimized_seen.add(minimized)
                    self.stats["minimized_unique_endpoints"] += 1
                    label = str(trial) + "-" + str(pass_index)
                    category, record, error = self.evaluate_unique(
                        minimized, initial_domains([int(c) if c != "." else 0
                                                   for c in minimized]), label)
                    record.update(search_trial=trial, search_seed=self.args.seed,
                                  unminimized_puzzle=puzzle,
                                  rejection_category=category, rejection_detail=error)
                    destination = self.args.output_dir / (category + "-" + label + ".json")
                    destination.write_text(json.dumps(record, indent=2) + "\n")
                    print("UNIQUE", category, "givens", sum(d != "." for d in minimized),
                          "path", destination, flush=True)
                    if category == "local-filter":
                        # Initial opposite-branch GAC failure is monotone in givens.
                        bad = [c for c in self.cells if minimized[c] != "."]
                        self.add({self.index[cell]: 1 for cell in bad}, -np.inf, len(bad)-1)
            self.log({"trial": trial, "result": category, "givens": len(selected),
                      "solutions_capped": len(completions), "puzzle": puzzle,
                      "enumeration_cutoff": counted["cutoff"],
                      "cuts": len(self.cuts), "milp_status": found.status})
            if trial % 10 == 0:
                print("progress", trial, self.stats, "seconds", round(time.monotonic()-self.started, 2), flush=True)
        self.finish()

    def finish(self):
        self.events.close()
        self.engine.close()
        self.enumerator.close()
        summary = {"baseline_commit": BASELINE, "seed": self.args.seed,
                   "script_sha256": self.script_sha256,
                   "outside_core_only": self.args.outside_core,
                   "stats": self.stats, "seconds": time.monotonic() - self.started,
                   "engine_calls": self.engine.calls,
                   "binary_sha256": hashlib.sha256(cached_engine().read_bytes()).hexdigest(),
                   "command": " ".join(sys.argv),
                   "limitations": ["Floating MILP/LP used for discovery only; exact positive certification is separate.",
                                   "No-good cuts reject exact tested sets only; time and heuristic bounds leave the search incomplete.",
                                   "Full dispatcher and nonforcing results are separate tiers."]}
        (self.args.output_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
        print(json.dumps(summary), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed", type=int, default=62841)
    parser.add_argument("--trials", type=int, default=1000)
    parser.add_argument("--seconds", type=float, default=900)
    parser.add_argument("--milp-seconds", type=float, default=2)
    parser.add_argument("--completions", type=int, default=8)
    parser.add_argument("--outside-core", action="store_true")
    parser.add_argument("--min-clues", type=int, default=28)
    parser.add_argument("--clue-span", type=int, default=6)
    parser.add_argument("--resume", type=Path, action="append", default=[])
    parser.add_argument("--minimize-passes", type=int, default=2)
    parser.add_argument("--full-steps", type=int, default=80)
    parser.add_argument("--qualify", type=Path, action="append", default=[])
    parser.add_argument("--output-dir", type=Path, default=ROOT / "target/expanded-cover-search")
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    search = Search(args)
    try:
        search.run()
    except KeyboardInterrupt:
        search.stats["interrupted"] = True
        search.finish()
