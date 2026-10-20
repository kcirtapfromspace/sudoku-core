"""Reproduce givens-derived three-house experiments without editing production Rust.

Standard-library Python independently counts completions. A temporary copy of
the current engine runs its existing detectors to a fixed point while preserving
eliminations across placements. The copy excludes forcing, backtracking, and the
arithmetic family; uniqueness is enabled only for independently unique puzzles.
No supplied solution is used to derive candidate domains.
"""

import argparse
import hashlib
import itertools
import json
from functools import lru_cache
from pathlib import Path
import shutil
import subprocess
import tempfile
import time
import urllib.request

if not __debug__:
    raise SystemExit("Run without -O: this research probe uses assertions.")

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
ALL = 511
def popcount(mask):
    return bin(mask).count("1")

ROWS = [[r * 9 + c for c in range(9)] for r in range(9)]
COLS = [[r * 9 + c for r in range(9)] for c in range(9)]
BOXES = [[(br + r) * 9 + bc + c for r in range(3) for c in range(3)]
         for br in (0, 3, 6) for bc in (0, 3, 6)]
HOUSES = ROWS + COLS + BOXES
SOURCES = {
    "top95": "https://www.norvig.com/top95.txt",
    "hardest": "https://www.norvig.com/hardest.txt",
}


def completions(puzzle, limit=2):
    """Independent MRV enumeration, using only givens and classic Sudoku rules."""
    board = [int(x) if x in "123456789" else 0 for x in puzzle]
    used = [0] * 27
    cell_houses = [[i // 9, 9 + i % 9, 18 + (i // 27) * 3 + i % 9 // 3]
                   for i in range(81)]
    for i, d in enumerate(board):
        if d:
            bit = 1 << (d - 1)
            if any(used[h] & bit for h in cell_houses[i]):
                return []
            for h in cell_houses[i]:
                used[h] |= bit
    answer = []

    def search():
        best, domain, size = -1, 0, 10
        for i, d in enumerate(board):
            if d:
                continue
            hs = cell_houses[i]
            mask = ALL & ~(used[hs[0]] | used[hs[1]] | used[hs[2]])
            n = popcount(mask)
            if n == 0:
                return
            if n < size:
                best, domain, size = i, mask, n
                if n == 1:
                    break
        if best < 0:
            answer.append("".join(map(str, board)))
            return
        while domain and len(answer) < limit:
            bit = domain & -domain
            domain ^= bit
            board[best] = bit.bit_length()
            for h in cell_houses[best]:
                used[h] |= bit
            search()
            for h in cell_houses[best]:
                used[h] ^= bit
            board[best] = 0

    search()
    return answer


def matching_support(domains, allowed):
    """All supported edges in injective assignments; independent of Rust code."""
    n = len(domains)
    order = sorted(range(n), key=lambda i: popcount(domains[i] & allowed))
    ordered = [domains[i] & allowed for i in order]

    @lru_cache(None)
    def feasible(depth, used):
        if depth == n:
            return True
        choices = ordered[depth] & ~used
        while choices:
            bit = choices & -choices
            choices ^= bit
            if feasible(depth + 1, used | bit):
                return True
        return False

    if not feasible(0, 0):
        return None
    supported, seen = [0] * n, set()

    def walk(depth, used):
        if depth == n or (depth, used) in seen:
            return
        seen.add((depth, used))
        choices = ordered[depth] & ~used
        while choices:
            bit = choices & -choices
            choices ^= bit
            if feasible(depth + 1, used | bit):
                supported[order[depth]] |= bit
                walk(depth + 1, used | bit)

    walk(0, 0)
    return supported


def scan_center(masks, center):
    r, c = divmod(center, 9)
    row, col, box = set(ROWS[r]), set(COLS[c]), set(BOXES[r // 3 * 3 + c // 3])
    atoms = [[center], sorted((row & box) - {center}),
             sorted((col & box) - {center}), sorted(row - box),
             sorted(col - box), sorted(box - row - col)]
    cache, supported, signatures = {}, {x: 0 for a in atoms for x in a}, 0

    def match(atom, mask):
        key = (atom, mask)
        if key not in cache:
            cache[key] = matching_support([masks[x] for x in atoms[atom]], mask)
        return cache[key]

    for d in range(9):
        dm = 1 << d
        if not masks[center] & dm:
            continue
        for av in itertools.combinations([k for k in range(9) if k != d], 2):
            am = sum(1 << k for k in av)
            sa, sr = match(1, am), match(3, ALL ^ (dm | am))
            if sa is None or sr is None:
                continue
            for bv in itertools.combinations([k for k in range(9)
                                               if not (dm | am) & (1 << k)], 2):
                bm = sum(1 << k for k in bv)
                sb = match(2, bm)
                sc, sp = match(4, ALL ^ (dm | bm)), match(5, ALL ^ (dm | am | bm))
                if sb is None or sc is None or sp is None:
                    continue
                signatures += 1
                for cells, support in zip(atoms, [[dm], sa, sb, sr, sc, sp]):
                    for x, values in zip(cells, support):
                        supported[x] |= values
    return {"center": center, "signatures": signatures,
            "eliminations": [{"cell": x, "mask": masks[x] & ~values}
                             for x, values in sorted(supported.items())
                             if masks[x] & ~values]}


RUST_EXAMPLE = r'''
use std::io::{self, BufRead};
use sudoku_core::{Grid, Solver, Position};
fn main() {
    let solver = Solver::new();
    for line in io::stdin().lock().lines() {
        let line = line.unwrap();
        let request: serde_json::Value = serde_json::from_str(&line).unwrap();
        let puzzle = request["puzzle"].as_str().unwrap();
        let grid = Grid::from_string(puzzle).unwrap();
        let started = std::time::Instant::now();
        let (grid, steps, fixed_point) = solver.research_close(&grid, request["unique"].as_bool().unwrap());
        let domains: Vec<u16> = (0..81).map(|i| {
            let p = Position::new(i / 9, i % 9);
            if let Some(d) = grid.get(p) {1u16 << (d-1)}
            else {grid.get_candidates(p).iter().fold(0, |m, d| m | (1u16 << (d-1)))}
        }).collect();
        let values: Vec<u8> = (0..81).map(|i|grid.get(Position::new(i/9,i%9)).unwrap_or(0)).collect();
        println!("{}", serde_json::json!({"domains":domains,"values":values,
            "steps":steps,"fixed_point":fixed_point,"complete":grid.is_complete(),
            "seconds":started.elapsed().as_secs_f64()}));
    }
}
'''


def build_probe(directory, extra_example=None):
    """Instrument a throwaway copy; the repository's src remains untouched."""
    copy = Path(directory)
    shutil.copytree(ROOT / "src", copy / "src")
    shutil.copy(ROOT / "Cargo.toml", copy / "Cargo.toml")
    shutil.copy(ROOT / "Cargo.lock", copy / "Cargo.lock")
    (copy / "examples").mkdir()
    (copy / "examples" / "basic.rs").write_text("fn main() {}\n")
    (copy / "examples" / "real_game_probe.rs").write_text(RUST_EXAMPLE)
    if extra_example:
        shutil.copy(HERE / (extra_example + ".rs"),
                    copy / "examples" / (extra_example + ".rs"))
    path = copy / "src" / "solver" / "mod.rs"
    source = path.read_text()
    start = source.index("    fn find_first_technique(&self, grid: &Grid)")
    end = source.index("    /// Solve the puzzle using human techniques", start)
    dispatch = source[start:end]
    dispatch = dispatch[:dispatch.index("        // Forcing chains need the Grid")]
    dispatch += "        None\n    }\n\n"
    dispatch = dispatch.replace("fn find_first_technique(&self, grid: &Grid)",
                                "fn research_find(&self, grid: &Grid, allow_uniqueness: bool)")
    import re
    dispatch = re.sub(r"        if let Some\(f\) = uniqueness::([^\n]+) \{\n            return Some\(f\);\n        \}",
                      r"        if allow_uniqueness { if let Some(f) = uniqueness::\1 { return Some(f); } }", dispatch)
    dispatch = re.sub(r"        if let Some\(f\) = arithmetic::find\(grid, &fab\) \{\n            return Some\(f\);\n        \}\n", "", dispatch)
    dispatch = re.sub(r"        if let Some\(f\) = three_house::find\(grid\) \{\n            return Some\(f\);\n        \}\n", "", dispatch)
    assert "three_house::" not in dispatch and "arithmetic::" not in dispatch
    closure = r'''
    pub fn research_baseline_finding(&self, grid: &Grid, unique: bool) -> Option<(Technique,usize,Vec<u8>,bool)> {
        self.research_find(grid, unique).map(|f|match f.inference {
            InferenceResult::Placement{cell,value} => (f.technique,cell,vec![value],true),
            InferenceResult::Elimination{cell,values} => (f.technique,cell,values,false),
        })
    }
    pub fn research_close(&self, grid: &Grid, unique: bool) -> (Grid, Vec<(Technique,usize,Vec<u8>,bool)>, bool) {
        let mut working = grid.deep_clone();
        working.recalculate_candidates();
        let mut trace = Vec::new();
        for _ in 0..2000 {
            let Some(f) = self.research_find(&working, unique) else {return (working,trace,true)};
            let old: [crate::BitSet;81] = std::array::from_fn(|i|working.get_candidates(idx_to_pos(i)));
            let step = match &f.inference {
                InferenceResult::Placement{cell,value} => (f.technique,*cell,vec![*value],true),
                InferenceResult::Elimination{cell,values} => (f.technique,*cell,values.clone(),false),
            };
            apply_finding(&mut working,&f);
            // Placement recalculation must not resurrect earlier eliminations.
            if step.3 {for i in 0..81 {
                let p=idx_to_pos(i);
                if working.get(p).is_none() {
                    let now=working.get_candidates(p);
                    let mut keep=crate::BitSet::empty();
                    for d in now.iter() {if old[i].contains(d) {keep.insert(d);}}
                    working.cell_mut(p).set_candidates(keep);
                }
            }}
            trace.push(step);
        }
        (working,trace,false)
    }
'''
    source = source[:start] + dispatch + closure + source[start:]
    path.write_text(source)
    example = extra_example or "real_game_probe"
    build = subprocess.run(["cargo", "build", "--release", "--example", example],
                           cwd=copy, text=True, capture_output=True)
    if build.returncode:
        raise RuntimeError(build.stderr)
    return copy / "target" / "release" / "examples" / example


def compact_results(result):
    """Keep exact positive states/traces and aggregate negative evidence in Git."""
    records = result["records"]
    positives = [record for record in records if record["three_house"]]
    cases = []
    for record in positives:
        cases.append({
            key: record[key] for key in ("id", "puzzle", "source_url",
                "source_index_1based", "source_sha256", "attribution", "solution",
                "solutions_capped_at_2", "three_house")
        })
        cases[-1]["state"] = dict(version=1, values=record["baseline"]["values"],
                                  domains=record["baseline"]["domains"])
        cases[-1]["baseline_trace"] = record["baseline"]["steps"]
    fixtures = {"schema_version": 1,
                "baseline_orchestrator_sha256": result["baseline_orchestrator_sha256"],
                "baseline_policy": result["baseline_policy"],
                "trace_note": "Each tuple is [existing technique, zero-based cell, digits, is_placement]. Detector choices may vary across runs because existing code uses randomized collections. These exact traces and states are retained.",
                "cases": cases}
    (HERE / "real_game_cases.json").write_text(json.dumps(fixtures, indent=2) + "\n")
    distinct = sum(len({(e["cell"], digit) for cluster in record["three_house"]
                        for e in cluster["eliminations"] for digit in range(1, 10)
                        if e["mask"] & (1 << (digit-1))}) for record in positives)
    summary = {
        "schema_version": 1,
        "baseline_orchestrator_sha256": result["baseline_orchestrator_sha256"],
        "baseline_policy": result["baseline_policy"],
        "puzzles_tested": len(records),
        "independently_unique": sum(r["solutions_capped_at_2"] == 1 for r in records),
        "baseline_solved": sum(r["baseline"]["complete"] for r in records),
        "baseline_stalled": sum(not r["baseline"]["complete"] for r in records),
        "positive_stalled_states": len(positives),
        "negative_stalled_states": sum(not r["baseline"]["complete"] and not r["three_house"] for r in records),
        "centers_checked": sum(r["centers_checked"] for r in records),
        "positive_centers": sum(len(r["three_house"]) for r in positives),
        "elimination_occurrences": sum(popcount(e["mask"]) for r in positives
                                       for c in r["three_house"] for e in c["eliminations"]),
        "distinct_puzzle_candidate_eliminations": distinct,
        "seconds": {"elapsed_including_build": result["elapsed_seconds"],
                    "independent_uniqueness": sum(r["uniqueness_seconds"] for r in records),
                    "existing_engine_closure": sum(r["baseline"]["seconds"] for r in records),
                    "python_three_house": sum(r["three_house_seconds"] for r in records)},
        "records": [{"id": r["id"], "unique": r["solutions_capped_at_2"] == 1,
                     "baseline_solved": r["baseline"]["complete"],
                     "baseline_steps": len(r["baseline"]["steps"]),
                     "positive_centers": len(r["three_house"])} for r in records],
    }
    (HERE / "real_game_summary.json").write_text(json.dumps(summary, indent=2) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--limit", type=int, default=106)
    parser.add_argument("--centers", type=int, default=81)
    parser.add_argument("--output", type=Path,
                        default=Path(tempfile.gettempdir()) / "sudoku-three-house-real-game-results.json")
    args = parser.parse_args()
    corpus_path = HERE / "real_game_corpus.json"
    if corpus_path.exists():
        corpus = json.loads(corpus_path.read_text())
    else:
        corpus = []
        for name, url in SOURCES.items():
            data = urllib.request.urlopen(url).read()
            for index, line in enumerate(data.decode().splitlines(), 1):
                if len(line.strip()) != 81:
                    continue
                corpus.append({"id": f"norvig-{name}-{index:03}", "source_url": url,
                               "source_index_1based": index, "puzzle": line.strip(),
                               "source_sha256": hashlib.sha256(data).hexdigest(),
                               "attribution": "Peter Norvig's author-hosted curated corpus; individual puzzle author unspecified"})
        corpus_path.write_text(json.dumps(corpus, indent=2) + "\n")
    started, records = time.perf_counter(), []
    source_digest = hashlib.sha256((ROOT / "src/solver/mod.rs").read_bytes()).hexdigest()
    with tempfile.TemporaryDirectory(prefix="sudoku-real-game-") as temp:
        binary = build_probe(temp)
        for item in corpus[:args.limit]:
            t = time.perf_counter()
            solutions = completions(item["puzzle"])
            uniqueness_seconds = time.perf_counter() - t
            assert solutions, item["id"]
            request = json.dumps({"puzzle": item["puzzle"], "unique": len(solutions) == 1})
            output = subprocess.check_output([binary], input=request + "\n", text=True)
            baseline = json.loads(output)
            assert baseline["fixed_point"], item["id"]
            solution = list(map(int, solutions[0]))
            for technique, cell, digits, placement in baseline["steps"]:
                assert (solution[cell] == digits[0] if placement else solution[cell] not in digits), (item["id"], technique, cell, digits)
            assert all(m & (1 << (d - 1)) for m, d in zip(baseline["domains"], solution))
            t = time.perf_counter()
            clusters = [scan_center(baseline["domains"], center)
                        for center in range(args.centers)] if not baseline["complete"] else []
            for cluster in clusters:
                assert cluster["signatures"] > 0
                for elim in cluster["eliminations"]:
                    assert not elim["mask"] & (1 << (solution[elim["cell"]] - 1))
            record = dict(item, solutions_capped_at_2=len(solutions), solution=solutions[0],
                          uniqueness_seconds=uniqueness_seconds, baseline=baseline,
                          three_house_seconds=time.perf_counter()-t,
                          centers_checked=len(clusters),
                          three_house=[c for c in clusters if c["eliminations"]])
            records.append(record)
            result = {"schema_version": 1, "baseline_orchestrator_sha256": source_digest,
                      "baseline_policy": "Known detectors excluding Three-House Matching, arithmetic, forcing and backtracking; independently justified uniqueness; preserve eliminations across placements",
                      "records": records, "elapsed_seconds": time.perf_counter()-started}
            args.output.write_text(json.dumps(result, indent=2) + "\n")
            count = sum(popcount(e["mask"]) for c in record["three_house"] for e in c["eliminations"])
            print(item["id"], "complete", baseline["complete"], "steps", len(baseline["steps"]),
                  "centers", len(clusters), "elimination_occurrences", count,
                  "baseline_seconds", round(baseline["seconds"], 3), flush=True)
    compact_results(result)


if __name__ == "__main__":
    main()
