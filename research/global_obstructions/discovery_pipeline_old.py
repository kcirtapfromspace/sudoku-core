"""Build the frozen older dispatcher with monotone hypothetical candidates.

This is an offline comparison instrument, not a production-engine patch or an
independent certificate checker. The original archive remains byte-for-byte
unchanged. Every solver-side candidate rebuild in its instrumented copy is
replaced by the intersection of the previous masks with the newly legal masks.
That includes rebuilds within assumptions, not merely the outer probe loop.
"""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess

from real_existing_check import build_engine, ROOT

BASELINE = "0532d63d9241a69c1fb57e9c34a1533effe1981d"
METHOD = """
impl Grid {
    /// Offline research: a solver assumption may only restrict candidates.
    pub fn research_restrict_candidates(&mut self) {
        let previous: Vec<_> = self.empty_positions().into_iter()
            .map(|pos| (pos, self.get_candidates(pos).as_raw())).collect();
        self.recalculate_candidates();
        for (pos, old) in previous {
            let keep = self.get_candidates(pos).as_raw() & old;
            self.cell_mut(pos).set_candidates(BitSet::from_raw(keep));
        }
    }
}
"""
TESTS = """
#[cfg(test)]
mod research_retained_candidate_tests {
    use super::*;
    #[test]
    fn research_monotone_rebuild_preserves_remote_elimination() {
        let mut grid = Grid::new_classic();
        let remote = Position::new(8, 8);
        grid.cell_mut(remote).remove_candidate(9);
        grid.set_cell_unchecked(Position::new(0, 0), Some(1));
        grid.research_restrict_candidates();
        assert!(!grid.get_candidates(remote).contains(9));
        assert!(!grid.get_candidates(Position::new(0, 1)).contains(1));
    }
    #[test]
    fn research_singles_assumption_preserves_remote_elimination() {
        let mut grid = Grid::new_classic();
        let remote = Position::new(8, 8);
        grid.cell_mut(remote).remove_candidate(9);
        let (branch, contradiction) = backtrack::propagate_singles(
            &grid, Position::new(0, 0), 1);
        assert!(!contradiction);
        assert!(!branch.get_candidates(remote).contains(9));
        assert!(!grid.get_candidates(remote).contains(9));
        assert_eq!(grid.get(Position::new(0, 0)), None);
    }
    #[test]
    fn research_full_assumption_preserves_remote_elimination() {
        let mut grid = Grid::new_classic();
        let remote = Position::new(8, 8);
        grid.cell_mut(remote).remove_candidate(9);
        let (branch, contradiction) = propagate_full(
            &grid, Position::new(0, 0), 1);
        assert!(!contradiction);
        assert!(!branch.get_candidates(remote).contains(9));
        assert_eq!(grid.get(Position::new(0, 0)), None);
    }
}
"""


def sha(data):
    return hashlib.sha256(data).hexdigest()


def checked(command, folder):
    result = subprocess.run(command, cwd=folder, text=True, capture_output=True)
    if result.returncode:
        raise RuntimeError(result.stdout[-3000:] + result.stderr[-6000:])
    return result.stdout


def build(folder):
    folder = folder.resolve()
    folder.mkdir(parents=True, exist_ok=False)
    binary = build_engine(folder)
    instrumented = folder / "instrumented"
    changes = []
    targets = sorted((instrumented / "src/solver").glob("*.rs"))
    targets.append(instrumented / "examples/real_existing_probe.rs")
    for source in targets:
        before = source.read_text()
        count = before.count(".recalculate_candidates()")
        if count:
            after = before.replace(".recalculate_candidates()", ".research_restrict_candidates()")
            source.write_text(after)
            changes.append({"path": str(source.relative_to(instrumented)),
                            "replacements": count, "before_sha256": sha(before.encode()),
                            "after_sha256": sha(after.encode())})
    grid = instrumented / "src/grid.rs"
    grid.write_text(grid.read_text() + METHOD)
    solver = instrumented / "src/solver/mod.rs"
    solver.write_text(solver.read_text() + TESTS)
    # A frozen unit test includes this fixture at compile time. Read the pinned
    # revision rather than the possibly changing current workspace copy.
    relative = "research/nonmodular/three-house-847.json"
    fixture = instrumented / relative
    fixture.parent.mkdir(parents=True, exist_ok=True)
    fixture.write_bytes(subprocess.check_output(["git", "show", BASELINE + ":" + relative], cwd=ROOT))
    checked(["cargo", "test", "--lib", "research_retained_candidate_tests"], instrumented)
    checked(["cargo", "build", "--release", "--example", "real_existing_probe"], instrumented)
    # Verify the historical source was never patched, including by a helper.
    for relative in ("src/grid.rs", "src/solver/mod.rs", "src/solver/backtrack.rs"):
        expected = subprocess.check_output(["git", "show", BASELINE + ":" + relative], cwd=ROOT)
        if (folder / "baseline" / relative).read_bytes() != expected:
            raise ValueError("The archived source changed: " + relative)
    result = {"schema_version": 1, "kind": "retained-frozen-dispatcher-v1",
              "baseline_commit": BASELINE, "binary_sha256": sha(binary.read_bytes()),
              "changes": changes, "retention_method_sha256": sha(METHOD.encode()),
              "targeted_rust_tests_passed": 3,
              "scope": "Frozen full older ladder including its original forcing, with monotone candidate rebuilds in root and hypothetical states. Its configured detector bounds and 200-round propagation cap remain. Positive endpoints must still be completion-validated; this instrument is not a proof certificate or an exhaustive known-method oracle."}
    (folder / "retained_manifest.json").write_text(json.dumps(result, indent=2) + "\n")
    return binary, result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--folder", type=Path, required=True)
    args = parser.parse_args()
    binary, result = build(args.folder)
    print(json.dumps({"binary": str(binary), **result}, sort_keys=True))


if __name__ == "__main__":
    main()
