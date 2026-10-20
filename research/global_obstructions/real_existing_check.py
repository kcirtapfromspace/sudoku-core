"""Reproduce old-method comparisons against the frozen pre-global engine.

The historical source is extracted into a temporary directory. Neither current
production sources nor the recorded logical premises are modified.
"""
import io
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tarfile
import tempfile

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
BASELINE = "0532d63"
sys.path.insert(0, str(HERE.parent / "nonmodular"))
import real_game_probe

def build_engine(temporary):
    folder = Path(temporary)
    baseline, instrumented = folder / "baseline", folder / "instrumented"
    baseline.mkdir()
    instrumented.mkdir()
    archive = subprocess.check_output(
        ["git", "archive", BASELINE, "src", "Cargo.toml", "Cargo.lock"], cwd=ROOT)
    with tarfile.open(fileobj=io.BytesIO(archive)) as files:
        for member in files.getmembers():
            destination = (baseline / member.name).resolve()
            if baseline.resolve() not in destination.parents:
                raise ValueError("Unexpected archive path")
            if member.isdir():
                destination.mkdir(parents=True, exist_ok=True)
            elif member.isfile():
                destination.parent.mkdir(parents=True, exist_ok=True)
                destination.write_bytes(files.extractfile(member).read())
            else:
                raise ValueError("Unexpected archive member type")
    real_game_probe.ROOT = baseline
    real_game_probe.build_probe(instrumented)
    source = instrumented / "src/solver/mod.rs"
    source.write_text(source.read_text() + """
impl Solver {
    pub fn research_full_finding(&self, grid: &Grid)
        -> Option<(Technique,usize,Vec<u8>,bool)> {
        self.find_first_technique(grid).map(|f| match f.inference {
            InferenceResult::Placement{cell,value} => (f.technique,cell,vec![value],true),
            InferenceResult::Elimination{cell,values} => (f.technique,cell,values,false),
        })
    }
}
""")
    shutil.copy(HERE / "real_existing_probe.rs",
                instrumented / "examples/real_existing_probe.rs")
    build = subprocess.run(["cargo", "build", "--release", "--example", "real_existing_probe"],
                           cwd=instrumented, capture_output=True, text=True)
    if build.returncode:
        raise SystemExit(build.stderr)
    binary = instrumented / "target/release/examples/real_existing_probe"
    return binary


def compare(binary, name):
    case = json.loads((HERE / ("real_" + name + ".json")).read_text())
    unique = case.get("solutions_capped_at_2", 1) == 1
    request = {"id": name, "puzzle": case["puzzle"], "state": case["state"],
               "full": unique, "unique": unique}
    output = subprocess.check_output([binary], input=json.dumps(request) + "\n",
                                     text=True, timeout=60)
    result = json.loads(output)
    result["baseline_commit"] = subprocess.check_output(
        ["git", "rev-parse", BASELINE], cwd=ROOT, text=True).strip()
    (HERE / ("real_" + name + "_existing.json")).write_text(
        json.dumps(result, indent=2) + "\n")
    print(name, "old structural steps", result["trace"],
          "arithmetic", result["arithmetic"]["hint"],
          "full dispatcher", result["full_current_finding"], flush=True)

if __name__ == "__main__":
    with tempfile.TemporaryDirectory(prefix="global-real-old-engine-") as temporary:
        binary = build_engine(temporary)
        for name in ("4662", "4401", "Loki", "four_digit_candidate", "four_digit_unique"):
            compare(binary, name)
