"""Compare exact synthetic q4 masks against explicitly bounded older methods.

Only a temporary source copy is instrumented. The new Global Coloring feature,
uniqueness rules, forcing chains/backtracking, and optional residue APIs are
excluded. A future older handler recognizing a step is reported, not suppressed.
"""

import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile

HERE = Path(__file__).resolve().parent
REPOSITORY = HERE.parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    options = parser.parse_args()
    digest = hashlib.sha256()
    for source in sorted((REPOSITORY / "src").rglob("*.rs")):
        digest.update(str(source.relative_to(REPOSITORY)).encode())
        digest.update(source.read_bytes())
    with tempfile.TemporaryDirectory(prefix="sudoku-structural-q4-") as folder:
        destination = Path(folder)
        for name in ("Cargo.toml", "Cargo.lock"):
            shutil.copy2(REPOSITORY / name, destination / name)
        shutil.copytree(REPOSITORY / "src", destination / "src")
        # Current production tests use relative include_str! research fixtures.
        shutil.copytree(REPOSITORY / "research", destination / "research",
                        ignore=shutil.ignore_patterns("__pycache__"))
        with (destination / "src/solver/mod.rs").open("a") as output:
            output.write("\n#[cfg(test)]\n#[path = " + json.dumps(str(HERE / "structural_engine_probe.rs"))
                         + "]\nmod structural_research_baseline;\n")
        environment = os.environ.copy()
        environment["CARGO_TARGET_DIR"] = str(REPOSITORY / "target")
        process = subprocess.run([
            "cargo", "test", "--manifest-path", str(destination / "Cargo.toml"),
            "--lib", "--profile", "coverage", "--locked",
            "structural_research_baseline", "--", "--nocapture",
        ], env=environment, capture_output=True, text=True, check=True)
        results = [json.loads(line.split("STRUCTURAL_BASELINE ", 1)[1])
                   for line in process.stdout.splitlines() if "STRUCTURAL_BASELINE " in line]
        assert len(results) == 2, process.stdout
    result = {"schema_version": 1, "source_sha256": digest.hexdigest(), "results": results,
              "excluded": ["new Global Coloring", "uniqueness", "forcing chains", "backtracking",
                           "optional partition-residue APIs", "external solvers and unimplemented named rules"]}
    encoded = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if options.output:
        options.output.write_text(encoded)
    print(encoded)


if __name__ == "__main__":
    main()
