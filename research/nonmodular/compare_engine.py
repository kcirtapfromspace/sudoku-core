"""Run the research comparison in a temporary source copy of this engine.

The public get_hint API recalculates candidates, so a private test probe is
needed to compare the exact recorded masks. The user's checkout is not edited.
"""
import json
from pathlib import Path
import shutil
import subprocess
import tempfile

here = Path(__file__).resolve().parent
repository = here.parents[1]
with tempfile.TemporaryDirectory(prefix="sudoku-nonmodular-") as folder:
    temporary = Path(folder)
    for name in ["Cargo.toml", "Cargo.lock"]:
        shutil.copy2(repository / name, temporary / name)
    shutil.copytree(repository / "src", temporary / "src")
    # Production unit tests include committed research fixtures at compile time.
    shutil.copytree(here, temporary / "research/nonmodular", ignore=shutil.ignore_patterns("__pycache__"))
    module = temporary / "src/solver/mod.rs"
    path_literal = json.dumps(str(here / "probe_existing_engine.rs"))
    with module.open("a") as output:
        output.write("\n#[cfg(test)]\n#[path = " + path_literal + "]\nmod nonmodular_research;\n")
    subprocess.run([
        "cargo", "test", "--manifest-path", str(temporary / "Cargo.toml"),
        "--lib", "--profile", "coverage", "--locked", "nonmodular_research", "--", "--nocapture",
    ], check=True)
