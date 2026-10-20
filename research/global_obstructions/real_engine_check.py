"""Check recorded real states through current public APIs in a temporary copy."""
import json
from pathlib import Path
import shutil
import subprocess
import tempfile

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
NAMES = ("4662", "4662_after_existing", "4401", "Loki",
         "four_digit_candidate", "four_digit_unique")


if __name__ == "__main__":
    with tempfile.TemporaryDirectory(prefix="global-real-public-api-") as temporary:
        engine = Path(temporary)
        shutil.copytree(ROOT / "src", engine / "src")
        for name in ("Cargo.toml", "Cargo.lock"):
            shutil.copy(ROOT / name, engine / name)
        (engine / "examples").mkdir()
        shutil.copy(HERE / "real_engine_probe.rs", engine / "examples/real_engine_probe.rs")
        build = subprocess.run(["cargo", "build", "--release", "--example", "real_engine_probe"],
                               cwd=engine, text=True, capture_output=True)
        if build.returncode:
            raise SystemExit(build.stderr)
        inputs = "".join((HERE / ("real_" + name + ".json")).read_text().replace("\n", "")
                         + "\n" for name in NAMES)
        output = subprocess.check_output([engine / "target/release/examples/real_engine_probe"],
                                         input=inputs, text=True)
        results = [json.loads(line) for line in output.splitlines()]
        (HERE / "real_engine_results.json").write_text(json.dumps(results, indent=2) + "\n")
        for name, result in zip(NAMES, results):
            print(name, "explicit proof nodes", result["explicit_nodes"],
                  "automatic nodes", result["automatic_nodes"],
                  "seconds", round(result["automatic_seconds"], 6))
