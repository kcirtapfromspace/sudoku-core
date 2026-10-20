"""Build the public-API fixture check in a temporary copy, leaving src untouched."""
import json
from pathlib import Path
import subprocess
import tempfile
from real_game_probe import build_probe

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
with tempfile.TemporaryDirectory(prefix="sudoku-real-game-api-") as temp:
    binary = build_probe(temp, "real_game_engine_probe")
    result = subprocess.check_output([binary, HERE / "real_game_cases.json",
                                      HERE / "three_house_example_replay.json"], text=True)
    json.loads(result)
    (HERE / "real_game_engine_checks.json").write_text(result)
    print(result)
