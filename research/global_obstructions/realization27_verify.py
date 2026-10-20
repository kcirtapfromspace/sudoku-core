"""Exact initial-givens and selected-graph negative controls for realization27.

Checks candidate provenance and a respecting opposite coloring independently.
These are graph-relaxation witnesses, not alternate full Sudoku solutions.
"""
import itertools
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent


def peer(a, b):
    return a != b and (a // 9 == b // 9 or a % 9 == b % 9
                       or (a // 27, a % 9 // 3) == (b // 27, b % 9 // 3))


def main():
    controls = json.loads((HERE / "realization27_selected_controls.json").read_text())
    assert controls["schema_version"] == 1 and len(controls["controls"]) == 2
    source = json.loads((HERE / "nonrank_graph_native27_candidate.json").read_text())
    for control in controls["controls"]:
        case = json.loads((HERE / control["fixture"]).read_text())
        values = [int(d) if d in "123456789" else 0 for d in case["puzzle"]]
        domains = []
        for cell, value in enumerate(values):
            blocked = {v for p, v in enumerate(values) if v and peer(cell, p)}
            domains.append(1 << (value - 1) if value else
                           sum(1 << (d - 1) for d in range(1, 10) if d not in blocked))
        assert case["initial_state"] == {"version": 1, "values": values, "domains": domains}
        assert domains == case["initial_domains"] and values[41] == 0
        assert sum(bool(v) for v in values) == case["givens"]
        cells, coloring = control["cells"], control["coloring"]
        assert cells == case["graph_cells"] == source["graph_cells"]
        assert len(coloring) == len(cells) == 27
        assert control["target"] == {"cell": 41, "digit": 1, "assumed": False}
        assert all(type(d) is int and 1 <= d <= 9 and domains[c] & (1 << (d - 1))
                   for c, d in zip(cells, coloring))
        mapped = dict(zip(cells, coloring))
        assert mapped[41] != 1
        assert all(mapped[a] != mapped[b] for a, b in itertools.combinations(cells, 2) if peer(a, b))
        assert case["local_domains"][41] == 1 and case["status"] == "solved"
        assert case["strict_gap"] is False
    request = json.loads((HERE / "realization27_proof_request.json").read_text())
    replay = json.loads((HERE / "replay_native27.json").read_text())
    assert request["state"] == replay["state"] == source["state"]
    assert request["cells"] == replay["proof"]["cells"] == source["graph_cells"]
    assert {k: request[k] for k in ("cell", "digit", "value")} == {
        k: replay["proof"][k] for k in ("cell", "digit", "value")}
    print("Two ordinary-givens snapshots and exact opposite selected-graph colorings verified; synthetic production request/replay binding: PASS.")


if __name__ == "__main__":
    if not __debug__:
        raise SystemExit("Run without -O: exact checks use assertions.")
    main()
