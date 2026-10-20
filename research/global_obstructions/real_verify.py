"""Independent exact checks for published tridagons and generated q4 controls.

No uniqueness premise is used by Hall, three-house, or coloring reasoning.
The optional exact-cover count establishes the separate uniqueness claim.
"""
import itertools
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "nonmodular"))
from real_game_probe import ALL, HOUSES, matching_support, scan_center

if not __debug__:
    raise SystemExit("Run without -O: this exact research checker uses assertions.")


def close(domains):
    domains, trace = domains.copy(), []
    for _ in range(100):
        changed = False
        for sector, cells in enumerate(HOUSES):
            support = matching_support([domains[c] for c in cells], ALL)
            if support is None:
                raise ValueError("Hall contradiction in sector " + str(sector))
            edits = [{"cell": c, "before": domains[c], "after": mask}
                     for c, mask in zip(cells, support) if domains[c] != mask]
            if edits:
                trace.append({"kind": "house-gac", "sector": sector, "edits": edits})
                for edit in edits:
                    domains[edit["cell"]] = edit["after"]
                changed = True
        for center in range(81):
            result = scan_center(domains, center)
            if result["signatures"] == 0:
                raise ValueError("Three-house contradiction at center " + str(center))
            if result["eliminations"]:
                edits = [{"cell": e["cell"], "before": domains[e["cell"]],
                          "after": domains[e["cell"]] & ~e["mask"]}
                         for e in result["eliminations"]]
                trace.append({"kind": "three-house-gac", "center": center, "edits": edits})
                for edit in edits:
                    domains[edit["cell"]] = edit["after"]
                changed = True
        if not changed:
            return domains, trace
    raise ValueError("Closure did not terminate within its monotone bound")


def exact_completions(puzzle, limit=2):
    """Independent exact-cover enumerator; no engine or published solution used."""
    rows, columns = {}, {j: set() for j in range(324)}
    for cell in range(81):
        r, c = divmod(cell, 9)
        box = r // 3 * 3 + c // 3
        for digit in range(9):
            row = cell * 9 + digit
            rows[row] = [cell, 81 + r * 9 + digit, 162 + c * 9 + digit,
                         243 + box * 9 + digit]
            for column in rows[row]:
                columns[column].add(row)

    def select(row):
        removed = []
        for column in rows[row]:
            conflicts = columns.pop(column)
            for other_row in conflicts:
                for other_column in rows[other_row]:
                    if other_column in columns:
                        columns[other_column].remove(other_row)
            removed.append((column, conflicts))
        return removed

    def restore(removed):
        for column, conflicts in reversed(removed):
            columns[column] = conflicts
            for other_row in conflicts:
                for other_column in rows[other_row]:
                    if other_column != column and other_column in columns:
                        columns[other_column].add(other_row)

    chosen, solutions = [], []
    for cell, char in enumerate(puzzle):
        if char in "123456789":
            row = cell * 9 + int(char) - 1
            if row not in columns.get(cell, set()):
                return []
            select(row)
            chosen.append(row)

    def search():
        if not columns:
            board = [0] * 81
            for row in chosen:
                board[row // 9] = row % 9 + 1
            solutions.append("".join(map(str, board)))
            return
        column = min(columns, key=lambda key: len(columns[key]))
        for row in list(columns[column]):
            if len(solutions) >= limit:
                return
            chosen.append(row)
            removed = select(row)
            search()
            restore(removed)
            chosen.pop()

    search()
    return solutions


def sees(a, b):
    ra, ca = divmod(a, 9)
    rb, cb = divmod(b, 9)
    return ra == rb or ca == cb or (ra // 3, ca // 3) == (rb // 3, cb // 3)


def graph_colorings(domains, cells, limit=1):
    adjacency = {c: {p for p in cells if p != c and sees(c, p)} for c in cells}
    assigned, count, nodes = {}, 0, 0

    def search():
        nonlocal count, nodes
        nodes += 1
        if count >= limit:
            return
        available = {}
        for cell in cells:
            if cell not in assigned:
                blocked = 0
                for peer in adjacency[cell]:
                    blocked |= assigned.get(peer, 0)
                available[cell] = domains[cell] & ~blocked
        if not available:
            count += 1
            return
        cell = min(available, key=lambda c: (bin(available[c]).count("1"), c))
        choices = available[cell]
        while choices:
            bit = choices & -choices
            choices ^= bit
            assigned[cell] = bit
            search()
            del assigned[cell]

    search()
    return count, nodes


def triple_feasible(domains, selected):
    """Direct CSP witness search for any three native houses, not atom DP."""
    cells = sorted(set().union(*(set(HOUSES[h]) for h in selected)))
    membership = {c: tuple(i for i, h in enumerate(selected) if c in HOUSES[h])
                  for c in cells}
    assigned, used = {}, [0, 0, 0]

    def search():
        available = {}
        for cell in cells:
            if cell not in assigned:
                mask = domains[cell]
                for house in membership[cell]:
                    mask &= ~used[house]
                if mask == 0:
                    return False
                available[cell] = mask
        if not available:
            return True
        cell = min(available, key=lambda c: (bin(available[c]).count("1"), c))
        choices = available[cell]
        while choices:
            bit = choices & -choices
            choices ^= bit
            assigned[cell] = bit
            for house in membership[cell]:
                used[house] |= bit
            if search():
                return True
            for house in membership[cell]:
                used[house] ^= bit
            del assigned[cell]
        return False

    return search()


def replay_trace(initial, trace, solution=None):
    domains = initial.copy()
    for event in trace:
        if event["kind"] == "house-gac":
            cells = HOUSES[event["sector"]]
            support = matching_support([domains[c] for c in cells], ALL)
            assert support is not None
            expected = [{"cell": c, "before": domains[c], "after": mask}
                        for c, mask in zip(cells, support) if domains[c] != mask]
        elif event["kind"] == "three-house-gac":
            result = scan_center(domains, event["center"])
            assert result["signatures"] > 0
            expected = [{"cell": e["cell"], "before": domains[e["cell"]],
                         "after": domains[e["cell"]] & ~e["mask"]}
                        for e in result["eliminations"]]
        elif event["kind"] == "existing-engine":
            assert solution is not None  # Independent uniqueness justifies this check.
            for technique, cell, digits, placement in event["steps"]:
                if placement:
                    assert digits == [solution[cell]]
                    domains[cell] = 1 << (digits[0] - 1)
                    for peer in range(81):
                        if peer != cell and sees(cell, peer):
                            domains[peer] &= ~(1 << (digits[0] - 1))
                else:
                    assert solution[cell] not in digits
                    for digit in digits:
                        assert domains[cell] & (1 << (digit - 1))
                        domains[cell] &= ~(1 << (digit - 1))
            assert domains == event["state_after"]["domains"]
            continue
        else:
            raise ValueError("Unknown trace step")
        assert expected == event["edits"]
        for edit in expected:
            domains[edit["cell"]] = edit["after"]
        if solution is not None:
            assert all(mask & (1 << (d - 1)) for mask, d in zip(domains, solution))
    return domains


def verify(path, all_triples=True):
    case = json.loads(path.read_text())
    assert len(case["puzzle"]) == len(case["solution"]) == 81
    givens = [int(c) if c in "123456789" else 0 for c in case["puzzle"]]
    initial = [1 << (digit - 1) if digit else ALL & ~sum({
        1 << (givens[peer] - 1) for house in HOUSES if cell in house
        for peer in house if givens[peer]}) for cell, digit in enumerate(givens)]
    assert initial == case["initial_domains"]
    solution = list(map(int, case["solution"]))
    completions = exact_completions(case["puzzle"])
    expected_count = case.get("solutions_capped_at_2", 1)
    assert len(completions) == expected_count
    if expected_count == 1:
        assert completions == [case["solution"]]
    for completion in completions:
        assert all(mask & (1 << (int(d) - 1))
                   for mask, d in zip(case["base_domains"], completion))
    trace = case["base_trace"]
    domains = replay_trace(initial, trace, solution)
    assert domains == case["base_domains"]
    assert all({solution[c] for c in house} == set(range(1, 10)) for house in HOUSES)
    assert all(not d or d == solution[c] for c, d in enumerate(givens))
    assert all(mask & (1 << (d - 1)) for mask, d in zip(domains, solution))
    target = case["target"]
    assert 0 <= target["cell"] < 81 and 1 <= target["digit"] <= 9
    assert target["value"] is False and solution[target["cell"]] == target["digit"]
    conclusion = case["conclusion"]
    assert conclusion["cell"] == target["cell"] and conclusion["digit"] == target["digit"]
    assert conclusion["value"] is True and conclusion["kind"] == "placement"
    assumed = domains.copy()
    assumed[target["cell"]] &= ~(1 << (target["digit"] - 1))
    branch_trace = case["assumption_trace"]
    branch = replay_trace(assumed, branch_trace)
    assert branch == case["assumed_domains"]
    for key, expected_domains in (("state", domains), ("assumed_state", branch)):
        state = case[key]
        assert state["version"] == 1 and state["domains"] == expected_domains
        assert len(state["values"]) == len(state["domains"]) == 81
        assert all(isinstance(mask, int) and 0 < mask <= ALL for mask in state["domains"])
        assert all(isinstance(value, int) and 0 <= value <= 9 and
                   (value == 0 or mask == 1 << (value - 1))
                   for value, mask in zip(state["values"], state["domains"]))
        assert all(not d or state["values"][cell] == d for cell, d in enumerate(givens))
    assert case["state"]["values"][target["cell"]] == 0
    assert case["graph_cells"] == sorted(set(case["graph_cells"]))
    assert target["cell"] in case["graph_cells"]
    assert all(0 <= cell < 81 for cell in case["graph_cells"])
    assert close(domains)[1] == [] and close(branch)[1] == []
    assert graph_colorings(branch, case["graph_cells"])[0] == 0
    expected_edges = [[a, b] for i, a in enumerate(case["graph_cells"])
                      for b in case["graph_cells"][i + 1:] if sees(a, b)]
    assert expected_edges == case["graph_edges"]
    if all_triples:
        assert all(triple_feasible(branch, selected)
                   for selected in itertools.combinations(range(27), 3))
    from linear_witness import verify as verify_linear
    if path.stem.startswith("real_four_digit"):
        witness_file = HERE / ("linear_" + path.stem + ".json")
    else:
        witness_name = "4662" if case["name"].startswith("4662") else case["name"]
        witness_file = HERE / ("linear_real_" + witness_name + ".json")
    lp = verify_linear(case, json.loads(witness_file.read_text()))
    return {"id": case["name"], "base_steps": len(trace),
            "assumption_steps": len(branch_trace), "coloring_nodes": graph_colorings(branch, case["graph_cells"])[1],
            "arbitrary_triples": 2925 if all_triples else 0, "exact_lp_equalities": lp["equalities"],
            "solutions_capped_at_2": len(completions)}


if __name__ == "__main__":
    import argparse
    import hashlib
    import urllib.parse
    import urllib.request
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fetch-sources", action="store_true",
                        help="Also fetch the immutable author-hosted sources and verify their hashes")
    args = parser.parse_args()
    for name in ("4662", "4662_after_existing", "4401", "Loki",
                 "four_digit_candidate", "four_digit_unique"):
        path = HERE / ("real_" + name + ".json")
        if args.fetch_sources:
            source = json.loads(path.read_text())["source"]
            if source.get("kind") == "generated-givens-puzzle":
                source = source["published_completion_origin"]
                url, expected_hash = source["source_url"], source["source_sha256"]
            else:
                url = "https://raw.githubusercontent.com/" + source["repository"] + "/" + source["commit"] + "/" + urllib.parse.quote(source["path"], safe="/")
                expected_hash = source["sha256"]
            data = urllib.request.urlopen(url, timeout=20).read()
            assert hashlib.sha256(data).hexdigest() == expected_hash
        print(json.dumps(verify(path)))
