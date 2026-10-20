"""Independently verify unique givens, exact trace replay, and Python atom support."""
import json
from real_game_probe import ALL, HERE, HOUSES, completions, scan_center

if not __debug__:
    raise SystemExit("Run without -O: this research checker uses assertions.")

fixture = json.loads((HERE / "real_game_cases.json").read_text())
total = 0
for case in fixture["cases"]:
    solutions = completions(case["puzzle"])
    assert solutions == [case["solution"]], case["id"]
    solution = list(map(int, solutions[0]))
    values = [int(x) if x in "123456789" else 0 for x in case["puzzle"]]
    peers = [set().union(*(set(h) for h in HOUSES if cell in h)) - {cell}
             for cell in range(81)]
    domains = [1 << (d-1) if d else ALL & ~sum({1 << (values[p]-1)
                   for p in peers[cell] if values[p]})
               for cell, d in enumerate(values)]
    for technique, cell, digits, placement in case["baseline_trace"]:
        assert not values[cell], (case["id"], technique)
        if placement:
            assert len(digits) == 1 and digits[0] == solution[cell]
            d = digits[0]
            assert domains[cell] & (1 << (d-1))
            values[cell], domains[cell] = d, 1 << (d-1)
            for peer in peers[cell]:
                if not values[peer]:
                    domains[peer] &= ~(1 << (d-1))
        else:
            assert solution[cell] not in digits
            for d in digits:
                assert domains[cell] & (1 << (d-1))
                domains[cell] &= ~(1 << (d-1))
        assert all(mask & (1 << (d-1)) for mask, d in zip(domains, solution))
    assert values == case["state"]["values"]
    assert domains == case["state"]["domains"]
    clusters = [scan_center(domains, center) for center in range(81)]
    assert all(cluster["signatures"] > 0 for cluster in clusters)
    assert [c for c in clusters if c["eliminations"]] == case["three_house"]
    total += len(case["baseline_trace"])
    print(case["id"], "independent unique completion, exact baseline trace and 81 centers verified")
print("Verified", len(fixture["cases"]), "real puzzle states and", total, "recorded detector steps.")
