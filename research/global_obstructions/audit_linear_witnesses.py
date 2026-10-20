"""Standard-library regression audit of exact native LP primal/dual witnesses.

No numerical optimizer is imported. This verifies recorded rational points and
negative controls, then checks that representative certificate tampering fails.
"""
import copy
import itertools
import json
from pathlib import Path

from linear_witness import binding, source_state, verify, verify_farkas

HERE = Path(__file__).resolve().parent

if not __debug__:
    raise SystemExit("Run without -O: assertions are part of the exact audit")


def load(name):
    return json.loads((HERE / name).read_text())


def reject(callback):
    try:
        callback()
    except (ValueError, KeyError, TypeError):
        return
    raise AssertionError("Malformed exact LP evidence was accepted")


def main():
    for name in ["4662", "4401", "Loki", "four_digit_candidate", "four_digit_unique"]:
        fixture = load("real_" + name + ".json")
        witness = load("linear_real_" + name + ".json")
        result = verify(fixture, witness)
        print(name, result)
        changed = copy.deepcopy(witness)
        changed["nonzero"].pop()
        reject(lambda: verify(fixture, changed))
        changed = copy.deepcopy(witness)
        changed["nonzero"][0]["numerator"] = -1
        reject(lambda: verify(fixture, changed))
        changed = copy.deepcopy(witness)
        changed["nonzero"][0]["denominator"] = 0
        reject(lambda: verify(fixture, changed))
        changed = copy.deepcopy(witness)
        changed["nonzero"][0]["numerator"] += 1
        reject(lambda: verify(fixture, changed))
        changed = copy.deepcopy(witness)
        changed["nonzero"].append(changed["nonzero"][0])
        reject(lambda: verify(fixture, changed))
        changed = copy.deepcopy(witness)
        changed["target"]["assumed"] = True
        reject(lambda: verify(fixture, changed))
        changed_fixture = copy.deepcopy(fixture)
        state = changed_fixture.get("assumed_state", changed_fixture.get("state", changed_fixture))
        domain_key = "assumed_domains" if "assumed_domains" in state else "domains"
        values = source_state(changed_fixture)[0]
        cell = next(c for c, (value, domain) in enumerate(zip(values, state[domain_key]))
                    if value == 0 and domain != 511)
        missing_bit = next(1 << d for d in range(9) if not state[domain_key][cell] & (1 << d))
        state[domain_key][cell] |= missing_bit
        reject(lambda: verify(changed_fixture, witness))

    fixture = load("structural_four_digit_13.json")
    witness = load("linear_four_digit_candidate.json")
    print("accepted thirteen-cell four-digit candidate", verify(fixture, witness))
    assert {entry["denominator"] for entry in witness["nonzero"]} == {4, 5}
    support = fixture["support"]
    def peer(a, b):
        return (a // 9 == b // 9 or a % 9 == b % 9
                or (a // 27, a % 9 // 3) == (b // 27, b % 9 // 3))
    assert len(support) == 13 and len(fixture["palette"]) == 4
    subsets = list(itertools.combinations(support, 4))
    assert len(subsets) == 715
    assert all(any(peer(a, b) for a, b in itertools.combinations(group, 2)) for group in subsets)
    independent_triples = [group for group in itertools.combinations(support, 3)
                           if all(not peer(a, b) for a, b in itertools.combinations(group, 2))]
    assert len(independent_triples) == 32
    palette_mask = sum(1 << (digit - 1) for digit in fixture["palette"])
    assert all(fixture["assumed_domains"][c] == palette_mask for c in support)
    print("Global capacity proof: alpha=3 (all 715 four-subsets contain an edge); four digits cover at most 12 < 13 cells.")

    fixture = load("linear_rejected_q4_input.json")
    dual = load("linear_rejected_q4_farkas.json")
    print("rejected four-digit candidate", verify_farkas(fixture, dual))
    changed = copy.deepcopy(dual)
    changed["terms"][0][1] *= -1
    reject(lambda: verify_farkas(fixture, changed))
    changed = copy.deepcopy(dual)
    changed["terms"].append(changed["terms"][0])
    reject(lambda: verify_farkas(fixture, changed))

    # Earlier fixture 847 is an exact negative control for confusing bounded
    # arithmetic failure with survival of the entire native continuous LP.
    fixture = json.loads((HERE.parent / "nonmodular" / "three-house-847.json").read_text())
    weights = {1: -1, 2: -1, 3: -1, 27: 1, 63: 1, 81: 1,
               83: 1, 162: -1, 163: -1, 167: -1, 244: 1}
    dual = {"schema_version": 1, "state_hash": binding(fixture),
            "target": source_state(fixture)[3], "weighted_rhs": -1,
            "terms": [[i, weight, 1] for i, weight in weights.items()]}
    print("trial847 negative control", verify_farkas(fixture, dual))
    print("All exact primal/dual and malformed-evidence checks passed.")


if __name__ == "__main__":
    main()
