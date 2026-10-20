"""Independent exact checker for exported GlobalColoringReplay certificates.

Reconstructs native Sudoku peers, the opposite target assumption, each branch's
complete current domain, and each empty-domain leaf. Does not import engine
search code or trust caller-supplied edges, propagation traces, or color-symmetry
shortcuts. This is finite list-coloring refutation, not a historical novelty claim.
"""
import argparse
import copy
import hashlib
import json
from pathlib import Path


def integer(value):
    return type(value) is int


def peers(a, b):
    r, c = divmod(a, 9)
    s, d = divmod(b, 9)
    return a != b and (r == s or c == d or (r // 3, c // 3) == (s // 3, d // 3))


def state_hash(state):
    h = hashlib.sha256(b"sudoku-global-coloring-v1-native")
    for value, domain in zip(state["values"], state["domains"]):
        h.update(bytes([value]))
        h.update(domain.to_bytes(2, "little"))
    return h.hexdigest()


def validate_state(state):
    if state.get("version") != 1 or len(state["values"]) != 81 or len(state["domains"]) != 81:
        raise ValueError("Invalid state version or size")
    for value, domain in zip(state["values"], state["domains"]):
        if not integer(value) or not integer(domain) or not 0 <= value <= 9 or not 0 < domain <= 511:
            raise ValueError("Invalid value/domain")
        if value and domain != 1 << (value - 1):
            raise ValueError("Placed value lacks singleton effective domain")
    units = ([[9 * r + c for c in range(9)] for r in range(9)] +
             [[9 * r + c for r in range(9)] for c in range(9)] +
             [[9 * r + c for r in range(3 * br, 3 * br + 3)
                         for c in range(3 * bc, 3 * bc + 3)]
              for br in range(3) for bc in range(3)])
    for unit in units:
        placed = [state["values"][c] for c in unit if state["values"][c]]
        if len(set(placed)) != len(placed):
            raise ValueError("Conflicting placed native values")
        union = 0
        for c in unit:
            union |= state["domains"][c]
        if union != 511:
            raise ValueError("A native house lacks some digit")


def verify(replay):
    if replay.get("version") != 1:
        raise ValueError("Unknown replay version")
    state, proof = replay["state"], replay["proof"]
    validate_state(state)
    if proof.get("version") != 1 or proof.get("state_hash") != state_hash(state):
        raise ValueError("Proof version/state hash mismatch")
    if "anchors" in proof:
        raise ValueError("Symmetry shortcuts are outside this checker")
    cells = proof["cells"]
    if not isinstance(cells, list) or not 1 <= len(cells) <= 24:
        raise ValueError("Invalid support size")
    if any(not integer(c) or not 0 <= c < 81 for c in cells) or cells != sorted(set(cells)):
        raise ValueError("Support cells must be sorted, unique native cells")
    cell, digit, result = proof["cell"], proof["digit"], proof["value"]
    if not integer(cell) or cell not in cells or not integer(digit) or not 1 <= digit <= 9 or type(result) is not bool:
        raise ValueError("Malformed target")
    if state["values"][cell] != 0 or not state["domains"][cell] & (1 << (digit - 1)):
        raise ValueError("Target must be available and unplaced")
    witness = proof["witness"]
    if len(witness) != len(cells):
        raise ValueError("Wrong local support witness length")
    for c, d in zip(cells, witness):
        if not integer(d) or not 1 <= d <= 9 or not state["domains"][c] & (1 << (d - 1)):
            raise ValueError("Local base support violates a recorded domain")
    for i, a in enumerate(cells):
        for j, b in enumerate(cells[:i]):
            if peers(a, b) and witness[i] == witness[j]:
                raise ValueError("Local base support has a native peer conflict")
    if (witness[cells.index(cell)] == digit) != result:
        raise ValueError("Base support does not witness the conclusion endpoint")
    nodes = proof["nodes"]
    if not isinstance(nodes, list) or not 1 <= len(nodes) <= 8192:
        raise ValueError("Invalid proof-node count")
    initial = {c: state["domains"][c] for c in cells}
    bit = 1 << (digit - 1)
    initial[cell] = initial[cell] & ~bit if result else bit
    visited = set()

    def descend(index, domains, assigned):
        if not integer(index) or not 0 <= index < len(nodes) or index in visited:
            raise ValueError("Invalid or reused proof node")
        visited.add(index)
        node = nodes[index]
        if not isinstance(node, dict) or len(node) != 1:
            raise ValueError("Malformed proof node")
        if "Conflict" in node:
            data = node["Conflict"]
            if set(data) != {"cell"} or not integer(data["cell"]) or data["cell"] not in domains or domains[data["cell"]] != 0:
                raise ValueError("Conflict leaf has no reconstructed empty domain")
            return
        if "Split" not in node:
            raise ValueError("Unknown proof-node tag")
        data = node["Split"]
        if set(data) != {"cell", "branches"}:
            raise ValueError("Malformed split")
        pivot = data["cell"]
        if not integer(pivot) or pivot not in domains or pivot in assigned or domains[pivot] == 0:
            raise ValueError("Invalid split pivot")
        available = {d for d in range(1, 10) if domains[pivot] & (1 << (d - 1))}
        branch_digits = []
        for branch in data["branches"]:
            if set(branch) != {"digit", "child"}:
                raise ValueError("Malformed branch")
            d, child = branch["digit"], branch["child"]
            if not integer(d) or d not in available or not integer(child) or not 0 <= child < index:
                raise ValueError("Invalid branch digit or non-postorder child")
            branch_digits.append(d)
            next_domains = domains.copy()
            next_domains[pivot] = 1 << (d - 1)
            for other in cells:
                if peers(pivot, other):
                    next_domains[other] &= ~(1 << (d - 1))
            descend(child, next_domains, assigned | {pivot})
        if len(set(branch_digits)) != len(branch_digits) or set(branch_digits) != available:
            raise ValueError("Split does not cover every current candidate exactly once")

    descend(len(nodes) - 1, initial, set())
    if len(visited) != len(nodes):
        raise ValueError("Unreachable proof evidence")
    return {"cells": len(cells), "nodes": len(nodes), "conclusion": "placement" if result else "elimination"}


def audit_mutations(replay):
    mutations = []
    changed = copy.deepcopy(replay)
    changed["proof"]["state_hash"] = "incorrect"
    mutations.append(changed)
    changed = copy.deepcopy(replay)
    changed["proof"]["value"] = not changed["proof"]["value"]
    mutations.append(changed)
    changed = copy.deepcopy(replay)
    changed["proof"]["witness"][0] = 0
    mutations.append(changed)
    changed = copy.deepcopy(replay)
    changed["proof"]["cells"][1] = changed["proof"]["cells"][0]
    mutations.append(changed)
    split_index = next(i for i, node in enumerate(replay["proof"]["nodes"])
                       if "Split" in node and len(node["Split"]["branches"]) > 1)
    changed = copy.deepcopy(replay)
    changed["proof"]["nodes"][split_index]["Split"]["branches"].pop()
    mutations.append(changed)
    changed = copy.deepcopy(replay)
    branches = changed["proof"]["nodes"][split_index]["Split"]["branches"]
    branches[1]["digit"] = branches[0]["digit"]
    mutations.append(changed)
    changed = copy.deepcopy(replay)
    changed["proof"]["nodes"][split_index]["Split"]["branches"][0]["child"] = split_index
    mutations.append(changed)
    changed = copy.deepcopy(replay)
    branches = changed["proof"]["nodes"][split_index]["Split"]["branches"]
    branches[1]["child"] = branches[0]["child"]
    mutations.append(changed)
    leaf_index = next(i for i, node in enumerate(replay["proof"]["nodes"]) if "Conflict" in node)
    changed = copy.deepcopy(replay)
    changed["proof"]["nodes"][leaf_index]["Conflict"]["cell"] = 81
    mutations.append(changed)
    changed = copy.deepcopy(replay)
    changed["proof"]["anchors"] = []
    mutations.append(changed)
    for changed in mutations:
        try:
            verify(changed)
        except (ValueError, TypeError, KeyError, IndexError):
            continue
        raise AssertionError("Checker accepted a tampered production certificate")
    return len(mutations)


def self_test():
    cells = [0, 3, 13, 27, 31]
    domains = [511] * 81
    for c in cells:
        domains[c] = 3
    domains[0] = 7
    state = {"version": 1, "values": [0] * 81, "domains": domains}
    nodes = []

    def tree(current, assigned):
        empty = next((c for c in cells if current[c] == 0), None)
        if empty is not None:
            nodes.append({"Conflict": {"cell": empty}})
            return len(nodes) - 1
        unassigned = [c for c in cells if c not in assigned]
        if not unassigned:
            raise AssertionError("Known C5 unexpectedly has a two-coloring")
        pivot = min(unassigned, key=lambda c: (bin(current[c]).count("1"), c))
        branches = []
        for digit in range(1, 10):
            if current[pivot] & (1 << (digit - 1)):
                nxt = current.copy()
                nxt[pivot] = 1 << (digit - 1)
                for c in cells:
                    if peers(c, pivot):
                        nxt[c] &= ~(1 << (digit - 1))
                branches.append({"digit": digit, "child": tree(nxt, assigned | {pivot})})
        nodes.append({"Split": {"cell": pivot, "branches": branches}})
        return len(nodes) - 1

    tree({c: 3 for c in cells}, set())
    proof = {"version": 1, "state_hash": state_hash(state), "cells": cells,
             "cell": 0, "digit": 3, "value": True, "witness": [3, 1, 2, 2, 1], "nodes": nodes}
    replay = {"version": 1, "state": state, "proof": proof}
    result = verify(replay)
    bad = copy.deepcopy(replay)
    bad["proof"]["nodes"][-1]["Split"]["branches"].pop()
    try:
        verify(bad)
    except ValueError:
        pass
    else:
        raise AssertionError("Checker accepted omitted branch")
    print("Known C5 sanity proof:", result)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("replay", nargs="*", type=Path)
    args = parser.parse_args()
    self_test()
    paths = args.replay or sorted(Path(__file__).resolve().parent.glob("replay_*.json"))
    for path in paths:
        replay = json.loads(path.read_text())
        print(path.name, verify(replay), "rejected mutations:", audit_mutations(replay))


if __name__ == "__main__":
    main()
