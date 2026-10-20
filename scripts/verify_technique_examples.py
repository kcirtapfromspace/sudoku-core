#!/usr/bin/env python3
"""Independent Python replay of uniqueness, setup and endpoint certificates.

This checker has no Rust dependency, calls no solver and does no proof search.
Version-two bundles also carry independently checked named pattern witnesses.
"""
import argparse
import hashlib
import json
from pathlib import Path


def require(condition, message):
    if not condition:
        raise ValueError(message)


def peers(a, b):
    return a != b and (
        a // 9 == b // 9
        or a % 9 == b % 9
        or (a // 27, a % 9 // 3) == (b // 27, b % 9 // 3)
    )


def from_givens(text):
    require(len(text) == 81 and set(text) <= set("0123456789."), "Invalid givens")
    values = [int(c) if c in "123456789" else 0 for c in text]
    domains = []
    for cell, value in enumerate(values):
        if value:
            require(
                all(
                    not peers(cell, other) or v != value
                    for other, v in enumerate(values)
                ),
                "Conflicting givens",
            )
            domains.append(1 << (value - 1))
        else:
            forbidden = {
                v for other, v in enumerate(values) if v and peers(cell, other)
            }
            domains.append(
                sum(1 << (d - 1) for d in range(1, 10) if d not in forbidden)
            )
    return {"version": 1, "values": values, "domains": domains}


def check(proof, state):
    require(
        state["version"] == 1 and len(state["values"]) == len(state["domains"]) == 81,
        "Invalid state shape",
    )
    for v, d in zip(state["values"], state["domains"]):
        require(
            type(v) is int and 0 <= v <= 9 and type(d) is int and 0 < d < 512,
            "Invalid state digit",
        )
        require(not v or d == 1 << (v - 1), "Placed domain is not a singleton")
    digest = hashlib.sha256(b"sudoku-exact-cover-v1-native")
    for v, d in zip(state["values"], state["domains"]):
        digest.update(bytes([v]))
        digest.update(d.to_bytes(2, "little"))
    require(
        proof["version"] == 1 and proof["state_hash"] == digest.hexdigest(),
        "Wrong version or premises",
    )
    witness = proof["witness"]
    require(
        len(witness) == 81 and all(type(d) is int and 1 <= d <= 9 for d in witness),
        "Invalid witness",
    )
    for c, d in enumerate(witness):
        require(
            (not state["values"][c] or state["values"][c] == d)
            and state["domains"][c] & (1 << (d - 1)),
            "Witness violates a premise",
        )
        require(
            all(not peers(c, o) or witness[o] != d for o in range(c)),
            "Witness violates Sudoku",
        )
    purpose = proof["purpose"]
    domains = state["domains"].copy()
    if purpose != "Unique":
        require(set(purpose) == {"Deduction"}, "Unknown proof purpose")
        p = purpose["Deduction"]
        c, ds, placement = p["cell"], p["digits"], p["placement"]
        require(
            type(c) is int and 0 <= c < 81 and not state["values"][c],
            "Invalid endpoint cell",
        )
        require(
            type(placement) is bool
            and ds
            and all(type(d) is int and 1 <= d <= 9 for d in ds)
            and len(set(ds)) == len(ds)
            and (not placement or len(ds) == 1),
            "Invalid endpoint digits",
        )
        mask = sum(1 << (d - 1) for d in ds)
        require(
            mask & ~domains[c] == 0 and (witness[c] in ds) == placement,
            "Endpoint lacks a base witness",
        )
        domains[c] &= ~mask if placement else mask
    nodes = proof["nodes"]
    require(0 < len(nodes) <= 2_000_000, "Unsupported tree size")
    seen, solutions = set(), []

    def support(requirement, board):
        require(
            type(requirement) is int and 0 <= requirement < 324, "Invalid requirement"
        )
        if requirement < 81:
            cells, digit = [requirement], None
            require(board[requirement] == 0, "Requirement already fulfilled")
        else:
            group, n, digit = (
                (requirement - 81) // 81,
                (requirement - 81) % 81 // 9,
                (requirement - 81) % 9 + 1,
            )
            cells = (
                [n * 9 + i for i in range(9)]
                if group == 0
                else (
                    [i * 9 + n for i in range(9)]
                    if group == 1
                    else [
                        (n // 3 * 3 + i // 3) * 9 + n % 3 * 3 + i % 3 for i in range(9)
                    ]
                )
            )
            require(
                all(board[c] != digit for c in cells), "Requirement already fulfilled"
            )
        return sorted(
            c * 9 + d - 1
            for c in cells
            for d in range(1, 10)
            if not board[c]
            and (digit is None or d == digit)
            and domains[c] & (1 << (d - 1))
            and all(not peers(c, o) or v != d for o, v in enumerate(board))
        )

    def visit(node_id, board, depth):
        require(
            type(node_id) is int
            and 0 <= node_id < len(nodes)
            and node_id not in seen
            and depth <= 81,
            "Invalid tree edge",
        )
        seen.add(node_id)
        node = nodes[node_id]
        if node == "Solution":
            require(board == witness, "Incorrect solution leaf")
            solutions.append(board.copy())
            return
        require(isinstance(node, dict) and len(node) == 1, "Invalid node")
        if "Conflict" in node:
            require(
                not support(node["Conflict"]["requirement"], board), "False conflict"
            )
            return
        require("Split" in node, "Unknown node type")
        split = node["Split"]
        candidates = support(split["requirement"], board)
        branches = split["branches"]
        require(
            candidates and [b["candidate"] for b in branches] == candidates,
            "Missing, repeated or wrong branch",
        )
        for b in branches:
            require(type(b["candidate"]) is int, "Invalid candidate type")
            require(
                type(b["child"]) is int and b["child"] < node_id, "Not a postorder tree"
            )
            c, d = b["candidate"] // 9, b["candidate"] % 9 + 1
            next_board = board.copy()
            next_board[c] = d
            visit(b["child"], next_board, depth + 1)

    visit(len(nodes) - 1, state["values"].copy(), 0)
    require(
        len(seen) == len(nodes) and len(solutions) == (1 if purpose == "Unique" else 0),
        "Incomplete proof",
    )


def verify(corpus, subset=False):
    registry = set(
        "NakedSingle HiddenSingle NakedPair HiddenPair NakedTriple HiddenTriple PointingPair BoxLineReduction XWing FinnedXWing Swordfish FinnedSwordfish Jellyfish FinnedJellyfish NakedQuad HiddenQuad EmptyRectangle AvoidableRectangle UniqueRectangle HiddenRectangle XYWing XYZWing WXYZWing WWing XChain ThreeDMedusa SueDeCoq AIC FrankenFish SiameseFish AlsXz ExtendedUniqueRectangle BivalueUniversalGrave AlsXyWing AlsChain MutantFish AlignedPairExclusion AlignedTripletExclusion DeathBlossom ThreeHouseMatching ArithmeticCounting GlobalColoring NishioForcingChain KrakenFish RegionForcingChain CellForcingChain DynamicForcingChain Backtracking".split()
    )
    require(
        corpus["version"] in (1, 2)
        and (subset or corpus["version"] == 2)
        and 0 < len(corpus["examples"]) <= 48,
        "Invalid corpus shape",
    )
    require(subset or len(corpus["examples"]) == 48, "Missing technique examples")
    starts = []
    for puzzle in corpus["puzzles"]:
        state = from_givens(puzzle["givens"])
        require(puzzle["uniqueness"]["purpose"] == "Unique", "No uniqueness proof")
        check(puzzle["uniqueness"], state)
        starts.append(state)
    methods, used = set(), set()
    for example in corpus["examples"]:
        index, state, hint = example["puzzle"], example["state"], example["hint"]
        require(
            type(index) is int and 0 <= index < len(starts), "Invalid puzzle reference"
        )
        used.add(index)
        method = hint["technique"]
        require(
            method in registry and method not in methods,
            "Unknown or duplicate technique",
        )
        methods.add(method)
        start, witness = (
            starts[index],
            corpus["puzzles"][index]["uniqueness"]["witness"],
        )
        require(
            len(state["values"]) == len(state["domains"]) == 81, "Invalid setup size"
        )
        for c in range(81):
            require(
                (not start["values"][c] or state["values"][c] == start["values"][c])
                and state["values"][c] in (0, witness[c])
                and state["domains"][c] & ~start["domains"][c] == 0
                and state["domains"][c] & (1 << (witness[c] - 1)),
                "Unauthenticated setup",
            )
        cert = example["certificate"]["Certified"]["deduction"]
        pattern = example["certificate"]["Certified"].get("pattern")
        require(
            not isinstance(pattern, dict) or "Certified" not in pattern,
            "Nested certification",
        )
        p = cert["purpose"]["Deduction"]
        action = hint["hint_type"]
        placement = "SetValue" in action
        action = action["SetValue" if placement else "EliminateCandidates"]
        require(
            all(
                type(action["pos"][k]) is int and 0 <= action["pos"][k] < 9
                for k in ("row", "col")
            ),
            "Invalid hint coordinates",
        )
        require(
            p
            == {
                "cell": action["pos"]["row"] * 9 + action["pos"]["col"],
                "digits": [action["value"]] if placement else action["values"],
                "placement": placement,
            },
            "Wrong hint endpoint",
        )
        check(cert, state)
        named = example["certificate"]["Certified"].get("named")
        require(
            corpus["version"] != 2 or named is not None,
            "Missing named pattern certificate",
        )
        if named is not None:
            from verify_named_patterns import verify_named

            verify_named(named, state, hint, corpus["puzzles"][index]["givens"], cert)
    require(len(used) == len(starts), "Unused source puzzle")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "path",
        nargs="?",
        default=Path(__file__).resolve().parent.parent
        / "tests/fixtures/technique_examples.json",
    )
    parser.add_argument("--subset", action="store_true")
    args = parser.parse_args()
    data = json.loads(Path(args.path).read_text())
    verify(data, args.subset)
    print(
        f"Independently verified {len(data['examples'])} examples and {len(data['puzzles'])} original puzzles"
    )
