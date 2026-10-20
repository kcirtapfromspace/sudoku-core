"""Independent three-house certificate-design audit and trial-847 counting proof.

Standard library only. This does not import or test the production Rust checker.
Matching search creates Hall blockers; the separate proof checker uses only
domain unions, cardinalities, deterministic signature coverage, and a witness.
The 11-source counting proof is fixed data, not an LP result trusted at runtime.
"""
import hashlib
import itertools
import json
from pathlib import Path

if not __debug__:
    raise SystemExit("Run without -O: this audit uses assertions.")

ALL = 511
FIXTURE = Path(__file__).with_name("three-house-847.json")


def popcount(mask):
    return bin(mask).count("1")


def atoms(center):
    r, c = divmod(center, 9)
    row = {9 * r + j for j in range(9)}
    col = {9 * i + c for i in range(9)}
    box = {9 * i + j for i in range(3 * (r // 3), 3 * (r // 3) + 3)
           for j in range(3 * (c // 3), 3 * (c // 3) + 3)}
    result = [[center], sorted((row & box) - {center}),
              sorted((col & box) - {center}), sorted(row - box),
              sorted(col - box), sorted(box - row - col)]
    assert list(map(len, result)) == [1, 2, 2, 6, 6, 4]
    assert len(set().union(*map(set, result))) == 21
    return result, [row, col, box]


def subsets(mask, size):
    digits = [1 << i for i in range(9) if mask & (1 << i)]
    for choice in itertools.combinations(digits, size):
        yield sum(choice)


def signatures(center_domain):
    for d in subsets(center_domain, 1):
        for a in subsets(ALL ^ d, 2):
            for b in subsets(ALL ^ (d | a), 2):
                yield [d, a, b, ALL ^ (d | a), ALL ^ (d | b), ALL ^ (d | a | b)]


def matching_exists(cell_domains, required):
    reachable = {0}
    for domain in cell_domains:
        reachable = {used | bit for used in reachable
                     for bit in subsets(domain & required & ~used, 1)}
    return required in reachable


def hall_block(cell_domains, required):
    for subset in range(1, 1 << len(cell_domains)):
        union = 0
        for i, domain in enumerate(cell_domains):
            if subset & (1 << i):
                union |= domain & required
        if popcount(union) < popcount(subset):
            return subset
    raise AssertionError("Matching/Hall disagreement")


def state_digest(domains):
    # Research format deliberately differs from production; no hash compatibility
    # claim. Fixture has no placed cells. Production must additionally bind values.
    return hashlib.sha256(bytes().join(m.to_bytes(2, "little") for m in domains)).hexdigest()


def check_proof(domains, center, target, blocks, base, digest):
    assert len(domains) == 81 and all(0 < m <= ALL for m in domains)
    assert 0 <= center < 81 and 0 <= target[0] < 81 and 1 <= target[1] <= 9
    groups, houses = atoms(center)
    cells = set().union(*map(set, groups))
    assert target[0] in cells
    target_bit = 1 << (target[1] - 1)
    assert domains[target[0]] & target_bit
    assert state_digest(domains) == digest
    assert set(base) == cells
    assert all(1 <= base[c] <= 9 and domains[c] & (1 << (base[c] - 1)) for c in cells)
    assert all({base[c] for c in h} == set(range(1, 10)) for h in houses)
    assert base[target[0]] != target[1]
    assumed = domains[:]
    assumed[target[0]] = target_bit
    indexed = {}
    assert len(blocks) <= 630
    for atom_id, required, subset in blocks:
        assert 1 <= atom_id <= 5
        assert 0 < required <= ALL and popcount(required) == len(groups[atom_id])
        assert 0 < subset < 1 << len(groups[atom_id])
        assert (atom_id, required) not in indexed
        union = 0
        for i, cell in enumerate(groups[atom_id]):
            if subset & (1 << i):
                union |= assumed[cell] & required
        assert popcount(union) < popcount(subset)
        indexed[atom_id, required] = subset
    covered = 0
    for signature in signatures(assumed[center]):
        assert any((a, signature[a]) in indexed for a in range(1, 6))
        covered += 1
    assert covered > 0
    return covered


def must_reject(callback):
    try:
        callback()
    except AssertionError:
        return
    raise AssertionError("Malformed research certificate was accepted")


def verify_counting_certificate(domains):
    # Native requirement numbering: 0..80 cells; 81..161 row/digit;
    # 162..242 column/digit; 243..323 box/digit, each in row-major order.
    weights = {1: -1, 2: -1, 3: -1, 27: 1, 63: 1,
               81: 1, 83: 1, 162: -1, 163: -1, 167: -1, 244: 1}
    assert len(weights) == 11 and sum(weights.values()) == -1
    assumed = domains[:]
    assumed[3] = 1
    negative_before_assumption = []
    positive = []
    for cell, domain in enumerate(domains):
        r, c = divmod(cell, 9)
        box = 3 * (r // 3) + c // 3
        for digit in range(1, 10):
            if not domain & (1 << (digit - 1)):
                continue
            sources = [cell, 81 + 9 * r + digit - 1,
                       162 + 9 * c + digit - 1, 243 + 9 * box + digit - 1]
            coefficient = sum(weights.get(source, 0) for source in sources)
            if assumed[cell] & (1 << (digit - 1)):
                assert coefficient >= 0
                if coefficient:
                    positive.append((cell, digit, coefficient))
            elif coefficient < 0:
                negative_before_assumption.append((cell, digit, coefficient))
    assert negative_before_assumption == [(3, 4, -1), (3, 6, -1)]
    print("11-source native integer counting certificate: nonnegative LHS = -1 under r1c4=1.")
    print("Negative coefficients before assumption:", negative_before_assumption)
    print("Therefore trial847 is not a separation from linear/rank/counting reasoning.")


def verify_failed_literal_gac(domains):
    """Known-method baseline: iterate individual-house GAC under one assumption."""
    restricted = domains[:]
    restricted[3] = 1
    _, houses = atoms(0)
    for pass_number in range(1, 10):
        old = restricted[:]
        for house_id, house in enumerate(houses):
            cells = sorted(house)
            prefix = [{0}]
            for cell in cells:
                prefix.append({used | bit for used in prefix[-1]
                               for bit in subsets(restricted[cell] & ~used, 1)})
            if ALL not in prefix[-1]:
                assert (pass_number, house_id) == (2, 1)
                print("Known failed-literal baseline: individual-house GAC refutes the target on pass 2 at column 1.")
                return
            suffix = [None] * 10
            suffix[9] = {0}
            for i in range(8, -1, -1):
                suffix[i] = {used | bit for used in suffix[i + 1]
                             for bit in subsets(restricted[cells[i]] & ~used, 1)}
            for i, cell in enumerate(cells):
                supported = 0
                for bit in subsets(restricted[cell], 1):
                    if any(not (used & bit) and (ALL ^ (used | bit)) in suffix[i + 1]
                           for used in prefix[i]):
                        supported |= bit
                restricted[cell] &= supported
        assert restricted != old, "Unexpected GAC fixed point without contradiction"
    raise AssertionError("Failed-literal baseline did not refute target")


def verify_exported_replay():
    """Check an actual Rust-produced JSON certificate without invoking Rust/DP."""
    replay = json.loads(FIXTURE.with_name("three_house_example_replay.json").read_text())
    assert set(replay) == {"version", "state", "proof"} and replay["version"] == 1
    state, proof = replay["state"], replay["proof"]
    assert set(state) == {"version", "values", "domains"} and state["version"] == 1
    values, domains = state["values"], state["domains"]
    assert len(values) == len(domains) == 81
    assert all(0 <= v <= 9 and 0 < m <= ALL and (v == 0 or m == 1 << (v - 1))
               for v, m in zip(values, domains))
    for center in range(81):
        for house in atoms(center)[1]:
            placed = [values[c] for c in house if values[c]]
            assert len(placed) == len(set(placed))
            union = 0
            for cell in house:
                union |= domains[cell]
            assert union == ALL
    assert set(proof) == {"version", "state_hash", "center", "cell", "digit", "witness", "blockers"}
    assert proof["version"] == 1 and 0 <= proof["cell"] < 81
    assert values[proof["cell"]] == 0 and len(proof["blockers"]) <= 366
    digest = hashlib.sha256(b"sudoku-three-house-v1-native")
    for value, domain in zip(values, domains):
        digest.update(bytes([value]))
        digest.update(domain.to_bytes(2, "little"))
    assert digest.hexdigest() == proof["state_hash"]
    groups, _ = atoms(proof["center"])
    cells = sorted(set().union(*map(set, groups)))
    assert len(proof["witness"]) == 21
    base = dict(zip(cells, proof["witness"]))
    blocks = []
    for block in proof["blockers"]:
        assert set(block) == {"atom", "digits", "cells"}
        blocks.append((block["atom"], block["digits"], block["cells"]))
    covered = check_proof(domains, proof["center"], (proof["cell"], proof["digit"]),
                          blocks, base, state_digest(domains))
    print("Actual Rust JSON replay independently verified:", covered, "signatures,", len(blocks), "Hall blockers.")


def main():
    fixture = json.loads(FIXTURE.read_text())
    domains = fixture["domains"]
    center, target = 0, (3, 1)
    groups, _ = atoms(center)
    assumed = domains[:]
    assumed[target[0]] = 1 << (target[1] - 1)
    blocked = {}
    for atom_id in range(1, 6):
        cell_domains = [assumed[c] for c in groups[atom_id]]
        for required in subsets(ALL, len(cell_domains)):
            if not matching_exists(cell_domains, required):
                blocked[atom_id, required] = hall_block(cell_domains, required)
    # Keep one valid blocker per signature, then greedily remove redundant ones.
    selected = {}
    signature_list = list(signatures(assumed[center]))
    for sig in signature_list:
        key = next((a, sig[a]) for a in range(1, 6) if (a, sig[a]) in blocked)
        selected[key] = blocked[key]
    for key in list(selected):
        if all(any((a, sig[a]) in selected and (a, sig[a]) != key for a in range(1, 6))
               for sig in signature_list):
            del selected[key]
    blocks = [(a, req, subset) for (a, req), subset in selected.items()]
    completion = list(map(int, fixture["completion"]))
    base = {c: completion[c] for group in groups for c in group}
    digest = state_digest(domains)
    covered = check_proof(domains, center, target, blocks, base, digest)
    assert all(len(list(signatures(ALL))) == 3780 for _ in [0])
    print("Independent coverage checker:", covered, "signatures blocked by", len(blocks), "Hall witnesses.")
    must_reject(lambda: check_proof(domains, center, target, blocks[1:], base, digest))
    bad = blocks[:]
    a, req, _ = bad[0]
    bad[0] = (a, req, 1 << len(groups[a]))
    must_reject(lambda: check_proof(domains, center, target, bad, base, digest))
    must_reject(lambda: check_proof(domains, center, (80, 1), blocks, base, digest))
    widened = domains[:]
    widened[0] |= 16
    must_reject(lambda: check_proof(widened, center, target, blocks, base, digest))
    vacuous = dict(base)
    vacuous[3] = 1
    must_reject(lambda: check_proof(domains, center, target, blocks, vacuous, digest))
    print("Rejected missing coverage, outside-atom bits, outside-cluster target, stale state, and invalid base witness.")
    verify_counting_certificate(domains)
    verify_failed_literal_gac(domains)
    verify_exported_replay()


if __name__ == "__main__":
    main()
