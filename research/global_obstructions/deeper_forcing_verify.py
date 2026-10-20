"""Stdlib verifier for a two-level established forcing explanation.

Leaves use explicit Hall deficiencies or compact Hall blockers covering EVERY
crossing-house boundary signature. No optimizer, SAT oracle, uniqueness premise,
or trusted contradiction string is used to accept the forcing tree.
"""
import argparse
from copy import deepcopy
import hashlib
import itertools
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from real_verify import HOUSES, close, replay_trace

DEFAULT_FIXTURE = HERE / "different_lead.json"
DEFAULT_CERTIFICATE = HERE / "deeper_forcing_certificate.json"


def count(mask):
    return bin(mask).count("1")


def integer(value, lower, upper):
    assert type(value) is int and lower <= value <= upper


def valid_domains(domains):
    assert type(domains) is list and len(domains) == 81
    for mask in domains:
        integer(mask,1,511)


def valid_trace(trace):
    assert type(trace) is list
    for event in trace:
        if event["kind"] == "house-gac":
            assert set(event) == {"kind","sector","edits"}
            integer(event["sector"],0,26)
        else:
            assert event["kind"] == "three-house-gac"
            assert set(event) == {"kind","center","edits"}
            integer(event["center"],0,80)
        assert type(event["edits"]) is list
        for edit in event["edits"]:
            assert set(edit) == {"cell","before","after"}
            integer(edit["cell"],0,80)
            integer(edit["before"],1,511)
            integer(edit["after"],1,511)


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def atoms(center):
    row, column, box = set(HOUSES[center//9]), set(HOUSES[9+center%9]), set(HOUSES[18+3*(center//27)+center%9//3])
    return [[center], sorted(row & box - {center}), sorted(column & box - {center}),
            sorted(row-box), sorted(column-box), sorted(box-row-column)]


def signatures(center_domain):
    for d in range(9):
        dm = 1 << d
        if not center_domain & dm:
            continue
        for aa in itertools.combinations([1 << k for k in range(9) if k != d], 2):
            a = aa[0] | aa[1]
            for bb in itertools.combinations([1 << k for k in range(9) if not (dm|a) & (1 << k)], 2):
                b = bb[0] | bb[1]
                yield dm, a, b, 511 ^ (dm|a), 511 ^ (dm|b), 511 ^ (dm|a|b)


def union(domains, cells):
    result = 0
    for c in cells:
        result |= domains[c]
    return result


def blockers(domains, center):
    """Untrusted generator: the separate acceptance loop rechecks the cover."""
    sigs, options = list(signatures(domains[center])), []
    for atom, cells in enumerate(atoms(center)):
        if atom == 0:
            continue
        for size in range(1, len(cells)+1):
            for selected in itertools.combinations(cells, size):
                neighbors = union(domains, selected)
                covered = sum(1 << i for i,sig in enumerate(sigs)
                              if count(neighbors & sig[atom]) < size)
                if covered:
                    options.append((atom, selected, neighbors, covered))
    remaining, chosen = (1 << len(sigs))-1, []
    while remaining:
        best = min(options, key=lambda item: (-count(item[3]&remaining),item[0],item[1]))
        atom, cells, neighbors, covered = best
        assert covered & remaining
        chosen.append({"atom":atom,"cells":list(cells),"neighbor_union":neighbors})
        remaining &= ~covered
    return {"signature_count":len(sigs),"blockers":chosen}


def verify_node(initial, node, stats):
    assert set(node) <= {"domains","closure_trace","kind","sector","witness","center","cover","cell","branches"}
    valid_domains(node["domains"])
    valid_trace(node["closure_trace"])
    domains = replay_trace(initial, node["closure_trace"])
    assert domains == node["domains"]
    stats["nodes"] += 1
    stats["closure_events"] += len(node["closure_trace"])
    if node["kind"] == "split":
        cell = node["cell"]
        integer(cell,0,80)
        assert count(domains[cell]) > 1
        assert close(domains) == (domains, [])
        expected = [d for d in range(1,10)if domains[cell]&(1 << (d-1))]
        assert [b["digit"]for b in node["branches"]] == expected
        depths = []
        for branch in node["branches"]:
            assert set(branch) == {"digit","node"}
            integer(branch["digit"],1,9)
            assumed = domains.copy()
            assumed[cell] = 1 << (branch["digit"]-1)
            depths.append(verify_node(assumed,branch["node"],stats))
        stats["splits"] += 1
        return 1 + max(depths)
    if node["kind"] == "hall-contradiction":
        integer(node["sector"],0,26)
        selected = node["witness"]["cells"]
        assert type(selected) is list
        for cell in selected:
            integer(cell,0,80)
        integer(node["witness"]["neighbor_union"],0,511)
        assert selected and len(selected) == len(set(selected))
        assert set(selected) <= set(HOUSES[node["sector"]])
        neighbors = union(domains, selected)
        assert neighbors == node["witness"]["neighbor_union"] and count(neighbors)<len(selected)
        stats["hall_leaves"] += 1
        return 0
    assert node["kind"] == "three-house-contradiction"
    center, cover = node["center"], node["cover"]
    integer(center,0,80)
    integer(cover["signature_count"],1,3780)
    partition = atoms(center)
    assert list(map(len,partition)) == [1,2,2,6,6,4]
    for blocker in cover["blockers"]:
        atom, cells = blocker["atom"], blocker["cells"]
        integer(atom,0,5)
        integer(blocker["neighbor_union"],0,511)
        assert type(cells) is list
        for cell in cells:
            integer(cell,0,80)
        assert cells and len(cells) == len(set(cells))
        assert set(cells) <= set(partition[atom])
        assert union(domains,cells) == blocker["neighbor_union"]
    checked = 0
    for signature in signatures(domains[center]):
        assert any(count(b["neighbor_union"] & signature[b["atom"]]) < len(b["cells"])
                   for b in cover["blockers"])
        checked += 1
    assert checked == cover["signature_count"] and checked > 0
    stats["three_house_leaves"] += 1
    stats["hall_blockers"] += len(cover["blockers"])
    stats["boundary_signatures_covered"] += checked
    return 0


def verify(fixture, cert):
    case = json.loads(fixture.read_text())
    integer(cert["schema_version"],1,1)
    assert cert["fixture_sha256"] == hashlib.sha256(fixture.read_bytes()).hexdigest()
    assert type(cert["target"]) is dict and set(cert["target"]) == {"cell","digit","value"}
    integer(cert["target"]["cell"],0,80)
    integer(cert["target"]["digit"],1,9)
    assert type(cert["target"]["value"]) is bool
    assert cert["target"] == case["target"] == {"cell":7,"digit":3,"value":True}
    assert case["conclusion"]["cell"] == 7 and case["conclusion"]["digit"] == 3 and case["conclusion"]["value"] is False
    base = case["base_domains"]
    valid_domains(base)
    assert base == case["state"]["domains"] and close(base) == (base, [])
    assert cert["base_state_hash"] == digest(base)
    raw = base.copy()
    raw[7] = 4
    opposite = replay_trace(raw,case["assumption_trace"])
    valid_domains(opposite)
    assert opposite == case["assumed_domains"] and close(opposite) == (opposite, [])
    assert cert["closed_opposite_state_hash"] == digest(opposite)
    stats = {k:0 for k in ["nodes","splits","hall_leaves","three_house_leaves","closure_events","hall_blockers","boundary_signatures_covered"]}
    depth = verify_node(opposite,cert["proof"],stats)
    assert depth == cert["additional_split_depth"] == 2
    audit = cert["local_only_depth1_survivors"]
    expected = [c for c,m in enumerate(opposite) if m&(m-1)]
    assert [r["pivot"]for r in audit] == expected
    for row in audit:
        cell,digit = row["pivot"],row["surviving_digit"]
        integer(cell,0,80)
        integer(digit,1,9)
        assert opposite[cell] & (1 << (digit-1))
        assumed = opposite.copy()
        assumed[cell] = 1 << (digit-1)
        final,_ = close(assumed)
        assert digest(final) == row["closed_state_hash"]
    assert len(audit) == cert["exhaustively_covered_depth1_pivots"] == 32
    checked = {**stats,"additional_split_depth":depth,"local_only_minimum_depth":2,
            "depth1_pivots_with_a_checked_survivor":len(audit),"exact_verified":True,
            "scope":"Minimum for this specified Hall/all81 local case-split calculus only; no historical novelty and no minimum over arbitrary known lemmas or rules claimed."}
    if "verification" in cert:
        assert cert["verification"] == checked
    return checked


def mutations(fixture, cert):
    """Reject malicious/accidental changes to all proof-critical components."""
    examples = []
    def changed(name, action):
        altered = deepcopy(cert)
        action(altered)
        examples.append((name,altered))
    changed("wrong fixture binding",lambda c:c.update(fixture_sha256="0"*64))
    changed("wrong target digit",lambda c:c["target"].update(digit=4))
    changed("wrong target polarity",lambda c:c["target"].update(value=False))
    changed("wrong target cell",lambda c:c["target"].update(cell=8))
    changed("wrong base state",lambda c:c.update(base_state_hash="0"*64))
    changed("wrong opposite state",lambda c:c.update(closed_opposite_state_hash="0"*64))
    changed("omitted branch",lambda c:c["proof"]["branches"].pop())
    changed("duplicated branch",lambda c:c["proof"]["branches"].append(deepcopy(c["proof"]["branches"][0])))
    changed("wrong branch digit",lambda c:c["proof"]["branches"][0].update(digit=1))
    changed("wrong split pivot",lambda c:c["proof"].update(cell=0))
    changed("omitted sound propagation prefix",lambda c:c["proof"]["branches"][0]["node"]["closure_trace"].pop())
    changed("forged propagation edit",lambda c:c["proof"]["branches"][0]["node"]["closure_trace"][0]["edits"][0].update(after=511))
    def leaves(node):
        if node["kind"] != "split":
            yield node
        for b in node.get("branches",[]):
            yield from leaves(b["node"])
    def hall(c):
        return next(n for n in leaves(c["proof"])if n["kind"]=="hall-contradiction")
    def triple(c):
        return next(n for n in leaves(c["proof"])if n["kind"]=="three-house-contradiction")
    changed("wrong Hall union",lambda c:hall(c)["witness"].update(neighbor_union=511))
    changed("empty Hall subset",lambda c:hall(c)["witness"].update(cells=[]))
    changed("wrong Hall house",lambda c:hall(c).update(sector=0))
    changed("omitted triple blocker",lambda c:triple(c)["cover"]["blockers"].pop())
    changed("wrong triple blocker union",lambda c:triple(c)["cover"]["blockers"][0].update(neighbor_union=511))
    changed("empty triple blocker cells",lambda c:triple(c)["cover"]["blockers"][0].update(cells=[]))
    changed("wrong triple signature count",lambda c:triple(c)["cover"].update(signature_count=0))
    changed("omitted minimum-depth pivot",lambda c:c["local_only_depth1_survivors"].pop())
    changed("wrong surviving digit",lambda c:c["local_only_depth1_survivors"][0].update(surviving_digit=1))
    changed("wrong survivor state",lambda c:c["local_only_depth1_survivors"][0].update(closed_state_hash="0"*64))
    changed("wrong claimed depth",lambda c:c.update(additional_split_depth=1))
    changed("wrong claimed node count",lambda c:c["verification"].update(nodes=1))
    for value in [81,-1,True]:
        changed("out-of-range/type triple center "+str(value),lambda c,v=value:triple(c).update(center=v))
    for value in [27,-26,True]:
        changed("out-of-range/type Hall sector "+str(value),lambda c,v=value:hall(c).update(sector=v))
    for value in [6,-1,True]:
        changed("out-of-range/type blocker atom "+str(value),lambda c,v=value:triple(c)["cover"]["blockers"][0].update(atom=v))
    for value in [81,-35,True]:
        changed("out-of-range/type split cell "+str(value),lambda c,v=value:c["proof"].update(cell=v))
    for value in [0,10,True]:
        changed("out-of-range/type branch digit "+str(value),lambda c,v=value:c["proof"]["branches"][0].update(digit=v))
    changed("boolean domain alias",lambda c:c["proof"]["domains"].__setitem__(0,True))
    changed("boolean propagation cell alias",lambda c:c["proof"]["branches"][0]["node"]["closure_trace"][0]["edits"][0].update(cell=True))
    changed("integer target polarity alias",lambda c:c["target"].update(value=1))
    for name, altered in examples:
        try:
            verify(fixture,altered)
        except (AssertionError,KeyError,ValueError,IndexError):
            continue
        raise AssertionError("Mutation accepted: "+name)
    return {"mutations_rejected":len(examples),"names":[name for name,_ in examples]}


def build(fixture, result_path, output):
    case = json.loads(fixture.read_text())
    result = json.loads(result_path.read_text())
    assert result["status"] == "refuted" and result["depth"] == 2 and result["known_lemmas"] is False
    proof = deepcopy(result["proof"])
    def add_cover(node):
        if node["kind"] == "three-house-contradiction":
            node["cover"] = blockers(node["domains"],node["center"])
        for branch in node.get("branches",[]):
            add_cover(branch["node"])
    add_cover(proof)
    opposite, survivors = case["assumed_domains"], []
    for cell,mask in enumerate(opposite):
        if mask & (mask-1) == 0:
            continue
        for digit in range(1,10):
            if not mask & (1 << (digit-1)):
                continue
            assumed = opposite.copy()
            assumed[cell] = 1 << (digit-1)
            try:
                final,_ = close(assumed)
            except ValueError:
                continue
            survivors.append({"pivot":cell,"surviving_digit":digit,"closed_state_hash":digest(final)})
            break
        else:
            raise AssertionError("A depth1 refutation exists; claimed minimum is false")
    cert = {"schema_version":1,"fixture_sha256":hashlib.sha256(fixture.read_bytes()).hexdigest(),
            "base_state_hash":digest(case["base_domains"]),"closed_opposite_state_hash":digest(opposite),
            "target":case["target"],"additional_split_depth":2,"proof":proof,
            "exhaustively_covered_depth1_pivots":len(survivors),"local_only_depth1_survivors":survivors,
            "calculus":"Native Hall GAC and all81 crossing-three-house GAC, closed to fixedpoint after every assumption; complete cell-value splits only. No uniqueness or named-pattern lemma used in the refutation."}
    cert["verification"] = verify(fixture,cert)
    output.write_text(json.dumps(cert,separators=(",",":"))+"\n")
    return cert["verification"]


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--fixture",type=Path,default=DEFAULT_FIXTURE)
    p.add_argument("--certificate",type=Path,default=DEFAULT_CERTIFICATE)
    p.add_argument("--build-from",type=Path)
    p.add_argument("--mutations",action="store_true")
    a = p.parse_args()
    checked = build(a.fixture,a.build_from,a.certificate) if a.build_from else verify(a.fixture,json.loads(a.certificate.read_text()))
    if a.mutations:
        checked = {**checked,**mutations(a.fixture,json.loads(a.certificate.read_text()))}
    print(json.dumps(checked,sort_keys=True))


if __name__ == "__main__":
    if not __debug__:
        raise SystemExit("Run without -O")
    main()
