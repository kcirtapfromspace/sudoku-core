"""Check a published EL12c38 OR1 step followed by the actual target refutation."""
import argparse
from copy import deepcopy
import hashlib
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
from catalog630_verify import native_map
from catalog630_match import proper_coloring
from deeper_forcing_verify import blockers,digest,integer,valid_domains,valid_trace,verify_node
from real_verify import close,replay_trace

FIXTURE = HERE / "different_lead.json"
CERTIFICATE = HERE / "deeper_forcing_catalogue_chain.json"


def check_trace(initial,trace,expected):
    valid_domains(initial)
    valid_domains(expected)
    valid_trace(trace)
    final = replay_trace(initial,trace)
    assert final == expected and close(final) == (final,[])
    return final


def verify(cert):
    case = json.loads(FIXTURE.read_text())
    integer(cert["schema_version"],1,1)
    assert cert["fixture_sha256"] == hashlib.sha256(FIXTURE.read_bytes()).hexdigest()
    assert cert["target"] == case["target"] == {"cell":7,"digit":3,"value":True}
    assert type(cert["target"]["cell"]) is int and type(cert["target"]["digit"]) is int and type(cert["target"]["value"]) is bool
    assert cert["conclusion"] == {"cell":7,"digit":3,"value":False}
    assert type(cert["conclusion"]["cell"]) is int and type(cert["conclusion"]["digit"]) is int and type(cert["conclusion"]["value"]) is bool
    base = case["base_domains"]
    valid_domains(base)
    assert digest(base) == cert["base_state_hash"] and close(base) == (base,[])
    path = HERE / "catalog630_catalogue.json"
    assert cert["catalogue_sha256"] == hashlib.sha256(path.read_bytes()).hexdigest()
    catalogue = json.loads(path.read_text())
    source = next(r for r in catalogue["patterns"]if r["id"]==cert["catalogue_id"])
    assert source["id"] == "EL12c38" and source["sha256"] == cert["rule_source_sha256"]
    occurrence = cert["known_occurrence"]
    assert occurrence["id"] == source["id"]
    native_map(source["cells"],base,occurrence)
    # Independent finite proof of the published OR lemma's graph premise.
    assert proper_coloring(source["cells"],3) is None
    assert occurrence["guardians"] == [[68,7]]
    assert cert["known_placement"] == {"cell":68,"digit":7}
    assert type(cert["known_placement"]["cell"]) is int and type(cert["known_placement"]["digit"]) is int
    # All possible extras on the mapped non-three-colorable support form one
    # literal, so this placement follows BEFORE assuming the target endpoint.
    known_raw = base.copy()
    known_raw[68] = 64
    known = check_trace(known_raw,cert["base_known_trace"],cert["base_known_domains"])
    assert known[7] & 4 and known[7] & ~4
    target_raw = known.copy()
    target_raw[7] = 4
    target = check_trace(target_raw,cert["target_trace_after_known"],cert["known_closed_opposite"])
    # Also bind the alternate route through the existing public opposite state.
    opposite_raw = base.copy()
    opposite_raw[7] = 4
    opposite = check_trace(opposite_raw,case["assumption_trace"],case["assumed_domains"])
    assert digest(opposite) == cert["original_closed_opposite_hash"]
    opposite[68] = 64
    assert check_trace(opposite,cert["known_trace_after_opposite"],target) == target
    assert digest(target) == cert["known_closed_opposite_hash"]
    stats = {k:0 for k in ["nodes","splits","hall_leaves","three_house_leaves","closure_events","hall_blockers","boundary_signatures_covered"]}
    depth = verify_node(target,cert["proof"],stats)
    assert depth == 1
    result = {**stats,"known_catalogue_id":"EL12c38","known_placement_cell":68,"known_placement_digit":7,
              "base_known_closure_events":len(cert["base_known_trace"]),
              "target_after_known_closure_events":len(cert["target_trace_after_known"]),
              "opposite_then_known_closure_events":len(cert["known_trace_after_opposite"]),
              "additional_split_depth_after_known_rule":1,"actual_target_refuted":True,
              "classification":"published catalogue OR1 plus one complete cell split and exact Hall/all81 propagation"}
    if "verification"in cert:
        assert result == cert["verification"]
    return result


def build():
    case = json.loads(FIXTURE.read_text())
    found = json.loads(Path("target/deeper-forcing-catalogue-chain.json").read_text())
    assert found["status"] == "refuted"
    proof = deepcopy(found["proof"])
    def cover(node):
        if node["kind"] == "three-house-contradiction":
            node["cover"] = blockers(node["domains"],node["center"])
        for b in node.get("branches",[]):cover(b["node"])
    cover(proof)
    base_raw = case["base_domains"].copy()
    base_raw[68] = 64
    known,base_trace = close(base_raw)
    target_raw = known.copy()
    target_raw[7] = 4
    target,target_trace = close(target_raw)
    occurrence = json.loads((HERE/"catalog630_match_result.json").read_text())["runs"]["full_base"]["occurrences"][0]
    catalogue_path = HERE/"catalog630_catalogue.json"
    source = next(r for r in json.loads(catalogue_path.read_text())["patterns"]if r["id"]=="EL12c38")
    cert = {"schema_version":1,"fixture_sha256":hashlib.sha256(FIXTURE.read_bytes()).hexdigest(),
            "base_state_hash":digest(case["base_domains"]),"target":case["target"],
            "conclusion":{"cell":7,"digit":3,"value":False},"catalogue_id":"EL12c38",
            "catalogue_sha256":hashlib.sha256(catalogue_path.read_bytes()).hexdigest(),
            "rule_source_sha256":source["sha256"],"known_occurrence":occurrence,
            "known_placement":{"cell":68,"digit":7},"base_known_domains":known,"base_known_trace":base_trace,
            "target_trace_after_known":target_trace,"known_closed_opposite":target,
            "known_closed_opposite_hash":digest(target),"original_closed_opposite_hash":digest(case["assumed_domains"]),
            "known_trace_after_opposite":found["known_propagation_trace"],"proof":proof}
    cert["verification"] = verify(cert)
    CERTIFICATE.write_text(json.dumps(cert,separators=(",",":"))+"\n")
    return cert


def mutations(cert):
    examples=[]
    def altered(name,action):
        c=deepcopy(cert);action(c);examples.append((name,c))
    altered("wrong rule ID",lambda c:c.update(catalogue_id="EL12c1"))
    altered("missing guardian",lambda c:c["known_occurrence"].update(guardians=[]))
    altered("wrong placement",lambda c:c["known_placement"].update(digit=8))
    altered("wrong palette",lambda c:c["known_occurrence"]["palette"].__setitem__(0,2))
    altered("wrong mapping",lambda c:c["known_occurrence"]["mapping"][0].__setitem__(1,0))
    altered("wrong target",lambda c:c["target"].update(digit=4))
    altered("wrong conclusion",lambda c:c["conclusion"].update(value=True))
    altered("missing base propagation",lambda c:c["base_known_trace"].pop())
    altered("missing target propagation",lambda c:c["target_trace_after_known"].pop())
    altered("missing alternate-route propagation",lambda c:c["known_trace_after_opposite"].pop())
    altered("wrong original opposite",lambda c:c.update(original_closed_opposite_hash="0"*64))
    altered("missing pivot branch",lambda c:c["proof"]["branches"].pop())
    altered("wrong Hall leaf",lambda c:c["proof"]["branches"][0]["node"]["witness"].update(neighbor_union=511))
    altered("missing triple blocker",lambda c:c["proof"]["branches"][1]["node"]["cover"].update(blockers=[]))
    for name,c in examples:
        try:verify(c)
        except(AssertionError,KeyError,ValueError,IndexError,StopIteration):continue
        raise AssertionError("Mutation accepted: "+name)
    return len(examples)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("--build",action="store_true");p.add_argument("--mutations",action="store_true");a=p.parse_args()
    cert=build()if a.build else json.loads(CERTIFICATE.read_text())
    result=verify(cert)
    if a.mutations:result["mutations_rejected"]=mutations(cert)
    print(json.dumps(result,sort_keys=True))


if __name__=="__main__":
    if not __debug__:raise SystemExit("Run without -O")
    main()
