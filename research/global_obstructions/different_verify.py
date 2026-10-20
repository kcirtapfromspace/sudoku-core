"""Exact acceptance of the unclassified lead and known-composition control.

The actual engine certificate is valid; historical novelty is explicitly
rejected by a separate known-composition audit. No SAT/LP optimizer is trusted
or required here. The historical dispatcher is freshly run by a separate job.
"""
import hashlib
import itertools
import json
from pathlib import Path

from audit_global_certificates import audit_mutations, verify as verify_tree
from nonrank_audit_template import verify as verify_template
from real_four_digit_mine import initial_domains
from real_verify import close, exact_completions, replay_trace, sees, triple_feasible

HERE = Path(__file__).resolve().parent


def verify_fixture(stem, cell, digit, value, template_file, novelty):
    path = HERE / (stem + ".json")
    case = json.loads(path.read_text())
    source = next(entry for entry in json.loads((HERE / "expanded_modern_corpus.json").read_text())
                  if entry["id"] == case["name"])
    assert case["source"] == source and case["puzzle"] == source["puzzle"]
    assert exact_completions(case["puzzle"], 2) == [case["solution"]]
    solution = list(map(int, case["solution"]))
    initial = initial_domains([int(c) if c in "123456789" else 0 for c in case["puzzle"]])
    assert initial == case["initial_domains"]
    domains = replay_trace(initial, case["base_trace"], solution)
    assert domains == case["base_domains"] == case["state"]["domains"]
    assert close(domains)[1] == []
    assert case["target"] == {"cell": cell, "digit": digit, "value": not value}
    assert case["conclusion"] == {"cell": cell, "digit": digit, "value": value,
                                  "kind": "placement" if value else "elimination"}
    assert (solution[cell] == digit) == value and case["state"]["values"][cell] == 0
    raw = domains.copy()
    bit = 1 << (digit - 1)
    raw[cell] = raw[cell] & ~bit if value else bit
    branch = replay_trace(raw, case["assumption_trace"])
    assert branch == case["assumed_domains"] == case["assumed_state"]["domains"]
    assert close(branch)[1] == []
    for key in ("state", "assumed_state"):
        state = case[key]
        assert state["version"] == 1 and len(state["values"]) == 81
        assert all(v == (m.bit_length() if m & (m-1) == 0 else 0)
                   for v, m in zip(state["values"], state["domains"]))
    assert case["graph_edges"] == [[a, b] for a, b in itertools.combinations(case["graph_cells"], 2)
                                    if sees(a, b)]
    assert all(triple_feasible(branch, hs) for hs in itertools.combinations(range(27), 3))
    template = verify_template(case, json.loads((HERE / template_file).read_text()))
    replay = json.loads((HERE / ("replay_" + stem + ".json")).read_text())
    assert replay["state"] == case["state"] and replay["proof"]["cells"] == case["graph_cells"]
    assert {k: replay["proof"][k] for k in ("cell", "digit", "value")} == {
        "cell": cell, "digit": digit, "value": value}
    tree = verify_tree(replay)
    mutations = audit_mutations(replay)
    audit = json.loads((HERE / (stem + "_old_audit.json")).read_text())
    assert audit["fixture_sha256"] == hashlib.sha256(path.read_bytes()).hexdigest()
    assert audit["baseline_commit"] == "0532d63d9241a69c1fb57e9c34a1533effe1981d"
    assert audit["state_sha256"] == hashlib.sha256(json.dumps(case["state"], sort_keys=True,
                                                            separators=(",", ":")).encode()).hexdigest()
    assert len(audit["full_older_engine_confirmations"]) == 2
    for old in audit["full_older_engine_confirmations"]:
        assert old["state"] == case["state"] and not old["trace"]
        assert old["full_dispatcher_executed"] is True and old["full_current_finding"] is None
        assert old["arithmetic"]["error"] is None and old["arithmetic"]["hint"] is None
        assert old["arithmetic"]["budget_exhausted"] is False
        assert old["three_house"]["error"] is None and old["three_house"]["hint"] is None
        assert old["three_house"]["centers_checked"] == 81
        assert old["three_house"]["exhausted"] is True
    return {"fixture": stem, "unique_completion": True, "exact_givens_trace": True,
                      "retained_forcing_inclusive_stall_calls": 2, "arbitrary_triples": 2925,
                      "template_check": template, "production_certificate": tree,
                      "rejected_certificate_mutations": mutations,
                      "novelty": novelty}


def main():
    for arguments in [
        ("different_control", 73, 5, True, "different_template_primal.json",
         "rejected; known conditional tridagon control"),
        ("different_lead", 7, 3, False, "different_lead_template_primal.json",
         "unclassified; bounded known-method no-hit is not originality evidence"),
    ]:
        print(json.dumps(verify_fixture(*arguments)))


if __name__ == "__main__":
    if not __debug__:
        raise SystemExit("Run without -O")
    main()
