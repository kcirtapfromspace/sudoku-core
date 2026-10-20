"""Export bounded corpus outcomes and a checked known-pattern classification."""
from collections import Counter
import hashlib
import json
from pathlib import Path

from real_verify import graph_colorings, sees

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent


def minimize(domains, cells):
    result = cells.copy()
    for cell in cells:
        without = [other for other in result if other != cell]
        if graph_colorings(domains, without)[0] == 0:
            result = without
    return result


def isomorphism(first, second):
    def graph(cells):
        return [[int(i != j and sees(a,b)) for j,b in enumerate(cells)]
                for i,a in enumerate(cells)]
    a,b = graph(first),graph(second)
    da,db = list(map(sum,a)),list(map(sum,b))
    mapping, used = {}, set()
    def visit():
        if len(mapping) == len(first):
            return True
        x = max((i for i in range(len(first)) if i not in mapping),
                key=lambda i:(sum(a[i][j] for j in mapping),da[i]))
        for y in range(len(second)):
            if y in used or da[x] != db[y]:
                continue
            if any(a[x][px] != b[y][py] for px,py in mapping.items()):
                continue
            mapping[x] = y; used.add(y)
            if visit():
                return True
            del mapping[x]; used.remove(y)
        return False
    assert visit()
    return [[first[x],second[y]] for x,y in sorted(mapping.items())]


def main():
    fixture = HERE / "nonrank_full_old_stall.json"
    case = json.loads(fixture.read_text())
    raw = case["base_domains"].copy()
    raw[35] &= ~256
    raw_core = minimize(raw,case["graph_cells"])
    closed_core = minimize(case["assumed_domains"],case["graph_cells"])
    known = json.loads((HERE / "real_Loki.json").read_text())
    mapping = isomorphism(closed_core,known["graph_cells"])
    classification = {"schema_version":1,
        "fixture_sha256":hashlib.sha256(fixture.read_bytes()).hexdigest(),
        "classification":"Four-digit list obstruction reducing, after permitted local branch closure, to an established tridagon/trivalue-oddagon graph.",
        "historical_novelty":"Not claimed. The twelve-cell graph is isomorphic to the published Loki tridagon control.",
        "direct_production_support":case["graph_cells"], "raw_list_core":raw_core,
        "raw_list_core_domains":[raw[c] for c in raw_core],
        "raw_refutation_nodes":graph_colorings(raw,raw_core)[1],
        "closed_tridagon_core":closed_core,
        "closed_tridagon_digits":[4,7,8],
        "closed_refutation_nodes":graph_colorings(case["assumed_domains"],closed_core)[1],
        "isomorphism_to_published_Loki_control":mapping,
        "known_control":"research/global_obstructions/real_Loki.json",
        "known_control_primary_source":known["source"],
        "primary_general_reference":"https://github.com/denis-berthier/CSP-Rules-Examples#the-tridagons-folder",
        "scope_note":"Uniform four-colorability is not enough to reject this candidate. Its actual domain lists matter; the old uniform chi5 corpus filter missed this type."}
    (HERE / "nonrank_puzzles_first_classification.json").write_text(json.dumps(classification,indent=2)+"\n")

    source = ROOT / "target/nonrank-puzzles-modern200/summary.json"
    summary = json.loads(source.read_text())
    status_codes = {"full-old-stall":0,"solved":1,"probe-timeout":2,
                    "invalid-or-error":3,"step-limit":4,"not-unique":5}
    rows = [[r["id"],r["puzzle"],status_codes[r["status"]],r.get("old_findings",0),
             r.get("unplaced",0),r.get("domains_sha256","")] for r in summary["records"]]
    data = "".join(json.dumps(row,separators=(",",":"))+"\n" for row in rows).encode()
    (HERE / "nonrank_puzzles_trials.jsonl").write_bytes(data)
    compact = {"schema_version":1,"baseline_commit":summary["baseline_commit"],
               "source_corpus":"research/global_obstructions/expanded_modern_corpus.json",
               "source_corpus_sha256":summary["corpus_sha256"] ,
               "source_order":"hard-first: Monhard TE3 then Mith TE3 then Monhard TE2",
               "completed_cases":len(rows),"status_counts":dict(Counter(r["status"] for r in summary["records"])),
               "seconds":summary["seconds"],"requested_limit":200,
               "stopping_reason":"Stopped after the first qualifying unique full-stall/template-LP candidate was independently accepted; no all200 result is claimed.",
               "command":summary["command"],"records":summary["records"],
               "manifest":{"file":"nonrank_puzzles_trials.jsonl","sha256":hashlib.sha256(data).hexdigest(),
                   "records":len(rows),"bytes":len(data),
                   "array_fields":["id","puzzle","status","old_findings","unplaced","domains_sha256"],
                   "status_codes":{str(value):key for key,value in status_codes.items()}},
               "classification_scope":"Bounded engine-comparison outcomes. Full traces retained under target; the promoted first fixture has separate exact acceptance and fresh-build audit."}
    (HERE / "nonrank_puzzles_results.json").write_text(json.dumps(compact,indent=2)+"\n")
    print(json.dumps({"completed_cases":len(rows),"statuses":compact["status_counts"],
                      "raw_core":len(raw_core),"closed_core":len(closed_core),"isomorphism":True}))


if __name__ == "__main__":
    if not __debug__:
        raise SystemExit("Run without -O")
    main()
