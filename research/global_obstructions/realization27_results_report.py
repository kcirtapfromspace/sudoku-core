"""Publish compact realization27 discovery records and exact rejected controls."""
import hashlib
import json
from pathlib import Path

from real_four_digit_mine import initial_domains
from real_verify import sees

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent.parent


def main():
    folder=ROOT/"target/realization27-completions"
    summary=json.loads((folder/"summary.json").read_text())
    trials=[json.loads(line) for line in (folder/"trials.jsonl").read_text().splitlines()]
    assert all(r["status"]=="maximal-outside-nonunique" and not r["count_cutoff"] for r in trials)
    rows=[[r["seed"],r["completion"],r["solutions_capped_at_2"]] for r in trials]
    data="".join(json.dumps(row,separators=(",",":"))+"\n" for row in rows).encode()
    (HERE/"realization27_trials.jsonl").write_bytes(data)
    source_path=HERE/"nonrank_graph_native27_candidate.json"
    source=json.loads(source_path.read_text())
    core=set(source["graph_cells"]);palette=set(source["palette"])
    guardian=source["target"]["cell"];extra=source["target"]["digit"]
    solution=source["completion"]
    allowed=[cell for cell in range(81) if cell not in core and all(
        int(solution[cell]) not in (palette|({extra} if peer==guardian else set()))
        for peer in core if sees(cell,peer))]
    puzzle="".join(solution[cell] if cell in allowed else "." for cell in range(81))
    initial=initial_domains([int(d) if d!="." else 0 for d in puzzle])
    swapped=solution.translate(str.maketrans({"2":"3","3":"2"}))
    stored=json.loads((folder/"stored-control.json").read_text())
    controls={"schema_version":1,"core":sorted(core),"palette":sorted(palette),
        "guardian":guardian,"extra":extra,
        "exact_initial_mask_realization":{"puzzle":puzzle,"givens":len(allowed),
            "initial_domains":initial,"support_domains":[initial[c] for c in sorted(core)],
            "respecting_completions":[solution,swapped],"solution_orbit_lower_bound":24,
            "status":"Exact full initial palette masks realized, but necessarily nonunique; no old-engine gap or LP acceptance claimed."},
        "maximal_outside_stored_completion":{"puzzle":stored["puzzle"],"givens":54,
            "respecting_completions":stored["independent_completions"],
            "reference_completion":solution,
            "status":"Two completions agree on every outside clue. Therefore every outside-clue subset of this completion is nonunique."}}
    controls_data=(json.dumps(controls,indent=2)+"\n").encode()
    (HERE/"realization27_results_controls.json").write_bytes(controls_data)
    generated_seeds={r[0] for r in rows if r[0]}
    summary.update(public_manifest={"file":"realization27_trials.jsonl","records":len(rows),
        "bytes":len(data),"sha256":hashlib.sha256(data).hexdigest(),
        "array_fields":["seed_zero_means_stored_completion","compatible_completion","solutions_capped_at_2"],
        "meaning":"Recorded discovery classifications. The checker validates completion legality and counts; it does not independently recount all4999 maximum-clue puzzles."},
        duplicate_generation_seeds=sorted(set(range(summary["seed_start"],summary["seed_start"]+summary["requested_samples"]))-generated_seeds),
        exact_controls={"file":"realization27_results_controls.json","sha256":hashlib.sha256(controls_data).hexdigest()},
        all_unclued_unique_tier="No candidate reached uniqueness; hence full-older-engine comparison was not eligible and no stall is claimed.",
        exact_full_palette_theorem="Every outside cell sees a support cell. A palette given would remove a required support candidate, and monotone earlier deletions cannot restore it. Thus all givens omit allfourpalette symbols; their24 relabelings preserve givens and disprove uniqueness for any completion.",
        stronger_support_trade_result="The separate structural exhaustive certificate proves no unique compatible completion exists for outside-only clues even when support lists shrink. Sampled negatives alone are not this universal proof.")
    (HERE/"realization27_results.json").write_text(json.dumps(summary,indent=2)+"\n")
    print(json.dumps({"records":len(rows),"manifest_bytes":len(data),
                      "exact_mask_control_givens":len(allowed),"duplicate_seeds":summary["duplicate_generation_seeds"]}))


if __name__=="__main__":
    main()
