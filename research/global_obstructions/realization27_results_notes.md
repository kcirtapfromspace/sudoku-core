# Givens realization of the exact native27 obstruction

The requested all27-unclued unique realization is impossible for this fixed
support, palette {2,3,5,9} and guardian r5c6=1. The separate
[`realization27_structure_report.md`](realization27_structure_report.md) and
[`realization27_structure_verify.py`](realization27_structure_verify.py) contain
the exhaustive universal support-trade proof: all1,320 labelled four-colorings
of the support minus the guardian fall into408 native-house-signature buckets,
each containing at least two colorings. Replacing one with another in its bucket
preserves every native house and every outside cell. Hence even all54 outside
clues cannot make any compatible completion unique. Native Sudoku symmetries
preserve this obstruction. Earlier sound candidate eliminations cannot discard
one of these actual completions.

There is also a simpler necessary obstruction to exact full-palette masks.
Every outside cell peers with a support cell. An outside palette given removes
one of that support cell's required four candidates, and later monotone domain
deletions cannot restore it. Thus exact full-palette masks permit no palette
givens. All24 permutations of the four palette symbols preserve the givens,
disproving uniqueness. The retained44-given control demonstrates that the exact
initial masks themselves are realizable; it is the uniqueness demand that fails.

## Bounded corroborating search

The independent randomized exact-cover generator made5,000 compatible-completion
attempts, plus the stored completion. Two samples repeated earlier completions;
the4,999 distinct full grids all had at least two solutions with every outside
cell fixed. There were no generation or uniqueness cutoffs. Runtime was16.085 s.
These finite observations corroborate the exhaustive structural proof; they are
not a substitute for it. No outside-only unique candidate reached the full-old
engine or LP comparison, so no old-engine stall is claimed for this tier.

`realization27_trials.jsonl` retains every distinct completion and seed.
`realization27_results.json` binds the public manifest hash and exact counts.
`realization27_results_controls.json` retains the44-given exact-mask control,
24-symbol-permutation orbit, and two explicit solutions of the maximal54-given
control. `realization27_results_verify.py` independently checks all completion
records for native legality and palette membership, plus the explicit rejected
controls. Its recorded discovery labels are not all independently recounted.

## Changed-premise internal-anchor controls

Permitting one selected non-guardian vertex to be a given changes the problem.
Of26 possible single anchors with all54 outside clues, four anchors—zero-based
cells66,69,74,75—make the stored completion unique. All four are solved by local
closure. Exact controls and completed old-dispatcher diagnostics are in
`realization27_results_internal.json`; none is classified as a live stall.

A separate bounded minimization made40 restarts (seeds92137–92176), testing2,160
clue removals and accepting1,232 while preserving uniqueness and the selected
palette-list restriction **after local closure**. There were zero count cutoffs
and zero live-guardian endpoints. Givens counts ranged21–26. Exact restart
records are in `realization27_results_anchored.json`.

`realization27_results_anchored_control.json` retains a24-given unique anchored
control and `realization27_results_anchored_minimal_control.json` a21-given one.
Their `initial_state` contains precisely the actual givens and ordinary native
candidates. Their original27 `graph_cells`, target and conclusion are retained.
Their initial lists include extra nonpalette candidates: containment in the
four-digit lists is obtained only after the exact recorded local closure, which
already places the guardian. Full older comparisons are solved outcomes, not
strict gaps. These controls neither realize the all27-live uniform state nor
justify replacing their raw candidate masks by the synthetic masks.

`realization27_results_raw_old.json` also runs the complete frozen older engine
directly from each ordinary givens-only initial state, without external Hall/81
closure or singleton promotion. Both solve in their first call. The24-given case
uses46 Naked Singles and11 Hidden Singles. The21-given case uses46 Naked Singles,
14 Hidden Singles,3 Naked Pairs,3 Naked Triples and4 Pointing Pairs. The full
dispatcher executes; its `None` finding is on the solved board, not a live
stall. The stdlib verifier binds inputs/calls and replays every endpoint against
the independently counted sole completion.

## Reproduce

```sh
# Stdlib exact controls and public-manifest validation.
python3 -S research/global_obstructions/realization27_results_verify.py --recount-sample 10

# Randomized compatible completions; all27 support cells excluded from givens.
python3 -S research/global_obstructions/realization27_search.py \
  --samples 5000 --seconds 360 --output-dir target/realization27-completions
python3 -S research/global_obstructions/realization27_results_report.py

# Clearly separate one-internal-given controls and bounded clue deletion.
python3 -S research/global_obstructions/realization27_search_internal.py
python3 -S research/global_obstructions/realization27_search_anchored.py --seconds 120 --restarts 40

# Direct full older-engine comparison; optional fresh CI binary/output paths.
python3 -S research/global_obstructions/realization27_search_raw_old.py \
  --binary target/nonrank-stdlib-baseline/instrumented/target/release/examples/real_existing_probe \
  --output target/realization27-raw-old-rerun.json
```

Completion discovery uses a small independent C++ exact-cover program. Numerical
library threads are capped at one; no floating LP acceptance is asserted here.
The selected anchored controls retain complete initial domains, local/older
traces and independent Python uniqueness checks.
