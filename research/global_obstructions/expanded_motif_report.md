# Expanded four-digit obstruction search

The strict requested result remains open: this search has not produced a
certified unique givens puzzle whose guardian remains live after the full older
dispatcher, including forcing. Synthetic graph discoveries and surviving
nonunique clue covers do not meet that requirement.

## Broader structural coverage

The graph search used the independently completed Norvig top95 puzzles 1–12,
rather than the original puzzle-20 completion. It retained 18 native support
embeddings in 14 graph-isomorphism classes, on 13–24 cells. Every graph differs
from the earlier 13-vertex/35-edge example. "Different" here concerns the graph
and embedding; it makes no historical originality claim.

The first fully completed search stage tested 209,276 support presentations on
eight completed grids. Of 89 non-four-colorable presentations, floating role-LP
screening rejected 75; graph-isomorphism filtering retained 12 embeddings. A
second interrupted stage added six embeddings from four further completions.
The second stage's counters are a saved checkpoint, not a completed universe.
Presentations can describe the same graph. No global classification or exhaustive
Sudoku search is claimed.

The exact classification subsequently proves that all 18 retained graphs have
chromatic number five and are vertex-critical: every one-vertex deletion is
four-colorable. Their respecting complete Sudokus give the five-color upper
bound. Native edges are reconstructed from row, column, and box coordinates.
All 18 synthetic candidate states also have exact rational witnesses satisfying
all 324 native linear equations after removing the guardian's extra candidate.
The SciPy proposal is accepted only after every equality, bound, domain
restriction, and target restriction passes `Fraction` checks.

Fourteen embeddings admit the uniform capacity argument `|S| > 4*alpha(S)`.
The other four initially looked stronger because their independence bounds do
not give that inequality. Exact fractional-coloring analysis rejects that
interpretation: **all four have fractional chromatic number 25/6**. They have
simple integer vertex weights, each either one or two, with total weight 25
and maximum independent-set weight six. Four digits can carry weight at most
24. These are classical weighted independent-set/rank arguments as well.

| Source completion | Cells / edges | Independence number | Cells with weight two; all others weight one |
| --- | --- | --- | --- |
| top95 2 | 24 / 66 | 6 | r3c4 |
| top95 3 | 22 / 61 | 6 | r4c1, r4c8, r5c8 |
| top95 10 | 23 / 67 | 6 | r1c1, r3c1 |
| top95 10 | 24 / 69 | 6 | r6c6 |

The compact fractional certificates include primal covers and dual weights.
The standard-library verifier reconstructs every maximal independent set, checks
all dual bounds, and checks primal coverage and equality of their exact totals.
Thus neither a floating optimum nor a supplied stable-set list is trusted.
This is useful falsification of a proposed stronger separation, not a new
mathematical principle. It is consistent with the older example's known
weighted odd-hole/hyperhole explanation in
[`prior_art_13_report.md`](prior_art_13_report.md).

## Givens and old-method search

For each retained support, the clue builder enumerates inclusion-minimal
per-digit peer covers for the five digits outside its palette. These create
the desired native candidate restrictions without imposing arbitrary masks.
Clue augmentation then tests three seeded orders: other digits first, palette
digits first, and a mixed order. The corrected run also permits palette clues
inside the selected support, while keeping the guardian unfilled. Such clues
may reduce a four-digit graph to a smaller list-coloring instance or a known
three-digit pattern; they would require that comparison before acceptance.

Every evaluation first screens the full native LP. The corrected protocol then
normalizes candidate domains by individual-house Hall propagation and all 81
crossing-three-house support propagators **before** converting singleton
domains into placed values for the frozen `0532d63` engine. It alternates this
local closure with older nonforcing handlers until both stop. Each reported
old step is checked against the respecting completion. The opposite guardian
assumption must then survive the complete local closure and full LP screen.
Only independently unique finalists can enable uniqueness handlers or the full
older dispatcher. A solved board or singleton guardian is always rejected as
a stalling endpoint.

Initial exploratory old-method counts in `expanded_motif_results.json` and
`expanded_motif_mixed_results.json` are explicitly superseded. The wrapper
initially converted raw singleton masks to placed values before propagating
their peer exclusions, violating assumptions of the older handlers. Apparent
bad Hidden Pair steps and probe panics from those inputs are invalid probe
premises, not demonstrated engine bugs or valid negative results. Graph
discovery and subsequent independent synthetic checks are unaffected.

The final normalized counters and strict-tier outcome are recorded in
[`expanded_motif_summary.json`](expanded_motif_summary.json). LP screening
failures remain floating discovery negatives, not exact infeasibility
certificates. Surviving clue states remain discovery-only until an independent
exact acceptance process validates every required claim.

The corrected run entered 16 of the 18 retained motifs within its 420-second
budget. It tested 2,153 clue presentations: 1,186 failed the raw LP screen;
967 passed that screen. After normalized old/local closure, 185 had their
guardian forced and 53 contradicted local propagation under the opposite
assumption. The remaining 729 discovery evaluations survived the local and LP
screens, but **none was unique**. Thus zero full older-dispatcher finalists
were tested. Ten saved endpoint examples have up to 57 givens and still have
at least two solutions. The corrected run had no execution failures or old
steps inconsistent with its respecting completion. These counts refer to
evaluations, which can revisit the same state in different clue orders.

## Reproduce

The fast exact synthetic proof audit needs only the standard library and takes
approximately 0.1 seconds:

```sh
python3 research/global_obstructions/expanded_motif_verify.py
```

Rebuild the synthetic classification, including exact checks of proposed role
witnesses; discovery dependencies are NumPy and SciPy:

```sh
python3 research/global_obstructions/expanded_motif_classify.py \
  research/global_obstructions/expanded_motif_results.json \
  research/global_obstructions/expanded_motif_mixed_results.json --exact-role-lp
```

Reproduce the bounded normalized clue search using the frozen old-engine probe:

```sh
python3 research/global_obstructions/expanded_motif_search.py \
  --seconds 420 --covers 64 --augment-seeds 3 --allow-core-clues \
  --motif-input research/global_obstructions/expanded_motif_results.json,research/global_obstructions/expanded_motif_mixed_results.json \
  --seed 2198343 --output research/global_obstructions/expanded_motif_normalized_results.json \
  --binary /tmp/global-real-engine-probe/target/release/examples/real_existing_probe
```

Build that probe with `real_existing_check.build_engine` if it is unavailable;
its source archive is pinned to commit `0532d63`. Runtime-bounded counters can
vary with machine speed, and old detector tie choices can vary. Bulk raw
candidate states and quarantines belong under `target/expanded-motif-raw`, not
the committed certificate corpus. The retained summaries record seeds, bounds,
support coordinates, source completions, and exact mathematical witnesses.

The graph search can be rerun without clue experiments using `--covers 0`:

```sh
python3 research/global_obstructions/expanded_motif_search.py \
  --seconds 720 --covers 0 --augment-seeds 0 --indices 1,2,3,4,5,6,7,8 \
  --seed 2198341 --output target/expanded-motif-raw/graph-reproduction.json
```

## References and next experiment

- [Norvig's author-hosted top95 corpus](https://www.norvig.com/top95.txt) supplies
  the original puzzles whose completions anchor this generated search; the
  newly generated clue covers are not attributed as authored Norvig puzzles.
- [Letchford, Rossi and Smriglio, *The stable set problem: Clique and nodal inequalities revisited*, §2.1](https://www.lancaster.ac.uk/staff/letchfoa/articles/2020-nodal.pdf)
  describes the established stable-set rank inequality family.

A focused next experiment should preserve three independent palette anchors
early, instead of greedily filling the five other digit classes first. Many
maximal nonunique survivors still leave two palette digits interchangeable.
Clue-set backtracking or counterexample-guided uniqueness cuts can search this
tradeoff more deliberately, always keeping exact local and LP acceptance ahead
of any novelty or engine-strength claim.
