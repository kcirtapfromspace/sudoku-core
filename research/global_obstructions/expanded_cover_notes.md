# Expanded four-palette clue-cover search

This experiment searches generated givens selections for the existing thirteen-cell
four-palette native-peer graph. Its respecting completion is the completion of
Norvig top95 #20 already recorded in `real_four_digit_candidate.json`. The new
puzzles are generated selections, not publicly authored Norvig puzzles, and no
historical mathematical novelty is asserted.

The previous experiment added clues greedily to one fixed nineteen-given seed.
Here a mixed-integer optimizer can remove or replace every earlier clue. It
chooses up to eighty clue selectors alongside all 729 continuous native Sudoku
variables. Native LP equations, guardian-false restriction and actual-given
coverage of every forbidden core candidate must hold together. An alternate
complete Sudoku contributes a necessary hitting-set constraint: a unique clue
set must include at least one position distinguishing that completion from the
reference completion. Random positive clue weights vary the minimum-cost cover.
Unique endpoints undergo random uniqueness-preserving clue deletion, while
retaining the givens-only core restrictions. The final run restricts selectors
to the 68 cells outside the core and changes the random seed.

An independent C++ exact-cover enumerator supplies alternate completions and
bounded uniqueness checks for discovery. It was checked against the five earlier
published/generated controls. A one-million-node enumeration cutoff is reported
explicitly and never counted as uniqueness. Positive acceptance uses the separate
independent Python enumeration and exact rational LP checker.

Before every historical baseline call, full Hall and all 81 crossing-three-house
GAC are closed to a fixed point. The wrapper asserts that a singleton converted
to a placement has already been removed from every native peer. All old candidate
eliminations are retained across placements. The historical binary is pinned to
`0532d63d9241a69c1fb57e9c34a1533effe1981d`, cached under `target`, and checked against
the known #4401 and #4662 baseline results before each run.

Candidate tiers are recorded separately:

- Initial or post-old closure makes the opposite guardian branch locally
  inconsistent: rejected by the named strict local filters.
- Old nonforcing closure places the guardian: rejected as an incremental gap.
- Old nonforcing closure stalls with guardian live, default Arithmetic returns
  no hint, and the opposite branch remains locally and LP feasible: a useful
  lesser tier requiring separate exact certification.
- Full older dispatcher, including forcing, also stalls while the guardian
  remains live: the stronger target. A solved board is never classified a stall.

For the lesser tier, the experiment additionally follows old forcing hints,
alternating them with complete local and old nonforcing closure. The complete
defeating trace for the certified twenty-nine-given puzzle is reproduced by
`expanded_cover_forcing_trace.py`; every endpoint is checked against its
independently counted sole completion. This does not reconstruct each named old
technique's internal proof independently.

## Reproduce

Run from the repository root. SciPy supplies MILP and floating LP discovery;
`c++` builds the independent discovery enumerator. All numerical-library and
HiGHS threads are capped at one and the search uses one worker.

```sh
# Fresh broad search, with no dependency on retained target logs.
python3 research/global_obstructions/expanded_cover_search.py \
  --seed 98247 --outside-core --min-clues 26 --clue-span 8 \
  --seconds 360 --trials 1000 --milp-seconds 1 \
  --output-dir target/expanded-cover-fresh

# Independently completion-validate the complete defeating old forcing trace.
python3 research/global_obstructions/expanded_cover_forcing_trace.py

# Verify both retained negative controls without running the engine or MILP.
python3 research/global_obstructions/expanded_cover_verify.py

# Regenerate aggregate exact counts from retained discovery logs.
python3 research/global_obstructions/expanded_cover_report.py
```

`expanded_cover_summary.json` retains exact per-run counts, seeds, commands,
rejection categories and deduplicated puzzle counts. Raw discovery states,
alternate-solution learning logs and source snapshots stay under `target`.
Per-MILP wall-time bounds can change the search trajectory between machines.
Floating infeasibility screens and the bounded search do not establish that a
stronger unique puzzle cannot exist. Separately verified fixed fixtures and
their exact proof/LP checks do not depend on discovery timing.

The public compact manifests preserve the historical count populations without
requiring the untracked `target` logs. Each `expanded_cover_trials.jsonl` line is
`[puzzle, status]`, where 0 means discovery-nonunique, 1 discovery-unique, and −1
an enumeration cutoff. Each `expanded_cover_endpoints.jsonl` line is
`[puzzle, category, forcing]`: categories 0–4 mean old-target-placed,
joint-target-placed, initial-local-filter, post-old-local-filter and
nonforcing-plus-Arithmetic gap respectively. Forcing 0 means not attempted and
1 means the later forcing closure placed the guardian. The summary binds their
SHA-256 hashes, byte lengths and schemas. `expanded_cover_verify.py` checks
deduplication and every aggregate count from these public files using only the
standard library. These are recorded discovery classifications; this check does
not independently certify all 438 raw uniqueness claims or all 611 endpoints.

## Recorded outcome

Excluding the seven-set preliminary debugging pilot, the runs fully logged 771
raw clue sets, or 767 distinct givens strings. Of these, 438 distinct raw sets
were unique. Uniqueness-preserving clue deletion and requalification produced
611 distinct evaluated unique endpoints.

| Run | Seed | Fully logged raw sets | Unique raw sets started | Distinct minimized endpoints within run | Elapsed seconds |
| --- | ---: | ---: | ---: | ---: | ---: |
| b, sparse cover pilot | 62841 | 17 | 0 | — | 61.626 |
| c, whole-set MILP | 62841 | 113 | 41 | — | 65.608 |
| d, clue deletion added | 62841 | 120 | 62 | 102 | 110.295 |
| e, full forcing closure added | 62841 | 170 | 82 | 127 | 130.548 |
| f, outside-core selectors | 98247 | 351 | 254 | 342 | 360.373 |

Runs b–e were interrupted at deliberate phase changes; their summaries can
include a started candidate not yet fully logged. The table populations differ:
a raw unique set can yield two distinct minimized endpoints, and resumed cuts
and requalified fixtures do not count as new generated puzzles.

The 611 distinct endpoint classifications were: 425 placed the guardian through
old nonforcing methods, 15 placed it during joint local closure, 16 failed the
initial local branch filters, 64 failed the filters after old nonforcing closure,
and 91 passed the nonforcing-plus-Arithmetic discovery tier. There were no
conflicting classifications for an identical givens string.

All 91 lesser-tier endpoints were defeated by subsequent old forcing closure;
**no forcing-inclusive stall was found**. The twenty-nine-given example promoted
to `expanded_unique_nonforcing.json` separately passed independent uniqueness,
all local/triple checks and an exact half-integral native LP witness. Other
discovery-screen positives are not represented as exact-LP-certified fixtures.

`expanded_cover_old_forcing_trace_29.json` records five Nishio eliminations—
r4c6≠4, r5c1≠6, r1c3≠1, r2c1≠7, r1c2≠8—with complete Hall/three-house closure
between them, after which guardian r2c8=1 is placed. The independent Python
enumerator verifies the sole completion, every older endpoint preserves it, and
the entire recorded domain trace replays exactly from the original givens.

`expanded_cover_rejected_nonunique.json` preserves a thirty-two-given endpoint
that retains the nonforcing/local gap but has two independently verified
respecting completions. It is rejected by the uniqueness requirement. Its two
completions and earlier local trace are checked by `expanded_cover_verify.py`;
no exact-LP separation acceptance is claimed for this negative control.
