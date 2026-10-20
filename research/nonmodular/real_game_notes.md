# Givens-derived three-house experiment

This is a finite experiment on actual publicly supplied givens, not synthetic
candidate masks. Sources are Peter Norvig's author-hosted curated
[95 hard puzzles](https://www.norvig.com/top95.txt) and
[11 hardest puzzles](https://www.norvig.com/hardest.txt), described in his
[Sudoku essay](https://www.norvig.com/sudoku.html). These are curated collections;
individual puzzle authors are not attributed in these files. The original
givens, source URL, one-based source index and source-byte SHA-256 are retained
in `real_game_corpus.json` and the positive fixtures.

The recorded run covered 106 distinct givens strings. A separate Python MRV
enumerator, written for this experiment, found exactly one solution for each
when counting up to two. No recorded solution was used to derive candidates.

The baseline uses the existing engine's non-forcing detectors, excluding the
new Three-House Matching technique and the arithmetic family. Uniqueness
techniques are permitted because uniqueness was independently established.
Every placement and elimination was checked against the independent completion.
Candidate eliminations are retained across placements: after the engine
recalculates candidates, its result is intersected with the previous masks.
This avoids accidentally resurrecting earlier eliminations. No forcing-chain
or backtracking deduction enters the baseline state.

| Recorded result | Count |
| --- | ---: |
| Independently unique puzzles | 106 |
| Solved completely by the baseline | 74 |
| Stalled at a baseline fixed point | 32 |
| Centers scanned in stalled states | 2,592 |
| Stalled states with a three-house deduction | 8 |
| Stalled states with no three-house deduction | 24 |
| Positive centers | 17 |
| Distinct puzzle/cell/digit eliminations | 15 |
| Center/cell/digit elimination occurrences | 25 |

| Source puzzle, one-based | Distinct eliminations |
| --- | --- |
| top95 #20 | r9c2 ≠ 3 |
| top95 #35 | r6c1 ≠ 8 |
| top95 #52 | r4c3 ≠ 8; r7c2 ≠ 9 |
| top95 #59 | r2c4 ≠ 7; r3c4 ≠ 7; r3c6 ≠ 8; r5c6 ≠ 7; r6c6 ≠ 7 |
| top95 #65 | r8c8 ≠ 8 |
| top95 #66 | r4c3 ≠ 1,9; r6c9 ≠ 3 |
| top95 #88 | r3c4 ≠ 4 |
| top95 #93 | r1c7 ≠ 1 |

All eight exact positive states and their 207 baseline steps are saved in
`real_game_cases.json`. The independent Python checker reconstructs candidates
from givens, replays every step to the exact stored masks and values, independently
confirms uniqueness again, and recomputes all 81 centers for each positive state.
Existing detector choices can vary between runs because some existing code uses
randomized collections; the saved traces make these particular states exact.
`real_game_summary.json` records the baseline orchestrator source hash and every puzzle's
solved/stalled/positive status. Full raw traces for the recorded 106-puzzle run
are intentionally outside the repository; running the probe regenerates them in
the system temporary directory.

The public Rust API was then tested on all eight exact snapshots. All 25
center/candidate deductions produced valid independent Hall certificates.
For every snapshot, the first hint preserved its input, its target disagreed
with the independently verified completion, and a JSON round trip through
`ThreeHouseReplay` verified and applied exactly the recorded elimination. A
second application of the stale replay was rejected. The current restricted
baseline was also rechecked on all eight states, with no finding.

Default Arithmetic Counting returned **no hint on all eight states**, after
7,190; 10,015; 7,531; 10,780; 7,316; 8,385; 7,503; and 8,046 combinations,
respectively. Every arithmetic search reported beam pruning, no combination
budget exhaustion, and no error. These are bounded-search results, not proof
that arithmetic or other known methods cannot make the deductions.

The full recorded run took 48.37 seconds including temporary compilation:
24.41 seconds for independent uniqueness enumeration, 5.25 seconds for existing
engine closure and 3.89 seconds for the Python three-house scans. On the eight
positive snapshots, the latest optimized Rust first-hint search took roughly
0.29–2.33 milliseconds. These measurements describe one local run, not a
portable benchmark or a claim that a puzzle is solved by the new technique.

Run from the repository root:

```sh
python3 research/nonmodular/real_game_verify.py
python3 research/nonmodular/real_game_engine_check.py
# Rebuild the corpus experiment; saves compact fixtures and summary in this folder.
python3 research/nonmodular/real_game_probe.py
```

Both Rust probes compile instrumented copies in temporary directories and leave
production sources untouched. The logical inference uses only native row,
column and box constraints and requires no uniqueness premise; the independently
established uniqueness premise belongs only to the selected baseline techniques.
This experiment demonstrates additional deductions beyond the tested bounded
engine implementations. It establishes neither historical novelty nor a
separation from all known Sudoku methods or general CSP reasoning.
