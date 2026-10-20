# A unique puzzle beyond the older logical engine and per-digit weighted rank

**A qualifying real benchmark was found. Its obstruction is an established
tridagon with guardian forcing; no new mathematical family is claimed.** The
existing Global Coloring implementation certifies **r4c9=9** where a fresh build
of the full older logical dispatcher, including forcing, returns no hint.
The opposite branch survives exact Hall and crossing-three-house closure,
all arbitrary three-house feasibility checks, and a stronger full-board
single-digit-template relaxation. No production algorithm change is needed.

The independently unique, 28-given puzzle is Monhard T&E3 source line 1:

```text
987......
6..95....
.4.......
3..21.65.
...5.3..2
....96.31
.3.1...2.
....6.59.
....2.1.3
```

Its [immutable author-maintained source](https://raw.githubusercontent.com/denis-berthier/Sudoku-classif/15c9d217b08d2d3ebf0f147d4f7c0beb5c2fcdb2/Monhard-until-2023-08-15-TE23/TE3/puzzles.txt)
is bound by the source hash and line number in the retained corpus. An independent
exact-cover enumeration, capped at two completions, finds precisely the stored
solution. The clue string, original domains, derivation trace, retained state,
opposite-assumption closure and endpoint are in
[`nonrank_full_old_stall.json`](nonrank_full_old_stall.json).

## What passed

| Obligation | Evidence |
| --- | --- |
| Valid premises from givens | Independent completion count and exact replay of 20 base events, including twelve older forcing findings; every old endpoint agrees with the sole solution. |
| Full older logical engine stalls | Two complete calls on the same live retained state, 50 cells unplaced, from a fresh archived `0532d63d9241a69c1fb57e9c34a1533effe1981d` build; both return no hint and no trace. |
| Complete local closure | Hall matching support in all 27 houses and all 81 crossing-three-house support checks reach unchanged fixed points before and after the opposite assumption. Six branch events are replayed. |
| Arbitrary triples | Every one of the 2,925 selections of three native houses has a candidate-respecting assignment in the closed opposite branch. This checks feasibility, not all arbitrary-triple GAC. |
| Beyond all per-digit weighted counting | An exact feasible mixture of full native digit templates: 33 positive rational terms, nine digit masses and 81 cell coverages exactly one, maximum denominator 112. |
| Native continuous LP | The same mixture projects to 169 nonzero indicators satisfying all 324 native equations and the opposite target exactly. |
| Integral contradiction | Actual default Rust export: 16 selected cells, 169 postorder proof nodes, r4c9=9. Independent reconstruction checks every branch; ten corruptions are rejected. |
| Engine integration | Regression tests run the complete current logical ladder on the retained state, default direct search, independent certificate verification, serialization and exact-state replay. |

The older comparison includes its actual configured forcing, arithmetic and
three-house behavior. Arithmetic records beam pruning. It does not claim
exhaustion of every possible arithmetic proof, unbounded forcing, or every
hybrid ordering of known methods. The older logical dispatcher excludes the
public hint API's backtracking fallback: an exact search solver can solve this
unique puzzle. The fresh archive/source hashes and both successful calls are
retained in [`nonrank_puzzles_first_audit.json`](nonrank_puzzles_first_audit.json).

The bounded corpus run completed six published puzzles: five full-old logical
stalls and one solved board. Only the first is promoted through this entire
acceptance suite. This is not a claim that all 200 modern corpus entries were
processed by the forcing-inclusive closure experiment. The public six-case
manifest, result counts and reproduction instructions are in
[`nonrank_puzzles_notes.md`](nonrank_puzzles_notes.md).

## Why weighted rank does not explain this snapshot

For each digit, distribute mass one over legal full placements with exactly one
cell in every row, column and box. Require total coverage one at each cell
across all digits. Every convex combination satisfies every linear inequality
valid for that digit's placement polytope. Therefore an exact feasible mixture
survives **all single-digit weighted stable-set, rank and hyperhole bounds**,
including bounds on subsets not selected by the graph detector.

The accepted opposite branch has such a mixture in
[`nonrank_template_primal.json`](nonrank_template_primal.json). Its independent
standard-library checker checks every listed template and all equations with
rational arithmetic. Numerical optimization proposes witnesses only. The
negative control is the previous 13-cell hyperhole: its exact Farkas certificate
has RHS -1 and is checked against every admissible full template. Thus this
stronger screen rejects the earlier rank example and accepts the new benchmark.

This is a statement about each digit separately. Joint multidigit cuts and
correlations between templates are outside the relaxation. Template methods,
fractional coloring and finite failed-literal refutations have established
prior art. The [independent mathematical audit](nonrank_audit_report.md) gives
the precise theorem, checker boundary, controls and primary references.

## The human explanation is a known tridagon composition

Assume r4c9 is not 9. In column 9, r2c9, r4c9 and r8c9 form a naked
{4,7,8} triple. The {6,7} cell r7c9 becomes 6, removing 6 from r9c8 in box 9.
The remaining twelve cells all have candidates {4,7,8}:

```text
r4c6 r4c9 r5c5 r5c8 r6c4 r6c7
r7c5 r7c7 r8c6 r8c9 r9c4 r9c8
```

Their peer graph is exactly isomorphic to three published tridagon controls;
all three explicit vertex maps are checked independently. This contradiction
proves r4c9=9. The unforced pattern has two extra candidates, 9r4c9 and 6r9c8.
The direct 16-cell certificate incorporates the surrounding restrictions
without trusting a named-pattern identity or a supplied graph. Its underlying
graph is four-colorable: incompatible candidate lists cause the contradiction.
It is not a new uniform four-color critical obstruction.

The separate [structural search](nonrank_graph_report.md) found a synthetic
27-cell graph requiring five colors with an exact fractional chromatic number
of four. Its live guardian is r5c6=1. Both local closures pass, and an independent
**49-term full-board template mixture** survives all per-digit-valid linear
inequalities on the closed opposite branch. All 2,925 selected-house triples
are separately feasible. An explicit
[`research refutation tree`](nonrank_graph_native27_tree.json) has 11,081 nodes;
its separate checker reconstructs every candidate branch and conflict.
This is a stronger synthetic
separation than the old 13-cell rank example. Another 27-cell graph fails Hall
and both LP screens and is retained with exact negative certificates.

The synthetic positive has no verified unique givens-only realization or
forcing-inclusive older-engine stall. Its support and exhaustive refutation
exceed the current production bounds. It remains research evidence with exact
independent checks; it is not added as a special production rule. Historical
novelty of that graph has not been established. No unique puzzle exhibiting a
verified new obstruction family has been found in this run.

## Reproduce and use the deduction

```sh
# Exact end-to-end acceptance, including all 2,925 triple feasibility checks.
python3 -S research/global_obstructions/nonrank_verify.py
# Premise provenance and retained full-old-dispatcher comparisons.
python3 -S research/global_obstructions/nonrank_puzzles_verify.py
# Independent stronger-relaxation controls and witness mutations.
python3 -S research/global_obstructions/nonrank_audit_template.py --self-test
# Rust automatic search, certificates and portable replays.
cargo test --test global_coloring --locked
cargo test full_candidate_ladder_advances_the_frozen_old_engine_stall --locked

# Freshly rebuild the old engine and rerun both full dispatcher calls.
# The destination must be new; an audit runs this in dedicated CI too.
python3 -S research/global_obstructions/nonrank_puzzles_audit.py \
  --fresh-folder target/nonrank-fresh-baseline \
  --output target/nonrank-fresh-audit.json
```

Import the fixture's `state` as a `CandidateState`, then call the direct
`get_global_coloring_hint` / `search_global_coloring` API, or apply the checked
[`replay_nonrank_full_old_stall.json`](replay_nonrank_full_old_stall.json) to
that exact state. The general `get_hint` API rebuilds candidates; it does not
preserve all prior eliminations merely by importing the derived value grid.
The regression exercises the actual internal full logical ladder with the
retained masks and original given flags. Replay premises, production search
bounds and proof trust limits are documented in
[`docs/global-coloring.md`](../../docs/global-coloring.md).
