# Published real-puzzle controls for global coloring proofs

These are **established tridagon (trivalue oddagon) patterns**, included as real
positive controls. They are not proposed new mathematical discoveries. Denis
Berthier's author-hosted repository supplies the original givens, intermediate
candidate states, explicit twelve-cell patterns and full solving traces:

| Puzzle | Primary source at the recorded immutable commit | Verified placement |
| --- | --- | --- |
| mith collection #4662 | [Published trace](https://github.com/denis-berthier/CSP-Rules-Examples/blob/fd9207f8dc1eb6e2e655048ba3b327d423900e1b/Sudoku-a/Tridagons/%234662-SFin%2BTrid.clp) | r5c4 = 2 |
| mith collection #4401 | [Published trace](https://github.com/denis-berthier/CSP-Rules-Examples/blob/fd9207f8dc1eb6e2e655048ba3b327d423900e1b/Sudoku-a/Tridagons/%234401-Trid%2BW8.clp) | r4c4 = 9 |
| Loki, created by mith | [Published trace](https://github.com/denis-berthier/CSP-Rules-Examples/blob/fd9207f8dc1eb6e2e655048ba3b327d423900e1b/Sudoku-a/Tridagons/Loki-11.9.clp) | r8c8 = 1 |

Each JSON fixture retains the exact givens, source URL and byte hash, a complete
respecting solution, all earlier local inference steps and the twelve-cell
native-peer graph. An independent exact-cover enumeration, using only the
givens and the 324 native requirements, found one solution when counting up to
two for each puzzle. The local and coloring proofs do not use uniqueness.

For each puzzle, start with the ordinary candidate domains derived from its
givens. Repeatedly enforce exact perfect-matching support on all 27 houses and
exact support on all 81 crossing row/column/box triples until no deletion
remains. The fixtures record every domain change. The target assumption removes
the guardian's extra digit. Repeat the same complete local closure on this
assumed branch.

All three assumed branches:

- Survive full individual-house matching GAC.
- Reach a fixed point for all 81 crossing three-house matching checks.
- Admit an assignment for **every one of the 2,925 arbitrary triples of native
  houses**. This stronger additional control tests feasibility; it does not
  enforce GAC for arbitrary triples.
- Admit an exactly verified rational solution to the **full 324-equation native
  continuous LP**, including all remaining domains, original givens and the
  guardian-false assumption. The separate `linear_real_*.json` witnesses use
  exact `Fraction` checking, not floating tolerances.
- Have no list coloring of the recorded twelve-cell native-peer graph, whose
  domains are the same three digits in every vertex. Exhaustive independent
  MRV checks refute them in 142, 154 and 142 visited nodes respectively.

Consequently, the guardian-false assumption is impossible and the recorded
placement follows. The global graph contradiction survives the named local
relaxations and full native LP. This does not establish a separation from
stronger local-marginal LPs, cutting planes, SAT, full CSP reasoning or known
forcing methods.

## Existing engine comparison

`real_existing_check.py` extracts the pre-global engine at commit
`0532d63d9241a69c1fb57e9c34a1533effe1981d` into a temporary directory and exposes
its existing detector dispatcher there. Current production sources are never
edited by the probe. The structural comparison excludes three-house matching,
arithmetic, forcing chains and backtracking, and preserves all stored candidate
eliminations across placements. Independently established uniqueness permits
the old uniqueness detectors in this separate baseline comparison.

| Puzzle | Old structural findings after local closure | Default Arithmetic Counting | First full old dispatcher finding, including forcing |
| --- | --- | --- | --- |
| #4662 | One Finned Swordfish: r4c2 ≠ 2; then fixed point | No hint | Dynamic Forcing Chain: r6c1 ≠ 9 |
| #4401 | None | No hint | Dynamic Forcing Chain: r5c3 ≠ 4 |
| Loki | None | No hint | Dynamic Forcing Chain: r8c1 ≠ 7 |

The forcing deductions are additional, unrelated findings. The experiment does
**not** show that the old engine has no possible hint or that forcing cannot
eventually solve these puzzles. Arithmetic search is bounded; its lack of a hint
is not an impossibility proof.

`real_4662_after_existing.json` adds the one old structural step before closing
the local filters again. Its guardian-false closed domains are identical to
those in the pure-local `real_4662.json`, so the same exact LP witness applies.
For a baseline with no remaining old structural hint, use this extra-step
fixture together with `real_4401.json` and `real_Loki.json`.

## Reproduce

```sh
# Standard-library exact checks: original givens, unique completion, all traces,
# both local fixed points, twelve-cell graph, all 2,925 house triples, full LP.
python3 research/global_obstructions/real_verify.py

# Additionally verify immutable primary-source byte hashes over the network.
python3 research/global_obstructions/real_verify.py --fetch-sources

# Compile a temporary historical engine and reproduce the comparison table.
python3 research/global_obstructions/real_existing_check.py

# Compile a temporary current engine; check public search/certificate/replay APIs.
python3 research/global_obstructions/real_engine_check.py
```

The exact checker depends on the earlier independent matching functions in
`research/nonmodular/real_game_probe.py` and the separate exact LP checker in
`linear_witness.py`. LP verification needs no NumPy or SciPy. Native graph
coloring and arbitrary-house-triple checks use a separate direct search rather
than the matching DP used for discovery.

## Generated four-digit givens controls

`real_four_digit_candidate.json` and `real_four_digit_unique.json` are new
**generated givens selections**, not puzzles published by Denis Berthier or
Peter Norvig. Their respecting full completion was independently obtained from
Norvig's top95 puzzle #20; only that completion has the cited public origin.
The fixture metadata distinguishes these sources explicitly.

The selected thirteen-cell native-peer graph cannot be colored with the four
digits {4,6,7,8}. The companion structural experiment proves that an independent
set has at most three vertices, so four digit classes can occupy at most twelve
of the thirteen cells. The guardian r2c8 therefore needs its extra digit 1.
This is an application of established graph-coloring and independent-set bounds;
historical mathematical novelty is not established.

To realize the candidate restrictions using only givens, the experiment selected
clues from the full completion outside the core and outside its four-digit
palette. It enumerated inclusion-minimal per-digit peer covers of all forbidden
core candidates: 2 × 5 × 8 × 7 × 9 = 5,040 possible combinations. Trial 132
produced a nineteen-given puzzle whose ordinary givens-derived domains, after
the recorded sound closure, suffice for the global proof. There are no arbitrary
candidate restrictions in this construction. `real_four_digit_mine.py` reproduces
this bounded discovery using SciPy only for its floating discovery screen;
acceptance uses the separate exact rational witness.

| Generated control | Givens | Completions counted up to two | Old nonforcing detector result |
| --- | ---: | ---: | --- |
| Sparse four-digit realization | 19 | 2, so uniqueness is not assumed | No finding; default arithmetic also returns no hint |
| Unique four-digit realization | 34 | 1 | Finned X-Wings and Pointing followed by singles solve the puzzle |

Both states and both guardian-false branches reach Hall GAC and all 81 crossing
three-house fixed points. Both assumed branches also pass all 2,925 arbitrary
house-triple feasibility checks and the exact full native 324-equation LP.
Independent thirteen-cell list-color searches refute their assumed branches in
1,193 and 14 visited nodes respectively. Every earlier local inference and every
givens-derived domain is checked by `real_verify.py`.

The nineteen-given case demonstrates a useful deduction on a nonunique puzzle
without a uniqueness premise. The thirty-four-given case demonstrates the
relaxation separation on a proper unique puzzle, but is a **negative control for
incremental engine value** because existing methods already solve it. Neither
case establishes that the pattern is absent from prior Sudoku literature.

An additional deterministic greedy search augmented givens while rejecting any
addition that let the old nonforcing baseline place the guardian. Thirty seeds
tested 1,470 clue additions and accepted 525 steps, ending with 30–42 givens.
No unique puzzle that also retained an old nonforcing gap was found in 376.8 s.
`real_four_digit_augmentation_summary.json` records the exact scope and
`real_four_digit_augment.py --seeds 30` reproduces it. The floating LP discovery
screens and restricted greedy clue order make this a bounded negative result,
not a proof that a qualifying unique puzzle cannot exist. Exploratory candidate
outputs go to `/tmp`, separate from the certified fixtures.
