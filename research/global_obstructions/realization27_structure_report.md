# The exact native27 support cannot have a unique outside-clue realization

There is an exact obstruction to the requested realization. **For this native
house incidence, every completed Sudoku respecting the 27-cell pattern has a
different completion agreeing at every cell outside the support.** Consequently
no givens puzzle with all 27 support cells unfilled can be unique while having
even one completion respecting this pattern. Changing the completion, adding
outside clues, shrinking candidate lists by sound reasoning, or applying a
native Sudoku symmetry cannot remove this obstruction.

This concerns the stored embedding. A graph-isomorphic embedding with a
different incidence of rows, columns and boxes remains a separate question.
Clues inside the support, additional variant constraints, or a different
support pattern are also different tasks.

## Exact finite proof

Let S be the 27 recorded cells, g=r5c6, and P={2,3,5,9}. The native graph on S
is not four-colorable. In a respecting completion, every cell of S except g
uses P, while g uses its extra digit 1. There are exactly **1,320 proper
labelled four-colorings of S minus g**.

The four cells r8c1,r8c3,r8c4,r8c7 form a native four-clique. Fixing their four
colors to a prescribed order gives exactly **55 normalized colorings**. Their
24 global color permutations are all distinct and recover the complete set of
1,320 labelled colorings. The independent checker additionally enumerates all
labelled colorings directly, without fixing this clique or pruning color
symmetries, and checks equality of the two complete sets.

For each labelled coloring, record its **27-house signature**: the set of
palette colors appearing in each native row, column and box inside S minus g.
Proper coloring ensures each appears at most once in a house, so these sets
also specify the exact digit multisets. The complete classification is:

| Colorings sharing one signature | Number of signatures |
| --- | --- |
| 2 | 192 |
| 3 | 96 |
| 5 | 72 |
| 6 | 48 |

There are 408 signatures and **none occurs only once**.

**Proposition.** Every respecting completion has a distinct completion agreeing
outside S.

**Proof.** Take its labelled coloring on S minus g. Choose a different coloring
with the same signature, which exists by the exhaustive classification.
Replace the selected cells by that coloring and leave g and all outside cells
unchanged. In each native house, the selected-cell digit multiset is unchanged;
the remaining cells are unchanged as well. Every house therefore still contains
digits 1 through 9 exactly once. The new Sudoku is different, and every given
outside S is preserved. Thus no such outside-clue puzzle is unique. ∎

The proof is invariant under native Sudoku automorphisms and global digit
renaming. It does not assert invariance under arbitrary graph-isomorphic
embeddings that realize edges using different houses.

The retained certificate
[`realization27_structure_proof.json`](realization27_structure_proof.json)
contains all 55 normalized colorings, exact counts, support, palette, anchor,
and an explicit completion pair. The standalone standard-library checker
[`realization27_structure_verify.py`](realization27_structure_verify.py)
reconstructs native geometry, independently enumerates all 1,320 labelled
colorings, checks all signatures, verifies the original graph is not
four-colorable, and rejects eight malformed certificate mutations. It does
not trust a numerical optimizer, a search count, or the stored representatives
as a completeness claim. A separate agent independently enumerated the same
counts in `realization27_audit_trades.py`.

## Explicit stored-completion trade

Even the maximal 54-given puzzle containing **every** outside value of the
stored completion has six completions. A smallest difference from the stored
completion changes these six zero-based cells:

```text
65,66,69,74,75,78
```

They are r8c3,r8c4,r8c7 and r9c3,r9c4,r9c7. The original rows contain
`(3,2,5)` and `(2,5,3)` in those columns. Swapping those triples gives another
completed Sudoku with all 54 outside clues preserved. The retained full pair
is checked directly. The universal signature proof does not assume this same
six-cell trade exists for every completion; it guarantees some internal trade.

The search agent independently sampled 4,999 compatible completions, all with
nonunique maximal outside puzzles. That finite sample suggested the theorem;
the exhaustive coloring classification, rather than sampling, proves it.

## A simpler obstruction for preserving all four live palette candidates

The 26 nonguardian support cells meet all nine boxes. Thus every possible
palette given anywhere outside S sees a full-palette support cell in its box
and removes that palette candidate. Preserving all four live candidates at
every support cell therefore forbids every palette given on the board. All
four palette digits remain globally interchangeable, giving at least 4!=24
completions from any one completion.

This simpler symmetry argument already rules out an exact uniform-mask unique
realization. The full signature theorem is stronger: it also rules out unique
realizations with outside palette clues and nonuniform subsets of P at the
support cells.

## What the LP checks establish

A scan of all 54 outside positions and four hypothetical palette digits gives
216 possible single clues. Of these, 172 are excluded by the original closed
guardian-false domains; the other 44 pass the native 324-equation LP screen.
They occupy cells `{19,21,22,23,26,28,46,55,57,59,62}`, with all four palette
choices at each. Floating screening results alone are not certificates.

There is also an exact positive control: adding the stored-completion clue
r3c2=3 preserves local Hall/all-81 closure and admits a full single-digit
template primal with 48 rational terms, verified exactly by the template
search checker. Its raw input and witness are retained under
`target/realization27-one-palette-*`. Therefore a stronger claim that *every*
palette clue destroys the relaxation is false. Template-LP feasibility does
not rescue uniqueness from the universal internal trade.

## Production implications

The underlying candidate deduction remains sound without a uniqueness premise.
Its inability to arise as a unique outside-clue puzzle does not invalidate the
explicit research proof. The automatic 24-cell bound should not silently grow
to manufacture a benchmark. A separate explicitly configured larger checker
can retain bounded cell/node counts, exact branch coverage, native geometry,
premise binding, and a respecting base witness. Production changes belong to
the parent implementation; this research made none.

A shorter certificate is possible in principle for this uniform four-color
case: every proper coloring assigns a native four-clique all four colors, and
a global palette permutation maps those assignments to one fixed ordering.
An independently checked symmetry rule could therefore verify one normalized
refutation instead of all 24 equivalent labelings. Such a rule must explicitly
check uniform domains after the assumption and the native clique. It is not
sound for arbitrary unequal lists without additional checks.

The next meaningful realization target must change at least one premise:
permit clues inside the unavoidable support, find another house incidence,
find a different obstruction, or add explicitly requested variant constraints.
Further completion or clue-subset sampling on this exact incidence cannot
produce the requested unique puzzle.

## Reproduce

```sh
python3 research/global_obstructions/realization27_structure_verify.py
python3 research/global_obstructions/realization27_structure_verify.py --self-test
```

The direct labelled search visits 16,217 nodes. The verifier uses only the
standard library and takes well under a second without mutation tests on the
development machine. Full raw completion lists and exploratory clue screens
are under `target/`; the public proof is compact and self-contained apart from
binding to the existing native27 source fixture.

This is a computer-assisted result for one specified support. The general
notion of an unavoidable set and the need for a unique clue set to intersect
it are established prior art; see the independent trade audit's primary-source
comparison. No historical priority claim is made for that general principle.
