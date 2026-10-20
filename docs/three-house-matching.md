# Three-House Matching

`Technique::ThreeHouseMatching` performs exact support filtering on the row,
column and box through one cell. It uses no modular arithmetic or uniqueness
assumption. The main hint and rating ladders try it after Death Blossom and
before Arithmetic Counting and forcing. Recursive forcing propagation omits it
to avoid multiplying its search cost. Its Extreme tier and numerical value
8.5 are uncalibrated engine estimates, not published Sudoku Explainer ratings.

This is an application of established simultaneous matching, Hall reasoning,
and constraint decomposition. Historical mathematical novelty is **not**
established. The research comparison even identifies standard ALS-XZ and linear
counting proofs for the original synthetic example; see the
[prior-art audit](../research/nonmodular/three_house_prior_art.md) and
[finite classification](../research/nonmodular/three_house_classification.md).

## Exact three-house decomposition

Select the row R, column C and box B through a center. Partition their union into
six disjoint membership atoms, ordering each atom by row-major cell index:

| Atom | Cells | Size | Required digits in a signature |
| --- | --- | ---: | --- |
| 0 | Center | 1 | {d} |
| 1 | (R intersect B) minus center | 2 | A |
| 2 | (C intersect B) minus center | 2 | E |
| 3 | R minus B | 6 | U minus ({d} union A) |
| 4 | C minus B | 6 | U minus ({d} union E) |
| 5 | B minus (R union C) | 4 | U minus ({d} union A union E) |

Here U={1,...,9}, A and E are two-digit sets, and {d}, A and E are pairwise
disjoint. There are at most 9*C(8,2)*C(6,2)=3,780 signatures per center.

**Support theorem.** A domain-respecting assignment to the three houses exists
if and only if some signature admits a perfect matching in every atom between
its cells and required digits. For the forward direction, take each atom's
assigned digit set from a respecting assignment. Conversely, combine the atom
matchings. Atoms partition the 21 cells, and their digit sets partition U within
each selected house, so each house contains exactly digits 1 through 9. Moreover,
a candidate has support precisely when an accepting signature has an atom
matching containing that edge. Search unions those matching-supported edges
using subset dynamic programming and proposes an unsupported live candidate.
This theorem concerns the selected-house relaxation, not a whole Sudoku solution.

## Independently checked certificates

The checker never invokes the search matching DP. `ThreeHouseProof` contains:

- Version 1, a SHA-256 digest binding all 81 placed values and effective domains,
  the selected center, and the target cell/digit to eliminate.
- An explicit assignment of all 21 union cells in ascending row-major order,
  respecting their domains and all three houses and avoiding the target. This
  rejects an already inconsistent three-house relaxation.
- At most 366 distinct `HallObstruction` records. Each identifies an atom, its
  exact required digit mask, and a nonempty subset of that atom's cells.

To check an obstruction, first restrict the target cell's domain to the assumed
value. Intersect each recorded subset cell's domain with the required digit set.
The union must have fewer digits than the number of subset cells. By the
pigeonhole principle, that atom cannot be filled injectively. The checker then
reconstructs **every** center-domain-permitted signature and requires a checked
obstruction with the matching atom and digit mask. Consequently every respecting
target assignment is impossible, so every full completion respecting the premises
avoids the target. No search budget or claimed count can substitute for coverage.

The checker rejects stale state hashes, unsupported versions, filled or absent
targets, targets outside the union, malformed masks/subsets, duplicate obstruction
keys, missing signature coverage and invalid base assignments. The digest is
state binding, not an authenticity signature.

[`verification/three_house.rs`](../verification/three_house.rs) machine-checks
the abstract injection/Hall and complete-signature-coverage theorems with Verus.
Native geometry, signature enumeration, bit masks, hashing, the production Rust
implementation, serialization and replay remain outside that formal proof.
Independent assignment enumeration and mutation tests check those interfaces.

## API and replay

```rust
use sudoku_core::{Grid, ProofCertificate, Solver, ThreeHouseReplay,
                  ThreeHouseSearchOptions};

let grid = Grid::new_classic();
let result = Solver::new().search_three_house(
    &grid, &ThreeHouseSearchOptions { max_centers: 81 });
if let Some(hint) = result.hint {
    if let Some(ProofCertificate::ThreeHouse(proof)) = &hint.proof {
        assert!(proof.verify(&grid));
        let replay = ThreeHouseReplay::capture(&grid, proof).unwrap();
        assert!(replay.verify());
        let mut next = grid.deep_clone();
        replay.apply(&mut next).unwrap();
    }
}
```

`get_three_house_hint` uses all 81 centers. Both direct search APIs preserve the
supplied candidates, do not mutate the grid, and require installed standard
Sudoku constraints. Additional restricting variants preserve the validity of
native deductions. Custom grids without native-rule provenance are rejected.
`max_centers` accepts 0..=81 and searches a deterministic prefix of centers.
`centers_checked`, `error`, and `exhausted` distinguish a truncated or invalid
search; `exhausted` is true only when all 81 triples have been checked without
a hint. Neither absence of a hint nor a local base witness establishes global
satisfiability or completeness for solving Sudoku.

`ThreeHouseProof::for_elimination` certifies a requested center/cell/digit directly.
`check` reports reconstructed signature and obstruction counts. Proofs and
`ThreeHouseReplay` derive Serde serialization; `Hint.proof` remains skipped in
ordinary hint JSON, so export the replay envelope explicitly.

A [concrete exported replay](../research/nonmodular/three_house_example_replay.json)
for Top95 puzzle #20 is independently checked in both Python and Rust. Regenerate
it with `python3 research/nonmodular/real_game_engine_check.py`.

Replay includes `CandidateState`, a versioned native projection containing exact
81-entry values and domain vectors. Placed values become singleton domains;
bit digit-1 permits a digit. `CandidateState` also backs the source-compatible
`ArithmeticState` alias; the existing arithmetic replay JSON layout is unchanged.
Replay restores real native constraints, verifies the certificate, requires the
identical current premises, and atomically removes just the certified candidate.
It preserves earlier deletions and refuses a stale second application. Added
variant constraints are checked on application but are not serialized in this
native projection. A local base witness does not prove variant satisfiability.

Candidate masks are premises. A proof does not authenticate how prior deletions
were obtained. General `get_hint`/rating methods rebuild or chain candidates as
before; their proof belongs to the internal state at the point of discovery.
Use the direct API with your stored state when exporting an exact replay.

## Evidence and limits

The [real-game experiment](../research/nonmodular/real_game_summary.json) tested
106 independently verified unique puzzles from Peter Norvig's author-hosted
Top95 and Hardest collections. The explicitly defined structural baseline solved
74 and stalled on 32. Exact three-house scans found 15 distinct puzzle/candidate
eliminations across eight stalled states and 17 centers (25 center/candidate
occurrences). All eight also produce no hint from the default bounded Arithmetic
Counting search, which reports beam pruning. This is a gap in those bounded
implementations, not a separation from general arithmetic or forcing methods.

All eight exact states, givens, provenance, unique completions and 207 baseline
steps are retained in [real_game_cases.json](../research/nonmodular/real_game_cases.json).
Independent Python checks replay those traces and enumerate each puzzle's
solutions. Rust integration tests reproduce the states, check all 25 certificates,
compare each deletion with the unique solution, and check portable atomic replay.
Observed release first-hint searches took roughly 0.4–4.3 ms on the local machine;
these small curated collections are not a general performance or prevalence study.

The finite synthetic classification examines 33 eliminations in 18 retained
states: 26 are already refutable using a house pair, while seven need all three
selected houses in their original assumed domains. All 33 are refuted by ordinary
single-house Hall propagation after assuming the candidate. No previously unknown
solving principle is established by this experiment.
