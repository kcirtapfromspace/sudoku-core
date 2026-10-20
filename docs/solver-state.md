# Retained solver state and checked forcing

`Grid` candidate masks are also used for player pencil marks. A solver must opt
in before treating those masks as logical premises. `Solver::get_hint` retains
its existing behavior: it rebuilds candidates from values. The new
`Solver::get_hint_with_candidates` runs the complete strategy ladder on the
supplied masks without rebuilding them.

For a sequence of deductions, use `SolverState`:

```rust
use sudoku_core::{Grid, Solver, SolverState};

let player_grid = Grid::from_string(
    "530070000600195000098000060800060003400803001700020006060000280000419005000080079"
).unwrap();
let solver = Solver::new();
let mut state = SolverState::new(&player_grid).unwrap();
while !state.grid().is_complete() {
    let hint = state.step(&solver).unwrap().expect("solvable puzzle");
    println!("{}", hint.explanation);
}
```

`new` clones the values and rebuilds solver masks independently of player notes.
`from_candidates` explicitly retains existing deductions. The source grid and
its notes remain independent of the solver. `grid()` provides read access;
`step` applies one hint, preserving earlier eliminations. `get_hint` on the
state returns a hint without applying it. Standard logical solving, difficulty
rating and technique profiles use the same retained step path.

Validation requires installed standard Sudoku rules, valid values, valid
nonempty domains, no candidates conflicting with placed peers, feasible native
house matchings, and no placed-value violation of additional installed rules.
This checks local consistency; it does not authenticate the derivation of
imported eliminations or establish global satisfiability. The full strategy
ladder retains its existing uniqueness assumptions. Preserve the actual game's
given flags when importing candidate masks: `CandidateState::to_grid` restores
a native proof projection and does not recover given/user provenance.

All solver placements, backtracking branches and assumption propagation now
remove candidates monotonically. Only entry points that start from player
values and explicit edit operations rebuild masks. `SolverState::edit` rebuilds
all masks and discards learned relations. `undo` restores the preceding candidate
snapshot and discards relations. A failed edit or proof application leaves the
state unchanged. These changes can affect difficulty ratings and generated
puzzle selection because deductions survive later placements.

## Shared dispatch and policies

Hinting, logical solving and forcing discovery call one ordered technique
dispatcher. `PropagationPolicy::Singles` permits naked and hidden singles;
`Logical` permits established techniques through Death Blossom; `Full` adds
three-house matching, arithmetic, global coloring and checked forcing. The
`get_hint_with_options` API accepts this policy and forcing budgets. Only `Full`
uses the existing backtracking fallback. `SolverState::step_with_options` applies
a step under the same policy and budgets. Forcing discovery uses the first two
policies, so it cannot recursively invoke the full forcing ladder.

`ForcingSearchOptions` separately selects the native proof calculus:
`ForcingPropagation::House` uses matching propagation in each house;
`ThreeHouse` additionally propagates all 81 crossing triples. The default permits
one additional cell-split level, 1,024 node visits and 16,384 rule evaluations.
Supported limits are eight additional levels, 4,096 visits and 32,768 evaluations.
These counts include abandoned search branches; rule evaluations include scans
that make no edits. The dispatcher shares its forcing budget between endpoint
discovery and certification. These limits do not bound other strategy engines
or the backtracking fallback.

## Forcing proofs

```rust,ignore
let result = solver.search_forcing(&grid, 7, 3, false, &options);
if let Some(proof) = result.proof {
    assert!(proof.verify(&grid));
    proof.apply(&mut grid).unwrap();
}
```

The endpoint above means eliminate digit 3 from cell 7 (zero-based r1c8).
`ForcingProof` contains an exact native-state hash, the conclusion, and a
postorder tree. The root assumes the opposite conclusion. Each node contains
recorded house or three-house propagation edits followed by either a complete
cell-value split or an explicit contradiction leaf. Each branch records its
assumed digit and child node. The configured depth counts additional splits
after the root endpoint assumption.

The checker replays the supplied tree; it never chooses pivots or calls the
producer's matching DP. It independently reconstructs propagation support by
Hall subsets and tight sets. It checks the exact before/after edits, complete
candidate coverage, unique node use and child ordering. Hall leaves contain a
native-house subset whose digit union is too small. Three-house leaves contain
Hall blockers covering every admissible boundary signature. Unknown fields,
invalid geometry, stale states and incomplete trees are rejected. `apply`
verifies first and commits atomically, preserving all other candidate masks and
checking the installed constraints.

Public ladder hints carry `ProofCertificate::ForcingProof`. The historical
`Forcing { source, branches }` variant remains available for compatibility with
existing evidence consumers; it is only an internal discovery summary and is
never accepted as a checked forcing hint. A proposed endpoint that cannot be
certified within the allowed calculus and budget is skipped.

`search_forcing` reports visits, evaluations, errors and budget exhaustion.
An incomplete search never emits a proof. A missing proof is not a claim that
no forcing deduction exists. A depth limit can prevent a proof even when the
visit and evaluation budgets remain available.

The committed `tests/fixtures/forcing_mith515.json` ports the verified research
example: seven nodes, three complete splits, one Hall leaf, three three-house
leaves, and 84 propagation events (the ten-event initial assumption prefix plus
74 tree events). Tests independently replay it, reject altered proofs, and
rediscover its endpoint with two additional split levels. Recreate the format
translation with:

```sh
python3 research/global_obstructions/export_forcing_replay.py
```

The original research verifiers and certificates remain available separately.

## Guardian OR relations

`GuardianProof::from_cells` proves that at least one supplied `GuardianLiteral`
is true. It refutes their simultaneous negation on a selected native peer graph,
with complete list-coloring branches and a respecting base witness. The proof
includes its candidate premises and can be serialized independently.

`SolverState::retain_guardian` verifies a relation against the exact current
state before storing it. `learn_guardian` performs a bounded scan for small
three-digit patterns with multiple guardian extras. The normal step path tries
this scan before falling back from the logical ladder; placement certificates
from Global Coloring also become retained single-literal relations.

Later deductions reevaluate all remaining literals. A single remaining literal
forces its value. A candidate conflicting with every live literal is removed.
These hints carry a `GuardianPropagation` certificate, which checks the original
relation and verifies that the current state descends monotonically from its
premises. Satisfied relations make no further deductions; losing every literal
reports a contradiction. Relations are invalidated on edits and undo, and the
state exposes no mutable grid access that could bypass this lifecycle.
