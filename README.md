# `sudoku-core` (Rust Engine)

`sudoku-core` is the shared Rust engine used by the TUI, WASM, and iOS targets. It provides:

- A `Grid` model with constraint-aware candidate tracking
- A solver with 45 human-style techniques plus backtracking
- A hint system (returns an explanation + affected cells)
- A puzzle generator that targets a requested difficulty with SE rating bounds
- Dual rating: technique-based difficulty tiers and Sudoku Explainer (SE) numerical ratings

## How The Engine Operates

### Data Model

- `Grid` (`src/grid.rs`) is a 9x9 array of `Cell`s.
- Each `Cell` stores:
  - `value: Option<u8>`
  - `given: bool` (part of the original puzzle)
  - `candidates: BitSet` (1-9)
- `Grid::recalculate_candidates()` recomputes candidates for every empty cell by:
  1. Resetting all empty cells to all candidates (1-9)
  2. Removing candidates that violate any constraint based on filled values

### Constraints (Variants)

The engine is constraint-driven (`src/constraint.rs` + `Grid::constraints()`):

- Classic: row/column/box constraints
- X-Sudoku: adds diagonal constraints (`Grid::new_x_sudoku()`)
- Killer: adds cage-sum constraints (`Grid::new_killer(...)`)
- Additional constraints exist (thermo, etc.) and are enforced through the same interface

The solver, validator, and candidate updates all operate through the active constraint set.

### Solver

`Solver` (`src/solver/`) is organized as a modular directory with dedicated engine files for each technique family (fish, ALS, AIC/chains, uniqueness, etc.), plus `types.rs` for the `Technique` enum and `fabric.rs` for the dual-indexed candidate state. It has two layers:

1. **Human-style techniques** (45 techniques, ordered by complexity)
   - The solver can search for a concrete next step (a `Hint`) by trying techniques in increasing complexity.
   - Techniques: Naked/Hidden Singles, Naked/Hidden Pairs/Triples, Pointing Pairs, Box-Line Reduction, X-Wing, Finned X-Wing, Swordfish, Finned Swordfish, Jellyfish, Finned Jellyfish, Naked/Hidden Quads, Empty Rectangle, Avoidable Rectangle, XY-Wing, XYZ-Wing, WXYZ-Wing, W-Wing, X-Chain, 3D Medusa, Sue de Coq, AIC, Franken Fish, Siamese Fish, ALS-XZ, ALS-XY-Wing, ALS Chain, Unique Rectangle (Types 1-4), Hidden Rectangle, Extended Unique Rectangle, Mutant Fish, Aligned Pair Exclusion, Aligned Triplet Exclusion, BUG+1, Death Blossom, Nishio Forcing Chain, Kraken Fish, Region Forcing Chain, Cell Forcing Chain, Dynamic Forcing Chain.
2. **Backtracking fallback**
   - `solve()` uses recursion with the MRV heuristic (choose the empty cell with the fewest candidates).
   - It tries each candidate, validates the grid, and recurses until solved.
   - `count_solutions(limit)` uses the same approach but counts solutions up to a cap (used for uniqueness checks).

### Hints

`get_hint()` returns the first available step from the technique ladder. Each hint includes a human-readable explanation, affected cells, and the technique's SE rating. If no human technique applies, it falls back to:

- Solve the puzzle fully
- Return a "Backtracking" hint for the next empty cell (a forced placement from the known solution)

## Difficulty Rating

The engine provides two complementary rating systems.

### Technique-Based Difficulty Tiers

`Solver::rate_difficulty()` assigns a difficulty tier by simulating human-style solving:

1. Clones the puzzle and recalculates candidates.
2. Repeatedly applies techniques in increasing difficulty until the grid is complete.
3. Tracks the **hardest technique actually required** (`max_technique`).
4. If it gets stuck (no technique applies and the grid is not complete), it classifies the puzzle as **Extreme** (backtracking required).

Final mapping is based on the hardest technique needed:

| Tier | Techniques Required |
|---|---|
| **Beginner** | Naked Singles only, <= 35 empties |
| **Easy** | Naked Singles only, > 35 empties |
| **Medium** | Hidden Singles |
| **Intermediate** | Naked/Hidden Pairs/Triples |
| **Hard** | Pointing Pair / Box-Line Reduction |
| **Expert** | Fish (X-Wing, Swordfish, Jellyfish, Finned variants), Naked/Hidden Quads, Empty Rectangle, Avoidable Rectangle, Unique Rectangle, Hidden Rectangle |
| **Master** | Wings (XY-Wing, XYZ-Wing, WXYZ-Wing, W-Wing), Chains (X-Chain, AIC), 3D Medusa, Sue de Coq, Franken/Siamese Fish, ALS-XZ, Extended UR, BUG+1 |
| **Extreme** | ALS (XY-Wing, Chain), Mutant Fish, Aligned Pair/Triplet Exclusion, Death Blossom, Forcing Chains (Nishio, Kraken, Region, Cell, Dynamic), Backtracking |

This is intentionally technique-based, not purely "clue count" based (although clue count is used as a generation constraint).

### Sudoku Explainer (SE) Numerical Rating

`Solver::rate_se()` returns a community-standard numerical rating (1.5 - 11.0 scale) based on the hardest technique used. Each technique maps to a fixed SE value via `Technique::se_rating()`:

| Technique | SE |
|---|---:|
| Hidden Single | 1.5 |
| Naked Single | 2.3 |
| Pointing Pair | 2.6 |
| Box-Line Reduction | 2.8 |
| Naked Pair | 3.0 |
| X-Wing | 3.2 |
| Finned X-Wing / Hidden Pair | 3.4 |
| Naked Triple | 3.6 |
| Swordfish / Hidden Triple | 3.8 |
| Finned Swordfish | 4.0 |
| XY-Wing | 4.2 |
| XYZ-Wing / W-Wing | 4.4 |
| X-Chain | 4.5 |
| Empty Rectangle / Avoidable Rectangle / Unique Rectangle / WXYZ-Wing | 4.6 |
| Hidden Rectangle | 4.7 |
| Naked Quad / 3D Medusa / Sue de Coq | 5.0 |
| Jellyfish | 5.2 |
| Finned Jellyfish / Hidden Quad | 5.4 |
| ALS-XZ / Franken Fish / Siamese Fish / Extended UR | 5.5 |
| BUG+1 | 5.6 |
| AIC | 6.0 |
| Aligned Pair Exclusion | 6.2 |
| Mutant Fish | 6.5 |
| ALS-XY-Wing | 7.0 |
| ALS Chain / Aligned Triplet Exclusion / Nishio Forcing Chain | 7.5 |
| Kraken Fish | 8.0 |
| Cell Forcing Chain | 8.3 |
| Death Blossom / Region Forcing Chain | 8.5 |
| Dynamic Forcing Chain | 9.3 |
| Backtracking | 11.0 |

**Note on SE vs Difficulty Tier ordering:** The SE system considers Hidden Singles (1.5) easier than Naked Singles (2.3), which inverts the Beginner/Easy vs Medium tier ordering. This is correct per the SE community standard — hidden singles are "last remaining in a house" (scanning), while naked singles require full candidate elimination. The two scales intentionally measure different axes: SE measures *technique complexity*, our tiers measure *pedagogical progression*.

## Puzzle Generation

Puzzle generation is in `src/generator.rs`.

Generation is "generate-and-test" with best-candidate tracking:

1. **Create a fully solved grid**
   - Fill the 3 diagonal 3x3 boxes randomly (they don't interact).
   - Use the solver to complete the grid.
2. **Remove givens while preserving uniqueness**
   - Iterate shuffled positions and remove values in symmetry pairs (configurable; `Extreme` disables symmetry).
   - After each removal, check the puzzle still has **exactly 1 solution** (`Solver::has_unique_solution`, which counts up to 2).
   - Stop removing once we're around the minimum givens threshold.
3. **Rate and filter**
   - Run `Solver::rate_difficulty()` on the candidate puzzle.
   - Accept puzzles whose rated difficulty is **equal to the target or one step easier**, whose number of givens falls within the configured range, and whose SE rating (if bounds are configured) falls within the specified range.
4. **Best-candidate fallback**
   - If no perfect match is found within the attempt budget, the generator returns the **closest candidate seen** (measured by distance from target difficulty) rather than a random last attempt.

### Generator Config Per Difficulty

Each difficulty has a preset configuration (`GeneratorConfig::{beginner,easy,...}`):

| Target | Symmetry | Givens | Attempts | Min SE |
|---|---|---:|---:|---:|
| Beginner | 180 deg | 45-55 | 30 | — |
| Easy | 180 deg | 36-45 | 50 | — |
| Medium | 180 deg | 32-38 | 100 | — |
| Intermediate | 180 deg | 28-34 | 150 | — |
| Hard | 180 deg | 24-30 | 200 | — |
| Expert | 180 deg | 22-26 | 500 | 3.0 |
| Master | 180 deg | 20-24 | 1000 | 4.5 |
| Extreme | none | 17-22 | 2000 | 6.0 |

The SE floor for Expert/Master/Extreme ensures meaningful differentiation between the upper tiers. The difficulty tier check already caps the upper bound (a puzzle rated "Expert" can only use fish-level techniques, so its SE is naturally <= 5.4).

### Requested vs Rated Difficulty

The generator accepts puzzles one tier easier than requested for generation speed. The FFI layer exposes both:

- `get_difficulty()` — the *requested* difficulty (what the user selected)
- `get_rated_difficulty()` — the *actual* rated difficulty of the generated puzzle

Both values are persisted through serialization. UIs can choose which to display — the iOS app currently uses the requested difficulty as the label with the SE numeric rating for precision.

### PuzzleId

The `PuzzleId` module (`src/puzzle_id.rs`) encodes puzzle parameters (seed, difficulty, variant) into short alphanumeric codes. This enables deterministic puzzle regeneration from a compact identifier, powering shareable puzzle links (e.g., on ukodus.com).

## Reproducibility

Use `Generator::with_seed(seed)` to generate deterministic sequences for tests and debugging.

## Minimal Usage

```rust
use sudoku_core::{Difficulty, Generator, Solver};

let mut gen = Generator::with_seed(42);
let puzzle = gen.generate(Difficulty::Expert);

let solver = Solver::new();
assert!(solver.has_unique_solution(&puzzle));

let rated = solver.rate_difficulty(&puzzle);
let se = solver.rate_se(&puzzle);
println!("Requested: Expert, Rated: {rated}, SE: {se}");
```

