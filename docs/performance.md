# Generation and uniqueness performance

Clue removal calls the solution counter repeatedly. Built-in classic grids now
use a stack-allocated counter with nine occupancy masks each for rows, columns,
and boxes. It chooses an empty cell with the fewest candidates, places a digit,
and restores the masks after each branch. It stops as soon as it reaches the
requested solution cap; uniqueness checks still require exactly one completion.
It reads placed values directly and ignores stored pencil marks and given flags,
without modifying the input grid.

The optimization applies only when the grid's installed constraints are known
to be exactly the built-in classic rules. X-Sudoku, Killer, custom constraints,
and deserialized grids awaiting constraint restoration retain the existing
constraint-aware implementation. Neither a Classic variant tag nor matching
constraint names enables the specialization.

Generation also avoids rebuilding a copy and rewriting all givens for every
removal. It tracks the remaining clue count, rejects symmetry groups that would
exceed the removal budget before searching, and obtains difficulty and SE ratings
from one `Solver::analyze` call. Symmetries, uniqueness requirements, rating
thresholds, attempt budgets, and candidate fallback rules are preserved.

## Local comparison

Measured on 2026-10-19 on Apple Silicon macOS with Rust 1.99.0, using release
builds. The baseline is commit `84696bec4003d129a2c11f12610cd1fb3cdb8c13` with
the same performance example. Each value below is the median of three alternating
before/after runs. Each generation run averages seeds 42, 99, and 314; counting
uses the fixed fixtures and iteration counts in `examples/performance.rs`.
Compilation and post-generation correctness checks are outside the timed region.

| Generation target | Before (ms/puzzle) | After (ms/puzzle) | Speedup |
| --- | ---: | ---: | ---: |
| Beginner | 5.68 | 2.82 | 2.0× |
| Easy | 2.44 | 0.84 | 2.9× |
| Medium | 3.61 | 0.89 | 4.1× |
| Hard | 67.67 | 37.08 | 1.8× |
| Expert | 292.19 | 84.23 | 3.5× |

| Counting workload (cap 2) | Before (µs/puzzle) | After (µs/puzzle) | Speedup |
| --- | ---: | ---: | ---: |
| Easy unique fixture | 27.96 | 2.36 | 11.8× |
| Hard unique fixture | 12,428.09 | 1,557.42 | 8.0× |
| Two-completion rectangle | 70.94 | 0.19 | 373× |
| Empty classic grid | 879.77 | 12.84 | 68.5× |

These are local measurements of the listed workloads, not latency guarantees for
every seed, difficulty, variant, or platform. Advanced logical rating searches
remain part of generation cost. All 15 generated benchmark puzzles have the same
fingerprints before and after the changes.

## Reproduce and validate

```sh
cargo run --release --locked --example performance
cargo run --release --locked --example performance -- count
cargo run --release --locked --example performance -- generate
cargo test --all-targets --all-features --profile coverage --locked
```

Use the same example, compiler, release profile, and otherwise idle machine when
comparing revisions. The benchmark checks reference counts, generated puzzle
validity, and uniqueness, and prints fingerprints for seeded output comparisons.

Regression tests compare the optimized counter against the general search on
96 deterministic boards at six clue densities with caps 0, 1, 2, and 3 (384
comparisons). They also check invalid values, row/column/box conflicts, impossible
boards, input preservation, variant/custom/deserialization dispatch, and the
15 pre-optimization seeded generation outputs. The full soundness suite remains
required; run `./scripts/test-coverage.sh` to enforce production coverage gates.
