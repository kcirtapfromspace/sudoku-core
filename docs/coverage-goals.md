# Production coverage goals

## Required thresholds

These targets measure the authored Rust engine in this repository. They include
every production implementation; callers such as Ukodus have separate coverage.

| Area | Minimum line coverage |
| --- | ---: |
| Entire production crate | 75% |
| Foundations | 90% |
| Solver orchestration, basic techniques and backtracking | 90% |
| Generator and diversity | 80% |
| Advanced technique engines, combined | 75% |
| Each advanced technique engine | 60% |

Foundations includes grid, cells, bitsets, positions, constraints, hashing, puzzle
IDs, candidate fabric, technique types and proof explanations. Advanced engines
are AIC, ALS, Arithmetic Counting, fish and uniqueness. `scripts/coverage.py`
lists the exact files for each group. Any new production source remains part of
the crate-wide gate and must have coverage mappings, even before it is assigned
a group.

The higher targets protect puzzle validity, solution counting, data preservation
and the public solver paths. The advanced engine targets require broad exercise
without pretending that executing every rare combinatorial branch establishes
correctness. These percentages are minimums; useful behavior and regression tests
remain required even when coverage already passes.

Rust branch coverage is not claimed. Technique soundness, solution preservation,
generation reproducibility and termination need explicit assertions in addition
to line coverage. Tests must not be weakened, production files excluded or test
code counted merely to satisfy a percentage.

## Measured results — 2026-09-15

The full local acceptance run passed every gate:

| Area | Original | Revised | Minimum |
| --- | ---: | ---: | ---: |
| Entire production crate | 59.09% | **95.66%** | 75% |
| Foundations | 60.36% | **96.04%** | 90% |
| Core solver | 68.74% | **96.94%** | 90% |
| Generator and diversity | 80.52% | **96.75%** | 80% |
| Advanced engines, combined | 52.10% | **94.97%** | 75% |
| AIC | 46.67% | **98.65%** | 60% |
| ALS | 39.70% | **98.97%** | 60% |
| Fish | 68.49% | **98.54%** | 60% |
| Uniqueness | 61.28% | **84.43%** | 60% |

The original baseline is commit
[`ad8f024`](https://github.com/kcirtapfromspace/sudoku-core/commit/ad8f024d507a52eff99fdd8b5173763487b30a31):
3,511 of 5,942 production lines covered, with all 70 original Rust tests passing.
The revised run covered 5,658 of 5,915 production lines and passed all 159 Rust
tests plus nine coverage-collector checks. Neither run ignored or filtered tests.
Formatting, strict Clippy, and the WebAssembly build check also passed locally.

Both measurements count unique production LCOV lines. The original binary and
source snapshot were preserved and verified against the baseline commit; its
trailing inline test modules were excluded from the line counts. The original
run used the `release` profile and the revised run uses `coverage`, which retains
debug assertions and overflow checks. Source changes and profile differences
can change the denominator. These results are a recorded local measurement;
current CI artifacts remain the source of truth for later revisions.

## Behavioral checks and regressions

The added tests check more than whether a search returns a result:

- Grid edits preserve valid candidate sets; cloning preserves custom constraints
  while keeping cell state independent. Bitsets, positions, hashes and puzzle IDs
  have boundary and round-trip checks.
- Reference puzzles have asserted solutions and solution counts. Hint chains must
  preserve those solutions and finish; difficulty and technique profiles agree.
- Fish, ALS, AIC and uniqueness fixtures check deductions and proof metadata.
  Independent Boolean, exact-cover or local-assignment checks support the relevant
  pattern assertions. Uniqueness techniques use genuinely unique puzzle fixtures.
  These technique fixtures exercise classic Sudoku; variant cloning and
  backtracking checks do not establish every variant's hint behavior.
- Generator checks cover all supported clue symmetries, exact clue-count bounds,
  seed reproducibility and bounded-search fallback. Diversity statistics are
  checked against the samples they summarize.

The tests exposed and fixed lost custom constraints during cloning, candidate
restoration errors, incomplete symmetry groups, minimum-clue boundary errors,
and unsound AIC, fish, ALS, BUG and Hidden Rectangle deductions. Empty Rectangle
detection and several proof descriptions were also corrected. The solver
specification records the assumptions behind the revised deductions.

Fixed-seed regression tests verify reproducibility for their fixture cases.
Corrected clue removal and difficulty deductions can change which puzzle a seed
produces compared with older revisions; seed-only puzzle IDs do not encode the
generator version.

Fish searches reuse sector masks and combinations. Forced-single propagation
keeps its initial candidate rebuild, then updates candidates affected by each new
placement. Differential comparisons against the previous implementations checked
finding metadata and complete grid state; these changes reduce the cost of
running the full suite.

The original soundness tests silently skipped an unsolvable fixture and used an
ambiguous expert fixture. Their replacements assert uniqueness and solvability;
stored-candidate chains now assert completion instead of silently reaching a step
limit. The generated-puzzle cases and public placement-hint checks remain required.

## Run the checks

Use Rust 1.95.0 and Python 3.11 or newer. Install coverage tooling once:

```sh
rustup component add llvm-tools-preview
cargo install cargo-llvm-cov --version 0.9.1 --locked
```

Run from the repository root:

```sh
cargo fmt --check
cargo clippy --all-features --locked -- -D warnings
./scripts/test-coverage.sh
```

The helper runs the entire test suite, including soundness tests, under the
`coverage` profile. This profile optimizes expensive solver/generator paths while
retaining the test profile's debug assertions and overflow checks. Tests run
sequentially to avoid contention on shared instrumentation counters during
solver searches. No tests are excluded by name. The dedicated
`target/coverage-build` directory keeps these
builds separate from ordinary development builds and saved baseline results.

To run tests without coverage instrumentation:

```sh
cargo test --all-targets --all-features --profile coverage --locked
```

## Measurement and artifacts

The collector counts each unique production LCOV source line once. A line is
covered if any instrumented execution hits it; zero-hit lines remain in the
denominator. The source inventory comes from the Cargo package/workspace and
`src/**/*.rs`, independently of LCOV. Missing function or file mappings fail the
gate. Files containing only declarations may have no executable mappings.

Existing inline test modules were moved into separate `*_tests.rs` files with
their private parent access preserved. Only test directories, `tests.rs` and
`*_tests.rs` files are excluded. Inline test modules are rejected so assertions
and fixtures cannot inflate the production percentage. Dependency implementation
and examples are outside the authored production library scope.

After an acceptance run:

- `target/coverage/summary-production.json` contains authoritative percentages,
  pass/fail status, missing mappings and uncovered source lines.
- `target/coverage/lcov.info` contains raw Rust coverage records.
- `target/coverage/rust/html/index.html` is LLVM's browsable report. Its aggregate
  percentages use different accounting from the production gate.

Test the collector itself with
`python3 -m unittest discover -s scripts -p test_coverage.py -v`.

## CI

The `CI` workflow runs formatting, Clippy, all tests including soundness, and every
production coverage gate on pull requests and pushes to `main`. It uploads reports
even when checks fail. The required acceptance statuses are `Lint & Format` and
`Tests and production coverage`. Repository branch protection must require those
statuses separately; a workflow definition alone does not configure protection.

Use current generated reports as the source of truth. Executable line counts can
change with compiler versions, formatting and production changes.
