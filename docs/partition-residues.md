# Partition residue certificates

This extension preserves native exactly-one choices while testing a weighted
Sudoku equation modulo a small integer. It adds a version-three proof terminal,
an opt-in bounded search, and portable replay through the existing
`ArithmeticReplay` API. It requires neither guessing a solution nor assuming
that the puzzle has a unique solution.

The code and the arguments below were independently authored for this engine.
**Historical mathematical novelty is not established.** Finite-domain modular
propagation has substantial prior art. The contribution here is an explicit
Sudoku certificate format, its reconstructed checker, partition-aware search,
and a reproducible comparison with this engine's earlier proof formats. The
small regression patterns below are known Sudoku deductions.

## Premises and notation

Let `V` be the available candidate variables in the recorded grid. A candidate
`x_(c,d)` is one when cell `c` contains digit `d`, and zero otherwise. Placed
cells contribute only their placed digit. The native cell and sector/digit rules
give 324 equations of the form

```text
Σ(j in E) x_j = 1.
```

Combine between one and six distinct native equations with nonzero integer
weights of magnitude at most eight. Reconstruct their sum as

```text
Σ(j in V) a_j x_j = β.
```

The version-three terminal also selects one of four partitions of `V` into 81
disjoint exactly-one groups:

| API variant | One group per | Choices within the group |
| --- | --- | --- |
| `Cells` | Cell | Available digits in that cell |
| `RowDigits` | Row and digit | Available cells in that row for that digit |
| `ColumnDigits` | Column and digit | Available cells in that column for that digit |
| `BoxDigits` | Box and digit | Available cells in that box for that digit |

These are four separate partitions, not one collection of overlapping groups.
Each available variable occurs in precisely one group of a selected partition.
The checker reconstructs these groups; certificates cannot supply custom groups.

Zero-coefficient alternatives **must remain in the groups**. They represent the
valid choice to contribute zero instead of a different candidate's coefficient.
The target may itself have coefficient zero or lie outside every weighted source.

For a target `t` and an assumed Boolean value `b`, restrict its group to `{t}`
when `b=1`, and delete `t` from its group when `b=0`. Leave other groups as they
are. Call the resulting groups `G_i(b)`. An empty group means no permitted choice.

## Proof 1: exact recurrence and sound deduction

Choose an integer modulus `q` between 2 and 16, inclusive, and define

```text
C_i(b) = { a_j mod q : j in G_i(b) }
R_0(b) = {0}
R_(i+1)(b) = { (r+c) mod q : r in R_i(b), c in C_i(b) }.
```

**Lemma.** `R_k(b)` is exactly the set of residues obtainable by choosing one
candidate from each of the first `k` restricted groups and adding its coefficient.

**Proof.** The empty selection has sum zero. If the statement holds for `k`,
every selection for `k+1` consists of a selection for the first `k` groups plus
one member of the next group. Its residue is therefore in the displayed
recurrence. Conversely, each residue in that recurrence has a witnessing
selection for its first summand by induction, and a witnessing next-group
candidate for its second summand. Combining those choices gives the residue.
An empty next group has no such choice, so the result is empty. This proves
both directions by induction. Duplicate coefficients may be collapsed within
a group, because only one member is selected. Different groups must still be
processed separately. ∎

**Soundness theorem.** If `β mod q` is absent from `R_81(b)`, every Sudoku
completion respecting the recorded values and candidates has `x_t != b`.

**Proof.** Suppose a respecting completion instead has `x_t=b`. Every weighted
source equation holds in that completion, so its full weighted sum is `β`.
It selects exactly one candidate in every native partition group, including the
restricted target group. The lemma therefore places `β mod q` in `R_81(b)`, a
contradiction. Thus `b=1` proves an elimination and `b=0` proves a placement. ∎

The implementation additionally requires `β mod q` to be present under the other
target value. This refuses a deduction when this relaxation refutes both values.
Presence in a residue set does **not** prove Sudoku satisfiability: the relaxation
does not enforce the other three partitions or every installed variant rule.

The candidate masks remain premises. A proof does not validate earlier deletions,
establish uniqueness, or certify that a completion exists. Grids must have the
standard constraints installed; additional variant constraints only restrict
the set of completions further.

## Proof 2: relation to the earlier Boolean residue check

For fixed `a`, `β`, `q`, `t`, and `b`, consider the earlier relaxation that allows
every variable other than the assumed target to be chosen independently.
Restore the target contribution to its residual-sum residues so both checks
describe full sums.

**Containment theorem.** The partition attainable set is a subset of that
independent Boolean attainable set.

**Proof.** Any selection of one candidate in each disjoint partition group
defines a Boolean assignment: set exactly those selected variables to one.
It obeys the target assumption and has the identical weighted sum. That
assignment is among those considered by the independent Boolean relaxation.
Thus every partition-attainable residue is Boolean-attainable. ∎

Consequently, absence in the Boolean relaxation implies absence in the
partition relaxation. On a satisfiable candidate state, any respecting
completion also supplies the accepted endpoint required by the version-three
checker. On inconsistent states, the two-endpoint guards can prevent comparison
of certificate acceptance even though the set containment still holds.

**Strictness witness, for the same source equations.** In row 1, allow digit 1
and digit 2 only in cells `A=r1c1` and `B=r1c4`. Give `A` candidates `{1,2,3}`
and `B` candidates `{1,2}`; other cells outside the stated row restrictions have
all candidates. Adding the row/digit equations gives

```text
x_A1 + x_A2 + x_B1 + x_B2 = 2.
```

For target `x_A3`, the coefficient is zero. The independent interval, gcd, and
Boolean residue tests cannot distinguish its endpoints, for any modulus.
In the cell partition, assuming `x_A3=1` makes `A` contribute zero and `B`
contribute exactly one. The full sum is odd, contradicting `2 mod 2 = 0`.
With `x_A3=0`, both cells contribute one and the required residue is reachable.
The certificate therefore excludes 3 from `A`.

This pattern is a hidden pair. It demonstrates a strictly stronger certificate
terminal for the same weighted equation, not a newly discovered Sudoku pattern,
nor a proof that the entire existing solver cannot find the step another way.
Independent exact-cover tests verify a respecting complete Sudoku and the
absence of any respecting completion with `A=3`.

Modulus choice matters: two source rows forming an X-Wing and a column/digit
partition can exclude a candidate using modulus 3, while modulus 2 aliases the
required sum 2 with the possible sum 0. The checker correctly refuses the
modulus-two certificate in that fixture.

## API and search

```rust
use sudoku_core::{ArithmeticPartition, ArithmeticSearchOptions, Grid, Solver};

let grid = Grid::new_classic(); // Supply the current grid with its stored masks.
let result = Solver::new().search_partition_arithmetic(
    &grid,
    &ArithmeticSearchOptions::default(),
    ArithmeticPartition::Cells,
    3,
);
assert!(result.error.is_none());
```

To check an explicit source combination, use
`ArithmeticProof::from_terms_with_terminal` with
`ArithmeticTerminal::PartitionResidue { modulus, partition }`. Accepted proofs
have `version=3`; `proof.check(&grid)` reports a `PartitionResidue` check with
the full right-hand side and the remainders under the rejected assumption.
`ArithmeticReplay::capture` exports the exact values and masks alongside the
proof. Replay independently reconstructs and verifies the state and can apply
the step atomically to that identical state.

Run the complete, deterministic strictness example and export its replay:

```sh
cargo run --release --example partition_residues > partition-proof.json
```

Search evaluates one selected partition and modulus per call. It shares the
existing limits: at most six weighted sources, search weights `±1,±2`, beam
width 1..512, and at most 1,000,000 evaluated combinations. The default budget
is 30,000. Checking explicit certificates supports weights through `±8`.
Search adjacency also follows candidates in the same partition group, so it
can connect disjoint source equations coupled by an exclusive choice.
Prefix/suffix residue sets accelerate candidate checks; the certificate checker
reconstructs the equation and recomputes both assumptions without those caches.

Search is deterministic and incomplete. Budget exhaustion and beam pruning are
reported; no hint is not a proof that no such deduction exists. The budget
counts combinations rather than elapsed time. This initial release is opt-in:
the automatic hint ladder, ratings, and generator retain their existing behavior.
No claim of a runtime improvement or complete Sudoku-solving method is made.

Versions one and two retain their existing JSON and checking semantics. The new
public enum variants require downstream exhaustive Rust matches to be updated;
older engine binaries cannot replay version-three certificates.

## Verification and research provenance

The tests exhaustively compare grouped residue addition with integer Cartesian
enumeration, including empty groups, equal alternatives, signed coefficients,
zero coefficients, and every supported modulus. They test all four partitions
against native requirements, cross-check concrete deductions with the existing
independent exact-cover oracle, verify deterministic bounded discovery, and
exercise JSON replay, stale-state rejection, and malformed proofs.

[`verification/partition_residues.rs`](../verification/partition_residues.rs)
also machine-checks the soundness invariant in Verus: any residue fold containing
zero initially and every permitted next-group transition contains the residue
of every allowed full selection. Absence of the required residue therefore
refutes every selection. The production bitset operations and the Sudoku-to-group
reconstruction are not formally verified by this file; the tests exercise those
remaining interfaces. The new proof runs with the existing `scripts/verify.sh`
and pinned CI verifier, alongside the earlier arithmetic proofs.

Historical novelty needs a separate research assessment. A finite web search
cannot establish “100% original,” and this review found related methods:

- [Dukes and Nimegeers, The linear system for Sudoku and a fractional completion
  threshold](https://alco.centre-mersenne.org/item/10.5802/alco.375.pdf)
  studies Sudoku incidence equations. Candidate-equation modeling is established.
- [Pesant, Meel, and Mohammadalitajrishi, On the Usefulness of Linear Modular
  Arithmetic in Constraint Programming](https://www.cs.toronto.edu/~meel/Papers/cpaior21-pmm.pdf)
  discusses domain filtering and dynamic programming for finite-domain modular
  constraints. Treating each group as a finite-domain coefficient choice fits
  that established setting; the general recurrence is not claimed as new.
- [Choco's finite-domain constraint documentation](https://choco-solver.org/docs/modeling/intconstraints/)
  documents sum, modulo, and table constraints, which are further related
  formulations rather than evidence of priority for this specific engine format.

The search and this comparison were performed on 19 October 2026. These sources
are conceptual references, not copied code or proof text. This extension is
released under the repository's existing MIT license. Any future historical
priority claim should compare the precise certificate construction against
earlier constraint-programming and Sudoku work and undergo independent review.

## Acceptance results — 19 October 2026

- The required `scripts/test-coverage.sh` passed all **211 Rust tests** and all
  **nine coverage-collector tests**, with no failed or ignored tests. All
  production coverage gates passed: **96.07%** overall and **98.65%** for the new
  partition module. This complete run used engine base `f59f7c8`, including the
  newly merged performance changes and the modulus-two alias regression.
- Grouped residue addition matched exhaustive integer enumeration on **166,665**
  group-list/modulus cases.
- All four files in `scripts/verify.sh` passed with the pinned CI Verus release
  `0.2026.09.13.671956e`: 45, 40, 26, and 3 verified items, respectively.
- `cargo fmt --check`, CI's strict Clippy command, strict example Clippy,
  `cargo check --target wasm32-unknown-unknown --locked`, and the runnable
  JSON-export/replay example passed. Local Rust checks used Rust 1.99.0; Verus
  used its required Rust 1.98.1 toolchain.

These are local acceptance measurements; the PR's CI checks provide their own
results on the pinned CI toolchains. No solver-speed or puzzle-difficulty
improvement is inferred from the coverage measurements.
