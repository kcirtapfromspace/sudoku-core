# Arithmetic Counting

`Technique::ArithmeticCounting` combines the engine's original exactly-one
requirements with small integer weights. It returns a placement or elimination
only when a separate certificate checker reconstructs the equations and verifies
a contradiction. It runs after Death Blossom and before forcing chains in the
normal hint and technique-profile paths.

## Soundness theorem

Give every allowed candidate a Boolean indicator `x_j`: one if it is selected in
a completion, zero otherwise. A placed cell is represented by its singleton
value, irrespective of its stored pencil marks. The 81 cell requirements and
243 row/column/box digit requirements are exactly-one equations:

```text
B x = 1,    x_j ∈ {0,1}.
```

Select distinct original equations and integer weights `w_i`. Their weighted sum
is `a·x = β`, where `a = wᵀB` and `β = Σw_i`. To prove that candidate `t` takes
value `v`, assume the opposite `b = 1−v`. The remaining variables must satisfy:

```text
Σ(j≠t) a_j x_j = R,       R = β − a_t b.
L = Σ(j≠t) min(0,a_j),    U = Σ(j≠t) max(0,a_j).
g = gcd(|a_j| : j≠t).
```

If `R` is outside `[L,U]`, the assumption is impossible: each Boolean variable
contributes either zero or its coefficient. If `g` does not divide `R`, the
assumption is also impossible: every integer combination of the remaining
coefficients is divisible by `g`. Zero divides only zero. These arguments prove
the reported placement (`v=1`) or elimination (`v=0`) in every completion
respecting the candidate state. They require no uniqueness assumption.

Version two adds a Boolean residue terminal. For a modulus `q` from 2 to 16,
start with the reachable set `{0}`. For each remaining candidate coefficient
`a_j`, replace the set with `R ∪ {(r+a_j) mod q : r∈R}` using the previous set.
If the required remainder is absent, the opposite candidate value is impossible.
Each variable contributes at most once: repeated coefficients still get separate
updates, and negative coefficients use nonnegative remainders. The implementation
uses a 16-bit set, with bit rotations for each coefficient. Version two also
requires the claimed endpoint's remainder to be reachable; it rejects a
combination that refutes both candidate values.

Overlapping requirements are allowed. A candidate present in two selected
equations gets both coefficients. Treating two overlapping sectors as independent
placements would be unsound; the checker always reconstructs the actual incidence.

## Examples exercised by the tests

### Five-equation parity (a familiar Guardian)

```text
a+b=1; b+c=1; c+d=1; d+e=1; e+a+t=1.
```

Adding the equations gives `2(a+b+c+d+e)+t=5`. If `t=0`, the remaining coefficients
have gcd two, which cannot divide five. Therefore `t=1`.

The fixture uses digit 1 at `a=r1c1`, `b=r1c4`, `c=r4c4`, `d=r4c2`,
`e=r2c2`, `t=r3c3`. Sources are row 1, column 4, row 4, column 2, and box 1.
Other positions for digit 1 in those sources are excluded in the fixture.

### Six equations with a coefficient of two

```text
t+a+d=1; a+b+e=1; b+c=1; c+d=1; d+f=1; e+f=1.
```

Weights `(2,−1,1,−1,−1,1)` give `2t+a=1`. If `t=1`, the remaining sum must be
`−1`, outside `[0,1]`. Therefore `t=0`. The test embeds these equations in actual
Sudoku sectors and independently checks compatible full completions.

### Six equations that require the new terminal

```text
g+p0+p4=1; p0+p1=1; p1+p2=1; p2+p3=1; p3+p4=1; g+t+u=1.
```

Multiply the first five equations by two and add the last:
`4(p0+p1+p2+p3+p4)+3g+t+u=11`. Assuming `t=1` requires the remaining
sum to be ten. Its possible remainders modulo four are `{0,1,3}`, while ten
has remainder two. Therefore `t=0`. The interval `[0,24]` and gcd one both
miss this contradiction. The native fixture uses digit 1 in column 4, row 4,
column 2, box 1, and row 2, followed by the cell rule for r1c4 `{1,2,3}`.
It removes 2 from r1c4. Tests validate the actual source incidence and a full
completion, and independently establish that the opposite has no completion
under those exact masks.

Run the [complete example](../examples/arithmetic_counting.rs) with
`cargo run --release --example arithmetic_counting` to check a Guardian certificate
and discover a deduction through the public search API.

## Using the engine

Ordinary `get_hint` and technique solving try this method in their existing
technique order. Like those APIs' other techniques, it receives the current
internal candidate state after the normal initial candidate rebuild.

For a caller that already has candidate eliminations, use the dedicated APIs.
They preserve the stored masks and never recalculate candidates:

```rust
use sudoku_core::{ArithmeticSearchOptions, Grid, ProofCertificate, Solver};

let grid = Grid::new_classic(); // Or a puzzle with a current candidate state.
let solver = Solver::new();
let result = solver.search_arithmetic(&grid, &ArithmeticSearchOptions::default());
if let Some(hint) = result.hint {
    if let Some(ProofCertificate::Arithmetic(proof)) = &hint.proof {
        assert!(proof.verify(&grid));
        println!("{}", hint.explanation);
    }
}
// solver.get_arithmetic_hint(&grid) is the convenience form without search stats.
```

`ArithmeticProof::from_terms` accepts a grid, source terms, cell index, digit,
and result value. It returns a proof only if it verifies. `proof.check(&grid)`
returns the reconstructed interval or divisibility contradiction. For residue
proofs, use `ArithmeticProof::from_residue_terms(grid, terms, cell, digit, value,
modulus)`, or `from_terms_with_terminal` with an explicit `ArithmeticTerminal`.
The resulting `ArithmeticCheck::Residue` includes the residual sum, modulus,
required remainder, and the independently recomputed reachable remainders. Cells are
zero-based row-major indices; sectors 0..8 are rows, 9..17 columns, and 18..26
boxes. Digits are 1..9.

Certificates carry a version and SHA-256 hash of all normalized domains and
placed values. Changing an unrelated candidate also invalidates the certificate.
The hash binds a proof to its state; it is not an authenticity signature.
A proof returned after the main solver has rebuilt candidates or chained earlier
eliminations belongs to that internal state. Replaying it requires the same
state; the dedicated arithmetic APIs make this explicit by using the caller's
stored masks.

`ArithmeticProof` can be serialized separately. As with existing certificates,
`Hint.proof` is skipped when a `Hint` is serialized.

Use `ArithmeticReplay::capture(&grid, &proof)` to export a one-step envelope
containing the full normalized candidate state and certificate. Serialize the
envelope explicitly, then call `check()` or `verify()` after deserialization.
`apply(&mut grid)` requires identical premises, independently checks the proof,
and preserves earlier candidate eliminations. Errors leave the grid unchanged.
`ArithmeticState::capture` and `to_grid` also expose the state representation
directly. The recorded rules are the native Sudoku projection: extra variant
constraints are retained when applying to an existing variant grid, but are not
serialized in the native proof envelope. Record engine provenance and search
options alongside exported envelopes at the application layer.

### Version compatibility and parity compilation

Old JSON has `version: 1` and no terminal field. It still deserializes to the
interval/gcd terminal, retains its original weight limit of two, and serializes
without an added field. Residue JSON uses `version: 2` and
`"terminal": {"Residue": {"modulus": 4}}`. Mismatched versions and terminals,
unsupported moduli, and missing version-two terminal fields fail verification.
The domain hash format remains the original one; proof semantics are determined
by the separately validated version and terminal.

`parity_proof.compile_parity_tail(grid, tail_requirement)` compiles an independently
checked unit-weight parity placement plus one distinct incident native rule into
direct residue exclusions. If parity forces `g=1` and the tail has `g` and `L`
other candidates, multiplying the parity combination by `L` and adding the tail
proves each other candidate false modulo `2L`. Native tails have at most nine
candidates, so `L≤8`, weights are at most eight, and the modulus is at most 16.
The combined proof must still fit six distinct sources. Unsupported inputs
return an empty vector, and every compiled proof is independently checked.

## Search limits and integration choices

- Defaults: at most six sources, weights `±1` and `±2`, beam width 64, and 30,000
  evaluated combinations. Options permit beam widths 1..512 and combination
  budgets 0..1,000,000. A zero budget evaluates nothing. At most 8,192 child
  combinations are retained at once, independently of the evaluation budget.
  Discarding children at this cap also sets `beam_pruned`.
- The deterministic search expands connected source groups, identifies equivalent
  global sign choices, and keeps promising partial combinations with few odd and
  nonzero coefficients. It reserves evaluation budget for deeper groups. Each
  combination tries interval/gcd checks first, then moduli 2..16. Remainder sets
  are cached by coefficient residue when testing multiple targets.
- A discovered parity placement can be compiled with an incident tail before
  returning its hint. These extra combinations count against the same budget
  and respect the configured source and weight limits. The public compiler can
  produce weights up to eight; the beam and its automatic compiler remain capped
  at the configured weight of one or two.
- `tested_combinations`, `budget_exhausted`, and `beam_pruned` describe the work
  performed. An absent hint is inconclusive even if the budget is not exhausted:
  this search is incomplete, and these bounded certificates do not express every
  Boolean consequence. `tested_combinations` counts source combinations, not
  individual target/modulus checks.
- Certificate verification is independent of the search budget and node caches.
  It accepts at most six distinct sources with nonzero weights of magnitude at
  most two for version one, or eight for version two. Invalid versions, terminal
  choices, moduli, IDs, weights, duplicates, targets, hashes, candidate
  masks, and locally invalid states are rejected.
- Standard constraints must actually be installed. Built-in Classic, X-Sudoku,
  and Killer grids qualify; added variant constraints only restrict completions
  further. Arbitrary custom constraint lists are conservatively unsupported.
  Deserialized grids need `restore_constraints()` before use.
- Candidate masks are logical premises. The checker does not establish that a
  previous elimination was correct or that the whole state is satisfiable.
- Arithmetic search is omitted from the propagation subset used recursively by
  forcing chains, to avoid running a beam search at every assumption.
- The technique is classified Extreme. Its `se_rating()` value of **8.5** is an
  engine-local, uncalibrated estimate; it is not an official Sudoku Explainer
  technique rating. Puzzle ratings and generated results can change when this
  technique replaces a harder step.

## Validation and attribution

Tests reconstruct X-Wing, Guardian, and six-source weighted certificates; reject
an overlapping-fish trap and stale/malformed proofs; and exercise placement,
serialization, deterministic search, state preservation, and bounded work. A
separate exact-cover oracle uses the original candidate masks to check both a
compatible full completion and the impossibility of the opposite result.
Residue tests additionally compare the production DP with independent Boolean
enumeration on 58,590 coefficient-list/modulus pairs, check all native tail
sizes 2..9, exercise signed and repeated coefficients, and preserve version-one
JSON and checker behavior. The default public search discovers a modulus-four
exclusion on the six-source fixture within its existing budget.

This implementation makes no historical-novelty claim. Exactly-one candidate
incidence, weighted equations, interval bounds, divisibility, and Guardian
reasoning have established antecedents. Relevant primary and specialist sources:

- [Dukes–Nimegeers, Sudoku incidence and fractional completion](https://alco.centre-mersenne.org/item/10.5802/alco.375.pdf).
- [The Lambda manuscript: linear combinations and interval tests](https://www.cs.purdue.edu/homes/ci/draft/Lambda.ps).
- [Wolfe–Tseng, the Power test](https://engineering.purdue.edu/~eigenman/ECE663/Handouts/The_Power_Test_for_data_dependence.pdf), including limitations of interval/gcd tests.
- [Andrew Stuart's Guardian explanation](https://www.sudokuwiki.org/Print_Guardians).
- [Pesant–Meel–Mohammadalitajrishi, modular subset-sum nonreachability](https://www.cs.toronto.edu/~meel/Papers/cpaior21-pmm.pdf).

The implementation theorem above is the elementary soundness argument used by
this engine. Broader source-count completeness and priority questions are separate
research claims, and are not required for accepting these certificates.

## Earlier version-one acceptance results — 15 September 2026

On the implementation branch based on `361c254`:

- Required `./scripts/test-coverage.sh`: **183 Rust tests and 9 collector tests
  passed**, with no failed or ignored tests and no missing production mappings.
- Production line coverage: **95.74%**; Arithmetic Counting: **97.77%**.
  Every configured coverage gate passed.
- `cargo fmt --check`, the CI strict Clippy command, WebAssembly target checking,
  and the runnable arithmetic example passed. The example discovers and verifies
  the Guardian placement after 11,358 evaluated combinations.
- An additional strict Clippy run with `--all-targets` reported six pre-existing
  test-code warnings: type complexity in `tests/variant_constraints.rs`, unit
  struct default construction in `tests/solver_behavior.rs`, three Boolean
  simplifications in `src/solver/technique_tests.rs`, and iterator flattening in
  `src/solver/behavior_tests.rs`. These are outside the required CI lint command;
  the new example also passed its own strict Clippy check.

Two indicative release comparisons of Extreme generation produced the same
unique 22-clue puzzles with and without arithmetic dispatch. Seed 42 took
10.15 s versus 9.13 s; seed 99 took 44.26 s versus 41.78 s. Concurrent coverage
limits timing precision; these two cases are not a general performance guarantee.
WebAssembly runtime latency has not been measured.
