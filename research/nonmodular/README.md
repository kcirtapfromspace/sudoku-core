# Non-modular research experiments — 19 October 2026

Three agents investigated structural matching, nonlinear certificates, and
trades/variant interactions. This directory preserves a concrete inference
fixture, exact verification scripts, and the limits of the engine comparison.
It is research code under the repository's MIT license, not an installed solver
technique. No historical mathematical novelty is established.

## Findings and priorities

| Priority | Research target | Precise question |
| --- | --- | --- |
| 1 | Three-house matching obstructions | Which minimal row/column/box contradictions survive the existing bounded, non-forcing detectors? Can recurring cases be given short certificates? |
| 2 | Four-digit and mixed-list obstructions | Which small native patterns require exact quadratic or quartic certificates beyond existing rules and known oddagon/tridagon families? |
| 3 | Mixed variant interactions | Which small contradictions or safe exchanges require simultaneous Killer, Thermo, diagonal, and classic constraints? |

These are hypotheses about useful narrow research gaps. Matching, polynomial
certificates, and constraint interchangeability themselves have substantial
prior art. A search that finds few papers does not establish that an area is
unexamined.

## A concrete three-house inference

[`three-house-847.json`](three-house-847.json) records all 81 candidate domains,
using bit `digit-1`, and a respecting full completion. It is a candidate-state
fixture, not a givens-only puzzle. Its conclusion is **r1c4 cannot be 1**.

An independent oracle enumerates every permutation of digits 1–9 for row 1,
column 1, and box 1 after assuming r1c4=1:

| Necessary condition | Result |
| --- | --- |
| Row 1 alone | 28 assignments |
| Column 1 alone | 300 assignments |
| Box 1 alone | 180 assignments |
| Row 1 and box 1 together | Feasible |
| Column 1 and box 1 together | Feasible |
| Row 1 and column 1 together | Feasible |
| All three together | **No compatible assignment** |

The checker joins actual ordered shared-cell values, not just marginal digit
counts. Every full solution must restrict to a compatible three-house assignment,
so the absent join refutes the assumed candidate. Two explicit full completions
respect the unforced state. No uniqueness premise is used.

### A short proof of the same result

The relevant original domains are:

```text
r1c1 = {1,2,3,4,8}    r1c2 = {1,2,3}    r1c3 = {1,2,3}
r2c1 = {2,3,4}        r3c1 = {2,4,7,8,9}
r5c1 = {3,4,5,8}      r6c1 = {4,5,7,8}
r7c1 = {3,4,9}        r9c1 = {3,4,5,7,9}
```

Assume r1c4=1. The two cells r1c2 and r1c3 must occupy 2 and 3. Both share
box 1 with r2c1, forcing r2c1=4. The row and box restrictions then force
r1c1=8. Therefore these five cells in column 1 have the following possibilities:

```text
r3c1: {7,9}
r5c1: {3,5}
r6c1: {5,7}
r7c1: {3,9}
r9c1: {3,5,7,9}
```

Five mutually distinct entries would have to use just four digits
{3,5,7,9}, contradicting Hall's necessary condition. Thus r1c4≠1. ∎

This argument uses ordinary matching and all-different constraints. It is an
independently written proof for a newly mined candidate state, not evidence
of a historically new solving principle.

### Engine comparison and its boundaries

At engine base `f59f7c839dd6a488c80c44010b6af624ed4580c2`, the stored state yields
no finding from the **34 explicitly listed non-forcing, uniqueness-independent
detectors** in [`probe_existing_engine.rs`](probe_existing_engine.rs).
The default Arithmetic Counting search also returns no hint after 19,215
combinations, with beam pruning reported. Search failure is not an impossibility
or completeness theorem.

Uniqueness techniques are excluded because the state has multiple completions.
Forcing chains, backtracking, optional partition-residue searches, stronger
external solvers, and unimplemented known techniques are outside this comparison.
In particular, this does not prove a separation from all known ALS, MSLS,
3D Sue de Coq, or general constraint-programming reasoning.

The public `get_hint()` rebuilds candidates, so it is unsuitable for comparing
this stored-mask fixture. The comparison script copies the source into a
temporary directory and adds a private test probe there. It leaves the checkout
unchanged and reports if a future engine version closes the gap.

The bounded discovery experiment used seed 83920 and 1,000 completion-preserving
random states, inspecting 1,080 row/column/box clusters. It found 33 eliminations
in 18 states beyond individual-house matching plus locked candidates. Against
the restricted existing detector baseline, 32 eliminations in 17 states remained.
These are biased synthetic samples, not prevalence estimates for generated games.
The fixture above is trial 847; the committed scripts reproduce its proof and
engine comparison, not the entire sampling run.

An earlier two-house fixture was rejected as a discovery: the existing ALS-XZ
handler already found the same elimination. This is why the comparison uses the
actual detector implementations rather than a singles/subsets-only baseline.

## Exact compression theorem for selected houses

Choose any family of native houses. Partition their cells into groups having
identical membership in those houses. Call these membership atoms. For each atom,
choose a digit set whose size equals its number of cells.

**Theorem.** The selected houses have a domain-respecting assignment if and only if:

1. Each atom has a candidate-respecting perfect matching to its chosen digit set.
2. Atoms in the same selected house have disjoint chosen digit sets.

**Proof.** Any valid assignment gives precisely these digit sets and matchings.
Conversely, independent atom matchings combine into one assignment because the
atoms partition the cells. Disjointness prevents repeats in every selected house.
Each house has nine cells over digits 1–9, so it contains every digit once. ∎

This is exact for the selected-house relaxation, not for the full Sudoku. It is
a general matching/constraint decomposition argument; priority is unverified.

For a cell's row, column, and box, the six atoms have sizes 1,2,2,6,6,4.
The boundary can be represented by a center digit and two disjoint two-digit
sets, giving at most `9 * C(8,2) * C(6,2) = 3,780` signatures. Failed matchings
can be certified by explicit Hall deficiencies. This is the proposed bounded
proof-search mechanism; a production checker is not included here.

## Nonlinear and geometric baselines

[`verify_quadratic_baseline.py`](verify_quadratic_baseline.py) checks a known
bivalue oddagon using an exact rational positive-semidefinite matrix. It verifies
a full completion before the target deletion and a fractional solution of all
324 native linear equations afterward. The fractional witness satisfies all
recorded domains. Thus the nonlinear proof excludes a state still feasible for
that full continuous linear relaxation. This is not a new oddagon pattern or a
comparison with every existing logical rule.

[`verify_graph_exclusions.py`](verify_graph_exclusions.py) checks why two tempting
graph templates do not supply new three-digit obstructions here:

- Among 2,376 diamond-pair Moser-spindle configurations at a representative cell,
  every native support contains four mutually seeing cells. Cell transitivity
  covers every choice of shared tip. Such a support restricted to three digits
  already has a Hall contradiction.
- Every cell's peer neighborhood is chordal: its six row peers outside its box
  are simplicial, followed by its six column peers outside its box; the remaining
  eight box peers form a clique. A wheel without a four-clique would require a
  triangle-free cycle in that neighborhood, which is impossible. The script
  verifies the complete elimination ordering as well as bounded cycle cases.

These independently derived geometric observations narrow the next search; no
historical novelty claim is made for them. They do not rule out the known
tridagon or other distributed obstructions.

## Run the checks

Python 3 and its standard library are sufficient for the mathematical checks:

```sh
python3 research/nonmodular/verify_three_house.py
python3 research/nonmodular/verify_quadratic_baseline.py
python3 research/nonmodular/verify_graph_exclusions.py
```

The comparison additionally requires Cargo and the engine's locked dependencies:

```sh
python3 research/nonmodular/compare_engine.py
```

Do not use Python's `-O` flag; the checkers reject it because they use assertions.
The artifacts are not a formal verification of a production implementation.
The three mathematical scripts run in the dedicated research CI workflow.
The comparison remains manual: a future detector closing the gap should be
reported, not treated as an engine correctness regression.

## Prior art and the next research question

- [Simonis, Sudoku as a Constraint Problem (2005)](https://ai.dmi.unibas.ch/_files/teaching/fs21/ai/material/ai26-simonis-cp2005ws.pdf)
  already treats matching, row/box value-set interactions, and flow propagation.
- [Bessière et al., Propagating Conjunctions of AllDifferent Constraints (2010)](https://arxiv.org/abs/1004.2626)
  studies simultaneous matchings and generalized Hall conditions.
- [Monfroy et al., Overlapping Alldiff constraints (2011)](https://academicjournals.org/journal/IJPS/article-abstract/7B3498125394)
  addresses overlap propagation and Sudoku encodings.
- [De Loera et al., Hilbert's Nullstellensatz and an Algorithm for Proving Combinatorial Infeasibility (2008)](https://arxiv.org/abs/0801.3788)
  establishes related low-degree polynomial proof methods.
- [Peyrl and Parrilo, Computing sum of squares decompositions with rational coefficients (2008)](https://www.mit.edu/~parrilo/pubs/files/PeyrlParrilo-ComputingSumOfSquaresDecompositionsWithRationalCoefficients.pdf)
  separates approximate numerical results from exact rational certificates.
- [Substitutability and interchangeability taxonomy (2010)](https://arxiv.org/abs/1010.4609)
  supplies generic theory behind safe local exchanges. A classical exchange can
  violate a Killer cage or Thermo, so every active intersecting constraint must
  be checked. Exchange-based eliminations additionally require justified uniqueness.
- [Efficiently Explaining CSPs with Unsatisfiable Subset Optimization](https://arxiv.org/abs/2303.11712)
  already studies minimal-cost explanations, including Sudoku.
- [Berthier's author-maintained examples](https://github.com/denis-berthier/CSP-Rules-Examples)
  include automated tridagon proofs and hard puzzle collections for comparison.

The next task is to minimize the three-house fixture, compare it against those
broader methods, and mine a recurring family with a short general proof. A second
track can seek exact quadratic/quartic certificates on mixed-list or four-digit
supports. Publish precise degree/support/rule budgets and negative searches.
Only a surviving, independently reviewed family should be described as a
potential mathematical discovery.
