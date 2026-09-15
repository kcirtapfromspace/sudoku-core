# Solver Engine Formal Specification

Derived from first principles of constraint satisfaction over the four
spaces of a standard 9x9 Sudoku grid. Every elimination theorem is proved
from the constraint axioms; named techniques are exhibited as parameter
instantiations of the technique engines. Arithmetic Counting additionally uses
weighted exactly-one equations and checked arithmetic certificates.

---

## 0. Foundational Axioms

### 0.1 The Four Spaces

**Cell Space** *C*.
81 positions arranged in a 9x9 matrix. Each cell *c* ∈ *C* is identified
by its linear index 0..80 (row-major). A cell is either *given* (fixed
value), *placed* (solver-assigned value), or *empty* (no value yet).

**Candidate Space** *X*.
The set of live (cell, digit) pairs: *X* ⊆ *C* × {1..9}. Initially
|*X*| ≤ 729. The element *(c, d)* ∈ *X* means "digit *d* is still
possible in cell *c*". We write *X(c)* for the candidate set of cell *c*
and *X_d(S)* for the set of cells in sector *S* that have candidate *d*.

**Sector Space** *S*.
27 units partitioning the grid: rows *r_0..r_8*, columns *c_0..c_8*,
boxes *b_0..b_8*. Index convention: 0..8 = rows, 9..17 = columns,
18..26 = boxes. Every cell belongs to exactly 3 sectors. Two cells *see*
each other iff they share at least one sector.

**Link Space** *L*.
Binary relations over nodes in *X*:

- **Strong link** *(p, q)*: exactly one of *p, q* is true in any
  solution. Arises when a digit has exactly 2 positions in a sector
  (conjugate pair) or a cell has exactly 2 candidates (bivalue cell).
  Formally: *p ∨ q* is a tautology and *¬(p ∧ q)* is a tautology.

- **Weak inference** *(p, q)*: at most one of *p, q* is true. Arises
  from the Sudoku constraint: two candidates of the same digit in the
  same sector, or two different candidates of the same cell.
  Formally: *¬(p ∧ q)* is a tautology.

  *Terminology note*: The community convention (StrmCkr, Sudopedia)
  prefers "weak inference" over "weak link" because the relationship
  is a logical NAND constraint, not a structural connection. A
  *strong link* is structural (XOR — the two candidates are bound
  by an exactly-one-true relationship), whereas a *weak inference* is
  derived from the constraint that at-most-one-true holds. In the
  implementation, the `LinkType::WeakInference` variant represents weak
  inferences. The terms are used interchangeably in older literature.

Every strong link is also a weak inference, but not conversely.

### 0.2 Constraint Axioms

For every solution *σ*: *C* → {1..9}:

1. **Uniqueness axiom**: ∀ sector *S*, ∀ digit *d* ∈ {1..9}:
   exactly one cell in *S* has value *d*.
2. **Completeness axiom**: ∀ cell *c*: *σ(c)* ∈ *X(c)*.
3. **Consistency axiom**: ∀ cells *c₁, c₂* sharing a sector:
   *σ(c₁) ≠ σ(c₂)*.

### 0.3 Inference Rules

Given the axioms, we may perform two operations:

- **Elimination**: Remove *(c, d)* from *X* when we can prove
  *σ(c) ≠ d* in every solution.
- **Placement**: Set *σ(c) = d* when *X(c) = {d}* (naked single) or
  when *d* has exactly one remaining position in some sector (hidden
  single).

---

## 1. Fish Engine: Sector-Rank Deficiency

### 1.1 Mathematical Objects

**Definition 1.1 (Fish Configuration).**
Fix a digit *d*. Let *B = {B₁, ..., Bₙ}* be a set of *n* sectors
(the *base set*) and *K = {K₁, ..., Kₙ}* be a set of *n* sectors
(the *cover set*) such that *B ∩ K = ∅*. The candidate sets within
each collection must be pairwise disjoint:
*X_d(Bᵢ) ∩ X_d(Bⱼ) = ∅* and *X_d(Kᵢ) ∩ X_d(Kⱼ) = ∅* for *i ≠ j*.
Distinct sector names alone do not establish independence: a row and
box can share a candidate cell. Define:

- *Base cells*:  *β = ⋃ᵢ X_d(Bᵢ)*  (cells with candidate *d* in
  base sectors)
- *Cover cells*: *κ = ⋃ⱼ X_d(Kⱼ)*  (cells with candidate *d* in
  cover sectors)
- *Fins*:        *φ = β \ κ*         (base cells not covered)
- *Eliminations*: *ε = κ \ β*        (cover cells not in base)

**Definition 1.2 (Sector Constraint).**
The classification of which sector types may appear in *B* and *K*:

| Constraint | Base types    | Cover types      | Name             |
|------------|--------------|------------------|------------------|
| Basic      | {row} or {col} | {col} or {row} | Standard fish    |
| Franken    | {row, col}   | {row, col, box}  | Franken fish     |
| Mutant     | {row, col, box} | {row, col, box} | Mutant fish   |

### 1.2 Core Theorem

**Theorem 1.1 (Basic Fish Elimination).**
If *φ = ∅* (no fins), then for every cell *c ∈ ε*, we may eliminate
*(c, d)* from *X*.

*Proof sketch.*
By the Uniqueness axiom, digit *d* must appear exactly once in each base
sector *Bᵢ*. Pairwise disjoint base candidate sets make these exactly
*n* distinct placements among *β*. Since *φ = ∅*, all of *β ⊆ κ*.
Each placement belongs to exactly one cover candidate set, and each
cover allows at most one placement. Thus all *n* covers are occupied
by the placements in *β*.
Therefore no cell in *ε = κ \ β* can contain *d*. ∎

**Theorem 1.2 (Finned Fish Elimination).**
If *φ ≠ ∅* and all fin cells lie in a single box *F*, then for every
cell *c ∈ ε* such that *c* shares box *F*, we may eliminate *(c, d)*.

*Proof sketch.*
Case analysis: either (a) no fin cell holds *d* in the solution, reducing
to Theorem 1.1, or (b) some fin cell in *F* holds *d*. In case (b), box
*F* already has *d*, so any cell in *ε ∩ F* cannot hold *d*. In both
cases, *d* is eliminated from *ε ∩ cells(F)*. ∎

### 1.3 Named Techniques as Special Cases

Every named fish technique is an instantiation of (digit *d*, size *n*,
constraint, fin status):

| n | Fins | Constraint | Technique            | SE   |
|---|------|-----------|----------------------|------|
| 2 | No   | Basic     | X-Wing               | 3.2  |
| 2 | Yes  | Basic     | Finned X-Wing        | 3.4  |
| 3 | No   | Basic     | Swordfish            | 3.8  |
| 3 | Yes  | Basic     | Finned Swordfish     | 4.0  |
| 4 | No   | Basic     | Jellyfish            | 5.2  |
| 4 | Yes  | Basic     | Finned Jellyfish     | 5.4  |
| * | *    | Franken   | Franken Fish         | 5.5  |
| * | *    | Mutant    | Mutant Fish          | 6.5  |

**Siamese Fish** are two overlapping finned fish configurations sharing
a fin box. Both fish independently restrict eliminations to the fin box;
the implementation returns eliminations common to both patterns. An
intersection cannot add eliminations beyond either individual pattern.

### 1.4 Degenerate Case: Pointing Pairs and Box/Line Reduction

**Theorem 1.3 (Intersection as Fish of Size 1).**
Pointing Pair and Box/Line Reduction are degenerate fish with *n = 1*.

*Proof.*
A **Pointing Pair** for digit *d*: the *d*-candidates in box *B* are
confined to row *R*. Set base *B = {box B}*, cover *K = {row R}*. Then:
- *β* = positions of *d* in box *B*
- *κ* = positions of *d* in row *R*
- *φ = β \ κ = ∅* (all *d*-cells in the box lie on row *R*)
- *ε = κ \ β* = cells in row *R* outside box *B*

By Theorem 1.1, eliminate *d* from *ε*. This is exactly the Pointing
Pair rule.

A **Box/Line Reduction** for digit *d*: the *d*-candidates in row *R*
are confined to box *B*. Set base *B = {row R}*, cover *K = {box B}*.
- *β* = positions of *d* in row *R*
- *κ* = positions of *d* in box *B*
- *φ = β \ κ = ∅* (all *d*-cells in the row lie in box *B*)
- *ε = κ \ β* = cells in box *B* outside row *R*

By Theorem 1.1, eliminate *d* from *ε*. This is exactly the Box/Line
Reduction rule. ∎

Therefore, intersections are a strict subset of the fish framework.
They need not be implemented separately.

### 1.5 Sector-Rank Interpretation

The implementation uses the sufficient counting conditions of Definition
1.1, checked by `find_independent_combination`. A candidate may occur in at
most one base group and at most one cover group. This makes the distinct
placement count in Theorem 1.1 valid for mixed row/column/box sectors too.

Full matrix rank alone does not justify that count: overlapping sectors
can share a placement. Configurations needing additional overlap or rank
reasoning are rejected by this implementation. This restriction avoids
false eliminations while retaining independent Franken and Mutant fish.

---

## 2. ALS Engine: Subset Degree-of-Freedom

### 2.1 Mathematical Objects

**Definition 2.1 (Almost Locked Set).**
A set *A* of *n* empty cells within a single sector, such that
|⋃_{c ∈ A} X(c)| = *n* + 1. That is, *n* cells collectively have
exactly *n* + 1 distinct candidates. We write *cands(A)* for this
union.

The "almost" refers to having one excess candidate beyond what a fully
locked set (naked subset) would have. This single degree of freedom is
what enables linking.

**Definition 2.2 (Restricted Common Candidate, RCC).**
Given two non-overlapping ALS *A* and *B*, a digit *x* is an RCC of
(*A, B*) if:
- *x* ∈ *cands(A)* ∩ *cands(B)*, and
- every *x*-cell in *A* sees every *x*-cell in *B*.

When *A* and *B* are linked by RCC *x*, they cannot both contain a
placed *x*. Thus at least one of them omits *x* in every solution;
both may omit it. An ALS of *n* cells that omits one of its *n* + 1
candidates must use all *n* remaining digits, since its cells share
a sector and therefore take distinct values.

### 2.2 Core Theorem

**Theorem 2.1 (ALS-XZ Elimination).**
Let *A* (size *n*) and *B* (size *m*) be non-overlapping ALS linked
by RCC *x*. Let *z* ∈ *cands(A)* ∩ *cands(B)*, *z ≠ x*. Then for
any cell *c* ∉ *A* ∪ *B* that sees every *z*-cell in *A* and every
*z*-cell in *B*, we may eliminate *(c, z)*.

*Proof.*
The RCC prevents *A* and *B* from both containing *x*. At least one
ALS therefore omits *x* and must use every remaining digit in its
candidate union, including *z* because *z ≠ x*. Cell *c* sees every
possible *z*-placement in either ALS, so it cannot hold *z*. ∎

This argument holds in every solution and does not assume the puzzle
has a unique solution.

### 2.3 Named Techniques in the Implementation

| A size | B size | Links | Technique     | SE   |
|--------|--------|-------|---------------|------|
| 1      | 1      | 1 RCC | XY-Wing       | 4.2  |
| 1      | 2      | 1 RCC | XYZ-Wing      | 4.4  |
| 2      | 1      | 1 RCC | XYZ-Wing      | 4.4  |
| any    | any    | 1 RCC, total=4 | WXYZ-Wing | 4.6 |
| any    | any    | 1 RCC, other sizes | ALS-XZ | 5.5 |

These are the size-based labels returned by `classify_als_pair` for
the shared ALS-XZ detector. Every returned elimination must satisfy
Theorem 2.1. The labels do not establish equivalence to every pattern
conventionally called a wing; in particular, the `(1, 1)` label here
describes two cells, not a three-cell pivot-and-wings construction.

#### 2.3.1 ALS Chains (Generalized)

**Theorem 2.2 (ALS Chain Elimination).**
Given a chain of distinct ALS *A₁ - A₂ - ... - Aₖ*, *k ≥ 2*, let
*xᵢ* be an RCC of *(Aᵢ, Aᵢ₊₁)*. Consecutive ALS are disjoint as
required by the RCC definition, and the endpoints are disjoint.
Require consecutive link digits to differ: *xᵢ ≠ xᵢ₊₁*.
Let *z* ∈ *cands(A₁)* ∩ *cands(Aₖ)* with
*z ≠ x₁* and *z ≠ x_{k-1}*.

For any cell *c* outside the two endpoints that sees all *z*-cells in
*A₁* and all *z*-cells in *Aₖ*, we may eliminate *(c, z)*. Membership
in an intermediate ALS does not invalidate the implication below.

*Proof.*
For *k = 2*, apply Theorem 2.1. For *k ≥ 3*, suppose *c = z*.
Both endpoints must omit *z*. Hence *A₁* must use
*x₁*, which makes *A₂* omit *x₁*. Since consecutive link digits
differ, *A₂* must then use *x₂*. Repeating this implication makes
*Aₖ* omit *x_{k-1}*. But *Aₖ* already omits the distinct digit *z*:
its *n* cells would have at most *n* − 1 available digits, contrary
to sector consistency. Therefore *c ≠ z*. ∎

For the three-ALS detector, *z* must differ from **both** link digits
*x* and *y*. Allowing *z = y* invalidates the last counting step and
can remove a candidate from a valid completion. For a four-ALS chain,
*z* may equal the middle link digit: only the two endpoint link digits
must differ from *z*. These deductions require no unique-solution
assumption.

| Chain length | Technique     | SE   |
|-------------|---------------|------|
| 2           | ALS-XZ        | 5.5  |
| 3           | ALS-XY-Wing   | 7.0  |
| 4           | ALS Chain     | 7.5  |

The implementation currently searches three- and four-ALS chains;
the theorem itself applies to longer chains meeting the same conditions.

### 2.4 Sue de Coq as ALS Decomposition

**Theorem 2.3 (Sue de Coq).**
Let *I* be the intersection of a box and a line (2-3 empty cells).
Let *cands(I)* be the union of candidates in *I*. If there exist
ALS *A* in (rest of box) and ALS *B* in (rest of line) such that:
- *cands(A) ∪ cands(B) = cands(I)*
- *cands(A) ∩ cands(B) = ∅*

Then:
- Eliminate *cands(A)* from rest-of-box cells not in *A*.
- Eliminate *cands(B)* from rest-of-line cells not in *B*.

*Proof.*
Write *U = cands(A)* and *V = cands(B)*. In any solution, the
*|U| − 1* cells of *A* use that many distinct digits of *U*. Since
every intersection cell sees all of *A*, at most one intersection
cell can use a digit of *U*. Similarly, at most one can use *V*.
The intersection has at least two cells and uses only *U ∪ V*, so a
solution requires exactly two intersection cells, one using the digit
of *U* omitted by *A* and the other the digit of *V* omitted by *B*.
Thus every digit of *U* occurs in *A ∪ I* within the box, and every
digit of *V* occurs in *B ∪ I* within the line. The stated eliminations
follow from sector consistency. A three-cell intersection with these
particular two disjoint ALS contributions admits no solution. ∎

The proof certificate records three descriptors: the intersection,
the box ALS, and the line ALS. Each descriptor's `sector` must contain
its cells: the intersection uses the line sector, the box ALS uses the
actual box sector (18..26), and the line ALS uses the actual row or
column sector (0..17). The intersection is a described region, not
necessarily an ALS under Definition 2.1. Here `rcc_values` stores the
two disjoint candidate contributions and `z_value` is absent; these
fields do not encode an ordinary RCC chain.

### 2.5 Death Blossom as ALS Star Graph

**Theorem 2.4 (Death Blossom).**
Let *s* be a stem cell with candidates *{d₁, ..., dₖ}*. For each
*dᵢ*, let *Pᵢ* be an ALS (a petal) such that:
- *s ∉ Pᵢ* and the petals are pairwise non-overlapping.
- *dᵢ* ∈ *cands(Pᵢ)* and **every** cell of *Pᵢ* containing candidate
  *dᵢ* sees *s*.

If *z* ∈ ⋂ᵢ *cands(Pᵢ)*, *z ∉ {d₁,...,dₖ}*, and cell *c* outside
the stem and petals sees all *z*-cells in every petal, then eliminate
*(c, z)*.

*Proof sketch.*
In any solution, *s* takes some value *dⱼ*. Then petal *Pⱼ* cannot
use *dⱼ* (since every *dⱼ*-cell of *Pⱼ* sees *s*), so *Pⱼ* becomes
fully locked on its remaining *n* candidates. The value *z* is
confined to *z*-cells of *Pⱼ*. Since *c* sees all of them, *c ≠ z*.
This holds for every possible *dⱼ*, so the elimination is valid. ∎

Checking only that some linking-digit occurrence sees the stem is
insufficient: an unseen occurrence may retain *dⱼ*, leaving the petal
unlocked. The detector checks all occurrences before selecting petals.
Its current search tries only the first available petal per stem digit;
it may miss another usable combination, but every returned combination
must meet the above conditions.

### 2.6 Distributed Disjoint Subset (DDS)

DDS is used here as an organizational description for deductions
involving a stem region and surrounding candidate sets. The implemented
Sue de Coq and Death Blossom searches have different proof obligations:
Sue de Coq uses disjoint candidate partitions and box/line capacity;
Death Blossom uses a restricted link for every possible stem value.
In the latter, the chosen stem value locks its corresponding petal;
it need not lock every other petal.

Candidate-union containment alone, or an unspecified partition of
"degrees of freedom", does not imply an elimination. There is no
separate generalized DDS detector or elimination theorem established
here. A finding must satisfy the concrete hypotheses of Theorem 2.3
or Theorem 2.4, including their visibility and cardinality requirements.

### 2.7 Aligned Pair/Triplet Exclusion (Legacy)

**Theorem 2.6 (Enumerated Assignment Elimination).**
Let *P* contain two or three mutually visible cells. Enumerate all
assignments choosing a candidate for each cell, with distinct digits
in every pair of cells. If this collection is nonempty and every
assignment uses digit *z*, eliminate *z* from any cell outside *P*
that sees every cell in *P*.

*Proof.*
Every solution restricts to one of the enumerated assignments, because
it respects the candidate sets and pairwise visibility. Hence some cell
in *P* contains *z* in every solution. A common peer cannot also contain
*z*. ∎

This is the direct proof of the legacy APE/ATE detectors in this engine.
It does not require that *P* or each individual cell be an ALS: the pair
search admits cells with 2–5 candidates, and the triplet search admits
2–4. Their `Als` certificate descriptors record these candidate sets;
the certificate variant does not itself prove Definition 2.1 holds.
The public functions remain deprecated for compatibility. No complete
subsumption by the bounded ALS-chain search is established here.

### 2.8 Almost Locked Candidates (ALC)

An Almost Locked Candidates (ALC) pattern is the intersection dual of an
ALS: *n* candidates each appearing in at most *n* + 1 cells within a
sector's intersection region. Where ALS reasons about cells having a
surplus candidate, ALC reasons about digits having a surplus position.

The implementation does not enumerate ALC as a separate technique.
Some intersection deductions have the direct sector-containment proof
given for Pointing Pairs and Box/Line Reductions in Section 1.4. A digit
having a small number of possible positions alone does not establish
such an elimination. This specification does not prove that the bounded
ALS and fish searches cover every more general ALC pattern.

---

## 3. AIC Engine: Link-Reachability Arguments

### 3.1 Mathematical Objects

**Definition 3.1 (Inference Graph).**
The AIC inference graph *G = (V, E_s, E_w)* where:
- *V* = *X* (the candidate space — each node is a (cell, digit) pair)
- *E_s* ⊆ *V × V*: strong links (exactly-one-true)
- *E_w* ⊆ *V × V*: weak inferences (at-most-one-true), *E_s ⊆ E_w*

**Definition 3.2 (Link Sources).**
Strong links arise from:
1. **Conjugate pairs**: digit *d* has exactly 2 positions in sector *S*
   → strong link between *(c₁, d)* and *(c₂, d)*.
2. **Bivalue cells**: cell *c* has exactly 2 candidates *{a, b}*
   → strong link between *(c, a)* and *(c, b)*.

Weak inferences arise from:
1. **Same digit, same sector**: *(c₁, d)* — *(c₂, d)* when *c₁, c₂*
   share a sector (at most one can be *d*).
2. **Same cell, different digit**: *(c, d₁)* — *(c, d₂)* (a cell holds
   at most one value).

**Definition 3.3 (Alternating Inference Chain).**
An AIC is a path *p₁ — p₂ — p₃ — ... — pₖ* in *G* where edges
alternate between strong and weak:

    p₁ ==[strong]== p₂ --[weak]-- p₃ ==[strong]== p₄ -- ... -- pₖ

The search starts with a strong link and traverses
strong → weak → strong → weak → ... . For the endpoint eliminations
below, the final edge must also be strong, so the chain has an even
number of nodes. The implementation checks these endpoints after
arriving by a strong link and requires at least four nodes.

Under the assumption *p₁* = false, propagation assigns OFF to
odd-indexed nodes and ON to even-indexed nodes (using the one-based
indices above). These are conditional truth values, not descriptions
of the edge type. The emitted general AIC/X-Chain certificates follow
this OFF → ON → OFF → ON order.

### 3.2 Core Theorem

**Theorem 3.1 (AIC Elimination, Type 1: Shared Digit, Different Cells).**
If an alternating chain starts and ends with strong links, with
endpoints *(c₁, d)* and *(c₂, d)* where *c₁ ≠ c₂*, then:

For any cell *c* ∉ {c₁, c₂} that sees both *c₁* and *c₂*, eliminate
*(c, d)*.

*Proof.*
Consider the first endpoint *p₁ = (c₁, d)*.

**Case A**: *p₁* is true. Since *c* sees *c₁*, *c* cannot hold *d*.

**Case B**: *p₁* is false. The alternating implications are:

- ¬*p₁* → *p₂* (strong link: one must be true)
- *p₂* → ¬*p₃* (weak inference: at most one is true)
- ¬*p₃* → *p₄* (strong link)
- ...continuing alternation...
- ¬*pₖ₋₁* → *pₖ* (the final strong link)

Thus *pₖ = (c₂, d)* is true, and *c* cannot hold *d* because it sees
*c₂*. In both cases at least one endpoint is true, so the elimination
holds. No additional loop-closing assumption is needed. ∎

A final **weak** inference would instead force the last endpoint OFF
under Case B. It does not establish that either endpoint is true and
cannot justify this elimination. For example, an S-W-S-W path admits
the assignment OFF, ON, OFF, ON, OFF: every XOR/NAND edge is satisfied
while both endpoints are false.

**Theorem 3.2 (AIC Elimination, Type 2: Same Cell, Different Digits).**
If an alternating chain starts and ends with strong links, with
endpoints *(c, d₁)* and *(c, d₂)* where *d₁ ≠ d₂*, then eliminate all
candidates of *c* except *d₁* and *d₂*.

*Proof.*
The same implication argument establishes that at least one endpoint
is true. Since a cell holds exactly one value, *σ(c) ∈ {d₁, d₂}*.
All other candidates of *c* are eliminated. ∎

### 3.3 Named Techniques as Special Cases

**Empty Rectangle** (grouped single-digit inference):
For digit *d*, all candidates in a box must lie in the union of a row
arm and a column arm, with a candidate on each arm outside their
intersection. The intersection itself may contain a candidate.
A conjugate pair in a column outside the box connects the row arm to
another row outside the box. Eliminate *d* at the intersection of that
other row and the box's column arm, if it is a candidate. The transposed
pattern exchanges rows and columns. This uses ordinary Sudoku
constraints and does not assume a unique puzzle solution. SE: 4.6.

For example, let the box-1 candidates for *d* be r1c2, r2c1 and r3c1,
and let column 5 have exactly two candidates for *d*: r1c5 and r5c5.
Then r5c1 cannot hold *d*. If it did, it would remove *d* from r2c1,
r3c1 and r5c5. Box 1 would force r1c2 = *d*, while column 5 would force
r1c5 = *d*, contradicting row 1.

The box arm is a **group** of candidates. It does not create an
exactly-one-true link between arbitrary individual cells in that arm.
The implementation records the complete grouped evidence using a
`Fish` certificate: the box and external conjugate column are the two
base sectors, the two crossing rows are cover sectors, and the box's
column-arm candidates outside the row cover are fins. In the example,
these are bases {box 1, column 5}, covers {row 1, row 5}, and fins
{r2c1, r3c1}; the eliminated cell sees every fin. The deduction remains
classified as `EmptyRectangle` in the AIC engine.

**X-Chain** (single-digit AIC):
An AIC where every node concerns the same digit *d*. All links are
conjugate-pair strong links and same-sector weak inferences. The single-
value restriction makes the chain search cheaper. SE: 4.5.

**W-Wing**:
Two bivalue cells *{a, b}* at positions *c₁, c₂* connected by a
conjugate pair (strong link) on digit *a* in some sector. The chain:

    (c₁, b) ==[strong, bivalue]== (c₁, a) --[weak, sees link]--
    (ℓ₁, a) ==[strong, conjugate]== (ℓ₂, a) --[weak, sees c₂]--
    (c₂, a) ==[strong, bivalue]== (c₂, b)

This 6-node AIC with endpoints *(c₁, b)* and *(c₂, b)* yields:
eliminate *b* from any cell seeing both *c₁* and *c₂*. SE: 4.4.

**Named Wings** (3-strong-link AIC taxonomy, StrmCkr):
Each strong link in a size-3 chain is classified as V (bivalue: same
cell, different digits) or L (bilocal: same digit, conjugate pair).
The 3-character signature identifies the wing type:

| Signature | Name    | Structure                             |
|-----------|---------|---------------------------------------|
| VVV       | XY-Wing | 3 bivalue cells                       |
| VLV       | W-Wing  | bivalue → conjugate pair → bivalue    |
| LVL       | S-Wing  | conjugate → bivalue → conjugate       |
| VLL       | M-Wing  | bivalue → conjugate → conjugate       |
| LLV       | M-Wing  | conjugate → conjugate → bivalue       |
| LLL       | L-Wing  | 3 conjugate pairs (single-digit)      |
| VVL       | H-Wing  | bivalue → bivalue → conjugate         |
| LVV       | H-Wing  | conjugate → bivalue → bivalue         |

The implementation classifies chains at output time and includes the
signature in the explanation variant string (e.g., "S-Wing (LVL)").
XY-Wing and W-Wing are found by dedicated functions before the
general AIC search, so the AIC engine primarily discovers S-Wing,
M-Wing, L-Wing, and H-Wing patterns.

**General AIC** (multi-digit):
Unrestricted alternating chains mixing conjugate-pair and bivalue
strong links, with weak inferences from sector/cell constraints. SE: 6.0.

### 3.4 Medusa as AIC Coloring (Legacy)

**Theorem 3.3 (3D Medusa is AIC Strong-Link Coloring).**
3D Medusa is the connected-component coloring of the strong-link
subgraph of the AIC inference graph.

*Proof.*
Medusa assigns two colors (0 and 1) to nodes by BFS traversal of
strong links only. This partitions each connected component of the
strong-link graph into two color classes. The coloring satisfies:
if node *p* has color *c*, then every strong neighbor has color
*1 − c*.

Since every strong link represents an exactly-one-true constraint:
exactly one color class is the "true" class (all its nodes hold in
the solution) and the other is the "false" class.

**Contradiction rules** (identify the false color class):
- Rule 1: Two same-color nodes share a digit *d* and a sector.
  Then that color class would require *d* twice in one sector. ✗
- Rule 2: Two same-color nodes are in the same cell.
  Then that color class would place two values in one cell. ✗

If a contradiction is found for color *c*, all nodes of color *c*
are false → eliminate those candidates.

**Rule 5** (uncolored elimination): If an uncolored candidate
*(cell, d)* sees *d*-nodes of both colors, eliminate it (one color
is true, so the *d*-node of the true color blocks *cell*).

This is precisely the AIC framework restricted to strong-link-only
BFS. Every Medusa elimination can be expressed as an AIC, but the
coloring approach finds them more efficiently by processing entire
connected components at once. ∎

*Legacy status*: 3D Medusa is considered retired by the community
(StrmCkr, Sudopedia). It is simply strong-link-only AIC coloring and
adds no eliminations beyond what the general AIC engine finds. The
implementation retains it for SE rating compatibility (SE 5.0) but it
is not an independent technique — it is a restricted AIC search
strategy.

### 3.5 Forcing Chains as Multi-Source AIC

**Definition 3.4 (Forcing Chain).**
A forcing chain explores all candidates of a source (cell or sector
position) via independent propagation. If all branches agree on an
outcome, that outcome is valid.

| Source        | Branches                  | Technique              | SE   |
|--------------|---------------------------|------------------------|------|
| Cell *c*     | Each candidate of *c*     | Cell Forcing Chain     | 8.3  |
| Sector/digit | Each position of *d* in *S* | Region Forcing Chain | 8.5  |
| Cell (full)  | With full technique propagation | Dynamic FC      | 9.3  |
| Single cand. | One branch contradicts    | Nishio FC              | 7.5  |

**Theorem 3.4 (Forcing Chain Correctness).**
Let *S = {s₁, ..., sₖ}* be an exhaustive set of hypotheses (candidates
of a cell, or positions of a digit in a sector). For each *sᵢ*,
let *Gᵢ* be the grid state after propagating the assumption *sᵢ = true*.

- **Common placement**: If ∀ *i*: *Gᵢ* places value *v* at cell *t*,
  then *σ(t) = v*.
- **Common elimination**: If ∀ *i*: *Gᵢ* eliminates *(t, v)*, then
  eliminate *(t, v)*.
- **Nishio**: If propagating *sᵢ* leads to contradiction, then
  *sᵢ* is false; eliminate it.

*Proof.*
By exhaustion: one of *{s₁, ..., sₖ}* must be true (completeness
axiom for cells; uniqueness axiom for sectors). Since every possible
truth leads to the same conclusion, that conclusion holds. ∎

**Kraken Fish**: A finned fish where the fin's effect is verified by
forcing-chain propagation. If propagating each fin cell with digit *d*
still eliminates *d* from the target, the elimination is valid — even
when the standard finned-fish box restriction would not apply.

### 3.6 Dynamic Forcing Chains

Dynamic Forcing Chains are the most powerful non-backtracking technique.
They differ from standard forcing chains in the propagation function:
instead of propagating only naked/hidden singles, they apply the full
technique repertoire (including AIC, ALS, and fish) within each branch.

This makes them strictly more powerful: they can discover implications
that simple propagation misses. The trade-off is computational cost.
SE: 9.3.

---

## 4. Completeness and Orthogonality

### 4.1 Three Orthogonal Arguments

The three engines cover three fundamentally different proof strategies:

| Engine | Core argument          | Operates primarily in |
|--------|------------------------|----------------------|
| Fish   | Sector-rank deficiency | Sector Space *S*     |
| ALS    | Subset degree-of-freedom | Candidate Space *X* |
| AIC    | Link reachability      | Link Space *L*       |

**Fish** reasons about digits in aggregate across sectors: "these *n*
base sectors consume exactly *n* cover sector slots for digit *d*."

**ALS** reasons about cell groups and their candidate surplus: "*n*
cells have *n* + 1 candidates, so linking two such groups forces
value confinement."

**AIC** reasons about individual (cell, digit) nodes and the logical
implications of truth propagation through strong links and weak inferences.

### 4.2 Subsumption Relationships

```
Pointing Pair  ──┐
Box/Line Reduc.──┼──► Fish Engine (n=1 degenerate case)
                 │
X-Wing ──────────┤
Swordfish ───────┼──► Fish Engine (n=2,3,4; Basic constraint)
Jellyfish ───────┤
                 │
Franken Fish ────┼──► Fish Engine (Franken constraint)
Siamese Fish ────┤──► Fish Engine (overlapping finned pair)
Mutant Fish ─────┘──► Fish Engine (Mutant constraint)

XY-Wing ─────────┐
XYZ-Wing ────────┼──► ALS Engine (size-1,2 ALS pairs)
WXYZ-Wing ───────┤
ALS-XZ ──────────┤
ALS-XY-Wing ─────┼──► ALS Engine (3-ALS chain)
ALS Chain ───────┤──► ALS Engine (4+ ALS chain)
Sue de Coq ──────┤──► ALS Engine (box/line decomposition)
Death Blossom ───┤──► ALS Engine (star graph topology)
APE/ATE ─────────┘──► ALS Engine (mutually-visible cell group)

Empty Rectangle ─┐──► AIC Engine (single-digit ERI chain)
W-Wing ──────────┤
X-Chain ─────────┤
3D Medusa ───────┼──► AIC Engine (strong-link coloring, legacy)
AIC ─────────────┤
Nishio FC ───────┤──► AIC Engine (single-branch contradiction)
Kraken Fish ─────┤──► AIC Engine (fish + forcing verification)
Cell FC ─────────┤──► AIC Engine (cell-source forcing)
Region FC ───────┤──► AIC Engine (sector-source forcing)
Dynamic FC ──────┘──► AIC Engine (full-propagation forcing)
```

### 4.3 Technique Coverage by Engine

The solver dispatches 46 technique variants, including backtracking. Their engine ownership:

| Engine        | Count | Technique variants |
|---------------|:-----:|---|
| Basic (direct) | 8 | NakedSingle, HiddenSingle, NakedPair/Triple/Quad, HiddenPair/Triple/Quad |
| Fish          | 11 | PointingPair, BoxLineReduction, X-Wing, Finned X-Wing, Swordfish, Finned Swordfish, Jellyfish, Finned Jellyfish, Franken Fish, Siamese Fish, Mutant Fish |
| ALS           | 10 | XY-Wing, XYZ-Wing, WXYZ-Wing, ALS-XZ, ALS-XY-Wing, ALS Chain, Sue de Coq, Death Blossom, Aligned Pair Exclusion, Aligned Triplet Exclusion |
| AIC           | 10 | Empty Rectangle, W-Wing, X-Chain, 3D Medusa, AIC, Nishio FC, Kraken Fish, Cell FC, Region FC, Dynamic FC |
| Uniqueness    | 5 | Avoidable Rectangle, Unique Rectangle, Hidden Rectangle, Extended UR, BUG |
| Arithmetic    | 1 | Arithmetic Counting |
| Backtracking  | 1 | Backtracking |

Note: Basic techniques (singles, subsets) are direct applications of the
axioms and do not require an engine. Uniqueness techniques rely on the
assumption that the puzzle has a unique solution (an additional axiom
beyond the standard constraint set).

### 4.4 Uniqueness Proof Requirements

For classic Sudoku, uniqueness deductions require a uniquely solvable puzzle
and a candidate state that retains its solution. The engine receives that
assumption from its caller; these detectors do not establish uniqueness by
counting solutions themselves.
The following conditions are checked in addition to that assumption.
The swap arguments below preserve rows, columns and boxes. Additional variant
constraints require a separate justification that the swap preserves them too.

**Hidden Rectangle.** Let four empty corners occupy two rows, two columns and
exactly two boxes, with candidates *a* and *b* present at every corner. For a
target corner *T*, the diagonally opposite corner must have exactly `{a, b}`.
Candidate *b* must form a conjugate pair between *T* and the adjacent corner in
**both** T's row and T's column. Then *a* can be removed from *T*.

Proof: if *T = a*, those two strong links force *b* into both adjacent corners.
The opposite bivalue corner must therefore contain *a*. Swapping *a* and *b*
at all four corners preserves every row, column and box and changes no given,
creating a second solution. This contradicts uniqueness. A single row or column
strong link, or an opposite corner with additional candidates, does not prove
this result and must not trigger this inference.

**BUG remainder.** For proposed extra candidates *E*, removing *E* must leave
exactly two candidates in every empty cell. In every sector, every digit not
already placed must occur in exactly two candidate cells; a placed digit must
occur in none. Merely finding odd candidate counts or mostly bivalue cells
does not establish this pattern. `is_bug_remainder` validates the entire state.

Proof: any solution confined to such a remainder has a distinct partner obtained
by swapping each empty cell to its other candidate. Each sector-digit pair
occurs at two cells, so the swap moves the digit between those cells and keeps
the sector valid. Given cells remain fixed. The remainder can therefore have
zero or multiple solutions, but cannot contain the puzzle's unique solution
alone. At least one candidate in *E* must be true.

- **BUG+1:** if *E* contains one candidate, place it in its cell.
- **BUG+n:** the implemented common-unit elimination requires every extra to
  be the same digit *d*. If all extra cells share a row, column or box, remove
  *d* from other candidate cells in that unit. Such a target conflicts with
  every possible extra, and at least one extra must hold. Extras of another
  digit or outside the common unit invalidate this particular deduction.

The implementation bounds BUG analysis to at most six extras. Odd counts are
used to propose extras, followed by the complete remainder check; they are not
the proof. Cases the proposal or common-unit rule cannot establish return no
BUG finding so subsequent techniques can continue.

---

## 5. Hint Delivery

The main solver exposes two hint paths. Dedicated arithmetic APIs also expose
checked deductions over the caller's exact candidate state (Section 8).

### 5.1 `get_hint()` — Display Hints

`get_hint(grid)` returns the first applicable technique as a `Hint`.
It may return either a `SetValue` (placement) or `EliminateCandidates`
(elimination). **It does not verify the result against the backtracking
solution.** This method returns a hint and its available proof evidence without mutating the
input. Every deduction must still be logically sound. Arithmetic findings are
checked against their source equations before entering this path.

### 5.2 `get_next_placement()` — Verified Placement Hints

`get_next_placement(grid)` is the safe hint path used by all frontends
when *applying* a hint to the game state. It works as follows:

1. Solve the grid once via backtracking to obtain the verified solution.
2. Find the first applicable technique.
3. If it is an **elimination**: verify that no eliminated candidate is
   the solution value. If sound, apply it internally and repeat from
   step 2 (up to 500 iterations).
4. If it is a **placement**: verify that the placed value matches the
   solution. If sound, return it as the hint.
5. If any step produces an unsound result (technique bug), **stop
   chaining** and fall back to a backtracking hint (always correct).

This loop chains eliminations internally so the caller always receives
a `SetValue` hint. The solution verification at each step ensures that
known technique bugs (e.g., Avoidable Rectangle on given cells, W-Wing
self-links, X-Chain conjugate errors) never produce wrong placements.

### 5.3 ProofCertificate

Each `Hint` carries an optional `ProofCertificate` providing structured
metadata for visualization:

| Variant        | Fields                                        |
|---------------|-----------------------------------------------|
| `Basic`       | involved cells                                |
| `Fish`        | base sectors, cover sectors, fin cells, digit |
| `Als`         | ALS chain (cells + candidates per ALS)        |
| `Aic`         | chain of (cell, digit, polarity) nodes and link types |
| `Uniqueness`  | floor cells, roof cells                       |
| `Forcing`     | source cell, branches                         |
| `Backtracking`| (no fields)                                   |

Frontends use these to render proof-detail overlays (e.g., base/cover
sector highlighting for fish, on/off coloring for AIC chains).
General AIC/X-Chain certificates start with the first candidate OFF and
alternate conditional truth values along the implication described in
Section 3.2. Empty Rectangle uses the `Fish` representation to preserve
all candidates in its grouped box arm (Section 3.3); a certificate's
variant need not have the same name as the hint technique.

---

## 6. Soundness Guarantees

### 6.1 Engine Soundness

Each engine's eliminations are sound if the axioms hold:

- **Fish**: Theorems 1.1 and 1.2 require the sector uniqueness axiom
  and the pairwise candidate-set independence in Definition 1.1.
- **ALS**: Sound by Theorems 2.1, 2.2, 2.3, 2.4. Depends on the
  Uniqueness axiom and Completeness axiom.
- **AIC**: Sound by Theorems 3.1, 3.2, 3.4. Depends on the logical
  semantics of strong links and weak inferences, which are derived from all
  three axioms.
- **Uniqueness**: additionally requires a unique solution. Hidden Rectangle
  and BUG must satisfy the complete pattern conditions in Section 4.4.

### 6.2 Implementation Verification

Three test suites verify soundness empirically:

- `test_hint_soundness`: For uniquely solvable fixtures, every finding
  in a stored-candidate technique chain is checked against the solution.
  Public `get_hint()` behavior is tested separately: it deliberately
  recalculates candidates for each independent request. A placement
  *(c, v)* must match the solution, and an elimination *(c, v)* must
  not remove the solution value.
- `test_hint_soundness_all_tiers`: Extends coverage across all
  difficulty tiers.
- `test_next_placement_soundness`: Verifies that `get_next_placement()`
  produces correct placements for every puzzle in the battery.

These are runtime verifications, not formal proofs, but they provide
high confidence that the implementation correctly realizes the
theorems above. The `get_next_placement()` verification loop
(Section 5.2) provides an additional per-call safety net in
production.

### 6.3 Termination

The solver terminates because:
1. Each technique application either places a value (reducing
   empty cells by 1) or eliminates at least one candidate.
2. Candidates are never re-added.
3. The candidate space *X* is finite and monotonically decreasing.
4. If no technique applies, backtracking terminates in finite time
   (brute-force search over a finite space).
5. `get_next_placement()` bounds its chaining loop to 500 iterations.

---

## 7. Difficulty Classification

### 7.1 Eight-Tier System

Puzzles are classified into eight difficulty tiers. Three independent
mechanisms use these tiers for different purposes.

| Tier         | SE Range    | Technique Hint          |
|-------------|-------------|-------------------------|
| Beginner    | 1.5 – 2.0   | Hidden singles          |
| Easy        | 2.0 – 2.5   | Naked singles           |
| Medium      | 2.5 – 3.4   | Pairs & triples         |
| Intermediate| 3.4 – 3.8   | Hidden triples          |
| Hard        | 3.8 – 4.5   | Box/line reduction      |
| Expert      | 4.5 – 5.5   | Fish & rectangles       |
| Master      | 5.5 – 7.0   | Wings & chains          |
| Extreme     | 7.0 – 11.0  | Advanced techniques     |

Master and Extreme are hidden tiers, unlocked by the player.

### 7.2 Three Difficulty Axes

**SE rating** (`rate_se`): The maximum SE rating among all techniques
used to solve the puzzle. This is a continuous numerical score that
uses the engine's technique rating table. Arithmetic Counting has an engine-local,
uncalibrated value of 8.5, not an official Sudoku Explainer rating.

**Technique-based classification** (`technique_to_difficulty`): Maps
the hardest technique used during solving to a difficulty tier. This
is a discrete classification based on the *kind* of reasoning required:

| Tier         | Techniques included |
|-------------|---|
| Beginner    | NakedSingle (≤35 empty cells) |
| Easy        | NakedSingle (>35 empty cells) |
| Medium      | HiddenSingle |
| Intermediate| NakedPair, HiddenPair, NakedTriple, HiddenTriple |
| Hard        | PointingPair, BoxLineReduction |
| Expert      | X-Wing (±fin), Swordfish (±fin), Jellyfish (±fin), NakedQuad, HiddenQuad, EmptyRectangle, AvoidableRectangle, UniqueRectangle, HiddenRectangle |
| Master      | XY/XYZ/WXYZ-Wing, W-Wing, X-Chain, 3D Medusa, SueDeCoq, AIC, FrankenFish, SiameseFish, ALS-XZ, ExtendedUR, BUG |
| Extreme     | ALS-XY-Wing, ALS Chain, MutantFish, APE, ATE, DeathBlossom, ArithmeticCounting, Nishio/Kraken/Cell/Region/Dynamic FC, Backtracking |

**Generation cap** (`max_technique`): Limits which techniques the
generator may require. This uses the `Technique` enum ordering (not SE
rating) to define an upper bound:

| Tier         | Max technique allowed     |
|-------------|---------------------------|
| Beginner    | NakedSingle               |
| Easy        | NakedSingle               |
| Medium      | HiddenSingle              |
| Intermediate| HiddenTriple              |
| Hard        | BoxLineReduction          |
| Expert      | HiddenRectangle           |
| Master      | BivalueUniversalGrave     |
| Extreme     | Backtracking              |

The generation cap and technique classification are intentionally
different: the generator uses a coarse enum-order gate to quickly
reject puzzles during generation, while classification uses the full
technique-to-tier mapping after solving. The SE range provides the
continuous scale used by `generate_for_se()`.

---

## 8. Arithmetic Counting

Select at most six distinct native cell or sector/digit exactly-one equations.
Weight each by an integer in `{-2,-1,1,2}`, giving `a·x=β`. Assume the opposite
Boolean value `b` for target `t` and form `R=β−a_t b`. Reject the assumption if:

1. `R` is outside `[Σ(j≠t) min(0,a_j), Σ(j≠t) max(0,a_j)]`, or
2. `gcd(|a_j| : j≠t)` does not divide `R` (zero divides only zero).

Every Boolean contribution is inside its coefficient interval and every integer
contribution is divisible by that gcd. Therefore either rejection proves the
reported placement or elimination for every completion respecting the masks.

The certificate records the source equations, weights, target, and an exact
candidate-state hash. Its checker reconstructs the equations independently.
Search uses a bounded connected beam and is incomplete. The technique runs after
Death Blossom in hint and profile dispatch; it is omitted from recursive forcing
propagation. See [Arithmetic Counting](../../docs/arithmetic-counting.md) for the
full proof, examples, API, validation, assumptions, and attribution.

---

## Appendix A: Notation Summary

| Symbol | Meaning |
|--------|---------|
| *C*    | Cell space (81 cells) |
| *X*    | Candidate space ⊆ *C* × {1..9} |
| *S*    | Sector space (27 units) |
| *L*    | Link space (strong + weak edges) |
| *X(c)* | Candidates of cell *c* |
| *X_d(S)* | Cells in sector *S* with candidate *d* |
| *β, κ* | Base cells, cover cells (fish) |
| *φ, ε* | Fin cells, elimination cells (fish) |
| *cands(A)* | Candidate union of ALS *A* |
| RCC    | Restricted Common Candidate |
| ALS    | Almost Locked Set |
| AIC    | Alternating Inference Chain |
| SE     | Sudoku Explainer rating |

## Appendix B: CandidateFabric Sector Convention

The implementation uses a flat sector index:

```
Sector  0.. 8  →  rows r₁..r₉
Sector  9..17  →  columns c₁..c₉
Sector 18..26  →  boxes b₁..b₉
```

`sector_digit_cells[s][d-1]` gives a `u16` bitmask of positions within
sector `s` that have candidate `d`. `cell_sectors[c]` gives the 3
sectors containing cell `c`: `[row, col, box]`.
