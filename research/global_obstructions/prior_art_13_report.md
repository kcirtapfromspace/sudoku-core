# Prior-art audit of the 13-cell four-digit obstruction

Research date: 20 October 2026 UTC (19 October in the working timezone).
Scope: the exact support `[9,15,16,17,25,27,34,36,43,46,47,51,53]`
in `structural_four_digit_13.json`, with ordinary 9x9 Sudoku peers.
This report separates an obstruction's mathematical explanation from historical
priority for a particular drawing or program.

## Finding

**This example does not establish a new mathematical solving principle.** Its
complete four-color obstruction already occurs in a spanning **7-hyperhole**:
a clique expansion of the seven-cycle, with clique sizes `(1,3,1,2,2,2,2)`.
The native graph adds six unnecessary edges. The published hyperhole coloring
formula applies directly to this explicitly mapped subgraph. After the guardian
is removed and box subsets are propagated, an ordinary grouped discontinuous
Nice Loop eliminates each of the four digits from one cell. The guardian
conclusion then follows by established forcing/dark-logic reasoning.

I did not establish an older publication of this exact 35-edge graph or its exact
Sudoku embedding. Their first appearances remain **unresolved**. Failure to find
an indexed match is not evidence of novelty. In particular, an independently
generated example, an independently written proof, and independent solver code
do not establish mathematical priority.

## Exact graph and geometry

For a lightweight standard-library CI check, run
`python3 research/global_obstructions/prior_art_13_verify.py`. It reconstructs
the spanning-family map, checks all 715 four-subsets, a concrete independent
triple, the native proper-five-color completion, actual strong-box/weak-house
premises, and all 128 seven-group Boolean truth assignments. It explicitly
rejects Python `-O`. This is mathematical verification, not historical proof.

Run `python3 research/global_obstructions/prior_art_13_invariants.py` to reproduce
the machine-readable `prior_art_13_invariants.json`. NetworkX 3.1 supplies graph
isomorphism/automorphism enumeration; graph construction, exact coloring,
independence enumeration, and the canonical-label convention are explicit.

| Invariant | Value |
| --- | --- |
| Vertices / native peer edges | 13 / 35 |
| Independence number / chromatic number / clique number | 3 / 5 / 4 |
| Independence polynomial | `1 + 13z + 43z^2 + 32z^3` |
| Degree multiset | eight 5s, five 6s |
| Triangles / graph automorphisms | 24 / 8 |
| Vertex critical | Every one-vertex deletion has an explicit 4-coloring |
| Edge critical | No; seven native edges can individually be deleted |
| Graph6, vertices in the stated support order | `L~]@[Gj@OIgVCn` |
| Canonical adjacency SHA-256, documented IR-v1 convention | `5423759c913cfaf0de99ffecc028f4ae60a204ba8282097f900abc14ac8358e3` |

The graph is not a line graph (NetworkX inverse-line-graph recognition). It is
not `C13^3` or the undirected Paley graph of order 13: those are 6-regular with
39 edges. These limited comparisons do not constitute a complete named-graph
catalogue search. Abstract graph isomorphism also does not identify a Sudoku
embedding, because it discards row/column/box labels.

Occupancies, in natural house order, are:

* rows: `[0,4,1,2,2,4,0,0,0]`;
* columns: `[3,1,1,0,0,0,2,4,2]`;
* boxes: `[1,0,4,4,0,4,0,0,0]`.

The lexicographically minimum sorted support under **all 3,359,232** standard
spatial Sudoku symmetries (bands/stacks, internal rows/columns, transposition) is
`[0,1,2,3,9,27,30,36,39,46,47,49,50]`. A witnessing map, applied without
transposition, is old-to-new row indices `[2,0,1,3,4,5,6,7,8]` and old-to-new
column indices `[3,4,5,6,7,8,1,0,2]`. Indices here are zero based. Digit
renaming does not affect this support key. The abstract canonical key was also
checked against 30 deterministic random vertex relabelings.

## A checkable match to an established graph family

Let `H` consist of cliques `X1,...,X7`, joined completely exactly when their
indices are successive modulo seven. The following is an explicit vertex map;
within a bag, list order can be used to label the abstract vertices `(i,j)`.

| Bag | Sudoku cells | Zero-based board indices | Size |
| --- | --- | --- | --- |
| X1 | r2c1 | 9 | 1 |
| X2 | r2c789 | 15,16,17 | 3 |
| X3 | r3c8 | 25 | 1 |
| X4 | r45c8 | 34,43 | 2 |
| X5 | r6c79 | 51,53 | 2 |
| X6 | r6c23 | 46,47 | 2 |
| X7 | r45c1 | 27,36 | 2 |

All within-bag edges and all complete joins of consecutive bags are native
Sudoku peer edges. Retain these 29 edges, and delete exactly these six extra
native edges (zero-based endpoints):
`(15,51),(16,34),(16,43),(17,53),(27,34),(36,43)`.
The resulting graph is **exactly** `H`, graph6 `L~]?SGI@OI_V?n` in the original
support order. It is a spanning subgraph, not an induced hyperhole in the
original 35-edge graph. The distinction matters because the cited paper uses
the word *contains* for induced subgraphs.

An independent set in `H` uses at most one vertex per bag and cannot use
successive bags, so it has at most three vertices. Three are possible. Thus
`alpha(H)=3`, while its largest consecutive pair of bags has size four, giving
`omega(H)=4`. Thirteen vertices cannot be partitioned into four independent
sets of size at most three. This proves `chi(H)>=5`; the artifact supplies an
explicit 5-coloring. It also supplies a 4-coloring for every one-edge deletion
of `H`, verifying that this 29-edge skeleton is edge-5-critical.

[Maffray, Penev and Vušković, *Coloring rings*, immutable arXiv v1,
27 July 2019](https://arxiv.org/pdf/1907.11905v1), page 3, defines hyperholes
and states Lemma 1.1:
`chi(H)=max{omega(H),ceil(|V(H)|/alpha(H))}`. Substituting the values above
gives five. This is an exact established-family match for the whole necessary
obstruction, not a superficial shared vertex count. The paper attributes this
formula to Narayanan and Shende, *Static frequency assignment in cellular
networks*, Algorithmica 29 (2001), 396–409,
[DOI 10.1007/s004530010067](https://doi.org/10.1007/s004530010067).
[Narayanan's author-hosted preprint](https://users.encs.concordia.ca/~lata/Joffline.pdf)
is dated 11 August 1998 and studies weighted-cycle multicoloring. Its broader
algorithmic claims are not needed for the elementary lower bound here.

## Exact reduction to known Sudoku grouped links

Let `P={4,6,7,8}`. The original guardian cell r2c8 also allows 1. Assume
r2c8 is not 1, so all thirteen selected cells have domains contained in `P`.
Each of boxes 3, 4 and 6 then has four selected cells confined to four digits:
they form naked quads. The ordinary subset consequences remove `P` from the
remaining cells in those boxes. This is a sound step under the assumption;
the strong grouped links below must not be asserted unconditionally in the
original guardian-present state.

The script verifies directly on the actual `assumed_domains` snapshot that,
for **each** `d in P`, the only positions for `d` in each relevant box are:

* box 3: `r2c789` or `r3c8`;
* box 6: `r45c8` or `r6c79`;
* box 4: `r6c23` or `r45c1`.

Consequently this grouped Nice Loop is valid for each of the four digits:

```text
(d)r2c1 - (d)r2c789 = (d)r3c8 - (d)r45c8
        = (d)r6c79 - (d)r6c23 = (d)r45c1 - (d)r2c1
```

Here `-` is a weak link and `=` a strong link. Every weak link is justified
by a single native row or column; every strong link is justified by the
exhaustive grouped positions in one of the three boxes. Starting with r2c1=d
forces, successively, the opposite group in each box and finally r2c1!=d.
Thus the loop eliminates d from r2c1. Repeating for four digits empties r2c1,
refuting the guardian-false assumption. This is also a direct role-by-role
explanation of the hyperhole obstruction.

[HoDoKu's author-maintained Chains and Loops guide](https://hodoku.sourceforge.net/en/tech_chains.php)
defines grouped strong/weak links, discontinuous loops and their contradiction
rule. [Its Last Resort guide](https://hodoku.sourceforge.net/en/tech_last.php)
describes combining implications under an assumption and contradiction. The
retrieved guides lack reliable paragraph-level publication dates; their current
contents establish the method correspondence, not the first appearance of this
exact mapped example.

[Andrew Stuart's Empty Rectangles page](https://www.sudokuwiki.org/Empty_Rectangles)
states: “Empty Rectangles can be re-expressed as Grouped X-Cycles.” It also
credits Klaus Brenner with Triple Empty Rectangle examples. Its displayed
creation date is 12 April 2008, but a 2023 change coexists with older modification
metadata, so those dates cannot date each paragraph. I did **not** prove that
its illustrated Triple Empty Rectangle is the same thirteen-cell drawing.

## Rank, MSLS, guardians and other candidate names

For any one digit `d`, its true selected cells form a stable set. Therefore
`sum(v in S) x[v,d] <= alpha(G[S]) = 3`. Summing four such established
stable-set rank inequalities gives capacity 12 against thirteen cells.
[Doignon, Fiorini and Joret, immutable arXiv v1 (22 September 2008)](https://arxiv.org/abs/0809.3724v1)
discuss rank inequalities and Chvátal's 1975 work. This identity, the
pigeonhole count, graph coloring and the use of an independence bound are
established mathematics.

The word *rank* has different conventions in graph polyhedra and Sudoku set
logic. They should not silently be equated. A plain rank-zero MSLS proof using
only linear combinations of the original native equalities/nonnegativity
cannot refute this assumed state: it has an independently checked rational
point satisfying all 324 native equalities and all candidate masks. That point
assigns `1/4` to each selected palette candidate and violates the **additional
integer-valid** stable-set inequality (`13/4>3`). This does not exclude richer
MSLS/general-logic analyses with Boolean case distinctions or imported illegal
substructures. Nor does it claim independence from the known chain/fish steps
found under the assumption.

[Allan Barker's author page on general Sudoku logic](https://sudoku.allanbarker.com/sweb/general.htm)
describes truths, links, rank, illegal logic and guardians, including reducing
illegal substructures to virtual truths. [His Dark Logic page](https://sudoku.allanbarker.com/sweb/gen2/blacklog.htm)
explicitly explains an impossible core with guardian candidates, and states:
“If there is one guardian, it can be assigned.” These pages provide primary
prior art for the general guardian rule. Their retrieved pages do not supply
a trustworthy original publication timestamp. This report uses the checked
loop map above, rather than merely assuming every named rule is implemented
by a particular general solver.

[Tom Collyer's *The MSLS: Revolutions*, 20 October 2020](https://tcollyer.blogspot.com/2020/10/the-msls-revolutions.html)
discusses MSLS and intersection-of-partitions deductions and recounts an
earlier Nikoli puzzle. It reinforces that multisector counting and alternative
descriptions are well established. It does not supply an exact map to this
motif. The primary MSLS forum thread was inaccessible during this audit.

The standard [SK Loop description](https://www.sudokuwiki.org/SK_Loops)
uses four boxes, eight linked cell pairs and a sixteen-cell loop, with smaller
variants when cells are already solved. This motif uses seven groups of sizes
shown above, including a three-cell group. I established no exact identity
with an older SK Loop template; that is not a proof that every generalized
SK/MSLS formulation is different. [Exocet's author-maintained documentation](https://www.sudokuwiki.org/Exocet)
describes base/target arrangements and credits the earlier Barker/Bird/Champagne
work. Merely using four digits is insufficient to identify an Exocet, and no
explicit Exocet map was established here.

The familiar twelve-cell three-digit tridagon is a 3-color obstruction; the
present uniform-palette core is 5-chromatic. That numerical difference does
not establish a new family: the hyperhole reduction already explains it.
Searches for *quadrivalue oddagon*, *quadagon* and *four-value oddagon* found
no relevant exact indexed match. Terminology and incomplete indexing make
such negative searches especially weak evidence.

## Catalogue and historical limits

[Brendan McKay's primary graph catalogue](https://users.cecs.anu.edu.au/~bdm/data/graphs.html)
and [Olivier Lalonde's author implementation](https://github.com/lalondeo/gencrit)
provide complete edge-critical graph lists through order 13, including
8,311,809 edge-5-critical graphs on thirteen vertices. The 35-edge native
graph is not edge critical, so direct membership in that file is not an
appropriate test. Its 29-edge hyperhole skeleton is edge critical, and is
already covered by the explicit established-family construction above. I did
not download and search the multi-million-graph file, enumerate every
vertex-critical supergraph, or complete a House of Graphs database lookup.

The original EnjoySudoku chromatic-patterns and MSLS threads returned errors;
a direct HTTPS request to the chromatic-patterns thread returned connection
refused. These inaccessible primary discussions are a material historical
coverage gap. Search-engine summaries or mutable tutorial wording cannot
establish that no earlier identical drawing exists.

The search log is `prior_art_13_queries.json`; URLs with versioned arXiv IDs
and DOIs provide stable scholarly references. Other linked sites are mutable
and are cited as observed on the research date. No full copyrighted article
has been copied into this repository.

## Claims justified for the engine

| Claim | Audit conclusion |
| --- | --- |
| A new general mathematical rule | No: stable-set capacity, hyperholes, grouped loops and guardians are established |
| A new abstract minimal obstruction | Not supported: its necessary 29-edge obstruction is an exactly mapped established hyperhole |
| First occurrence of this exact 35-edge native graph | Unresolved |
| First occurrence of this exact Sudoku embedding | Unresolved |
| Thirteen cells is globally smallest for this separation | Not established by this audit |
| A reproducible separation from specified Hall/three-house/native-LP baselines | Supported by the separate exact fixtures/checks; not a historical-priority claim |
| Original implementation and checked certificate packaging in this engine | May be accurately described as this project's implementation; the underlying proof ideas require attribution |

Suggested description: **a generated thirteen-cell four-digit benchmark whose
obstruction is a spanning seven-hyperhole, with exact baseline separation and
replayable native-peer certificates**. A worthwhile follow-up is a bounded
canonical census of *baseline-separating embeddings*, explicitly quotienting
known hyperhole/chaining constructions rather than assigning new names to
their rediscoveries.
