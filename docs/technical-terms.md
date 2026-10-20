# Technical terms

This glossary defines the project terms used in the documentation.
API identifiers and published technique names retain their original spelling.

| Term | Meaning |
| --- | --- |
| candidate | A digit that a cell can still contain under the current premises. |
| candidate mask | A bitset that records the candidates of a cell. |
| candidate state | The placed values and effective candidate masks of all 81 cells. |
| cell | One position in the Sudoku board. |
| given | A fixed value supplied by the original puzzle. |
| placement | A deduction that sets a cell to a digit. |
| elimination | A deduction that removes one or more candidates. |
| house | One native row, column, or box. |
| sector | The house identifier used by the solver implementation. |
| peer | A different cell that shares a native house with the specified cell. |
| band | Three adjacent rows aligned with the boxes. |
| stack | Three adjacent columns aligned with the boxes. |
| domain | The allowed values of a variable. A cell domain is its candidate set. |
| completion | A complete board that satisfies the applicable rules and premises. |
| unique puzzle | A givens-only puzzle with exactly one completion. |
| premise | A condition assumed by a theorem or a proof certificate. |
| witness | Explicit data that a checker uses to confirm a stated property. |
| certificate | Finite evidence that a checker can verify against stated premises. |
| replay | A candidate state and certificate used to check or apply a deduction. |
| state hash | A digest that binds evidence to specified values and candidate masks. |
| contradiction | A condition that no assignment can satisfy under the stated premises. |
| refutation | A proof that the specified premises have no respecting assignment. |
| fixed point | A state unchanged by a completed pass of the specified rules. |
| cutoff | A limit that stops a search before the specified work is complete. |
| baseline | The exact rule set, implementation, limits, and input used for comparison. |
| native rules | The standard Sudoku row, column, box, and cell requirements. |
| relaxation | A problem that retains only specified constraints of the complete problem. |
| Hall obstruction | A cell subset whose allowed digit union is smaller than the subset. |
| matching | An assignment of distinct allowed digits to the cells in a specified set. |
| atom | A disjoint cell group defined by membership in selected houses. |
| boundary signature | The digit sets assigned to the atoms of a selected-house decomposition. |
| guardian | An additional candidate that can prevent a specified contradiction. |
| guardian relation | A proved OR statement that requires at least one guardian literal to be true. |
| literal | A candidate assertion or its negation. |
| forcing | Reasoning that tests an assumption and follows its consequences. |
| exhaustive split | Branches that cover every allowed value of the selected variable. |
| template | A full placement of one digit with one cell in each row, column, and box. |
| template compatibility | The requirement that templates for different digits have disjoint cell sets. |
| support | An assignment that respects the specified constraints and contains the specified candidate or template. |
| arc consistency | The property that every retained value has support across each specified binary constraint. |
| linear relaxation | The specified linear constraints with continuous variables instead of Boolean variables. |
| LP | Linear programming. |
| rank inequality | A bound on how many selected variables can be true together. |
| independent set | A vertex set with no edge between any two selected vertices. |
| independence number | The largest number of vertices in an independent set of a specified graph. |
| clique | A vertex set in which every pair has an edge. |
| coloring | An assignment of labels to graph vertices. |
| proper coloring | A coloring in which adjacent vertices have different labels. |
| chromatic number | The smallest number of labels needed for a proper coloring. |
| fractional coloring | Nonnegative weights on independent sets that provide the specified vertex coverage. |
| hyperhole | A cycle of at least four nonempty cliques, with complete edges only between adjacent cliques. |
| chordal graph | A graph without an induced cycle of length four or more. |
| simplicial vertex | A vertex whose neighbors form a clique. |
| primal witness | An explicit feasible assignment for the stated optimization constraints. |
| Farkas certificate | Weighted constraints that prove the specified linear system has no feasible point. |
| trade | A replacement of selected values that preserves the required external conditions. |
| unavoidable set | A cell set that every uniquely solving clue set must intersect under the specified completed-board model. |
| historical novelty | A claim about whether earlier work contains the stated method or result. |
| soundness | The property that every accepted deduction follows from its stated premises. |
| completeness | The property that a method finds every result in its explicitly stated problem class. |
| atomic application | An operation that commits the complete verified change or leaves the input unchanged. |
| serialization | Conversion of structured data to a stored or transmitted format. |
| scene | Portable display data for a grid explanation and its evidence status. |
| frame | One board snapshot with related descriptions, highlights, and links. |
| highlight role | The meaning of a cell or candidate mark, independent of its display color. |
| provenance | The recorded source and derivation of a state or fact. |
| prior art | Earlier work relevant to a mathematical method, pattern, or implementation claim. |
| GAC | Generalized arc consistency: every retained candidate has support in the specified constraint. |
| dynamic program | An algorithm that combines stored results of smaller subproblems. |
| exact cover | A selection of sets in which each required item occurs exactly once. |
| canonical fingerprint | A deterministic representation under the explicitly specified equivalence operations. |
| iff | The mathematical abbreviation for *if and only if*. |

The mathematical verbs *enumerate*, *refute*, and *propagate* describe specified finite operations.
The software verbs *serialize*, *deserialize*, and *replay* describe operations on stored data.
The related document defines the exact operation and its limits.

## Cell and sector indices

Board coordinates use `r1c1` through `r9c9`.
API cell indices use zero-based row-major order: `cell = 9 * row + column`.
The candidate bit for digit `d` is `1 << (d - 1)` in `CandidateState`.
The internal `BitSet` uses bit `d` instead.

`CandidateFabric` sector indices use rows 0–8, columns 9–17, and boxes 18–26.
The solver specification records any other convention used by an internal module.
