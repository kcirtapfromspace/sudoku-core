# Different-obstruction campaign: independent classification and limits

**Follow-up:** the [expanded audit](novelty20_report.md) verifies a published
EL12c38 tridagon followed by one cell split as a complete explanation of
r1c8≠3. A separate [seven-node local-only forcing certificate](deeper_forcing_report.md)
uses two additional split levels. The earlier one-pivot no-hit here does not
exclude either stronger explanation.

Audit date: 2026-10-20. This is a bounded search and classification report, not
a claim of a new mathematical principle or 100% originality. The previously
studied uniform native-27 incidence is excluded from this campaign.

## Outcomes

| Discovery arm | Recorded core/endpoint cases | Checked one-pivot known compositions | No checked selected-core-pivot composition | Short-budget inconclusive |
| --- | ---: | ---: | ---: | ---: |
| First two SAT batches | 327 | 316 | 11 | 0 |
| Supports deliberately omitting standard Loki12 geometry | 20 | 20 | 0 | 0 |
| Final diverse Mith batch | 135 | 92 | 41 | 2 |

These are case records, not 482 independent puzzles or necessarily 482 distinct
conclusions. Different selected supports can describe the same endpoint.
The 327 cases all contain a standard four-box Loki12 embedding. In the 20
geometry-omitting cases the *selected support* omits that embedding, but 18
already expose a direct Loki12 contradiction in the closed **full opposite
state**. The other two reduce by an additional case split and a two-color odd
cycle. Geometry omission alone therefore failed as a novelty criterion.

The Mith batch has 123 standard Loki12 overlaps and 12 supports without one.
Among the latter, ten were not classified by the initial pivot test. A targeted
follow-up tested every effective non-singleton pivot on the **whole board** for
all ten, without timeouts; none produced a complete reduction in the tested
lemma menu. Their exact template primals were independently checked. This is
limited negative evidence, not exclusion of every known solving method.

Compact coordinate mappings, status counts and semantic proof hashes are in
`different_audit_manifest.json`, `different_audit_geometryfree_manifest.json`
and `different_audit_mith_manifest.json`. Bulky discovery fixtures and branch
traces stay under `target/`. Default CI verifies the committed mappings and
aggregate metadata, while the two retained public cases receive full exact
audits. It does not independently replay all absent bulk traces from the
compact manifests. Key full-state mappings for geometry-omitting cases 00/12
are retained in `different_audit_geometryfree_direct.json`.

## A valid deduction that is a known composition

`different_control.json` is the source-attributed Monhard113 unique puzzle.
The conclusion is `r9c2=5`. Its 13-cell Rust proof has 137 nodes, and the exact
single-digit-template primal has 12 positive rational terms. Independent
uniqueness, givens/base/opposite traces, all 2925 conjunctions of three native
houses, template equalities and Rust replay pass in `different_audit_verify.py`.
The three-house checks use only the three selected house constraints, not
additional peer edges from unselected houses.

Despite those strong separation checks, the deduction has a short familiar
explanation. Assume `r9c2!=5`, then split `r8c8` over `{5,9}`. The value 5 gives
a three-house contradiction at center 9; the value 9 leaves the published
Loki12 graph on cells
`[29,35,36,43,46,51,56,62,63,69,73,79]` with a common three-digit palette.
Every branch is covered, and the stored coordinate isomorphism is checked.
`different_audit_control_classification.json` and
`different_audit_control_result.json` retain this negative-novelty control.
Four corruptions and positive/negative recognizer controls also pass.

Two additional explicit Mith515 cases, under `target/different-explicit-corpus`
as cases 17/18, are separate from the 135-case Mith batch. They eliminate 3/4
at `r4c2`. Both have an even shorter established explanation: after the
opposite target, split `r2c1` over `{3,4}`; both alternatives give local
Hall/three-house contradictions. No tridagon lemma is needed.

## Retained research lead, not an originality result

`different_lead.json` is the source-attributed Mith515 unique puzzle, with the
conclusion `r1c8!=3`. The raw selected native graph has 20 vertices and 51
edges. Its ordinary chromatic number is exactly four: an exhaustive uniform
three-color refutation and a directly checked four-color witness are recorded.
Actual opposite-state lists use seven labels `{1,2,3,4,7,8,9}`. It is a
**list-coloring obstruction**, not a uniform four-color graph requiring five
colors.

The public replay `replay_different_lead.json` has 85 nodes and binds the exact
base snapshot and target polarity. The full single-digit-template relaxation
has an exact 18-term positive primal. The independent audit also checks the
unique completion, givens/base/opposite traces, all 2925 three-house
conjunctions, and the recorded whole-board pivot follow-up.

The standard four-box geometry recognizer completes without a hit. More
strongly, a separate exhaustive edge-preserving-injection search of the
abstract published Loki12 graph into the selected 20-cell graph completes
after 1020 search visits with no embedding, even allowing extra native edges.
This is not a timed-out isomorphism search.

Every effective non-singleton pivot on the full board was tested after the
opposite target. For each pivot at least one branch remains unclassified by
the specific menu: exact Hall/all81 crossing closure, a direct local
contradiction, a published Loki12 with actual lists contained in three digits,
or a uniform two-color odd cycle. The negative record is
`different_audit_lead_classification.json`; its snapshot hash and counterbranch
closures are independently rechecked. This **does not** exclude conditional
hyperholes, all 630 published impossible patterns, other list-coloring
patterns, Exocet/MSLS derivations, generalized chains, two additional pivots,
or arbitrary forcing. The exact incidence's historical priority is unresolved.

Reproduce the complete independent audit with
`python3 research/global_obstructions/different_audit_verify.py`.
To rerun the lead's 32 eligible full-board pivots, use
`python3 research/global_obstructions/different_audit_conditional.py research/global_obstructions/different_lead.json --all-board-pivots --max-pivots 81 --seconds 12 --output target/different-lead-pivot-replay.json`.
The stored negative record, rather than elapsed time, establishes that all
32 effective non-singleton pivots completed in the retained run.

## What the primary sources already establish

* Andrew Stuart's author-maintained
  [Tridagons](https://www.sudokuwiki.org/Tridagons) page was created
  2025-09-19 and reports earlier 2022–2025 development. Its sparse-list
  examples allow the triple to be distributed as `{12},{13},{23}` rather than
  identical three-digit lists. It documents multiple guardians, pointing-pair
  guardians and pseudo-subset combinations. Consequently sparse lists or
  several guardians alone do not establish a new family. The page contains
  historical uniqueness terminology; the non-colorability/guardian deduction
  used here needs no uniqueness premise.
* The immutable author-source Loki example is
  [Berthier's `Loki-11.9.clp` at fd9207f8](https://github.com/denis-berthier/CSP-Rules-Examples/blob/fd9207f8dc1eb6e2e655048ba3b327d423900e1b/Sudoku-a/Tridagons/Loki-11.9.clp).
  The local `real_Loki.json` records its source hash. Every claimed Loki match
  here supplies an actual coordinate map and verifies native adjacency.
* Berthier's author-maintained
  [Sudoku-classif](https://github.com/denis-berthier/Sudoku-classif) describes
  630 published impossible patterns, generated OR-relations, and combinations
  with OR-chains. This audit does not load that complete catalogue, so failing
  the two specific recognizers cannot establish a new pattern.
* Allan Barker's
  [Dark Logic](https://sudoku.allanbarker.com/sweb/gen2/blacklog.htm) already
  treats impossible cores, guardians and virtual truths. No original date for
  that mutable author page was recovered. A guardian deduction or a
  case-split presentation is not a new principle.
* Exact list-coloring search and stable-set-column relaxations are established;
  see Lucci, Nasini and Severin,
  [arXiv:1812.00040v1](https://arxiv.org/abs/1812.00040v1), submitted
  2018-11-30. Eppstein's
  [Solving Single-digit Sudoku Subproblems](https://arxiv.org/abs/1202.5074v2),
  revised 2012-03-04, studies single-digit Nishio/template problems. Feasibility
  of our template mixture is a meaningful relaxation separation, not a new
  list-coloring principle or a separation from full multi-digit POM/forcing.

Queries in this phase included `"chromatic-patterns-t39885"`,
`"Sudoku" "list coloring" impossible patterns`, and
`"Sudoku" "tridagon" "four" "guardian"`, plus direct primary-author page
and repository reads. The primary chromatic-patterns forum thread remained
inaccessible; the 2010 Ivanyi–Nemeth PDF host also failed to serve its file.
Unindexed or inaccessible prior art remains a material uncertainty. No newly
named rule, new production algorithm, or 100% historical originality is
claimed by this campaign.
