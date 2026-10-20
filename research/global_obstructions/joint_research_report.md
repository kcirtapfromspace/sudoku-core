# Stronger discovery baseline and joint digit placements

The stronger offline baseline resolves **all eight** previously recorded full
older-engine stalls. These puzzles are retained as known-method controls. No new
solving family qualifies for production promotion in this campaign.

The twenty-cell Mith515 certificate, replay and regression tests remain valid.
Its classification is a known-method benchmark. The previous audit already
gave an established catalogue-and-forcing explanation; the revised comparison
also solves it using the older dispatcher after hypothetical candidate retention
is corrected.

## What changed in the comparison

| Comparison stage | Solved former stalls | Live configured stalls | Inconclusive |
| --- | ---: | ---: | ---: |
| Original frozen ladder, recorded earlier | 0 | 8 | 0 |
| Same frozen ladder with candidate-preserving hypothetical rebuilds | 4 | 4 | 0 |
| Catalogue/chain/declared pattern closure on those retained states, with complete supplied givens-domain history scans | 8 | 0 | 0 after larger-budget reruns |

The initial historical catalogue profile had two node cutoffs, including one
already solved input. Those outcomes are preserved in
`strong_baseline_campaign_results.json`. Separate larger-budget reruns in
`strong_baseline_expanded_results.json` completed both scans. A cutoff was never
counted as a stall. Completion of a solved board does not assert that every
possible forcing query or named detector was run after the board was solved.

`discovery_pipeline_old.py` extracts immutable revision
`0532d63d9241a69c1fb57e9c34a1533effe1981d` into a temporary build directory. In its
instrumented solver, a candidate rebuild intersects the previous masks with
newly legal masks. This retains eliminations inside singles and full logical
assumptions as well as outside them. Three targeted Rust regressions exercise
the rebuild and both propagation paths. The original extracted source is
checked against Git; the binary and source transformations have recorded hashes.

The eight public input puzzles are independently counted unique and their
givens-to-state traces replayed. The corrected comparison adds 251 outer
endpoints, all checked against the sole completion. These are method-comparison
traces, not independent certificates for every older detector's internal proof.
The eight givens patterns have different symmetry-invariant signatures, so they
belong to at least eight different native Sudoku symmetry orbits.

## Precisely declared offline profile

`strong_baseline_core.py` repeatedly applies:

1. Strict native occurrences of the 592 public applicable graphs in the pinned
   630-pattern catalogue, over all 84 three-digit palettes and both orientations.
2. Full single-house Hall support and all 81 crossing three-house support checks.
3. Persistent guardian OR clauses, unit consequences and signed binary native
   implication paths, with explicit path reasons.
4. The previously specified basic and relaxed Junior Exocet variants and
   rectangular rank-zero MSLS variants.
5. Bounded failed literals with up to two additional cell or guardian splits,
   exact exhaustive branch coverage and retained candidates and guardian facts.

Each clause binds the catalogue hash, native embedding, palette, guardian list,
origin domains and puzzle context. It survives monotone candidate restrictions,
even if its strict occurrence is no longer recognizable. Branch-derived facts
remain local to that branch. A complete strict scan on a looser ancestor covers
all later strict occurrences; this does not cover sparse variants excluded by
the declared strict model. The 38 no-free-cell catalogue cases retain the
specific exclusion explained in `catalog630_report.md`.

Default forcing includes fresh catalogue discovery in assumptions when no
complete ancestor scan covers them. Junior/MSLS scans run at the root; enabling
them inside hypotheses is a separate profile option. Guardian splits have a
declared width bound of six. A focused endpoint list establishes only a focused
fixed point. A phase that changes domains restarts the menu. Time, node, clause,
round or branch limits produce an inconclusive outcome.

This is an explicitly bounded subset of known methods. Its unit propagation,
binary reachability and guardian branching do **not** implement the full family
of CSP-Rules ORk-whips or ORk-forcing-whips described in the
[catalogue author's pinned README](https://github.com/denis-berthier/CSP-Rules-V2.1/blob/7eefd71566b37e55c15f330d474140bfc4f4371d/README.md).
An eventual survivor must be compared with those additional implementations or
carry a precise limitation; this baseline is not an oracle for all known methods.

## Joint-placement search and its theorem

Every digit has 46,656 full native placement templates. Retain those compatible
with the candidate state and its singleton values. Every Sudoku completion
selects one template per digit, and templates selected for distinct digits have
disjoint cell sets. Thus a template lacking a disjoint support in another digit
can be removed; so can one lacking a mutually disjoint pair of supports in two
other digits. Induction over these removals preserves every completion. Emptying
one digit's template domain refutes the state. Apply this to a conclusion's
opposite literal to certify the conclusion. Uniqueness is unnecessary for this
theorem; it is required separately for research-host acceptance.

`joint_templates_search.py` implements exact pair support closure and ternary
support filtering across the entire board, including every band and stack.
`joint_audit_verify.py` independently enumerates templates using recursive
row/column/box masks and replays all claimed unsupported-template deletions.
The searcher uses nine-column permutations. A ternary separation certificate
requires a nonempty, fully checked nine-digit pairwise arc-consistent checkpoint
before the first ternary deletion. This is a precise separation of two declared
template relaxations; it is not a proof of novelty or separation from forcing.

The two stable controls include a seven-event pair refutation of the known
twenty-cell opposite and an eleven-event refutation of the known thirteen-cell
control after a nonempty pairwise checkpoint. Their original per-digit template
LP witnesses remain valid. Feasibility of an LP restricted to the templates
remaining after pairwise filtering is a separate, unchecked assertion.

| Screening epoch | Unique hosts | Opposites | Local refutations | Pair refutations | Refutations after pair checkpoint using ternary support | Ternary support fixed points |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Original older-engine states | 8 | 1,200 | 15 | 458 | 639 | 88 |
| Candidate-preserving older-engine states | 4 | 631 | 1 | 199 | 343 | 88 |

Both bounded sweeps completed without query cutoffs. These epochs overlap and
must not be presented as independent puzzle samples. The 88 surviving queries
are inconclusive under ternary support filtering; they are not Sudoku solutions
or stronger-baseline survivors. Endpoint manifests expose the complete candidate
menus, source states, raw and closed opposite hashes and per-query outcomes.
Sixteen retained controls from eight unique hosts received independent exact
proof and premise replay. Screening counts alone do not certify every reported
refutation.

The next bounded campaign used twelve previously untested puzzles, alternating
six Monhard and six Mith sources from the pinned corpus. All twelve are unique
and were solved by the candidate-preserving older ladder before a stronger
phase was needed. Their public stage logs replay 790 additional outer endpoints;
the run took 146.5 seconds. This is further negative discovery evidence, not an
exhaustive search of the corpus or evidence that all joint obstructions reduce
to known methods. See `discovery_pipeline_expanded_results.json` and its exact
premise/stage checker.

## Acceptance and historical classification

The [independent audit](joint_audit_report.md) identifies established Pattern
Overlay cross-digit filtering, Brouwer's digit-pattern jigsaw model and classical
constraint consistency as prior art. Joint compatibility is a useful search
model; its general principle is already known. Historical priority for a future
specific family needs a separate primary-source and pattern audit.

A candidate can advance only with an explicit theorem and independently checked
exact proof, an unchanged completed union baseline at the same retained state,
exact local and LP witnesses for any claimed separations, and examples in at
least three unique givens-only puzzles proved to occupy different native
symmetry orbits. The givens history, all enabled baseline variants and resource
limits must be bound to those states. Literature status is recorded separately
from mathematical validity. An older-engine stall, known template consistency
separation or a cutoff satisfies none of the missing obligations.

## Reproduction

All acceptance checks use Python's standard library:

```sh
python3 -S research/global_obstructions/discovery_pipeline_verify.py
python3 -S research/global_obstructions/discovery_pipeline_expand_verify.py
python3 -S research/global_obstructions/joint_audit_verify.py --self-test
python3 -S research/global_obstructions/joint_audit_relations.py --self-test
python3 -S research/global_obstructions/joint_audit_patterns.py --self-test
python3 -S research/global_obstructions/strong_baseline_verify.py --self-test
python3 -S research/global_obstructions/joint_templates_manifest.py verify --count --replay-closures 8
```

The offline retained engine requires Rust to rebuild:

```sh
python3 -S research/global_obstructions/discovery_pipeline_old.py --folder target/retained-reproduction
python3 -S research/global_obstructions/discovery_pipeline_compare.py --build-folder target/retained-reproduction --output target/retained-reproduction-results.json
python3 -S research/global_obstructions/strong_baseline_campaign.py --retained target/retained-reproduction-results.json --output target/strong-reproduction
```

Use the declared larger catalogue profile to revisit any first-pass cutoff.
`discovery_pipeline_expand.py` continues into previously untested pinned author
puzzles, alternating both baseline parts and retaining their separate trace
grammars. Its solved cases, fixed points and cutoffs are separate outcomes; a
new positive deduction remains unaccepted until its full stage log is replayed.
