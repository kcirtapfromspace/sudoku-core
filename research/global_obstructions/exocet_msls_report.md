# Exocet and MSLS audit of Mith515

**No direct basic Junior Exocet or native rank-zero MSLS explanation of r1c8≠3 was found.** There is an exact exclusion for the specified native MSLS class and for basic Junior Exocet's target restriction. These are scoped mathematical statements, not a claim that every technique called Exocet or MSLS has been excluded. The separate [deeper forcing audit](deeper_forcing_report.md) gives an established two-level forcing explanation of this same deduction; nothing here reclassifies that proof as Exocet/MSLS or as a new family.

The audited source is [`different_lead.json`](different_lead.json), the unforced retained state of the unique Mith source-line-515 puzzle. Its opposite assumption is **r1c8=3**, followed by the recorded exact Hall/all-81 crossing-three-house closure. Both snapshots were inspected on the full board and with pattern support restricted to the selected 20 cells. A move somewhere else on the board would not, by itself, explain the target elimination.

## Definitions and primary sources

[Allan Barker's own *Basic Principles*](https://sudoku.allanbarker.com/sweb/general.htm), Part 1, defines native cell/house-digit truths and cover links, disjoint bases, and equal-count rank-zero eliminations. Parts 2–3 separately introduce overlapping-triplet case analysis and dark/derived logic. The theorem below covers the native counting formulation, including arbitrary shapes and mixed native house types; it does not silently include those additional mechanisms.

[Tom Collyer's author-written MSLS discussion](https://tcollyer.blogspot.com/2020/10/the-msls-revolutions.html) explicitly notes disagreement about a universal formal definition. We therefore specify the exact formal class instead of treating “MSLS” as a uniquely standardized algorithm.

[Andrew Stuart's author-written Exocet implementation notes](https://www.sudokuwiki.org/Exocet) credit David P Bird and specify the familiar two aligned bases, two targets, companions and S-cell covers. The pinned [kyoyama-kazusa/Sudoku source](https://github.com/kyoyama-kazusa/Sudoku/blob/89667e24c67384a8234f40972d46776ce791ce6d/src/Sudoku.Analytics/Analytics/StepSearchers/Exocets/ExocetStepSearcher.cs) and [its author's tutorial](https://github.com/kyoyama-kazusa/Sudoku/blob/89667e24c67384a8234f40972d46776ce791ce6d/docs/tutorial/exocet/01-reasoning-of-exocet.md) provide a second primary implementation account, using exact digit-packing capacity in the cross-lines. They also document many broader variants.

Bird's original [JExocet Compendium](https://forum.enjoysudoku.com/jexocet-compendium-t32370.html) and [MSLS thread](https://forum.enjoysudoku.com/using-multi-sector-locked-sets-t31222.html) were unavailable: the web tool returned errors and the direct public client connection was refused. We do not claim to have read those documents. [`exocet_msls_sources.json`](exocet_msls_sources.json) records retrieval limits, the pinned implementation commit and source hashes. The current external library targets .NET11; no `dotnet` executable was installed here, so **the full external Exocet implementation was not run**. Our scan is a separately written, bounded implementation of the stated premises.

## Exact native MSLS exclusion

Let nonnegative indicators \(x_{c,d}\) satisfy each cell equation and each native row/column/box digit equation. Let \(S\) be any set of cells and \(C\) a collection of native house-digit covers that collectively cover every allowed candidate in \(S\). If \(|C|=|S|\), then

\[
 |S|=\sum_{c\in S,d}x_{c,d}
 \leq\sum_{(h,d)\in C}\sum_{c\in S\cap h}x_{c,d}
 \leq |C|=|S|.
\]

Both inequalities must be equalities. Thus every cover is fully occupied inside \(S\), so its candidate mass outside \(S\) is zero. A candidate counted by multiple covers inside \(S\) must also have zero mass. This proves the ordinary external and overlapping-cover/cannibal rank-zero conclusions using only native linear constraints. Fewer covers than cells would contradict the same feasible point outright.

The lead's existing **18-term exact full-board digit-template primal** projects to such an indicator point and has \(x_{r1c8,3}=1\). The checker verifies all 324 native equations, every domain restriction, the exact target assumption, and the original givens/closure traces. The closed opposite masks are subsets of the unforced masks, so this same point is feasible in both snapshots. Its restriction also satisfies every selected-cell truth and native weak-link capacity on the selected20.

Consequently **no direct native rank-zero cell-truth/house-digit-cover MSLS of any shape or size can eliminate this target or refute this opposite state**. This is stronger than a finite rectangle scan. A sequence using only these same linear MSLS eliminations also preserves the point, because each removed candidate has zero mass. This does not exclude imported almost-locked/derived truths, joint multidigit constraints, Boolean case distinctions, dark logic, or a larger conditional proof described informally as MSLS. It does not assert that MSLS cannot yield another valid deduction whose candidate already has zero mass.

## Basic Junior Exocet target restriction

There is a separate precise exclusion for the basic target restriction, often called ER1. Let two peer base cells have candidates only in digit set \(D\). There are three cross-lines, a chosen band/stack, two target cells, and no base digit in companions; remaining inner cross-line cells see both bases. Every base digit has native S-cell packing capacity at most two.

For one digit \(d\), examine any legal full native \(d\)-template. It contains at most one base cell. If it contains a base cell, base visibility and the companion restriction exclude the other inner positions, while the three cross-lines require three occurrences and S permits at most two. Hence at least one occurrence lies in the targets. In indicator notation every admissible one-digit template satisfies

\[
 \sum_{t\in T}z_{t,d}\geq\sum_{b\in B}z_{b,d}.
\]

The inequality is preserved by mixtures of templates. Summing over \(d\in D\), the two base-cell equations require target base-digit mass at least two. Two target-cell equations supply total mass exactly two. Therefore all non-base target candidates have zero mass.

The retained full-board template point gives the opposite target candidate mass one, so **no direct basic Junior Exocet ER1 restriction can explain its removal**. This conclusion relies on the full one-digit-template relaxation, not merely the 324-equation native LP. It does not exclude compatible-pair checks or other joint-digit, locked-member, Senior, Weak, Double, Franken/Mutant or conditional Exocet extensions.

## Finite configuration scans

The basic scan uses two live bases in one box-line, a union of three/four digits, three cross-lines with one in each perpendicular stack/band, two live nonseeing targets, no base digits in companions, and exact S-cell native packing capacity. Solved S-cell values are included. The capacity test is broader than requiring an explicit two-house S-cover, so a failure with three mutually nonseeing occurrences also invalidates the cover condition.

| Full-board snapshot | Eligible base pairs | Cross-line geometries | Structural target layouts | Capacity-valid patterns |
| --- | ---: | ---: | ---: | ---: |
| Unforced | 27 | 243 | 14 | 0 |
| Closed r1c8=3 opposite | 18 | 162 | 5 | 0 |

A relaxed additional scan permits targets that see each other and removes the three/four-digit bound on the base union. It still uses two bases, three cross-lines and no locked-member extension. Its counts are **44/396/37/0** unforced and **27/243/21/0** opposite. Explicit three-cell native independent-set witnesses refute the S-capacity premise in all **58** relaxed structural layouts. The basic19 are included in these58. The independent checker uses a different geometry iterator and checks every witness's digit availability and native nonvisibility. With all relevant bases, targets and base-digit S-occurrences required to lie within selected20, neither scan reaches an eligible structural layout.

The rectangular MSLS scan follows the pinned author's [pattern catalogue](https://github.com/kyoyama-kazusa/Sudoku/blob/89667e24c67384a8234f40972d46776ce791ce6d/src/Sudoku.Analytics/Analytics/Construction/Patterns/MultisectorLockedSetPattern.cs) and [cover-count algorithm](https://github.com/kyoyama-kazusa/Sudoku/blob/89667e24c67384a8234f40972d46776ce791ce6d/src/Sudoku.Analytics/Analytics/StepSearchers/LockedSets/MultisectorLockedSetsStepSearcher.cs): sizes3×3,3×4,4×3,4×4,4×5,5×4; neither row/column selection lies wholly in one band/stack. Each digit independently chooses its smallest row, column or box cover family; ties are considered for eliminations. A 3×3 rectangle permits at most one placed hole; larger rectangles permit at most two.

There are **74,601** rectangles per snapshot, of which **207 unforced** and **41 opposite** meet those vacancy limits. None has rank zero, so this exact variant finds **no other move**, as well as no explanation of the target. Restricting its unplaced pattern cells to selected20 leaves971/3504 rectangles but no vacancy-eligible configuration. This restriction is a support scan on the real state, not a replacement of all outside masks by unrestricted candidates.

## Reproduce and acceptance limits

```sh
python3 -S research/global_obstructions/exocet_msls_scan.py
python3 -S research/global_obstructions/exocet_msls_verify.py
```

The standard-library verifier recounts the unique completion, replays the givens and opposite traces, verifies the existing exact template certificate and12 malformed-witness controls, reconstructs all58 relaxed geometries independently, checks every explicit three-packing, and reruns the bounded scans. An artificial, explicitly completed positive Junior control produces12 expected eliminations; adding one S-cell occurrence breaks its capacity premise and produces no pattern. These are detector controls, not additional puzzle discoveries. It runs in approximately one second. Raw source copies stay under `target/`; no external library or numerical solver is needed for acceptance.

There is **no positive Exocet/MSLS explanation to certify from these scans**. The separate exact forcing certificate does explain the target. The exclusions here apply only to the stated native-cover class and basic Junior ER1 rule; they supply no “all known methods” or historical-originality conclusion.
