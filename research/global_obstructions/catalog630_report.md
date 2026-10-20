# Imp630 audit of the Mith515 twenty-cell deduction

The audit found a known **EL12c38 tridagon with one guardian in the full base state**. It proves `r8c6=7`. This is another valid deduction in the puzzle; the match alone does **not** prove the requested elimination `r1c8≠3`. No strict catalogue occurrence remains after the opposite assumption and local closure. That later absence does not invalidate an OR relation established earlier. This audit provides no basis for claiming an original mathematical family.

The input is `different_lead.json`, SHA256 `45d6f8e79f9dfb859ba6827dadb3934a1cc5c1b1ada822537b5346b926941358`. Its research target assumes zero-based cell7 is3; the engine conclusion eliminates that candidate.

## Primary source and the38 omitted geometries

The implementation is [Denis Berthier's CSP-Rules-V2.1 Imp630 directory](https://github.com/denis-berthier/CSP-Rules-V2.1/tree/7eefd71566b37e55c15f330d474140bfc4f4371d/SudoRules-V20.1/EXOTIC/Imp630), pinned at commit `7eefd71566b37e55c15f330d474140bfc4f4371d` (2026-06-20). The patterns are attributed to **eleven**. `catalog630_catalogue.json` records source paths and SHA256 hashes for every parsed file; original CLIPS source stays under `target/catalog630-rules`.

The author's [manual, sections15.2–15.4, printed pages274–279](https://github.com/denis-berthier/CSP-Rules-V2.1/blob/7eefd71566b37e55c15f330d474140bfc4f4371d/Docs/2023-UMRN-r1.pdf) explicitly describes **592 public geometries and38 excluded shapes**, rather than630 public geometry files. The original family sizes for10,12,13,14,15,16 cells are31,38,290,159,102,10. The omitted counts are0,0,5,7,24,2; their exact IDs are retained in the manifest. The manual deliberately excludes missing-candidate variants. Its update history dates the generated-rule publication to2023-03-16. Manual SHA256: `a9876d670c668bdfff3ee40714d5c09add75f96ca6c7579a6a6c04a14786a09a`.

The38 shapes have no free cell, according to that primary source. Under strict premises, every pattern cell contains all three palette candidates. A given of a palette digit therefore cannot share any pattern row, column or box. A shape with no free cell forbids givens of all three labels. This fixture has eight distinct given labels, so such a strict occurrence is impossible. More generally, any uniquely solvable standard Sudoku must have at least eight given labels: exchanging two labels absent from the givens gives a second solution. This exclusion applies under **native Sudoku geometry symmetries**. It does not transfer automatically to arbitrary abstract graph embeddings or sparse guarded variants.

Accordingly, all630 IDs are accounted for in the strict author model; **only592 geometries were actually retrieved and parsed**. This is not a claim to possess or parse the missing38 geometries. The earlier CSP-Rules-Examples lists contain380 shapes and are not a complete substitute. The original forum catalogue was inaccessible during this audit.

## Actual list and OR semantics

For each public rule, the parser checks the normal and transposed81-position geometry and all three candidate facts at every pattern cell. Extra candidates become guardians. Independently, the verifier exhaustively proves that each of the592 native pattern graphs is not three-colorable; hence at least one guardian must be true in every valid completion.

The exact scan uses all84 label triples and all injections preserving equality and inequality of rows, columns, bands and stacks, with optional transpose. These partial mappings extend to native Sudoku symmetries. There is no search budget or early cutoff. A map must satisfy all three actual candidate premises at each mapped cell, and the stored guardian list must include every extra candidate. This is stricter than finding an arbitrary non-three-colorable subgraph.

| Snapshot | Matched strict occurrences |
|---|---:|
| Selected20, raw opposite assumption |0|
| Selected20, closed opposite assumption |0|
| Full base |1|
| Full raw opposite assumption |1|
| Full closed opposite assumption |0|

The [EL12c38 source](https://github.com/denis-berthier/CSP-Rules-V2.1/blob/7eefd71566b37e55c15f330d474140bfc4f4371d/SudoRules-V20.1/EXOTIC/Imp630/12c/EL12c38.clp) maps source cells `[5,8,13,16,21,24,32,35,40,42,48,52]` respectively to board cells `[19,23,2,3,9,13,64,68,74,76,54,57]`, palette `{3,4,8}`, with the sole guardian `(68,7)`. `catalog630_match_result.json` retains the exact map in both snapshots. The independent coordinate checker validates it and rejects four corruptions.

The separate `novelty20_triplets.py/json` witnesses color every actual-list subproblem contained in each three-label palette. They exclude every **direct**, unguarded three-digit obstruction in the closed opposite snapshot, including sparse geometries absent from this public implementation. They do not exclude guarded OR relations, chains, conditional branches, or earlier snapshots. This catalogue scan likewise does not implement the author's OR-chain closure. The independently checked deeper-forcing certificate is the relevant known-method explanation of the requested target; the OR1 catalogue hit must not be presented as that explanation without additional proof.

That additional proof is now retained in [`deeper_forcing_catalogue_chain.json`](deeper_forcing_catalogue_chain.json): apply the checked EL12c38 OR1 placement on the base, assume r1c8=3, close locally, then split r9c5={4,8}. Both alternatives contradict native local constraints. The [expanded audit](novelty20_report.md) distinguishes the separately verified catalogue lemma from the three-node local part of this complete target explanation.

Reproduce offline acceptance and the complete rescan:

```sh
python3 research/global_obstructions/catalog630_verify.py --rescan
python3 research/global_obstructions/novelty20_triplets.py
```

Re-fetch pinned source and recheck all592 source hashes and candidate premises:

```sh
python3 research/global_obstructions/catalog630_sources.py
python3 research/global_obstructions/catalog630_verify.py --rescan --source-cache target/catalog630-rules
```

The verifier rejects Python `-O`. The published derived catalogue and exact maps are reviewable mathematical data, with attribution and source provenance; they do not vendor the full author implementation.
