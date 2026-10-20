# Global Sudoku obstructions

This experiment searches beyond complete single-house Hall propagation, the
engine's 81 crossing three-house support checks, and the full continuous native
324-equation relaxation. All three filters are checked after the opposite target
assumption and local closure. All 2,925 arbitrary selected-house triples are
additionally checked for feasibility, a separate weaker obligation than complete
support filtering for those arbitrary triples.

Verified deductions are integrated as Global Coloring, with independently
checked finite list-coloring refutation trees. The [production theorem and API](../../docs/global-coloring.md)
explain the exact state premises, resource bounds, proof grammar and replay.
No historical mathematical novelty is established.

| Evidence | Meaning |
| --- | --- |
| [13-cell four-digit support](structural_report.md) | A native graph with independence number 3 and chromatic number 5, exact local/LP separation, and a short four-box rank proof. Synthetic original state. |
| [Published tridagon cases](real_notes.md) | Independently verified unique givens from immutable author sources, exact traces, local filters and rational LP witnesses. Established patterns. |
| Generated four-digit givens | A nonunique 19-given board reaches the obstruction beyond the checked old structural baseline. A unique 34-given control is already solved by old fish/singles. Neither is a published author-source puzzle. |
| [Expanded unique nonforcing gap](expanded_search_report.md) | A certified 29-given unique board passes the local/native-LP filters and stalls older nonforcing plus default arithmetic. Full old Nishio still finds a step. |
| [Full older logical engine stall](nonrank_report.md) | A unique published 28-given puzzle stalls the complete frozen logical dispatcher including forcing; current Global Coloring certifies r4c9=9 with 169 nodes. Established tridagon composition. |
| [Beyond per-digit weighted rank](nonrank_audit_report.md) | Exact 33-term full-board digit-template mixture survives every single-digit-valid linear inequality. Earlier 13-cell hyperhole is rejected by an exact template Farkas control. |
| [Broader graph search](nonrank_graph_report.md) | Synthetic 27-cell five-chromatic graph with an exact fractional four-color cover and full-board template witness; unique-givens realization remains unverified. Includes exact rejected controls and bounded negative results. |
| [Exact motif prior art](prior_art_13_report.md) | An explicit spanning seven-hyperhole and grouped-loop map explain the 13-cell obstruction by established methods. First appearance of the exact Sudoku drawing remains unresolved. |
| [LP and prior-art audit](audit_prior_art.md) | Exact rational primal witnesses, rejected Farkas controls, certificate mutations and precise limits on method and novelty claims. |

The positive continuous-LP result is an exact rational witness, not a numerical
solver status. It does not imply survival of every stronger LP encoding,
integer cut, AIC or forcing procedure. The rejected 21-cell four-digit candidate
has a checked 35-term linear refutation and was excluded from the positive
claims. Earlier three-house trial 847 has a checked 11-term linear refutation.

## Exact checks

Python 3 standard library is sufficient for all acceptance checks:

```sh
python3 research/global_obstructions/structural_verify.py
python3 research/global_obstructions/real_verify.py
python3 research/global_obstructions/audit_linear_witnesses.py
python3 research/global_obstructions/audit_global_certificates.py
python3 research/global_obstructions/expanded_unique_verify.py
python3 research/global_obstructions/prior_art_13_verify.py
python3 -S research/global_obstructions/nonrank_verify.py
python3 -S research/global_obstructions/nonrank_puzzles_verify.py
python3 -S research/global_obstructions/nonrank_audit_template.py --self-test
python3 -S research/global_obstructions/nonrank_graph_verify.py
python3 -S research/global_obstructions/nonrank_graph_certificate.py --self-test
```

These run in dedicated CI. Numerical LP discovery can optionally use SciPy;
acceptance uses `fractions.Fraction`. Scripts reject Python `-O` where assertions
are part of validation. Committed proofs and witnesses use only native rules.
Rust tests additionally check actual API searches, direct proof construction,
serialization, wrong-state rejection, both result polarities and atomic replay.

The manually invoked comparisons instrument temporary source copies:

```sh
python3 research/global_obstructions/structural_compare_engine.py
python3 research/global_obstructions/real_existing_check.py
python3 research/global_obstructions/real_engine_check.py
```

Their precisely enumerated older detectors exclude the new feature; recorded
search limits and exclusions accompany results. An older method finding a future
step is a research result, not grounds to suppress it. A failing comparison is
not automatically a Sudoku correctness regression.

The Rust integration suite checks proof construction, automatic search, JSON
round trips, replay application, stale rejection and input preservation for all
seven givens-derived snapshots. A separate CI job freshly rebuilds the frozen
older engine and reconfirms the full logical fixed point twice.

Optional discovery reruns and source-fetch commands are documented in the
structural and real-puzzle reports. Generated clues, primary-source provenance,
unique versus multiple completions, and prior-method-solvable controls are
retained separately. The general coloring and stable-set principles are known;
the historical status of this particular four-digit Sudoku motif remains
unverified.
