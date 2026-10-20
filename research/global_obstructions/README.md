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

The current-engine probe checks proof construction, automatic search, JSON
round trips, replay application, stale rejection and input preservation for all
six givens-derived snapshots.

Optional discovery reruns and source-fetch commands are documented in the
structural and real-puzzle reports. Generated clues, primary-source provenance,
unique versus multiple completions, and prior-method-solvable controls are
retained separately. The general coloring and stable-set principles are known;
the historical status of this particular four-digit Sudoku motif remains
unverified.
