# A unique real state beyond the full frozen older engine

The bounded corpus experiment completed six author-hosted Monhard T&E3 puzzles
from the pinned `expanded_modern_corpus.json`, in its recorded hard-first order.
Every puzzle was independently unique. Five reached a live fixed point of the
full older dispatcher, including forcing; one was solved. The search stopped
when the first full-stall state yielded the separately accepted list-coloring
deduction. This is not a result for all 200 corpus entries.

`nonrank_puzzles_trials.jsonl` publicly retains all six clue strings, identifiers,
statuses, old finding counts, remaining unplaced cells and domain hashes.
`nonrank_puzzles_results.json` binds its hash, schema and counts. Full traces for
the six outcomes remain under `target/nonrank-puzzles-modern200/`; the accepted
first premise trace is also retained in `nonrank_full_old_stall.json`.

For each case the experiment independently counts completions, starts from
ordinary givens-derived candidates, and alternates exact Hall support and all 81
crossing-three-house checks with every finding from the complete older engine.
All candidate eliminations survive placements. Every older endpoint is checked
against the sole completion and the entire domain trace is replayed. A stall
requires a live board and two successful complete old-dispatcher calls that
return no finding and leave the state unchanged. Timeouts, invalid deductions,
search limits and solved boards have separate statuses.

The first accepted case is Monhard T&E3 source line 1. Twelve older forcing
findings precede the final state, which still has 50 unplaced cells. A fresh
build from archived commit `0532d63d9241a69c1fb57e9c34a1533effe1981d` confirms the
same live fixed point twice. The archive is checked byte-for-byte against the
Git source and contains no GlobalColoring module. The saved successful calls,
binary/source hashes and stable promoted-fixture hash are in
`nonrank_puzzles_first_audit.json`. The older Arithmetic search reports beam
pruning and no hint; its configured bounded behavior is retained explicitly.
No claim of exhaustive separation from every possible arithmetic or forcing
proof is made.

## What the obstruction is

The engine's direct certificate uses sixteen cells and the palette {4,6,7,8} to
place r4c9=9. An independent selected-graph reduction leaves a fourteen-cell
raw-assumption list obstruction. Its underlying graph **is four-colorable**;
the actual candidate lists cause the contradiction. This explains why the older
uniform non-four-colorability corpus filter missed it.

After the six sound local steps on the r4c9≠9 branch, a twelve-cell core has
uniform domains {4,7,8}. The saved bijection in
`nonrank_puzzles_first_classification.json` preserves every native peer edge and
maps this graph exactly to the published Loki tridagon control. The established
tridagon/trivalue-oddagon interpretation is therefore confirmed rather than
merely inferred from its size. Denis Berthier's
[primary Loki solving trace](https://raw.githubusercontent.com/denis-berthier/CSP-Rules-Examples/fd9207f8dc1eb6e2e655048ba3b327d423900e1b/Sudoku-a/Tridagons/Loki-11.9.clp)
identifies the twelve-cell pattern and its contradiction. This is a verified
new engine capability and relaxation-separation example; historical
mathematical novelty is not claimed.

The separate root acceptance retains the exact full single-digit-template LP
witness and production graph certificate. Those files establish the requested
survival of per-digit weighted-rank relaxation; the present files independently
audit the older-engine comparison, givens premises and known-pattern identity.

## Reproduce

```sh
# CI-safe stdlib verification: no compiler, historical engine or network.
python3 -S research/global_obstructions/nonrank_puzzles_verify.py

# Manual fresh historical build and two complete fixedpoint calls.
# Use a new target folder; this command deliberately does not reuse one.
python3 -S research/global_obstructions/nonrank_puzzles_audit.py \
  --fresh-folder target/nonrank-rebuilt-0532d63

# Optional reuse of a manually fresh-built historical binary and archive.
python3 -S research/global_obstructions/nonrank_puzzles_audit.py \
  --binary target/nonrank-rebuilt-0532d63/instrumented/target/release/examples/real_existing_probe \
  --archive-source target/nonrank-rebuilt-0532d63/baseline/src/solver/mod.rs

# Bounded corpus closure discovery; one worker, numerical-library threads=1.
python3 research/global_obstructions/nonrank_puzzles_search.py \
  --limit 6 --order hard-first --seconds 960 \
  --output-dir target/nonrank-puzzles-modern200

# Re-extract compact records and the graph-isomorphism classification.
python3 research/global_obstructions/nonrank_puzzles_report.py
```

The stdlib check verifies the retained historical calls and their exact bound
premises. Reexecuting the historical engine is intentionally a separate manual
audit, since fresh compilation and forcing execution do not belong in the
stdlib certificate checker. The fresh audit itself needs only standard-library
Python plus Rust/Cargo; NumPy and SciPy are not imported on this path.
