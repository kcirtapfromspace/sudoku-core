# Exact deeper forcing audit of Mith515

The deduction **r1c8≠3 has two checked established explanations**: a two-level
local forcing tree, and a published catalogue rule followed by one local split.
This audit supplies a small exact certificate; it supplies no evidence of a new
mathematical obstruction family. No tridagon, odd-cycle, uniqueness, SAT answer,
or optimizer answer is used in the refutation.

Start with the independently audited candidate state in
[`different_lead.json`](different_lead.json), assume r1c8=3, and replay its ten
recorded Hall/crossing-house propagation events. At the resulting fixed point:

| First split | Second split after propagation | Exact contradiction |
| --- | --- | --- |
| r6c2=6 | r2c1=3 | No compatible assignment to row8, column1, box7 |
| r6c2=6 | r2c1=8 | No compatible assignment to row8, column1, box7 |
| r6c2=9 | r9c5=4 | r2c1 and r2c5 both require digit3 in row2 |
| r6c2=9 | r9c5=8 | No compatible assignment to row1, column4, box2 |

The first split covers exactly r6c2={6,9}. Propagation after its value6 leaves
r2c1={3,8}; propagation after value9 leaves r9c5={4,8}. Thus every branch value
is covered, proving the initial target assumption impossible.

The [15,417-byte certificate](deeper_forcing_certificate.json) has seven tree
nodes, three complete splits and four contradiction leaves. It records 74 exact
propagation events after the already closed opposite state. One leaf contains
an explicit Hall deficiency. Three leaves contain seven Hall-subset blockers
covering all 1,260 admissible boundary signatures: 420 per leaf, after excluding
center digits already absent from that branch's domain. The checker recomputes
each propagation edit and every blocker union and signature; it never accepts
an exception message as a contradiction proof.

For a precise lower bound, every one of the 32 live pivot cells in the closed
opposite state has a recorded candidate whose deterministic local closure
survives. The checker independently closes each of those 32 branches. Therefore
no one-level complete cell-value case split can refute this state using only
this local propagator. The minimum additional split depth in this specified
calculus is exactly two.

Counting the initial target assumption, each longest branch contains three
nested value hypotheses from the base state. This is not an unqualified claim
of classic T&E(2).

Here the propagator is perfect-matching GAC in each native house, alternating
with all 81 crossing-three-house GAC checks until no candidate changes. The
calculus permits complete cell-value splits and these local contradiction
leaves. This lower bound does not allow additional named lemmas, failed-literal
learning between sibling branches, arbitrary rule applications, or additional
proof shortcuts. It is not a minimum classic Sudoku T&E rating and is not a
minimum over all known mathematical methods.

The seven-node refutation was found in 0.203 seconds, with no cutoff. The
exhaustive depth-one local search took 2.710 seconds. Both used one worker and
Python's standard library. The original 85-node production list-coloring tree
and this tree are different certificates of the same target inference; their
node counts use different propagation rules and should not be compared as a
common proof-complexity measure.

The successful explanation is consistent with the older engine's live stall.
In frozen commit `0532d63`, the forcing helper `propagate_full` omits arithmetic,
crossing-three-house propagation and recursive forcing inside assumptions, and
recalculates candidates after pinning a value. The new proof instead preserves
its candidate premises, applies Hall/all81 propagation under every assumption,
and uses two nested split levels after the target assumption. These are
different configured inference procedures. This audit establishes no engine
bug or unsound historical result.

Verify the exact proof and reject 42 altered targets, states, branches, edits,
leaves, lower-bound witnesses and claimed counts:

```sh
python3 -S research/global_obstructions/deeper_forcing_verify.py --mutations
```

Reproduce discovery independently of the retained certificate:

```sh
python3 -S research/global_obstructions/deeper_forcing_search.py \
  research/global_obstructions/different_lead.json \
  --depth 2 --seconds 200 --output target/deeper-forcing-local-depth2.json
python3 -S research/global_obstructions/deeper_forcing_search.py \
  research/global_obstructions/different_lead.json \
  --depth 1 --seconds 200 --output target/deeper-forcing-local-depth1.json
```

Discovery time/node bounds return an explicit cutoff status; they never produce
a no-refutation or minimum-depth claim. The retained lower-bound witnesses and
all positive proof leaves are checked without such search-oracle acceptance.
The candidate-state premise's derivation from the unique puzzle's givens is
audited separately by the existing `different_verify.py` and historical-engine
audit; the new certificate binds that exact public fixture and both its base
and closed opposite domain states.

## Published catalogue rule followed by the actual target proof

The complete-base catalogue audit additionally finds the published rule
**EL12c38** on cells
`[2,3,9,13,19,23,54,57,64,68,74,76]`, using palette {3,4,8}. These are zero-based
cell indices. The only candidate outside that palette on this non-three-colorable
support is **r8c6=7**. The known guardian OR therefore has one literal, proving
that placement unconditionally from the base candidate state.

This placement alone does not eliminate r1c8=3. Six Hall/all81 events after the
placement leave that target candidate available. Assume r1c8=3 and replay eleven
more events. Now split r9c5={4,8}: value4 gives the row2 Hall contradiction and
value8 gives the row1/column4/box2 contradiction described above.

The [chained certificate](deeper_forcing_catalogue_chain.json) records **three
local forcing nodes plus one separately verified published catalogue lemma**.
Its 26 branch-propagation events end at one explicit Hall deficiency and one
Hall blocker covering all 420 admissible triple boundary signatures. The checker
validates the source rule's native map, graph non-three-colorability, actual
palette candidates and exact sole guardian. It also verifies that applying the
known rule to the original closed opposite, then closing seven events, reaches
the identical state. Thus it explains the requested target, rather than merely
finding a different hint.

The local-only depth-two minimum above remains valid for its specified calculus.
The catalogue lemma is an additional allowed inference, so this shorter chain
does not contradict that lower bound. Neither certificate provides a new
obstruction-family or classic T&E-minimum claim.

```sh
python3 -S research/global_obstructions/deeper_forcing_catalogue_chain.py --mutations
```

The chained verifier rejects fourteen altered rule, guardian, target, mapping,
propagation, branch and contradiction certificates. The catalogue's broader
source coverage is checked separately by `catalog630_verify.py`.
