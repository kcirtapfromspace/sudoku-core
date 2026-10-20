//! Abstract Hall-cardinality core for three-house obstruction certificates.
//!
//! Distinct assigned digits in a Hall subset inject into its reconstructed
//! neighborhood. Consequently a smaller neighborhood refutes that assignment;
//! covering every possible boundary signature refutes every respecting target
//! assignment. This verifies the mathematical core, not the production Rust
//! checker, bit masks, native topology, hashing, serialization, or replay.
//! The interface obligation is that every respecting completion induces an
//! enumerated signature and an injective digit sequence in the certified union.
//! In the pinned vstd, `Set<T>` is a finite-set type (unlike `ISet<T>`), and
//! `Set::finite()` is deprecated and identically true. Thus `neighbors` and
//! `signatures` are finite by type; no additional finiteness premise is needed.
//!
//! Run with scripts/verify.sh and the pinned Verus release.

use vstd::prelude::*;
use vstd::set_lib::lemma_len_subset;

verus! {

/// Any injective assignment into a neighborhood needs at least as many digits
/// as assigned cells. `assigned` contains the digits of a Hall cell subset.
pub proof fn injective_assignment_needs_enough_neighbors(
    assigned: Seq<int>, neighbors: Set<int>,
)
    requires
        assigned.no_duplicates(),
        forall|i: int| 0 <= i < assigned.len()
            ==> #[trigger] neighbors.contains(assigned[i]),
    ensures
        assigned.len() <= neighbors.len(),
{
    assigned.unique_seq_to_set();
    assigned.to_set_ensures();
    assert forall|digit: int| assigned.to_set().contains(digit)
        implies #[trigger] neighbors.contains(digit) by {
        assert(assigned.contains(digit));
        let i = choose|i: int| 0 <= i < assigned.len() && assigned[i] == digit;
        assert(neighbors.contains(assigned[i]));
    }
    assert(assigned.to_set().subset_of(neighbors));
    lemma_len_subset(assigned.to_set(), neighbors);
}

/// A valid Hall deficiency rules out every injective assignment into its
/// neighborhood. No search algorithm or maximum-matching result is assumed.
pub proof fn hall_deficiency_refutes_assignment(
    assigned: Seq<int>, neighbors: Set<int>,
)
    requires
        neighbors.len() < assigned.len(),
        assigned.no_duplicates(),
        forall|i: int| 0 <= i < assigned.len()
            ==> #[trigger] neighbors.contains(assigned[i]),
    ensures
        false,
{
    injective_assignment_needs_enough_neighbors(assigned, neighbors);
}

/// Abstract complete-coverage theorem. Each enumerated signature has a chosen
/// blocker with its reconstructed union and cell count. A purported respecting
/// target assignment supplies one of those signatures and its actual assigned
/// digits in that blocker. The checked deficiency yields a contradiction.
pub proof fn complete_signature_coverage_refutes_assignment(
    signatures: Set<int>,
    neighborhood: spec_fn(int) -> Set<int>,
    blocker_size: spec_fn(int) -> nat,
    actual_signature: int,
    assigned: Seq<int>,
)
    requires
        signatures.contains(actual_signature),
        forall|signature: int| #[trigger] signatures.contains(signature)
            ==> neighborhood(signature).len() < blocker_size(signature),
        assigned.len() == blocker_size(actual_signature),
        assigned.no_duplicates(),
        forall|i: int| 0 <= i < assigned.len()
            ==> #[trigger] neighborhood(actual_signature).contains(assigned[i]),
    ensures
        false,
{
    assert(neighborhood(actual_signature).len() < blocker_size(actual_signature));
    hall_deficiency_refutes_assignment(assigned, neighborhood(actual_signature));
}

} // verus!
