//! Abstract soundness core of exhaustive native list-coloring refutation.
//!
//! A respecting proper coloring survives its own branch's pinning and peer
//! exclusions. Empty domains reject such a coloring. A split that covers the
//! current domain and refutes every child therefore refutes the parent.
//! Applying the split theorem inductively gives finite refutation-tree
//! soundness. This does NOT verify the production flat-tree checker, postorder
//! indexing, bit masks, Sudoku geometry, hashing, serde, search, or replay.
//! The trusted interface must identify actual native edges, correct effective
//! domains, the opposite target assumption, and exhaustive valid child proofs.
//! Pinned vstd `Set<T>` is finite by type; `Set::finite()` is deprecated/true.

use vstd::prelude::*;

verus! {

pub open spec fn respecting_coloring(
    domains: Seq<Set<int>>, assigned: Seq<int>, edge: spec_fn(int, int) -> bool,
) -> bool {
    &&& assigned.len() == domains.len()
    &&& forall|i: int| 0 <= i < domains.len()
        ==> #[trigger] domains[i].contains(assigned[i])
    &&& forall|i: int, j: int| 0 <= i < domains.len() && 0 <= j < domains.len()
        && #[trigger] edge(i, j) ==> assigned[i] != assigned[j]
}

/// Pin the pivot and remove its digit from exactly its recorded neighbors.
pub open spec fn branch_domains(
    domains: Seq<Set<int>>, edge: spec_fn(int, int) -> bool, pivot: int, digit: int,
) -> Seq<Set<int>> {
    Seq::new(domains.len(), |j: int|
        if j == pivot { Set::<int>::empty().insert(digit) }
        else if edge(pivot, j) { domains[j].remove(digit) }
        else { domains[j] })
}

/// The actual digit of any respecting coloring is never lost by forward checking.
pub proof fn own_branch_preserves_coloring(
    domains: Seq<Set<int>>, assigned: Seq<int>, edge: spec_fn(int, int) -> bool,
    pivot: int, digit: int,
)
    requires
        respecting_coloring(domains, assigned, edge),
        0 <= pivot < domains.len(),
        assigned[pivot] == digit,
    ensures
        respecting_coloring(branch_domains(domains, edge, pivot, digit), assigned, edge),
{
    let next = branch_domains(domains, edge, pivot, digit);
    assert forall|j: int| 0 <= j < domains.len()
        implies #[trigger] next[j].contains(assigned[j]) by {
        assert(domains[j].contains(assigned[j]));
        if j == pivot {
            assert(assigned[j] == digit);
        } else if edge(pivot, j) {
            assert(assigned[pivot] != assigned[j]);
            assert(assigned[j] != digit);
        }
    }
}

/// The only production leaf terminal is a reconstructed empty candidate domain.
pub proof fn empty_domain_refutes_coloring(
    domains: Seq<Set<int>>, assigned: Seq<int>, edge: spec_fn(int, int) -> bool,
    empty_cell: int,
)
    requires
        0 <= empty_cell < domains.len(),
        domains[empty_cell] == Set::<int>::empty(),
    ensures
        !respecting_coloring(domains, assigned, edge),
{
    if respecting_coloring(domains, assigned, edge) {
        assert(domains[empty_cell].contains(assigned[empty_cell]));
    }
}

/// Induction step: exact branch coverage plus sound child refutations excludes
/// every purported respecting coloring at this parent domain state.
pub proof fn exhaustive_split_refutes_coloring(
    domains: Seq<Set<int>>, assigned: Seq<int>, edge: spec_fn(int, int) -> bool,
    pivot: int, branches: Set<int>,
)
    requires
        0 <= pivot < domains.len(),
        branches == domains[pivot],
        forall|digit: int| #[trigger] branches.contains(digit)
            ==> !respecting_coloring(branch_domains(domains, edge, pivot, digit), assigned, edge),
    ensures
        !respecting_coloring(domains, assigned, edge),
{
    if respecting_coloring(domains, assigned, edge) {
        let digit = assigned[pivot];
        assert(branches.contains(digit));
        own_branch_preserves_coloring(domains, assigned, edge, pivot, digit);
    }
}

} // verus!
