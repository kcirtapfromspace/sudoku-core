//! Abstract exactly-one residue theorem for docs/partition-residues.md.
//!
//! Every selection of one coefficient from each group reaches its sum modulo q.
//! Therefore absence of the required residue refutes every such selection.
//! This verifies the recurrence's soundness invariant, not the production u16
//! bitset operations or the reconstruction of Sudoku equations/groups. Those
//! interfaces remain trusted and are exercised by independent tests.
//!
//! Run with scripts/verify.sh, using the pinned CI Verus release.

use vstd::arithmetic::div_mod::{lemma_add_mod_noop_right, lemma_small_mod};
use vstd::prelude::*;

verus! {

pub open spec fn selected_sum(selection: Seq<int>, n: int) -> int
    decreases n,
{
    if n <= 0 { 0 } else { selected_sum(selection, n - 1) + selection[n - 1] }
}

/// The recurrence may overapproximate attainable sums, but must include zero
/// initially and every transition formed by one next-group choice.
pub open spec fn residue_fold(domains: Seq<Set<int>>, states: Seq<Set<int>>, q: int) -> bool {
    &&& q > 0
    &&& states.len() == domains.len() + 1
    &&& states[0].contains(0)
    &&& forall|i: int, prior: int, coefficient: int|
        #![trigger states[i].contains(prior), domains[i].contains(coefficient)]
        0 <= i < domains.len()
        && states[i].contains(prior)
        && domains[i].contains(coefficient)
        ==> states[i + 1].contains((prior + coefficient) % q)
}

pub proof fn selection_reaches_its_sum(
    domains: Seq<Set<int>>, states: Seq<Set<int>>, selection: Seq<int>, q: int, n: int,
)
    requires
        residue_fold(domains, states, q),
        0 <= n <= domains.len(),
        selection.len() == domains.len(),
        forall|i: int| 0 <= i < n ==> #[trigger] domains[i].contains(selection[i]),
    ensures
        states[n].contains(selected_sum(selection, n) % q),
    decreases n,
{
    if n == 0 {
        lemma_small_mod(0, q as nat);
    } else {
        selection_reaches_its_sum(domains, states, selection, q, n - 1);
        let previous_sum = selected_sum(selection, n - 1);
        let coefficient = selection[n - 1];
        lemma_add_mod_noop_right(coefficient, previous_sum, q);
        assert(domains[n - 1].contains(coefficient));
        assert(states[n - 1].contains(previous_sum % q));
        assert(states[n].contains((previous_sum % q + coefficient) % q));
    }
}

/// Applied to target-restricted domains, this refutes the Boolean assumption.
pub proof fn partition_residue_soundness(
    domains: Seq<Set<int>>, states: Seq<Set<int>>, selection: Seq<int>, q: int, rhs: int,
)
    requires
        residue_fold(domains, states, q),
        selection.len() == domains.len(),
        forall|i: int| 0 <= i < domains.len() ==> #[trigger] domains[i].contains(selection[i]),
        !states[domains.len() as int].contains(rhs % q),
    ensures
        selected_sum(selection, selection.len() as int) != rhs,
{
    selection_reaches_its_sum(domains, states, selection, q, domains.len() as int);
}

} // verus!
