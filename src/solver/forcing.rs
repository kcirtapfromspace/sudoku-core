//! Bounded native-rule forcing with a replay-only, Hall-based proof checker.
//! Proofs are conditional on exact candidate premises, not on player notes or
//! uniqueness. The checker never calls pivot search or the matching DP producer.
use super::candidate_state::CandidateState;
use super::explain::{Finding, InferenceResult, ProofCertificate};
use super::fabric::sector_cells;
use super::three_house::{self, HallObstruction, Layout};
use crate::Grid;
use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};

const ALL: u16 = 511;
const MAX_NODES: usize = 4096;
const MAX_STEPS: usize = 32768;
const MAX_DEPTH: usize = 8;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum ForcingPropagation {
    House,
    ThreeHouse,
}

#[derive(Debug, Clone)]
pub struct ForcingSearchOptions {
    /// Additional cell splits after the endpoint assumption (0..=8).
    pub max_depth: usize,
    /// Shared visits across every attempted pivot and candidate branch.
    pub max_nodes: usize,
    /// Shared propagation rule evaluations, including evaluations with no edits.
    pub max_steps: usize,
    pub propagation: ForcingPropagation,
}
impl Default for ForcingSearchOptions {
    fn default() -> Self {
        Self {
            max_depth: 1,
            max_nodes: 1024,
            max_steps: 16384,
            propagation: ForcingPropagation::House,
        }
    }
}
impl ForcingSearchOptions {
    pub(super) fn validate(&self) -> Result<(), String> {
        if self.max_depth > MAX_DEPTH || self.max_nodes > MAX_NODES || self.max_steps > MAX_STEPS {
            return Err("Forcing search limits exceed supported bounds".into());
        }
        Ok(())
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct DomainEdit {
    pub cell: usize,
    pub before: u16,
    pub after: u16,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub enum ForcingStep {
    House {
        sector: usize,
        edits: Vec<DomainEdit>,
    },
    ThreeHouse {
        center: usize,
        edits: Vec<DomainEdit>,
    },
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub enum ForcingLeaf {
    Hall {
        sector: usize,
        cells: Vec<usize>,
    },
    ThreeHouse {
        center: usize,
        blockers: Vec<HallObstruction>,
    },
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ForcingBranch {
    pub digit: u8,
    pub child: usize,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ForcingNode {
    pub steps: Vec<ForcingStep>,
    /// A node has exactly one of `leaf` and `cell`; splits cover every live digit.
    pub leaf: Option<ForcingLeaf>,
    pub cell: Option<usize>,
    pub branches: Vec<ForcingBranch>,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ForcingProof {
    pub version: u8,
    pub state_hash: String,
    pub cell: usize,
    pub digit: u8,
    /// The conclusion. The tree assumes the opposite endpoint.
    pub value: bool,
    /// Postorder tree, root last; every node must be used exactly once.
    pub nodes: Vec<ForcingNode>,
}
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct ForcingCheck {
    pub nodes: usize,
    pub splits: usize,
    pub hall_leaves: usize,
    pub three_house_leaves: usize,
    pub steps: usize,
    pub depth: usize,
}
#[derive(Debug, Clone)]
pub struct ForcingSearchResult {
    pub proof: Option<ForcingProof>,
    pub nodes_visited: usize,
    pub steps_evaluated: usize,
    pub budget_exhausted: bool,
    pub error: Option<String>,
}

pub(super) fn state_hash(state: &CandidateState) -> String {
    let mut h = Sha256::new();
    h.update(b"sudoku-forcing-v1-native");
    for (&v, &d) in state.values.iter().zip(&state.domains) {
        h.update([v]);
        h.update(d.to_le_bytes());
    }
    format!("{:x}", h.finalize())
}
fn assume(domains: &mut [u16], cell: usize, digit: u8, value: bool) {
    let bit = 1 << (digit - 1);
    domains[cell] &= if value { bit } else { ALL ^ bit };
}
fn edits(before: &[u16], after: &[u16]) -> Vec<DomainEdit> {
    (0..81)
        .filter(|&c| before[c] != after[c])
        .map(|cell| DomainEdit {
            cell,
            before: before[cell],
            after: after[cell],
        })
        .collect()
}

// Independent Hall theorem implementation: deficient sets refute a matching;
// tight sets reserve their digits and eliminate them from every outside cell.
fn hall_support(cells: &[usize], domains: &[u16], digits: u16) -> Option<Vec<u16>> {
    let mut support: Vec<_> = cells.iter().map(|&c| domains[c] & digits).collect();
    for subset in 1usize..(1 << cells.len()) {
        let union = cells.iter().enumerate().fold(0u16, |m, (i, &c)| {
            m | if subset & (1 << i) != 0 {
                domains[c] & digits
            } else {
                0
            }
        });
        let size = subset.count_ones();
        if union.count_ones() < size {
            return None;
        }
        if union.count_ones() == size {
            for (i, mask) in support.iter_mut().enumerate() {
                if subset & (1 << i) == 0 {
                    *mask &= !union;
                }
            }
        }
    }
    Some(support)
}
fn house_support(domains: &[u16], sector: usize) -> Option<Vec<u16>> {
    if sector >= 27 {
        return None;
    }
    let cells = sector_cells(sector);
    let mut next = domains.to_vec();
    for (&c, mask) in cells.iter().zip(hall_support(&cells, domains, ALL)?) {
        next[c] = mask;
    }
    Some(next)
}
fn checked_three_support(domains: &[u16], center: usize) -> Option<Vec<u16>> {
    let layout = Layout::new(center)?;
    let mut cache = vec![vec![None; 512]; 6];
    let mut result = domains.to_vec();
    for &c in &layout.cells {
        result[c] = 0;
    }
    let mut feasible = false;
    three_house::signatures(domains[center], |sets| {
        for atom in 0..6 {
            cache[atom][sets[atom] as usize]
                .get_or_insert_with(|| hall_support(&layout.atoms[atom], domains, sets[atom]));
        }
        if (0..6).any(|a| cache[a][sets[a] as usize].as_ref().unwrap().is_none()) {
            return true;
        }
        feasible = true;
        for a in 0..6 {
            for (&c, &m) in layout.atoms[a].iter().zip(
                cache[a][sets[a] as usize]
                    .as_ref()
                    .unwrap()
                    .as_ref()
                    .unwrap(),
            ) {
                result[c] |= m;
            }
        }
        true
    });
    feasible.then_some(result)
}
fn hall_leaf(domains: &[u16], sector: usize) -> Option<ForcingLeaf> {
    let cells = sector_cells(sector);
    for subset in 1usize..512 {
        let selected: Vec<_> = cells
            .iter()
            .enumerate()
            .filter(|(i, _)| subset & (1 << i) != 0)
            .map(|(_, &c)| c)
            .collect();
        let union = selected.iter().fold(0u16, |m, &c| m | domains[c]);
        if (union.count_ones() as usize) < selected.len() {
            return Some(ForcingLeaf::Hall {
                sector,
                cells: selected,
            });
        }
    }
    None
}
fn check_leaf(domains: &[u16], leaf: &ForcingLeaf) -> Option<()> {
    match leaf {
        ForcingLeaf::Hall { sector, cells } => {
            if *sector >= 27
                || cells.is_empty()
                || cells.len() > 9
                || cells
                    .iter()
                    .any(|&c| c >= 81 || !sector_cells(*sector).contains(&c))
            {
                return None;
            }
            let mut sorted = cells.clone();
            sorted.sort_unstable();
            sorted.dedup();
            if sorted.len() != cells.len() {
                return None;
            }
            let union = cells.iter().fold(0u16, |m, &c| m | domains[c]);
            ((union.count_ones() as usize) < cells.len()).then_some(())
        }
        ForcingLeaf::ThreeHouse { center, blockers } => {
            let layout = Layout::new(*center)?;
            if blockers.len() > 366 || domains[*center] == 0 {
                return None;
            }
            let mut blocked = [[false; 512]; 6];
            for b in blockers {
                let a = b.atom as usize;
                if !(1..6).contains(&a)
                    || b.digits & !ALL != 0
                    || b.digits.count_ones() as usize != layout.atoms[a].len()
                    || b.cells == 0
                    || b.cells as usize >= (1 << layout.atoms[a].len())
                    || blocked[a][b.digits as usize]
                {
                    return None;
                }
                let union = layout.atoms[a].iter().enumerate().fold(0u16, |m, (i, &c)| {
                    m | if b.cells & (1 << i) != 0 {
                        domains[c] & b.digits
                    } else {
                        0
                    }
                });
                if union.count_ones() >= b.cells.count_ones() {
                    return None;
                }
                blocked[a][b.digits as usize] = true;
            }
            three_house::signatures(domains[*center], |sets| {
                (1..6).any(|a| blocked[a][sets[a] as usize])
            })
            .then_some(())
        }
    }
}

impl ForcingProof {
    pub fn verify(&self, grid: &Grid) -> bool {
        self.check(grid).is_some()
    }
    pub fn check(&self, grid: &Grid) -> Option<ForcingCheck> {
        let state = CandidateState::capture(grid).ok()?;
        // Check local premise feasibility by Hall subsets, independently of the
        // producer's matching DP and the retained-state validation fast path.
        for sector in 0..27 {
            house_support(&state.domains, sector)?;
        }
        if self.version != 1
            || self.state_hash != state_hash(&state)
            || self.cell >= 81
            || !(1..=9).contains(&self.digit)
            || state.values[self.cell] != 0
            || state.domains[self.cell] & (1 << (self.digit - 1)) == 0
            || self.nodes.is_empty()
            || self.nodes.len() > MAX_NODES
        {
            return None;
        }
        if self.nodes.iter().map(|n| n.steps.len()).sum::<usize>() > MAX_STEPS {
            return None;
        }
        let mut domains = state.domains;
        assume(&mut domains, self.cell, self.digit, !self.value);
        let mut stats = ForcingCheck::default();
        let mut used = vec![false; self.nodes.len()];
        self.check_node(self.nodes.len() - 1, domains, 0, &mut used, &mut stats)?;
        used.iter().all(|&v| v).then_some(stats)
    }
    fn check_node(
        &self,
        id: usize,
        mut domains: Vec<u16>,
        depth: usize,
        used: &mut [bool],
        stats: &mut ForcingCheck,
    ) -> Option<()> {
        if id >= self.nodes.len() || used[id] || depth > MAX_DEPTH {
            return None;
        }
        used[id] = true;
        stats.nodes += 1;
        stats.depth = stats.depth.max(depth);
        let node = &self.nodes[id];
        for step in &node.steps {
            let (next, supplied) = match step {
                ForcingStep::House { sector, edits } => (house_support(&domains, *sector)?, edits),
                ForcingStep::ThreeHouse { center, edits } => {
                    (checked_three_support(&domains, *center)?, edits)
                }
            };
            if supplied.is_empty() || *supplied != edits(&domains, &next) {
                return None;
            }
            domains = next;
            stats.steps += 1;
        }
        if let Some(leaf) = &node.leaf {
            if node.cell.is_some() || !node.branches.is_empty() {
                return None;
            }
            check_leaf(&domains, leaf)?;
            match leaf {
                ForcingLeaf::Hall { .. } => stats.hall_leaves += 1,
                ForcingLeaf::ThreeHouse { .. } => stats.three_house_leaves += 1,
            }
        } else {
            let cell = node.cell?;
            if cell >= 81 || domains[cell].count_ones() < 2 || node.branches.len() > 9 {
                return None;
            }
            let mut covered = 0u16;
            for branch in &node.branches {
                if !(1..=9).contains(&branch.digit) || branch.child >= id {
                    return None;
                }
                let bit = 1 << (branch.digit - 1);
                if domains[cell] & bit == 0 || covered & bit != 0 {
                    return None;
                }
                covered |= bit;
                let mut next = domains.clone();
                next[cell] = bit;
                self.check_node(branch.child, next, depth + 1, used, stats)?;
            }
            if covered != domains[cell] {
                return None;
            }
            stats.splits += 1;
        }
        Some(())
    }
    /// Verify and apply atomically to the exact bound state.
    pub fn apply(&self, grid: &mut Grid) -> Result<(), String> {
        if !self.verify(grid) {
            return Err("Forcing certificate does not match or verify".into());
        }
        let mut next = grid.deep_clone();
        let pos = super::fabric::idx_to_pos(self.cell);
        if self.value {
            next.place_solver_value(pos, self.digit);
        } else {
            next.cell_mut(pos).remove_candidate(self.digit);
        }
        super::solver_state::validate(&next)?;
        *grid = next;
        Ok(())
    }
}

pub(super) struct SearchBudget {
    nodes: usize,
    steps: usize,
    max_nodes: usize,
    max_steps: usize,
    cutoff: bool,
    depth_limited: bool,
}
impl SearchBudget {
    pub(super) fn new(o: &ForcingSearchOptions) -> Self {
        Self {
            nodes: 0,
            steps: 0,
            max_nodes: o.max_nodes,
            max_steps: o.max_steps,
            cutoff: false,
            depth_limited: false,
        }
    }
    pub(super) fn exhausted(&self) -> bool {
        self.cutoff || self.nodes >= self.max_nodes || self.steps >= self.max_steps
    }
}
// Shared rule scheduler. Policies choose scopes, never a second technique list.
fn close(
    domains: &mut Vec<u16>,
    o: &ForcingSearchOptions,
    budget: &mut SearchBudget,
) -> Option<(Vec<ForcingStep>, Option<ForcingLeaf>)> {
    let mut steps = Vec::new();
    loop {
        let mut changed = false;
        for scope in 0..if o.propagation == ForcingPropagation::ThreeHouse {
            108
        } else {
            27
        } {
            if budget.steps == budget.max_steps {
                budget.cutoff = true;
                return None;
            }
            budget.steps += 1;
            let next = if scope < 27 {
                three_house::house_domain_support(domains, scope)
            } else {
                three_house::domain_support(domains, scope - 27)
            };
            let Some(next) = next else {
                let leaf = if scope < 27 {
                    hall_leaf(domains, scope)?
                } else {
                    ForcingLeaf::ThreeHouse {
                        center: scope - 27,
                        blockers: three_house::contradiction_cover(domains, scope - 27)?,
                    }
                };
                return Some((steps, Some(leaf)));
            };
            let changes = edits(domains, &next);
            if !changes.is_empty() {
                steps.push(if scope < 27 {
                    ForcingStep::House {
                        sector: scope,
                        edits: changes,
                    }
                } else {
                    ForcingStep::ThreeHouse {
                        center: scope - 27,
                        edits: changes,
                    }
                });
                *domains = next;
                changed = true;
            }
        }
        if !changed {
            return Some((steps, None));
        }
    }
}
fn refute(
    mut domains: Vec<u16>,
    depth: usize,
    o: &ForcingSearchOptions,
    budget: &mut SearchBudget,
    nodes: &mut Vec<ForcingNode>,
) -> Option<usize> {
    if budget.nodes == budget.max_nodes {
        budget.cutoff = true;
        return None;
    }
    budget.nodes += 1;
    let (steps, leaf) = close(&mut domains, o, budget)?;
    if let Some(leaf) = leaf {
        let id = nodes.len();
        nodes.push(ForcingNode {
            steps,
            leaf: Some(leaf),
            cell: None,
            branches: vec![],
        });
        return Some(id);
    }
    if depth == 0 {
        budget.depth_limited |= domains.iter().any(|d| d.count_ones() > 1);
        return None;
    }
    let mut pivots: Vec<_> = (0..81).filter(|&c| domains[c].count_ones() > 1).collect();
    pivots.sort_by_key(|&c| (domains[c].count_ones(), c));
    for cell in pivots {
        let start = nodes.len();
        let mut branches = Vec::new();
        let mut complete = true;
        for digit in 1..=9 {
            if domains[cell] & (1 << (digit - 1)) == 0 {
                continue;
            }
            let mut next = domains.clone();
            assume(&mut next, cell, digit, true);
            if let Some(child) = refute(next, depth - 1, o, budget, nodes) {
                branches.push(ForcingBranch { digit, child });
            } else {
                complete = false;
                break;
            }
        }
        if complete {
            let id = nodes.len();
            nodes.push(ForcingNode {
                steps,
                leaf: None,
                cell: Some(cell),
                branches,
            });
            return Some(id);
        }
        nodes.truncate(start);
        if budget.exhausted() {
            return None;
        }
    }
    None
}
fn produce(
    grid: &Grid,
    cell: usize,
    digit: u8,
    value: bool,
    o: &ForcingSearchOptions,
    budget: &mut SearchBudget,
) -> Option<ForcingProof> {
    let state = CandidateState::capture(grid).ok()?;
    if cell >= 81
        || !(1..=9).contains(&digit)
        || state.values[cell] != 0
        || state.domains[cell] & (1 << (digit - 1)) == 0
    {
        return None;
    }
    let hash = state_hash(&state);
    let mut domains = state.domains;
    assume(&mut domains, cell, digit, !value);
    let mut nodes = Vec::new();
    refute(domains, o.max_depth, o, budget, &mut nodes)?;
    let proof = ForcingProof {
        version: 1,
        state_hash: hash,
        cell,
        digit,
        value,
        nodes,
    };
    proof.verify(grid).then_some(proof)
}
pub(super) fn search(
    grid: &Grid,
    cell: usize,
    digit: u8,
    value: bool,
    o: &ForcingSearchOptions,
) -> ForcingSearchResult {
    let mut budget = SearchBudget::new(o);
    let error = o
        .validate()
        .and_then(|_| super::solver_state::validate(grid))
        .and_then(|_| {
            if cell >= 81
                || !(1..=9).contains(&digit)
                || grid.get(super::fabric::idx_to_pos(cell)).is_some()
                || !grid
                    .get_candidates(super::fabric::idx_to_pos(cell))
                    .contains(digit)
            {
                Err("Invalid forcing endpoint".into())
            } else {
                Ok(())
            }
        })
        .err();
    let proof = if error.is_none() {
        produce(grid, cell, digit, value, o, &mut budget)
    } else {
        None
    };
    let budget_exhausted = proof.is_none() && (budget.cutoff || budget.depth_limited);
    ForcingSearchResult {
        proof,
        nodes_visited: budget.nodes,
        steps_evaluated: budget.steps,
        budget_exhausted,
        error,
    }
}
pub(super) fn certify(
    grid: &Grid,
    mut f: Finding,
    o: &ForcingSearchOptions,
    budget: &mut SearchBudget,
) -> Option<Finding> {
    let (cell, digit, value) = match &f.inference {
        InferenceResult::Placement { cell, value } => (*cell, *value, true),
        InferenceResult::Elimination { cell, values } if values.len() == 1 => {
            (*cell, values[0], false)
        }
        _ => return None,
    };
    f.proof = Some(ProofCertificate::ForcingProof(produce(
        grid, cell, digit, value, o, budget,
    )?));
    Some(f)
}

/// Legacy endpoint discovery shares the same node/step budget as certification.
/// Exhaustion is an unresolved branch, never a contradiction.
pub(super) fn propagate_budgeted(
    grid: &Grid,
    pos: crate::Position,
    val: u8,
    policy: super::PropagationPolicy,
    budget: &mut SearchBudget,
) -> (Grid, bool) {
    let mut g = grid.deep_clone();
    if budget.exhausted() {
        budget.cutoff = true;
        return (g, false);
    }
    budget.nodes += 1;
    g.place_solver_value(pos, val);
    while !g.is_complete() {
        if super::backtrack::has_contradiction(&g) {
            return (g, true);
        }
        if budget.steps == budget.max_steps {
            budget.cutoff = true;
            break;
        }
        budget.steps += 1;
        match super::Solver.dispatch(&g, policy, &ForcingSearchOptions::default()) {
            Some(f) => super::apply_finding(&mut g, &f),
            None => break,
        }
    }
    let conflict = super::backtrack::has_contradiction(&g);
    (g, conflict)
}
