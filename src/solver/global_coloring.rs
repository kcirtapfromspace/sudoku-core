//! Bounded native list-coloring obstructions with explicit refutation trees.
//!
//! Search chooses small cell sets; the checker reconstructs native peer edges
//! and checks every branch and empty-domain leaf. Neither search exhaustion nor
//! a heuristic failure is a proof. Candidate domains are conditional premises.
use super::candidate_state::CandidateState;
use super::explain::{ExplanationData, Finding, InferenceResult, ProofCertificate};
use super::fabric::idx_to_pos;
use super::{Hint, Technique};
use crate::Grid;
use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};

const VERSION: u8 = 1;
// Explicit supplied supports may be larger than the automatic hint search.
// Both paths retain fixed bounds; automatic hint cost remains unchanged.
const MAX_CELLS: usize = 32;
const MAX_PROOF_NODES: usize = 32768;
const SEARCH_MAX_CELLS: usize = 24;
const SEARCH_MAX_PROOF_NODES: usize = 8192;
const MAX_PATTERNS: usize = 16384;
const MAX_SEARCH_NODES: usize = 1_000_000;

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ColoringBranch {
    pub digit: u8,
    /// Postorder child index, strictly smaller than the parent index.
    pub child: usize,
}

/// An explicit case split or an independently reconstructed contradiction.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub enum ColoringNode {
    /// This selected cell has no digit left under the preceding assignments.
    Conflict { cell: usize },
    /// Exactly one child for every currently available digit of this cell.
    Split {
        cell: usize,
        branches: Vec<ColoringBranch>,
    },
}

/// Conditional proof of one placement or elimination using native peer edges.
///
/// Cells are sorted, unique board indices (at most 32), including the target.
/// Nodes form a postorder tree (at most 32,768), rooted at the last node. The
/// witness is a proper coloring of the unforced selected graph supporting the
/// RESULT endpoint. It rejects a contradictory selected graph, but does not
/// establish global or variant satisfiability or authenticate previous deletions.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct GlobalColoringProof {
    pub version: u8,
    pub state_hash: String,
    pub cells: Vec<usize>,
    pub cell: usize,
    pub digit: u8,
    /// True means placement; false means elimination. The tree refutes the opposite.
    pub value: bool,
    pub witness: Vec<u8>,
    pub nodes: Vec<ColoringNode>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GlobalColoringCheck {
    pub cells: usize,
    pub nodes: usize,
    pub splits: usize,
    pub conflicts: usize,
}

/// Bounds for guardian searches over three- and four-digit palettes.
#[derive(Debug, Clone)]
pub struct GlobalColoringSearchOptions {
    pub max_patterns: usize,
    /// Shared across all refutation and opposite-witness searches.
    pub max_nodes: usize,
}
impl Default for GlobalColoringSearchOptions {
    fn default() -> Self {
        Self {
            max_patterns: 512,
            max_nodes: 50_000,
        }
    }
}

#[derive(Debug, Clone)]
pub struct GlobalColoringSearchResult {
    pub hint: Option<Hint>,
    pub patterns_tested: usize,
    pub nodes_visited: usize,
    /// A global or per-pattern bound prevented a complete attempted search.
    pub budget_exhausted: bool,
    /// Components exceeding the certificate cell limit are not searched.
    pub oversized_components: usize,
    pub error: Option<String>,
}

fn peers(a: usize, b: usize) -> bool {
    a != b && (a / 9 == b / 9 || a % 9 == b % 9 || (a / 27 == b / 27 && a % 9 / 3 == b % 9 / 3))
}
fn state_hash(state: &CandidateState) -> String {
    let mut hash = Sha256::new();
    hash.update(b"sudoku-global-coloring-v1-native");
    for (&value, &domain) in state.values.iter().zip(&state.domains) {
        hash.update([value]);
        hash.update(domain.to_le_bytes());
    }
    format!("{:x}", hash.finalize())
}
fn restrict(domain: &mut u16, digit: u8, value: bool) {
    let bit = 1 << (digit - 1);
    *domain = if value { *domain & bit } else { *domain & !bit };
}

impl GlobalColoringProof {
    /// Search an explicit selected cell set. None includes bounded-search failure.
    /// This constructor neither checks LP feasibility nor establishes originality.
    /// Supports at most 32 cells and 32,768 proof nodes, with 50,000 total
    /// refutation/witness visits. Automatic hint searches keep smaller bounds.
    pub fn from_cells(
        grid: &Grid,
        mut cells: Vec<usize>,
        cell: usize,
        digit: u8,
        value: bool,
    ) -> Option<Self> {
        let state = CandidateState::capture(grid).ok()?;
        cells.sort_unstable();
        let mut budget = Budget::new(50_000);
        budget.pattern_limit = MAX_PROOF_NODES;
        build_proof(&state, &cells, cell, digit, value, &mut budget)
    }
    pub fn verify(&self, grid: &Grid) -> bool {
        self.check(grid).is_some()
    }
    pub fn check(&self, grid: &Grid) -> Option<GlobalColoringCheck> {
        self.check_state(&CandidateState::capture(grid).ok()?)
    }
    fn check_state(&self, state: &CandidateState) -> Option<GlobalColoringCheck> {
        if self.version != VERSION
            || self.state_hash != state_hash(state)
            || !valid_target(state, &self.cells, self.cell, self.digit)
            || self.witness.len() != self.cells.len()
            || self.nodes.is_empty()
            || self.nodes.len() > MAX_PROOF_NODES
        {
            return None;
        }
        let target = self.cells.binary_search(&self.cell).ok()?;
        for (i, (&cell, &digit)) in self.cells.iter().zip(&self.witness).enumerate() {
            if !(1..=9).contains(&digit) || state.domains[cell] & (1 << (digit - 1)) == 0 {
                return None;
            }
            for j in 0..i {
                if peers(cell, self.cells[j]) && digit == self.witness[j] {
                    return None;
                }
            }
        }
        if (self.witness[target] == self.digit) != self.value {
            return None;
        }
        let mut domains: Vec<_> = self.cells.iter().map(|&cell| state.domains[cell]).collect();
        restrict(&mut domains[target], self.digit, !self.value);
        let mut visited = vec![false; self.nodes.len()];
        let mut check = GlobalColoringCheck {
            cells: self.cells.len(),
            nodes: 0,
            splits: 0,
            conflicts: 0,
        };
        self.check_node(self.nodes.len() - 1, &domains, 0, &mut visited, &mut check)?;
        if !visited.into_iter().all(|used| used) {
            return None;
        }
        Some(check)
    }

    // Check the supplied tree; never choose a pivot or call the search DFS.
    fn check_node(
        &self,
        node: usize,
        domains: &[u16],
        assigned: u32,
        visited: &mut [bool],
        check: &mut GlobalColoringCheck,
    ) -> Option<()> {
        if node >= self.nodes.len() || visited[node] {
            return None;
        }
        visited[node] = true;
        check.nodes += 1;
        match &self.nodes[node] {
            ColoringNode::Conflict { cell } => {
                let index = self.cells.binary_search(cell).ok()?;
                if domains[index] != 0 {
                    return None;
                }
                check.conflicts += 1;
            }
            ColoringNode::Split { cell, branches } => {
                let index = self.cells.binary_search(cell).ok()?;
                if assigned & (1 << index) != 0
                    || domains[index] == 0
                    || branches.is_empty()
                    || branches.len() > 9
                {
                    return None;
                }
                let mut covered = 0u16;
                for branch in branches {
                    if !(1..=9).contains(&branch.digit) || branch.child >= node {
                        return None;
                    }
                    let bit = 1 << (branch.digit - 1);
                    if domains[index] & bit == 0 || covered & bit != 0 {
                        return None;
                    }
                    covered |= bit;
                    let mut next = domains.to_vec();
                    next[index] = bit;
                    // Direct native geometry reconstruction, independent of search adjacency.
                    for (j, &other) in self.cells.iter().enumerate() {
                        if peers(*cell, other) {
                            next[j] &= !bit;
                        }
                    }
                    self.check_node(branch.child, &next, assigned | (1 << index), visited, check)?;
                }
                if covered != domains[index] {
                    return None;
                }
                check.splits += 1;
            }
        }
        Some(())
    }
}

fn valid_target(state: &CandidateState, cells: &[usize], cell: usize, digit: u8) -> bool {
    !cells.is_empty()
        && cells.len() <= MAX_CELLS
        && cells.iter().all(|&c| c < 81)
        && cells.windows(2).all(|pair| pair[0] < pair[1])
        && cell < 81
        && (1..=9).contains(&digit)
        && state.values[cell] == 0
        && state.domains[cell] & (1 << (digit - 1)) != 0
        && cells.binary_search(&cell).is_ok()
}

struct Graph<'a> {
    cells: &'a [usize],
    edges: Vec<u32>,
}
impl<'a> Graph<'a> {
    fn new(cells: &'a [usize]) -> Self {
        let edges = cells
            .iter()
            .map(|&cell| {
                cells.iter().enumerate().fold(0, |mask, (i, &other)| {
                    mask | if peers(cell, other) { 1 << i } else { 0 }
                })
            })
            .collect();
        Self { cells, edges }
    }
    fn assign(&self, domains: &mut [u16], index: usize, bit: u16) {
        domains[index] = bit;
        let mut neighbors = self.edges[index];
        while neighbors != 0 {
            let j = neighbors.trailing_zeros() as usize;
            neighbors &= neighbors - 1;
            domains[j] &= !bit;
        }
    }
}
struct Budget {
    remaining: usize,
    visited: usize,
    exhausted: bool,
    pattern_limit: usize,
}
impl Budget {
    fn new(nodes: usize) -> Self {
        Self {
            remaining: nodes,
            visited: 0,
            exhausted: false,
            pattern_limit: SEARCH_MAX_PROOF_NODES,
        }
    }
}
enum Outcome {
    Refuted(usize),
    Supported(Vec<u8>),
    Cutoff,
}

// Search only: MRV/degree pivot selection and explicit forward-checking tree.
fn solve(
    graph: &Graph<'_>,
    domains: &[u16],
    assigned: u32,
    nodes: &mut Vec<ColoringNode>,
    budget: &mut Budget,
    calls: &mut usize,
) -> Outcome {
    if budget.remaining == 0 || *calls == budget.pattern_limit {
        budget.exhausted = true;
        return Outcome::Cutoff;
    }
    budget.remaining -= 1;
    budget.visited += 1;
    *calls += 1;
    if let Some(i) = domains.iter().position(|&domain| domain == 0) {
        let id = nodes.len();
        nodes.push(ColoringNode::Conflict {
            cell: graph.cells[i],
        });
        return Outcome::Refuted(id);
    }
    let pivot = (0..graph.cells.len())
        .filter(|&i| assigned & (1 << i) == 0)
        .min_by_key(|&i| {
            (
                domains[i].count_ones(),
                std::cmp::Reverse((graph.edges[i] & !assigned).count_ones()),
                i,
            )
        });
    let Some(index) = pivot else {
        return Outcome::Supported(
            domains
                .iter()
                .map(|m| m.trailing_zeros() as u8 + 1)
                .collect(),
        );
    };
    let mut choices = domains[index];
    let mut branches = Vec::new();
    while choices != 0 {
        let bit = 1 << choices.trailing_zeros();
        choices ^= bit;
        let mut next = domains.to_vec();
        graph.assign(&mut next, index, bit);
        match solve(graph, &next, assigned | (1 << index), nodes, budget, calls) {
            Outcome::Refuted(child) => branches.push(ColoringBranch {
                digit: bit.trailing_zeros() as u8 + 1,
                child,
            }),
            other => return other,
        }
    }
    let id = nodes.len();
    nodes.push(ColoringNode::Split {
        cell: graph.cells[index],
        branches,
    });
    Outcome::Refuted(id)
}

fn build_proof(
    state: &CandidateState,
    cells: &[usize],
    cell: usize,
    digit: u8,
    value: bool,
    budget: &mut Budget,
) -> Option<GlobalColoringProof> {
    if !valid_target(state, cells, cell, digit) {
        return None;
    }
    let target = cells.binary_search(&cell).ok()?;
    let graph = Graph::new(cells);
    let original: Vec<_> = cells.iter().map(|&c| state.domains[c]).collect();
    let mut assumed = original.clone();
    restrict(&mut assumed[target], digit, !value);
    let mut nodes = Vec::new();
    if !matches!(
        solve(&graph, &assumed, 0, &mut nodes, budget, &mut 0),
        Outcome::Refuted(_)
    ) {
        return None;
    }
    let mut opposite = original;
    restrict(&mut opposite[target], digit, value);
    let witness = match solve(&graph, &opposite, 0, &mut Vec::new(), budget, &mut 0) {
        Outcome::Supported(witness) => witness,
        _ => return None,
    };
    let proof = GlobalColoringProof {
        version: VERSION,
        state_hash: state_hash(state),
        cells: cells.to_vec(),
        cell,
        digit,
        value,
        witness,
        nodes,
    };
    proof.check_state(state)?;
    Some(proof)
}

fn finding(proof: GlobalColoringProof) -> Finding {
    let (cell, digit, value) = (proof.cell, proof.digit, proof.value);
    Finding {
        technique: Technique::GlobalColoring,
        inference: if value { InferenceResult::Placement { cell, value: digit } }
            else { InferenceResult::Elimination { cell, values: vec![digit] } },
        involved_cells: proof.cells.clone(),
        explanation: ExplanationData::Raw(format!(
            "{} {} at r{}c{}: the opposite assumption leaves a selected native peer graph with no candidate-respecting coloring. Every case is covered by the checked refutation tree.",
            if value { "Place" } else { "Remove" }, digit, cell / 9 + 1, cell % 9 + 1)),
        proof: Some(ProofCertificate::GlobalColoring(proof)),
    }
}

fn component(pool: &[usize], seed: usize) -> Vec<usize> {
    let mut seen = [false; 81];
    seen[seed] = true;
    let mut stack = vec![seed];
    while let Some(cell) = stack.pop() {
        for &other in pool {
            if !seen[other] && peers(cell, other) {
                seen[other] = true;
                stack.push(other);
            }
        }
    }
    pool.iter().copied().filter(|&cell| seen[cell]).collect()
}

pub(super) fn search(
    grid: &Grid,
    options: &GlobalColoringSearchOptions,
) -> GlobalColoringSearchResult {
    let mut result = GlobalColoringSearchResult {
        hint: None,
        patterns_tested: 0,
        nodes_visited: 0,
        budget_exhausted: false,
        oversized_components: 0,
        error: None,
    };
    if options.max_patterns > MAX_PATTERNS || options.max_nodes > MAX_SEARCH_NODES {
        result.error = Some("Global coloring limits exceed the supported bounds".into());
        return result;
    }
    let state = match CandidateState::capture(grid) {
        Ok(state) => state,
        Err(error) => {
            result.error = Some(error);
            return result;
        }
    };
    let mut budget = Budget::new(options.max_nodes);
    'palettes: for size in [3, 4] {
        for palette in 1u16..=0x1ff {
            if palette.count_ones() != size {
                continue;
            }
            let pure: Vec<_> = (0..81)
                .filter(|&c| state.values[c] == 0 && state.domains[c] & !palette == 0)
                .collect();
            if pure.len() < 3 {
                continue;
            }
            for guardian in 0..81 {
                if state.values[guardian] != 0 {
                    continue;
                }
                let domain = state.domains[guardian];
                let extra = domain & !palette;
                if extra.count_ones() != 1 || (domain & palette).count_ones() < 2 {
                    continue;
                }
                let mut pool = pure.clone();
                pool.push(guardian);
                pool.sort_unstable();
                let cells = component(&pool, guardian);
                if cells.len() < 4 {
                    continue;
                }
                if cells.len() > SEARCH_MAX_CELLS {
                    result.oversized_components += 1;
                    continue;
                }
                if result.patterns_tested == options.max_patterns || budget.remaining == 0 {
                    budget.exhausted = true;
                    break 'palettes;
                }
                result.patterns_tested += 1;
                let digit = extra.trailing_zeros() as u8 + 1;
                if let Some(proof) = build_proof(&state, &cells, guardian, digit, true, &mut budget)
                {
                    result.hint = Some(finding(proof).to_hint());
                    break 'palettes;
                }
            }
        }
    }
    result.nodes_visited = budget.visited;
    result.budget_exhausted = budget.exhausted;
    result
}

pub(super) fn find(grid: &Grid) -> Option<Finding> {
    search(grid, &GlobalColoringSearchOptions::default())
        .hint
        .and_then(|hint| match hint.proof? {
            ProofCertificate::GlobalColoring(proof) => Some(finding(proof)),
            _ => None,
        })
}

/// Portable one-step native projection; certificates are conditional on domains.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct GlobalColoringReplay {
    pub version: u8,
    pub state: CandidateState,
    pub proof: GlobalColoringProof,
}
impl GlobalColoringReplay {
    pub fn capture(grid: &Grid, proof: &GlobalColoringProof) -> Result<Self, String> {
        let replay = Self {
            version: VERSION,
            state: CandidateState::capture(grid)?,
            proof: proof.clone(),
        };
        replay.check()?;
        Ok(replay)
    }
    pub fn check(&self) -> Result<GlobalColoringCheck, String> {
        if self.version != VERSION {
            return Err("Unsupported global coloring replay version".into());
        }
        self.state.to_grid()?;
        self.proof
            .check_state(&self.state)
            .ok_or_else(|| "Global coloring certificate does not verify".into())
    }
    pub fn verify(&self) -> bool {
        self.check().is_ok()
    }
    /// Apply atomically to identical native premises; never rebuild candidates.
    pub fn apply(&self, grid: &mut Grid) -> Result<(), String> {
        self.check()?;
        if CandidateState::capture(grid)? != self.state {
            return Err("Global coloring replay does not match current state".into());
        }
        let mut next = grid.deep_clone();
        let pos = idx_to_pos(self.proof.cell);
        if self.proof.value {
            next.set_cell_unchecked(pos, Some(self.proof.digit));
            next.update_candidates_after_move(pos, self.proof.digit);
        } else {
            next.cell_mut(pos).remove_candidate(self.proof.digit);
        }
        if !next.validate().is_valid {
            return Err("Global coloring replay would violate an installed constraint".into());
        }
        CandidateState::capture(&next)?;
        *grid = next;
        Ok(())
    }
}

#[cfg(test)]
#[path = "global_coloring_tests.rs"]
mod tests;
