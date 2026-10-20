//! Algorithm X proof producer. No named detector or production solver is called.
use super::*;

fn columns(candidate: usize) -> [usize; 4] {
    let cell = candidate / 9;
    let digit = candidate % 9;
    [
        cell,
        81 + cell / 9 * 9 + digit,
        162 + cell % 9 * 9 + digit,
        243 + (cell / 27 * 3 + cell % 9 / 3) * 9 + digit,
    ]
}

struct Search {
    rows: Vec<Vec<usize>>,
    covered: [bool; 324],
    assignments: Vec<u8>,
    nodes: Vec<ExactNode>,
    solutions: Vec<Vec<u8>>,
    visits: usize,
    limit: usize,
    cutoff: bool,
    stop_after: usize,
}
impl Search {
    fn new(state: &CandidateState, domains: &[u16], limit: usize, stop_after: usize) -> Self {
        let mut s = Self {
            rows: vec![vec![]; 324],
            covered: [false; 324],
            assignments: state.values.clone(),
            nodes: vec![],
            solutions: vec![],
            visits: 0,
            limit,
            cutoff: false,
            stop_after,
        };
        for (cell, &domain) in domains.iter().enumerate() {
            if state.values[cell] != 0 {
                for c in columns(cell * 9 + state.values[cell] as usize - 1) {
                    s.covered[c] = true;
                }
            } else {
                for digit in 0..9 {
                    if domain & (1 << digit) == 0 {
                        continue;
                    }
                    let candidate = cell * 9 + digit;
                    for c in columns(candidate) {
                        s.rows[c].push(candidate);
                    }
                }
            }
        }
        s
    }
    fn walk(&mut self) -> Option<usize> {
        if self.visits >= self.limit {
            self.cutoff = true;
            return None;
        }
        self.visits += 1;
        let mut chosen = None;
        let mut candidates = Vec::new();
        for column in 0..324 {
            if self.covered[column] {
                continue;
            }
            let possible: Vec<_> = self.rows[column]
                .iter()
                .copied()
                .filter(|&r| columns(r).iter().all(|&c| !self.covered[c]))
                .collect();
            if chosen.is_none() || possible.len() < candidates.len() {
                chosen = Some(column);
                candidates = possible;
                if candidates.is_empty() {
                    break;
                }
            }
        }
        let node = if let Some(requirement) = chosen {
            if candidates.is_empty() {
                ExactNode::Conflict {
                    requirement: requirement as u16,
                }
            } else {
                let mut branches = vec![];
                for candidate in candidates {
                    let cols = columns(candidate);
                    for c in cols {
                        self.covered[c] = true;
                    }
                    self.assignments[candidate / 9] = (candidate % 9 + 1) as u8;
                    let child = self.walk()?;
                    for c in cols {
                        self.covered[c] = false;
                    }
                    self.assignments[candidate / 9] = 0;
                    branches.push(ExactBranch {
                        candidate: candidate as u16,
                        child,
                    });
                    if self.solutions.len() >= self.stop_after {
                        return None;
                    }
                }
                ExactNode::Split {
                    requirement: requirement as u16,
                    branches,
                }
            }
        } else {
            self.solutions.push(self.assignments.clone());
            ExactNode::Solution
        };
        let id = self.nodes.len();
        self.nodes.push(node);
        Some(id)
    }
}

pub(super) fn search(
    grid: &Grid,
    purpose: ExactPurpose,
    options: &ExactCoverSearchOptions,
) -> ExactCoverSearchResult {
    let mut result = ExactCoverSearchResult::default();
    if options.max_nodes > MAX_NODES || !grid.has_classic_sudoku_constraints() {
        result.error = Some(
            "Independent exact-cover proofs require built-in classic rules and supported limits"
                .into(),
        );
        return result;
    }
    let state = match CandidateState::capture(grid) {
        Ok(s) => s,
        Err(e) => {
            result.error = Some(e);
            return result;
        }
    };
    let Some(domains) = opposite(&state, &purpose) else {
        result.error = Some("Invalid deduction endpoint".into());
        return result;
    };
    let mut witness_search = Search::new(&state, &state.domains, options.max_nodes, 1);
    witness_search.walk();
    result.nodes_visited = witness_search.visits;
    let Some(witness) = witness_search.solutions.first().cloned() else {
        result.budget_exhausted = witness_search.cutoff;
        return result;
    };
    if let ExactPurpose::Deduction {
        cell,
        digits,
        placement,
    } = &purpose
    {
        if digits.contains(&witness[*cell]) != *placement {
            result.counterexample = Some(witness);
            return result;
        }
    }
    let mut tree = Search::new(
        &state,
        &domains,
        options.max_nodes - result.nodes_visited,
        if purpose == ExactPurpose::Unique {
            2
        } else {
            1
        },
    );
    let root = tree.walk();
    result.nodes_visited += tree.visits;
    result.budget_exhausted = tree.cutoff;
    let expected = if purpose == ExactPurpose::Unique {
        1
    } else {
        0
    };
    if tree.solutions.len() != expected || root.is_none() {
        result.counterexample = tree.solutions.last().cloned();
        return result;
    }
    let proof = ExactCoverProof {
        version: 1,
        state_hash: hash(&state),
        purpose,
        witness,
        nodes: tree.nodes,
    };
    if proof.verify(grid) {
        result.proof = Some(proof);
    } else {
        result.error = Some("Produced exact-cover tree failed independent replay".into());
    }
    result
}
