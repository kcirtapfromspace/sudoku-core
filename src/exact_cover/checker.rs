//! Replay-only checker. Geometry and support are reconstructed from board cells,
//! independently of the producer's column incidence lists and pivot heuristic.
use super::*;

fn unit(requirement: usize) -> Option<(Vec<usize>, Option<u8>)> {
    if requirement < 81 {
        return Some((vec![requirement], None));
    }
    if requirement >= 324 {
        return None;
    }
    let group = (requirement - 81) / 81;
    let number = (requirement - 81) % 81 / 9;
    let digit = ((requirement - 81) % 9 + 1) as u8;
    let cells = (0..9)
        .map(|n| match group {
            0 => number * 9 + n,
            1 => n * 9 + number,
            _ => (number / 3 * 3 + n / 3) * 9 + number % 3 * 3 + n % 3,
        })
        .collect();
    Some((cells, Some(digit)))
}

fn allowed(board: &[u8], cell: usize, digit: u8) -> bool {
    board[cell] == 0
        && (0..81).all(|other| {
            board[other] != digit
                || (other / 9 != cell / 9
                    && other % 9 != cell % 9
                    && (other / 27 != cell / 27 || other % 9 / 3 != cell % 9 / 3))
        })
}

fn support(requirement: usize, board: &[u8], domains: &[u16]) -> Option<Vec<u16>> {
    let (cells, digit) = unit(requirement)?;
    if digit.is_some_and(|d| cells.iter().any(|&c| board[c] == d))
        || (digit.is_none() && board[cells[0]] != 0)
    {
        return None;
    }
    let mut result = vec![];
    for c in cells {
        for d in 1..=9 {
            if digit.is_none_or(|v| d == v)
                && domains[c] & (1 << (d - 1)) != 0
                && allowed(board, c, d)
            {
                result.push((c * 9 + usize::from(d) - 1) as u16);
            }
        }
    }
    result.sort_unstable();
    Some(result)
}

fn witness_valid(witness: &[u8], state: &CandidateState) -> bool {
    if witness.len() != 81 || witness.iter().any(|d| !(1..=9).contains(d)) {
        return false;
    }
    for (c, &digit) in witness.iter().enumerate() {
        if (state.values[c] != 0 && digit != state.values[c])
            || state.domains[c] & (1 << (digit - 1)) == 0
        {
            return false;
        }
    }
    for req in 81..324 {
        let (cells, digit) = unit(req).unwrap();
        if cells
            .iter()
            .filter(|&&c| witness[c] == digit.unwrap())
            .count()
            != 1
        {
            return false;
        }
    }
    true
}

struct Replay<'a> {
    proof: &'a ExactCoverProof,
    domains: Vec<u16>,
    used: Vec<bool>,
    solutions: usize,
}
impl Replay<'_> {
    fn visit(&mut self, id: usize, board: &mut [u8], depth: usize) -> Option<()> {
        if id >= self.used.len() || self.used[id] || depth > 81 {
            return None;
        }
        self.used[id] = true;
        match &self.proof.nodes[id] {
            ExactNode::Conflict { requirement } => {
                if !support(usize::from(*requirement), board, &self.domains)?.is_empty() {
                    return None;
                }
            }
            ExactNode::Split {
                requirement,
                branches,
            } => {
                let expected = support(usize::from(*requirement), board, &self.domains)?;
                if expected.is_empty() || branches.len() != expected.len() {
                    return None;
                }
                for (&candidate, branch) in expected.iter().zip(branches) {
                    if branch.candidate != candidate || branch.child >= id {
                        return None;
                    }
                    let cell = usize::from(candidate) / 9;
                    board[cell] = (candidate % 9 + 1) as u8;
                    self.visit(branch.child, board, depth + 1)?;
                    board[cell] = 0;
                }
            }
            ExactNode::Solution => {
                if board.contains(&0) || board != self.proof.witness {
                    return None;
                }
                self.solutions += 1;
            }
        }
        Some(())
    }
}

pub(super) fn check(proof: &ExactCoverProof, grid: &Grid) -> Option<ExactCoverCheck> {
    if !grid.has_classic_sudoku_constraints()
        || proof.version != 1
        || proof.nodes.is_empty()
        || proof.nodes.len() > MAX_NODES
    {
        return None;
    }
    let state = CandidateState::capture(grid).ok()?;
    if proof.state_hash != hash(&state) || !witness_valid(&proof.witness, &state) {
        return None;
    }
    if let ExactPurpose::Deduction {
        cell,
        digits,
        placement,
    } = &proof.purpose
    {
        // Validate shape before indexing or shifting any supplied digit.
        opposite(&state, &proof.purpose)?;
        if digits.contains(&proof.witness[*cell]) != *placement {
            return None;
        }
    }
    let domains = opposite(&state, &proof.purpose)?;
    let mut replay = Replay {
        proof,
        domains,
        used: vec![false; proof.nodes.len()],
        solutions: 0,
    };
    replay.visit(proof.nodes.len() - 1, &mut state.values.clone(), 0)?;
    let expected = if proof.purpose == ExactPurpose::Unique {
        1
    } else {
        0
    };
    if replay.solutions != expected || !replay.used.iter().all(|&v| v) {
        return None;
    }
    Some(ExactCoverCheck {
        nodes: proof.nodes.len(),
        solutions: replay.solutions,
    })
}
