//! Portable exhaustive proofs, independent of named technique discovery.
//! The producer uses exact-cover columns. The checker reconstructs Sudoku rules
//! directly from cells and digits and never invokes the producer or the solver.
mod checker;
mod producer;

use crate::{CandidateState, Grid, Hint, HintType, ProofCertificate};
use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};

const MAX_NODES: usize = 2_000_000;

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub enum ExactPurpose {
    Unique,
    Deduction {
        cell: usize,
        digits: Vec<u8>,
        placement: bool,
    },
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ExactBranch {
    /// Row-major candidate index: 9 * cell + digit - 1.
    pub candidate: u16,
    /// Postorder child index; strictly smaller than the parent index.
    pub child: usize,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub enum ExactNode {
    /// No candidate can satisfy this unfulfilled native requirement.
    Conflict { requirement: u16 },
    /// Every available candidate for the selected requirement has one child.
    Split {
        requirement: u16,
        branches: Vec<ExactBranch>,
    },
    /// The current assignments complete every native requirement.
    Solution,
}

/// A complete native exact-cover tree, plus a checked feasible base witness.
/// A deduction tree refutes the opposite conclusion. A uniqueness tree contains
/// exactly one solution. Neither purpose certifies a named method's classification.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ExactCoverProof {
    pub version: u8,
    pub state_hash: String,
    pub purpose: ExactPurpose,
    pub witness: Vec<u8>,
    pub nodes: Vec<ExactNode>,
}

#[derive(Debug, Clone, Copy)]
pub struct ExactCoverSearchOptions {
    /// Shared witness and exhaustive-tree visits; at most two million.
    pub max_nodes: usize,
}
impl Default for ExactCoverSearchOptions {
    fn default() -> Self {
        Self { max_nodes: 100_000 }
    }
}

#[derive(Debug, Clone, Default)]
pub struct ExactCoverSearchResult {
    pub proof: Option<ExactCoverProof>,
    pub nodes_visited: usize,
    pub budget_exhausted: bool,
    pub counterexample: Option<Vec<u8>>,
    pub error: Option<String>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ExactCoverCheck {
    pub nodes: usize,
    pub solutions: usize,
}

pub(crate) fn hash(state: &CandidateState) -> String {
    let mut h = Sha256::new();
    h.update(b"sudoku-exact-cover-v1-native");
    for (&v, &d) in state.values.iter().zip(&state.domains) {
        h.update([v]);
        h.update(d.to_le_bytes());
    }
    format!("{:x}", h.finalize())
}

pub(crate) fn purpose(hint: &Hint) -> Option<ExactPurpose> {
    let (p, digits, placement) = match &hint.hint_type {
        HintType::SetValue { pos, value } => (*pos, vec![*value], true),
        HintType::EliminateCandidates { pos, values } => (*pos, values.clone(), false),
    };
    if p.row >= 9 || p.col >= 9 {
        return None;
    }
    Some(ExactPurpose::Deduction {
        cell: p.row * 9 + p.col,
        digits,
        placement,
    })
}

fn opposite(state: &CandidateState, purpose: &ExactPurpose) -> Option<Vec<u16>> {
    let mut domains = state.domains.clone();
    if let ExactPurpose::Deduction {
        cell,
        digits,
        placement,
    } = purpose
    {
        if *cell >= 81
            || state.values[*cell] != 0
            || digits.is_empty()
            || (*placement && digits.len() != 1)
            || digits
                .iter()
                .enumerate()
                .any(|(i, d)| !(1..=9).contains(d) || digits[..i].contains(d))
        {
            return None;
        }
        let mask = digits.iter().fold(0u16, |m, d| m | (1 << (d - 1)));
        if mask & !domains[*cell] != 0 {
            return None;
        }
        domains[*cell] &= if *placement { !mask } else { mask };
    }
    Some(domains)
}

impl ExactCoverProof {
    pub fn search_unique(grid: &Grid, options: &ExactCoverSearchOptions) -> ExactCoverSearchResult {
        producer::search(grid, ExactPurpose::Unique, options)
    }
    pub fn search_hint(
        grid: &Grid,
        hint: &Hint,
        options: &ExactCoverSearchOptions,
    ) -> ExactCoverSearchResult {
        match purpose(hint) {
            Some(p) => producer::search(grid, p, options),
            None => ExactCoverSearchResult {
                error: Some("Invalid hint coordinates".into()),
                ..Default::default()
            },
        }
    }
    pub fn check(&self, grid: &Grid) -> Option<ExactCoverCheck> {
        checker::check(self, grid)
    }
    pub fn verify(&self, grid: &Grid) -> bool {
        self.check(grid).is_some()
    }
    pub fn matches_hint(&self, hint: &Hint) -> bool {
        purpose(hint).as_ref() == Some(&self.purpose)
    }

    /// Add independent endpoint validation without replacing detector metadata.
    /// Search cutoffs and counterexamples return errors, never verified badges.
    pub fn certify_hint(
        grid: &Grid,
        hint: &Hint,
        options: &ExactCoverSearchOptions,
    ) -> Result<Hint, String> {
        let result = Self::search_hint(grid, hint, options);
        let proof = result.proof.ok_or_else(|| {
            result.error.unwrap_or_else(|| {
                if result.budget_exhausted {
                    "Independent proof search exhausted its node budget".into()
                } else {
                    "The independent checker could not establish this deduction".into()
                }
            })
        })?;
        let mut checked = hint.clone();
        let pattern = checked.proof.take().map(Box::new);
        checked.proof = Some(ProofCertificate::Certified {
            named: None,
            pattern,
            deduction: proof,
        });
        Ok(checked)
    }
}
