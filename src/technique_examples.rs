//! Reproducible real-puzzle examples with independently checked setup premises.
use crate::{
    CandidateState, ExactCoverProof, ExactPurpose, Grid, Hint, Position, ProofCertificate,
    Technique,
};
use serde::{Deserialize, Serialize};
use std::collections::HashSet;

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct TechniqueExamplePuzzle {
    pub id: String,
    pub givens: String,
    pub source_url: String,
    pub attribution: String,
    pub source_sha256: Option<String>,
    pub source_index_1based: Option<usize>,
    /// Checked against a fresh givens-only board, never against imported notes.
    pub uniqueness: ExactCoverProof,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct TechniqueExample {
    pub puzzle: usize,
    pub state: CandidateState,
    pub hint: Hint,
    /// Hint JSON omits its proof; retain the complete certificate separately.
    pub certificate: ProofCertificate,
}

impl TechniqueExample {
    pub fn checked_hint(&self) -> Hint {
        let mut h = self.hint.clone();
        h.proof = Some(self.certificate.clone());
        h
    }

    pub fn grid(&self, puzzle: &TechniqueExamplePuzzle) -> Result<Grid, String> {
        let givens = Grid::from_string(&puzzle.givens).ok_or("Invalid original puzzle")?;
        let mut grid = self.state.to_grid()?;
        for c in 0..81 {
            let p = Position::new(c / 9, c % 9);
            grid.cell_mut(p).set_given(givens.get(p).is_some());
        }
        Ok(grid)
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct TechniqueExampleCorpus {
    pub version: u8,
    pub puzzles: Vec<TechniqueExamplePuzzle>,
    pub examples: Vec<TechniqueExample>,
}

impl TechniqueExampleCorpus {
    /// Check the complete registry without running a named detector or proof search.
    /// Uniqueness of the original givens authenticates every setup placement and
    /// candidate removal: each retains the original puzzle's checked sole witness.
    pub fn verify(&self) -> Result<(), String> {
        self.check(true)
    }

    /// Validate a nonempty subset without claiming complete registry coverage.
    pub fn verify_subset(&self) -> Result<(), String> {
        self.check(false)
    }

    fn check(&self, complete: bool) -> Result<(), String> {
        if !matches!(self.version, 1 | 2)
            || (complete && self.version != 2)
            || self.examples.is_empty()
            || self.examples.len() > Technique::all().len()
            || (complete && self.examples.len() != Technique::all().len())
            || self.puzzles.is_empty()
            || self.puzzles.len() > 48
        {
            return Err("Wrong example corpus version or coverage".into());
        }
        let mut starts = vec![];
        for puzzle in &self.puzzles {
            if puzzle.givens.len() != 81
                || !puzzle
                    .givens
                    .bytes()
                    .all(|c| c == b'.' || c.is_ascii_digit())
            {
                return Err("Malformed givens string".into());
            }
            let grid = Grid::from_string(&puzzle.givens).ok_or("Invalid example puzzle")?;
            if puzzle.uniqueness.purpose != ExactPurpose::Unique || !puzzle.uniqueness.verify(&grid)
            {
                return Err(format!(
                    "Original-givens uniqueness proof rejected: {}",
                    puzzle.id
                ));
            }
            starts.push(CandidateState::capture(&grid)?);
        }
        let mut methods = HashSet::new();
        let mut used = vec![false; self.puzzles.len()];
        for e in &self.examples {
            if !methods.insert(e.hint.technique) || e.puzzle >= self.puzzles.len() {
                return Err("Duplicate technique or invalid puzzle reference".into());
            }
            used[e.puzzle] = true;
            let original = &self.puzzles[e.puzzle];
            let start = &starts[e.puzzle];
            let grid = e.grid(original)?;
            let witness = &original.uniqueness.witness;
            for (c, &digit) in witness.iter().enumerate() {
                if (start.values[c] != 0 && e.state.values[c] != start.values[c])
                    || (e.state.values[c] != 0 && e.state.values[c] != digit)
                    || e.state.domains[c] & !start.domains[c] != 0
                    || e.state.domains[c] & (1 << (digit - 1)) == 0
                {
                    return Err(format!(
                        "Unauthenticated setup premise for {:?} at cell {c}",
                        e.hint.technique
                    ));
                }
            }
            let ProofCertificate::Certified {
                pattern,
                deduction,
                named,
            } = &e.certificate
            else {
                return Err("Example lacks independent validation".into());
            };
            if pattern
                .as_deref()
                .is_some_and(|p| matches!(p, ProofCertificate::Certified { .. }))
                || !deduction.matches_hint(&e.hint)
                || !deduction.verify(&grid)
            {
                return Err(format!(
                    "Independent deduction proof rejected: {:?}",
                    e.hint.technique
                ));
            }
            if self.version == 2 && named.is_none() {
                return Err("Example lacks a named pattern certificate".into());
            }
            if let Some(p) = named {
                p.check(&grid, &e.hint)?;
            }
            if e.hint
                .involved_cells
                .iter()
                .any(|p| p.row >= 9 || p.col >= 9)
            {
                return Err("Invalid pattern coordinate".into());
            }
        }
        if (complete && Technique::all().iter().any(|t| !methods.contains(t)))
            || used.iter().any(|u| !u)
        {
            return Err("Incomplete technique or puzzle coverage".into());
        }
        Ok(())
    }
}
