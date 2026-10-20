//! Portable native-rule premises and exact-state arithmetic proof replay.
//!
//! This envelope certifies one step conditional on the recorded domains. It
//! does not certify the puzzle's uniqueness or the history of earlier deletions.
pub use super::candidate_state::CandidateState as ArithmeticState;
use super::{ArithmeticCheck, ArithmeticProof};
use crate::{Grid, Position};
use serde::{Deserialize, Serialize};

const VERSION: u8 = 1;

/// Exportable one-step certificate with all of its candidate premises.
///
/// The envelope version describes state serialization. The proof has its own
/// independent version identifying the arithmetic grammar. Engine source
/// provenance should accompany exported artifacts at the application layer.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ArithmeticReplay {
    pub version: u8,
    pub state: ArithmeticState,
    pub proof: ArithmeticProof,
}

impl ArithmeticReplay {
    pub fn capture(grid: &Grid, proof: &ArithmeticProof) -> Result<Self, String> {
        let replay = Self {
            version: VERSION,
            state: ArithmeticState::capture(grid)?,
            proof: proof.clone(),
        };
        replay.check()?;
        Ok(replay)
    }

    /// Reconstruct premises and independently check the proof's claimed step.
    pub fn check(&self) -> Result<ArithmeticCheck, String> {
        if self.version != VERSION {
            return Err("Unsupported arithmetic replay version".into());
        }
        let grid = self.state.to_grid()?;
        self.proof
            .check(&grid)
            .ok_or_else(|| "Arithmetic certificate does not verify for its recorded state".into())
    }

    pub fn verify(&self) -> bool {
        self.check().is_ok()
    }

    /// Apply only to the identical native state, preserving earlier eliminations.
    /// The operation is atomic: errors leave the supplied grid untouched.
    pub fn apply(&self, grid: &mut Grid) -> Result<(), String> {
        self.check()?;
        if ArithmeticState::capture(grid)? != self.state {
            return Err("Arithmetic replay does not match the current candidate state".into());
        }
        let mut next = grid.deep_clone();
        let pos = position(self.proof.cell);
        if self.proof.value {
            next.place_solver_value(pos, self.proof.digit);
        } else {
            next.cell_mut(pos).remove_candidate(self.proof.digit);
        }
        if !next.validate().is_valid {
            return Err("Arithmetic replay would violate an installed constraint".into());
        }
        // Reject locally contradictory outcomes without partially changing state.
        ArithmeticState::capture(&next)?;
        *grid = next;
        Ok(())
    }
}

fn position(cell: usize) -> Position {
    Position::new(cell / 9, cell % 9)
}

#[cfg(test)]
#[path = "arithmetic_replay_tests.rs"]
mod tests;
