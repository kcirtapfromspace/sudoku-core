//! Portable native-rule premises and exact-state arithmetic proof replay.
//!
//! This envelope certifies one step conditional on the recorded domains. It
//! does not certify the puzzle's uniqueness or the history of earlier deletions.
use super::{ArithmeticCheck, ArithmeticProof};
use crate::{BitSet, Grid, Position};
use serde::{Deserialize, Serialize};

const VERSION: u8 = 1;
const ALL: u16 = 0x1ff;

/// The native classic-Sudoku projection of a candidate state.
///
/// Vectors contain exactly 81 entries in row-major order. Values use zero for
/// empty cells; domain bit `digit - 1` permits a digit. Placed cells have singleton
/// domains. Extra variant constraints and given/user distinction are not needed
/// for proofs that use only the 324 native rules, and are not serialized here.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ArithmeticState {
    pub version: u8,
    pub values: Vec<u8>,
    pub domains: Vec<u16>,
}

impl ArithmeticState {
    /// Capture every value and candidate, with no candidate recalculation.
    pub fn capture(grid: &Grid) -> Result<Self, String> {
        if !grid.has_standard_sudoku_constraints() {
            return Err("Arithmetic replay requires installed standard Sudoku constraints".into());
        }
        let mut state = Self {
            version: VERSION,
            values: Vec::with_capacity(81),
            domains: Vec::with_capacity(81),
        };
        for cell in 0..81 {
            let pos = position(cell);
            let value = match grid.get(pos) {
                Some(value) if (1..=9).contains(&value) => value,
                Some(_) => return Err("Placed digit is outside 1..9".into()),
                None => 0,
            };
            state.values.push(value);
            if value != 0 {
                state.domains.push(1 << (value - 1));
            } else {
                let raw = grid.get_candidates(pos).as_raw();
                if raw == 0 || raw & !(ALL << 1) != 0 {
                    return Err("Empty cell has an empty or invalid candidate mask".into());
                }
                state.domains.push(raw >> 1);
            }
        }
        state.to_grid()?;
        Ok(state)
    }

    /// Restore actual native constraints and the exact recorded candidate masks.
    pub fn to_grid(&self) -> Result<Grid, String> {
        if self.version != VERSION || self.values.len() != 81 || self.domains.len() != 81 {
            return Err("Unsupported arithmetic state version or size".into());
        }
        let mut grid = Grid::new_classic();
        for cell in 0..81 {
            let (value, domain) = (self.values[cell], self.domains[cell]);
            if value > 9 || domain == 0 || domain & !ALL != 0 {
                return Err("Invalid arithmetic state value or domain".into());
            }
            if value != 0 {
                if domain != 1 << (value - 1) {
                    return Err("Placed value must have its singleton domain".into());
                }
                grid.set_cell_unchecked(position(cell), Some(value));
            }
        }
        // Setting values can affect peers. Install all recorded masks afterwards.
        for cell in 0..81 {
            if self.values[cell] == 0 {
                grid.cell_mut(position(cell))
                    .set_candidates(BitSet::from_raw(self.domains[cell] << 1));
            }
        }
        if !grid.validate().is_valid {
            return Err("Conflicting placed values in arithmetic state".into());
        }
        for sector in 0..27 {
            let possible = super::fabric::sector_cells(sector)
                .iter()
                .fold(0, |mask, &cell| mask | self.domains[cell]);
            if possible != ALL {
                return Err("A native sector has no remaining position for a digit".into());
            }
        }
        Ok(grid)
    }
}

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
            next.set_cell_unchecked(pos, Some(self.proof.digit));
            next.update_candidates_after_move(pos, self.proof.digit);
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
