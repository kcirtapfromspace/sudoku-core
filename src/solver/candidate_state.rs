//! Portable native-rule premises shared by checked deduction engines.
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
pub struct CandidateState {
    pub version: u8,
    pub values: Vec<u8>,
    pub domains: Vec<u16>,
}

impl CandidateState {
    /// Capture every value and candidate, with no candidate recalculation.
    pub fn capture(grid: &Grid) -> Result<Self, String> {
        if !grid.has_standard_sudoku_constraints() {
            return Err("Candidate capture requires installed standard Sudoku constraints".into());
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
            return Err("Unsupported candidate state version or size".into());
        }
        let mut grid = Grid::new_classic();
        for cell in 0..81 {
            let (value, domain) = (self.values[cell], self.domains[cell]);
            if value > 9 || domain == 0 || domain & !ALL != 0 {
                return Err("Invalid candidate state value or domain".into());
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
            return Err("Conflicting placed values in candidate state".into());
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

fn position(cell: usize) -> Position {
    Position::new(cell / 9, cell % 9)
}
