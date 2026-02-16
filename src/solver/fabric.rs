//! CandidateFabric: dual-indexed candidate state built from Grid.
//!
//! Provides O(1) lookups for "which cells in sector S have candidate d?" and
//! "which sectors does cell C belong to?". All engines operate on `&CandidateFabric`.

use crate::{BitSet, Grid, Position};

/// Sector index convention: 0..8 = rows, 9..17 = columns, 18..26 = boxes.
pub const SECTOR_ROW_BASE: usize = 0;
pub const SECTOR_COL_BASE: usize = 9;
pub const SECTOR_BOX_BASE: usize = 18;

/// Dual-indexed candidate state, built once per solve step from Grid.
pub struct CandidateFabric {
    /// Per-cell candidates (indexed by linear cell index 0..80)
    pub cell_cands: [BitSet; 81],
    /// Placed values (None if empty)
    pub values: [Option<u8>; 81],
    /// sector_digit_cells[sector][digit-1] = bitmask of cell positions within that sector
    /// that have `digit` as a candidate. Uses a u16 bitmask (bits 0..8 for the 9 cells in sector).
    pub sector_digit_cells: [[u16; 9]; 27],
    /// Candidate count per sector per digit
    pub sector_digit_count: [[u8; 9]; 27],
    /// Which 3 sectors each cell belongs to: [row_sector, col_sector, box_sector]
    pub cell_sectors: [[usize; 3]; 81],
    /// Precomputed 20 peers per cell (cells in same row/col/box, excluding self)
    pub peers: [[u8; 20]; 81],
    /// Number of empty cells
    pub empty_count: usize,
    /// Whether cell is a given
    pub is_given: [bool; 81],
}

/// Convert (row, col) to linear cell index
#[inline]
pub fn cell_index(row: usize, col: usize) -> usize {
    row * 9 + col
}

/// Convert linear cell index back to (row, col)
#[inline]
pub fn cell_pos(idx: usize) -> (usize, usize) {
    (idx / 9, idx % 9)
}

/// Convert linear cell index to Position
#[inline]
pub fn idx_to_pos(idx: usize) -> Position {
    let (r, c) = cell_pos(idx);
    Position::new(r, c)
}

/// Convert Position to linear cell index
#[inline]
#[allow(dead_code)]
pub fn pos_to_idx(pos: Position) -> usize {
    cell_index(pos.row, pos.col)
}

/// Get the 9 cell indices belonging to a sector
pub fn sector_cells(sector: usize) -> [usize; 9] {
    if sector < 9 {
        // Row
        let row = sector;
        std::array::from_fn(|col| cell_index(row, col))
    } else if sector < 18 {
        // Column
        let col = sector - 9;
        std::array::from_fn(|row| cell_index(row, col))
    } else {
        // Box
        let box_idx = sector - 18;
        let box_row = (box_idx / 3) * 3;
        let box_col = (box_idx % 3) * 3;
        let mut cells = [0usize; 9];
        for (i, cell) in cells.iter_mut().enumerate() {
            *cell = cell_index(box_row + i / 3, box_col + i % 3);
        }
        cells
    }
}

/// Compute the 20 peers of a cell (same row/col/box, excluding self)
fn compute_peers(idx: usize) -> [u8; 20] {
    let (row, col) = cell_pos(idx);
    let box_row = (row / 3) * 3;
    let box_col = (col / 3) * 3;
    let mut peers = [0u8; 20];
    let mut count = 0;

    // Row peers
    for c in 0..9 {
        if c != col {
            peers[count] = cell_index(row, c) as u8;
            count += 1;
        }
    }
    // Column peers (excluding row already counted)
    for r in 0..9 {
        if r != row {
            peers[count] = cell_index(r, col) as u8;
            count += 1;
        }
    }
    // Box peers not already counted
    for dr in 0..3 {
        for dc in 0..3 {
            let r = box_row + dr;
            let c = box_col + dc;
            if r != row && c != col {
                peers[count] = cell_index(r, c) as u8;
                count += 1;
            }
        }
    }
    debug_assert_eq!(count, 20);
    peers
}

/// Precompute which 3 sectors each cell belongs to
fn compute_cell_sectors(idx: usize) -> [usize; 3] {
    let (row, col) = cell_pos(idx);
    let box_idx = (row / 3) * 3 + col / 3;
    [
        SECTOR_ROW_BASE + row,
        SECTOR_COL_BASE + col,
        SECTOR_BOX_BASE + box_idx,
    ]
}

impl CandidateFabric {
    /// Build the fabric from a Grid snapshot. Call once per solve step.
    pub fn from_grid(grid: &Grid) -> Self {
        let mut fab = CandidateFabric {
            cell_cands: [BitSet::empty(); 81],
            values: [None; 81],
            sector_digit_cells: [[0u16; 9]; 27],
            sector_digit_count: [[0u8; 9]; 27],
            cell_sectors: [[0; 3]; 81],
            peers: [[0; 20]; 81],
            empty_count: 0,
            is_given: [false; 81],
        };

        // Precompute static topology
        for idx in 0..81 {
            fab.cell_sectors[idx] = compute_cell_sectors(idx);
            fab.peers[idx] = compute_peers(idx);
        }

        // Populate cell data from grid
        for idx in 0..81 {
            let pos = idx_to_pos(idx);
            let cell = grid.cell(pos);
            fab.is_given[idx] = cell.is_given();

            if let Some(v) = cell.value() {
                fab.values[idx] = Some(v);
            } else {
                fab.empty_count += 1;
                let cands = grid.get_candidates(pos);
                fab.cell_cands[idx] = cands;

                // Update sector-digit indexes
                let sectors = fab.cell_sectors[idx];
                for d in cands.iter() {
                    let di = (d - 1) as usize;
                    for &sec in &sectors {
                        // Find position of this cell within the sector
                        let pos_in_sector = sector_cell_position(sec, idx);
                        fab.sector_digit_cells[sec][di] |= 1u16 << pos_in_sector;
                        fab.sector_digit_count[sec][di] += 1;
                    }
                }
            }
        }

        fab
    }

    /// Check if two cells see each other (same row, col, or box)
    #[inline]
    pub fn sees(&self, a: usize, b: usize) -> bool {
        self.cell_sectors[a][0] == self.cell_sectors[b][0]  // same row
            || self.cell_sectors[a][1] == self.cell_sectors[b][1]  // same col
            || self.cell_sectors[a][2] == self.cell_sectors[b][2] // same box
    }

    /// Get all empty cell indices
    #[allow(dead_code)]
    pub fn empty_cells(&self) -> Vec<usize> {
        (0..81).filter(|&i| self.values[i].is_none()).collect()
    }

    /// Get the number of candidates for a cell
    #[inline]
    #[allow(dead_code)]
    pub fn cand_count(&self, idx: usize) -> u32 {
        self.cell_cands[idx].count()
    }

    /// Check if cell has candidate
    #[inline]
    #[allow(dead_code)]
    pub fn has_cand(&self, idx: usize, digit: u8) -> bool {
        self.cell_cands[idx].contains(digit)
    }

    /// Get cells in a sector that have a given candidate, as a list of cell indices
    pub fn sector_cells_with_candidate(&self, sector: usize, digit: u8) -> Vec<usize> {
        let di = (digit - 1) as usize;
        let mask = self.sector_digit_cells[sector][di];
        let cells = sector_cells(sector);
        let mut result = Vec::new();
        for (i, &cell) in cells.iter().enumerate() {
            if mask & (1u16 << i) != 0 {
                result.push(cell);
            }
        }
        result
    }

    /// Get candidate count for digit in sector
    #[inline]
    pub fn sector_cand_count(&self, sector: usize, digit: u8) -> u8 {
        self.sector_digit_count[sector][(digit - 1) as usize]
    }

    /// Get all cells that see both a and b (intersection of their peer sets)
    #[allow(dead_code)]
    pub fn common_peers(&self, a: usize, b: usize) -> Vec<usize> {
        let mut result = Vec::new();
        for &peer_a in &self.peers[a] {
            let pa = peer_a as usize;
            if pa == b {
                continue;
            }
            // Check if pa is also a peer of b
            if self.sees(pa, b) {
                result.push(pa);
            }
        }
        result
    }
}

/// Find the position (0..8) of cell `idx` within the given sector.
fn sector_cell_position(sector: usize, idx: usize) -> usize {
    let cells = sector_cells(sector);
    cells
        .iter()
        .position(|&c| c == idx)
        .expect("cell not in sector")
}

// ==================== Four Grid Spaces ====================
//
// The Sudoku constraint system decomposes into four mathematical spaces.
// These view types provide typed access into the corresponding partition
// of CandidateFabric, documenting which space each engine operates in.
//
// 1. Cell Space:      81 positions → placed values / emptiness
// 2. Candidate Space: ≤729 (cell, digit) pairs → Boolean candidate membership
// 3. Sector Space:    27 houses → cell memberships, peer relationships
// 4. Link Space:      Binary relations between candidates (strong links / weak inferences)
//
// Link Space is represented by `LinkGraph` in aic_engine.rs (built on demand).

/// Cell Space: maps the 81 grid positions to their placed values.
///
/// This is the "what is known" space — each cell is either filled (given or
/// deduced) or empty. Singles operate purely in this space.
#[allow(dead_code)]
pub struct CellSpaceView<'a> {
    pub values: &'a [Option<u8>; 81],
    pub is_given: &'a [bool; 81],
    pub empty_count: usize,
}

/// Candidate Space: the set of (cell, digit) pairs that remain possible.
///
/// Each empty cell has a BitSet of candidate digits. The sector-digit index
/// provides O(1) lookup of "which cells in sector S can hold digit d?".
/// Fish and ALS engines primarily operate in this space.
#[allow(dead_code)]
pub struct CandidateSpaceView<'a> {
    pub cell_cands: &'a [BitSet; 81],
    pub sector_digit_cells: &'a [[u16; 9]; 27],
    pub sector_digit_count: &'a [[u8; 9]; 27],
}

/// Sector Space: the 27 houses (9 rows, 9 columns, 9 boxes) and their topology.
///
/// Each cell belongs to exactly 3 sectors. The peer array gives the 20 cells
/// visible from each position. This is the structural backbone connecting
/// Cell Space to Candidate Space.
#[allow(dead_code)]
pub struct SectorSpaceView<'a> {
    pub cell_sectors: &'a [[usize; 3]; 81],
    pub peers: &'a [[u8; 20]; 81],
}

#[allow(dead_code)]
impl CandidateFabric {
    /// View into Cell Space: placed values and emptiness.
    pub fn cell_space(&self) -> CellSpaceView<'_> {
        CellSpaceView {
            values: &self.values,
            is_given: &self.is_given,
            empty_count: self.empty_count,
        }
    }

    /// View into Candidate Space: per-cell candidates and sector-digit indexes.
    pub fn candidate_space(&self) -> CandidateSpaceView<'_> {
        CandidateSpaceView {
            cell_cands: &self.cell_cands,
            sector_digit_cells: &self.sector_digit_cells,
            sector_digit_count: &self.sector_digit_count,
        }
    }

    /// View into Sector Space: house memberships and peer relationships.
    pub fn sector_space(&self) -> SectorSpaceView<'_> {
        SectorSpaceView {
            cell_sectors: &self.cell_sectors,
            peers: &self.peers,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_cell_index_roundtrip() {
        for row in 0..9 {
            for col in 0..9 {
                let idx = cell_index(row, col);
                assert_eq!(cell_pos(idx), (row, col));
            }
        }
    }

    #[test]
    fn test_sector_cells() {
        // Row 0
        let row0 = sector_cells(0);
        assert_eq!(row0, [0, 1, 2, 3, 4, 5, 6, 7, 8]);

        // Col 0
        let col0 = sector_cells(9);
        assert_eq!(col0, [0, 9, 18, 27, 36, 45, 54, 63, 72]);

        // Box 0
        let box0 = sector_cells(18);
        assert_eq!(box0, [0, 1, 2, 9, 10, 11, 18, 19, 20]);
    }

    #[test]
    fn test_peers() {
        let peers = compute_peers(0); // cell (0,0)
        assert_eq!(peers.len(), 20);
        // Should include (0,1)..(0,8), (1,0)..(8,0), and box peers not in row/col
        assert!(peers.contains(&1)); // (0,1)
        assert!(peers.contains(&9)); // (1,0)
        assert!(peers.contains(&10)); // (1,1) - box peer
    }

    #[test]
    fn test_fabric_from_grid() {
        let puzzle =
            "530070000600195000098000060800060003400803001700020006060000280000419005000080079";
        let grid = Grid::from_string(puzzle).unwrap();
        let fab = CandidateFabric::from_grid(&grid);

        // Cell (0,0) = 5, should be a given
        assert_eq!(fab.values[0], Some(5));
        assert!(fab.is_given[0]);

        // Cell (0,2) is empty, should have candidates
        let idx = cell_index(0, 2);
        assert!(fab.values[idx].is_none());
        assert!(!fab.cell_cands[idx].is_empty());

        // Verify candidate doesn't include 5 (same row as (0,0))
        assert!(!fab.cell_cands[idx].contains(5));
    }

    #[test]
    fn test_sees() {
        let puzzle =
            "530070000600195000098000060800060003400803001700020006060000280000419005000080079";
        let grid = Grid::from_string(puzzle).unwrap();
        let fab = CandidateFabric::from_grid(&grid);

        // Same row
        assert!(fab.sees(0, 5));
        // Same col
        assert!(fab.sees(0, 9));
        // Same box
        assert!(fab.sees(0, 10));
        // Not seeing each other
        assert!(!fab.sees(0, 40)); // (0,0) and (4,4)
    }
}
