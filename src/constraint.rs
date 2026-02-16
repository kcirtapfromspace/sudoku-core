use crate::Position;
use serde::{Deserialize, Serialize};

/// A constraint that must be satisfied by the Sudoku grid
pub trait Constraint: Send + Sync {
    /// Check if placing a value at a position is valid according to this constraint
    fn validate(&self, cells: &[[Option<u8>; 9]; 9], pos: Position, value: u8) -> bool;

    /// Get all cells affected by this constraint from a given position
    fn affected_cells(&self, pos: Position) -> Vec<Position>;

    /// Get a name for this constraint type
    fn name(&self) -> &'static str;
}

/// A boxed constraint for dynamic dispatch
pub type ConstraintBox = Box<dyn Constraint>;

/// Standard row uniqueness constraint
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct RowConstraint;

impl Constraint for RowConstraint {
    #[allow(clippy::needless_range_loop)]
    fn validate(&self, cells: &[[Option<u8>; 9]; 9], pos: Position, value: u8) -> bool {
        for col in 0..9 {
            if col != pos.col && cells[pos.row][col] == Some(value) {
                return false;
            }
        }
        true
    }

    fn affected_cells(&self, pos: Position) -> Vec<Position> {
        (0..9)
            .filter(|&col| col != pos.col)
            .map(|col| Position::new(pos.row, col))
            .collect()
    }

    fn name(&self) -> &'static str {
        "Row"
    }
}

/// Standard column uniqueness constraint
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct ColumnConstraint;

impl Constraint for ColumnConstraint {
    #[allow(clippy::needless_range_loop)]
    fn validate(&self, cells: &[[Option<u8>; 9]; 9], pos: Position, value: u8) -> bool {
        for row in 0..9 {
            if row != pos.row && cells[row][pos.col] == Some(value) {
                return false;
            }
        }
        true
    }

    fn affected_cells(&self, pos: Position) -> Vec<Position> {
        (0..9)
            .filter(|&row| row != pos.row)
            .map(|row| Position::new(row, pos.col))
            .collect()
    }

    fn name(&self) -> &'static str {
        "Column"
    }
}

/// Standard 3x3 box uniqueness constraint
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct BoxConstraint;

impl Constraint for BoxConstraint {
    #[allow(clippy::needless_range_loop)]
    fn validate(&self, cells: &[[Option<u8>; 9]; 9], pos: Position, value: u8) -> bool {
        let box_row = (pos.row / 3) * 3;
        let box_col = (pos.col / 3) * 3;

        for row in box_row..box_row + 3 {
            for col in box_col..box_col + 3 {
                if (row != pos.row || col != pos.col) && cells[row][col] == Some(value) {
                    return false;
                }
            }
        }
        true
    }

    fn affected_cells(&self, pos: Position) -> Vec<Position> {
        let box_row = (pos.row / 3) * 3;
        let box_col = (pos.col / 3) * 3;

        let mut cells = Vec::with_capacity(8);
        for row in box_row..box_row + 3 {
            for col in box_col..box_col + 3 {
                if row != pos.row || col != pos.col {
                    cells.push(Position::new(row, col));
                }
            }
        }
        cells
    }

    fn name(&self) -> &'static str {
        "Box"
    }
}

/// Diagonal constraint for X-Sudoku (both diagonals must have unique values)
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct DiagonalConstraint;

impl Constraint for DiagonalConstraint {
    #[allow(clippy::needless_range_loop)]
    fn validate(&self, cells: &[[Option<u8>; 9]; 9], pos: Position, value: u8) -> bool {
        // Check main diagonal (top-left to bottom-right)
        if pos.is_on_main_diagonal(9) {
            for i in 0..9 {
                if i != pos.row && cells[i][i] == Some(value) {
                    return false;
                }
            }
        }

        // Check anti-diagonal (top-right to bottom-left)
        if pos.is_on_anti_diagonal(9) {
            for i in 0..9 {
                if i != pos.row && cells[i][8 - i] == Some(value) {
                    return false;
                }
            }
        }

        true
    }

    fn affected_cells(&self, pos: Position) -> Vec<Position> {
        let mut cells = Vec::new();

        if pos.is_on_main_diagonal(9) {
            for i in 0..9 {
                if i != pos.row {
                    cells.push(Position::new(i, i));
                }
            }
        }

        if pos.is_on_anti_diagonal(9) {
            for i in 0..9 {
                if i != pos.row {
                    cells.push(Position::new(i, 8 - i));
                }
            }
        }

        cells
    }

    fn name(&self) -> &'static str {
        "Diagonal"
    }
}

/// Killer Sudoku cage constraint (cells must sum to a target and be unique within the cage)
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct KillerCageConstraint {
    /// The cells that make up this cage
    pub cells: Vec<Position>,
    /// The target sum for the cage
    pub sum: u8,
}

impl KillerCageConstraint {
    pub fn new(cells: Vec<Position>, sum: u8) -> Self {
        Self { cells, sum }
    }

    /// Check if a position is in this cage
    pub fn contains(&self, pos: Position) -> bool {
        self.cells.contains(&pos)
    }
}

impl Constraint for KillerCageConstraint {
    fn validate(&self, cells: &[[Option<u8>; 9]; 9], pos: Position, value: u8) -> bool {
        if !self.contains(pos) {
            return true; // Not in this cage, doesn't apply
        }

        // Check uniqueness within cage
        for &cage_pos in &self.cells {
            if cage_pos != pos && cells[cage_pos.row][cage_pos.col] == Some(value) {
                return false;
            }
        }

        // Calculate current sum and count empty cells
        let mut current_sum: u16 = value as u16;
        let mut empty_count = 0;

        for &cage_pos in &self.cells {
            if cage_pos != pos {
                if let Some(v) = cells[cage_pos.row][cage_pos.col] {
                    current_sum += v as u16;
                } else {
                    empty_count += 1;
                }
            }
        }

        // If all cells would be filled, sum must equal target
        if empty_count == 0 {
            return current_sum == self.sum as u16;
        }

        // Otherwise, sum must not exceed target
        current_sum <= self.sum as u16
    }

    fn affected_cells(&self, pos: Position) -> Vec<Position> {
        if !self.contains(pos) {
            return Vec::new();
        }

        self.cells.iter().filter(|&&p| p != pos).copied().collect()
    }

    fn name(&self) -> &'static str {
        "KillerCage"
    }
}

/// Thermo constraint (values must increase along a path)
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ThermoConstraint {
    /// The cells in order (values must increase from first to last)
    pub path: Vec<Position>,
}

impl ThermoConstraint {
    pub fn new(path: Vec<Position>) -> Self {
        Self { path }
    }

    /// Get the index of a position in the path
    fn position_index(&self, pos: Position) -> Option<usize> {
        self.path.iter().position(|&p| p == pos)
    }
}

impl Constraint for ThermoConstraint {
    fn validate(&self, cells: &[[Option<u8>; 9]; 9], pos: Position, value: u8) -> bool {
        let Some(idx) = self.position_index(pos) else {
            return true; // Not in this thermo
        };

        // Check all cells before this one in the path (must be smaller)
        for i in 0..idx {
            let prev_pos = self.path[i];
            if let Some(prev_value) = cells[prev_pos.row][prev_pos.col] {
                if prev_value >= value {
                    return false;
                }
            }
        }

        // Check all cells after this one in the path (must be larger)
        for i in (idx + 1)..self.path.len() {
            let next_pos = self.path[i];
            if let Some(next_value) = cells[next_pos.row][next_pos.col] {
                if next_value <= value {
                    return false;
                }
            }
        }

        true
    }

    fn affected_cells(&self, pos: Position) -> Vec<Position> {
        if self.position_index(pos).is_none() {
            return Vec::new();
        }

        self.path.iter().filter(|&&p| p != pos).copied().collect()
    }

    fn name(&self) -> &'static str {
        "Thermo"
    }
}

/// Get the standard constraints for a classic 9x9 Sudoku
pub fn classic_constraints() -> Vec<ConstraintBox> {
    vec![
        Box::new(RowConstraint),
        Box::new(ColumnConstraint),
        Box::new(BoxConstraint),
    ]
}

/// Get constraints for X-Sudoku (classic + diagonals)
pub fn x_sudoku_constraints() -> Vec<ConstraintBox> {
    vec![
        Box::new(RowConstraint),
        Box::new(ColumnConstraint),
        Box::new(BoxConstraint),
        Box::new(DiagonalConstraint),
    ]
}

#[cfg(test)]
mod tests {
    use super::*;

    fn empty_grid() -> [[Option<u8>; 9]; 9] {
        [[None; 9]; 9]
    }

    #[test]
    fn test_row_constraint() {
        let mut grid = empty_grid();
        grid[0][0] = Some(5);

        let constraint = RowConstraint;

        // Same value in same row should fail
        assert!(!constraint.validate(&grid, Position::new(0, 5), 5));

        // Same value in different row should pass
        assert!(constraint.validate(&grid, Position::new(1, 5), 5));

        // Different value in same row should pass
        assert!(constraint.validate(&grid, Position::new(0, 5), 3));
    }

    #[test]
    fn test_column_constraint() {
        let mut grid = empty_grid();
        grid[0][0] = Some(5);

        let constraint = ColumnConstraint;

        // Same value in same column should fail
        assert!(!constraint.validate(&grid, Position::new(5, 0), 5));

        // Same value in different column should pass
        assert!(constraint.validate(&grid, Position::new(5, 1), 5));
    }

    #[test]
    fn test_box_constraint() {
        let mut grid = empty_grid();
        grid[0][0] = Some(5);

        let constraint = BoxConstraint;

        // Same value in same box should fail
        assert!(!constraint.validate(&grid, Position::new(2, 2), 5));

        // Same value in different box should pass
        assert!(constraint.validate(&grid, Position::new(3, 3), 5));
    }

    #[test]
    fn test_diagonal_constraint() {
        let mut grid = empty_grid();
        grid[0][0] = Some(5);
        grid[0][8] = Some(3);

        let constraint = DiagonalConstraint;

        // Same value on main diagonal should fail
        assert!(!constraint.validate(&grid, Position::new(4, 4), 5));

        // Same value on anti-diagonal should fail
        assert!(!constraint.validate(&grid, Position::new(4, 4), 3));

        // Different value on diagonals should pass
        assert!(constraint.validate(&grid, Position::new(4, 4), 7));

        // Cell not on diagonal shouldn't care about diagonal values
        assert!(constraint.validate(&grid, Position::new(0, 1), 5));
    }

    #[test]
    fn test_killer_cage_constraint() {
        let mut grid = empty_grid();
        let cage = KillerCageConstraint::new(
            vec![
                Position::new(0, 0),
                Position::new(0, 1),
                Position::new(1, 0),
            ],
            15, // Target sum is 15
        );

        grid[0][0] = Some(3);
        grid[0][1] = Some(5);

        // Sum would be 3 + 5 + 7 = 15, and all unique
        assert!(cage.validate(&grid, Position::new(1, 0), 7));

        // Sum would be 3 + 5 + 8 = 16 > 15
        assert!(!cage.validate(&grid, Position::new(1, 0), 8));

        // Duplicate value in cage (3 already exists)
        assert!(!cage.validate(&grid, Position::new(1, 0), 3));
    }

    #[test]
    fn test_thermo_constraint() {
        let mut grid = empty_grid();
        let thermo = ThermoConstraint::new(vec![
            Position::new(0, 0),
            Position::new(0, 1),
            Position::new(0, 2),
        ]);

        grid[0][0] = Some(2);
        grid[0][2] = Some(6);

        // Value must be > 2 and < 6
        assert!(thermo.validate(&grid, Position::new(0, 1), 4));
        assert!(!thermo.validate(&grid, Position::new(0, 1), 2));
        assert!(!thermo.validate(&grid, Position::new(0, 1), 6));
        assert!(!thermo.validate(&grid, Position::new(0, 1), 1));
        assert!(!thermo.validate(&grid, Position::new(0, 1), 7));
    }
}
