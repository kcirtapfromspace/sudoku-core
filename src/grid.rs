use crate::{
    constraint::{classic_constraints, x_sudoku_constraints},
    BitSet, Cell, ConstraintBox, KillerCageConstraint, Position,
};
use serde::{Deserialize, Serialize};
use std::sync::Arc;

/// Error that can occur when making a move
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub enum MoveError {
    /// The cell is a given (part of the original puzzle)
    CellIsGiven,
    /// The value is out of range (must be 1-9)
    ValueOutOfRange,
    /// The position is out of bounds
    PositionOutOfBounds,
    /// The move violates a constraint
    ConstraintViolation(String),
}

impl std::fmt::Display for MoveError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            MoveError::CellIsGiven => write!(f, "Cannot modify a given cell"),
            MoveError::ValueOutOfRange => write!(f, "Value must be between 1 and 9"),
            MoveError::PositionOutOfBounds => write!(f, "Position is out of bounds"),
            MoveError::ConstraintViolation(msg) => write!(f, "Constraint violation: {}", msg),
        }
    }
}

impl std::error::Error for MoveError {}

/// Result of validating the entire grid
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct ValidationResult {
    /// Whether the grid is valid (no constraint violations)
    pub is_valid: bool,
    /// Positions of cells that violate constraints
    pub invalid_cells: Vec<Position>,
    /// Description of violations
    pub violations: Vec<String>,
}

impl ValidationResult {
    pub fn valid() -> Self {
        Self {
            is_valid: true,
            invalid_cells: Vec::new(),
            violations: Vec::new(),
        }
    }

    pub fn invalid(cells: Vec<Position>, violations: Vec<String>) -> Self {
        Self {
            is_valid: false,
            invalid_cells: cells,
            violations,
        }
    }
}

/// A 9x9 Sudoku grid
#[derive(Serialize, Deserialize)]
pub struct Grid {
    cells: [[Cell; 9]; 9],
    #[serde(skip)]
    constraints: Arc<Vec<ConstraintBox>>,
    /// Provenance of the installed rules, independent of names and variant tags.
    #[serde(skip)]
    standard_sudoku_constraints: bool,
    /// Variant type for serialization
    variant: GridVariant,
    /// Killer cages (serialized separately)
    #[serde(default)]
    killer_cages: Vec<(Vec<Position>, u8)>,
}

impl Clone for Grid {
    fn clone(&self) -> Self {
        self.deep_clone()
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, Default)]
pub enum GridVariant {
    #[default]
    Classic,
    XSudoku,
    Killer,
}

impl Grid {
    /// Create a new classic 9x9 Sudoku grid
    pub fn new_classic() -> Self {
        Self {
            cells: std::array::from_fn(|_| std::array::from_fn(|_| Cell::new_empty())),
            constraints: Arc::new(classic_constraints()),
            standard_sudoku_constraints: true,
            variant: GridVariant::Classic,
            killer_cages: Vec::new(),
        }
    }

    /// Create an X-Sudoku grid (with diagonal constraints)
    pub fn new_x_sudoku() -> Self {
        Self {
            cells: std::array::from_fn(|_| std::array::from_fn(|_| Cell::new_empty())),
            constraints: Arc::new(x_sudoku_constraints()),
            standard_sudoku_constraints: true,
            variant: GridVariant::XSudoku,
            killer_cages: Vec::new(),
        }
    }

    /// Create a grid with custom constraints
    pub fn new_with_constraints(constraints: Vec<ConstraintBox>) -> Self {
        Self {
            cells: std::array::from_fn(|_| std::array::from_fn(|_| Cell::new_empty())),
            constraints: Arc::new(constraints),
            standard_sudoku_constraints: false,
            variant: GridVariant::Classic,
            killer_cages: Vec::new(),
        }
    }

    /// Create a Killer Sudoku grid with the given cages
    pub fn new_killer(cages: Vec<KillerCageConstraint>) -> Self {
        let mut constraints: Vec<ConstraintBox> = classic_constraints();
        let killer_cages: Vec<(Vec<Position>, u8)> =
            cages.iter().map(|c| (c.cells.clone(), c.sum)).collect();

        for cage in cages {
            constraints.push(Box::new(cage));
        }

        Self {
            cells: std::array::from_fn(|_| std::array::from_fn(|_| Cell::new_empty())),
            constraints: Arc::new(constraints),
            standard_sudoku_constraints: true,
            variant: GridVariant::Killer,
            killer_cages,
        }
    }

    /// Restore constraints after deserialization
    pub fn restore_constraints(&mut self) {
        self.constraints = Arc::new(match self.variant {
            GridVariant::Classic => classic_constraints(),
            GridVariant::XSudoku => x_sudoku_constraints(),
            GridVariant::Killer => {
                let mut constraints: Vec<ConstraintBox> = classic_constraints();
                for (cells, sum) in &self.killer_cages {
                    constraints.push(Box::new(KillerCageConstraint::new(cells.clone(), *sum)));
                }
                constraints
            }
        });
        self.standard_sudoku_constraints = true;
    }

    /// Whether the installed rules are known to include standard Sudoku rules.
    ///
    /// Custom constraint lists are conservatively unverified, regardless of their
    /// names. Deserialized grids remain unverified until constraints are restored.
    pub(crate) fn has_standard_sudoku_constraints(&self) -> bool {
        self.standard_sudoku_constraints
    }

    /// Get the grid variant
    pub fn variant(&self) -> GridVariant {
        self.variant
    }

    /// Get a reference to a cell
    pub fn cell(&self, pos: Position) -> &Cell {
        &self.cells[pos.row][pos.col]
    }

    /// Get a mutable reference to a cell
    pub fn cell_mut(&mut self, pos: Position) -> &mut Cell {
        &mut self.cells[pos.row][pos.col]
    }

    /// Get the value at a position
    pub fn get(&self, pos: Position) -> Option<u8> {
        self.cells[pos.row][pos.col].value()
    }

    /// Get the values as a 2D array (for constraint checking)
    #[allow(clippy::needless_range_loop)]
    pub fn values(&self) -> [[Option<u8>; 9]; 9] {
        let mut values = [[None; 9]; 9];
        for row in 0..9 {
            for col in 0..9 {
                values[row][col] = self.cells[row][col].value();
            }
        }
        values
    }

    /// Set a cell's value with validation
    pub fn set_cell(&mut self, pos: Position, value: u8) -> Result<(), MoveError> {
        if pos.row >= 9 || pos.col >= 9 {
            return Err(MoveError::PositionOutOfBounds);
        }

        if !(1..=9).contains(&value) {
            return Err(MoveError::ValueOutOfRange);
        }

        if self.cells[pos.row][pos.col].is_given() {
            return Err(MoveError::CellIsGiven);
        }

        // Check constraints
        let values = self.values();
        for constraint in self.constraints.iter() {
            if !constraint.validate(&values, pos, value) {
                return Err(MoveError::ConstraintViolation(
                    constraint.name().to_string(),
                ));
            }
        }

        let previous = self.cells[pos.row][pos.col].value();
        self.cells[pos.row][pos.col].set_value(Some(value));
        self.update_candidates_after_move(pos, value);
        if let Some(previous) = previous.filter(|&old| old != value) {
            self.restore_candidates_after_clear(pos, previous);
        }
        Ok(())
    }

    /// Set a cell's value without validation (for solver/generator use)
    pub fn set_cell_unchecked(&mut self, pos: Position, value: Option<u8>) {
        self.cells[pos.row][pos.col].set_value(value);
    }

    /// Set a cell as a given (part of the puzzle)
    pub fn set_given(&mut self, pos: Position, value: u8) {
        self.cells[pos.row][pos.col] = Cell::new_given(value);
        self.update_candidates_after_move(pos, value);
    }

    /// Clear a cell (if not given)
    pub fn clear_cell(&mut self, pos: Position) -> Result<(), MoveError> {
        if pos.row >= 9 || pos.col >= 9 {
            return Err(MoveError::PositionOutOfBounds);
        }

        if self.cells[pos.row][pos.col].is_given() {
            return Err(MoveError::CellIsGiven);
        }

        let old_value = self.cells[pos.row][pos.col].value();
        self.cells[pos.row][pos.col].clear();
        self.cells[pos.row][pos.col].set_candidates(self.compute_candidates(pos));

        // Restore candidates that may have been removed by this cell
        if let Some(value) = old_value {
            self.restore_candidates_after_clear(pos, value);
        }

        Ok(())
    }

    /// Update candidates in affected cells after a move
    pub fn update_candidates_after_move(&mut self, pos: Position, value: u8) {
        for constraint in self.constraints.iter() {
            for affected_pos in constraint.affected_cells(pos) {
                self.cells[affected_pos.row][affected_pos.col].remove_candidate(value);
            }
        }
    }

    /// Restore candidates in affected cells after clearing a cell
    fn restore_candidates_after_clear(&mut self, pos: Position, value: u8) {
        let values = self.values();

        for constraint in self.constraints.iter() {
            for affected_pos in constraint.affected_cells(pos) {
                // A cleared row/column can overlap another constraint that still
                // forbids this value. Restored candidates must satisfy every rule.
                let can_restore = self
                    .constraints
                    .iter()
                    .all(|rule| rule.validate(&values, affected_pos, value));
                if can_restore && self.cells[affected_pos.row][affected_pos.col].is_empty() {
                    self.cells[affected_pos.row][affected_pos.col].add_candidate(value);
                }
            }
        }
    }

    /// Compute valid candidates for a cell based on current grid state
    /// (does not modify stored candidates)
    pub fn compute_candidates(&self, pos: Position) -> BitSet {
        if self.cells[pos.row][pos.col].is_filled() {
            return BitSet::empty();
        }
        let mut candidates = BitSet::all_9();
        let values = self.values();
        for constraint in self.constraints.iter() {
            for affected_pos in constraint.affected_cells(pos) {
                if let Some(v) = values[affected_pos.row][affected_pos.col] {
                    candidates.remove(v);
                }
            }
        }
        candidates
    }

    /// Get candidates for a cell
    pub fn get_candidates(&self, pos: Position) -> BitSet {
        self.cells[pos.row][pos.col].candidates()
    }

    /// Recalculate all candidates based on current values
    pub fn recalculate_candidates(&mut self) {
        // First, reset all empty cells to have all candidates
        for row in 0..9 {
            for col in 0..9 {
                if self.cells[row][col].is_empty() {
                    self.cells[row][col].set_candidates(BitSet::all_9());
                }
            }
        }

        // Then remove candidates based on filled cells
        let values = self.values();
        #[allow(clippy::needless_range_loop)]
        for row in 0..9 {
            for col in 0..9 {
                if let Some(value) = values[row][col] {
                    let pos = Position::new(row, col);
                    for constraint in self.constraints.iter() {
                        for affected_pos in constraint.affected_cells(pos) {
                            self.cells[affected_pos.row][affected_pos.col].remove_candidate(value);
                        }
                    }
                }
            }
        }
    }

    /// Clear all candidates from empty cells (for manual note-taking)
    pub fn clear_all_candidates(&mut self) {
        for row in 0..9 {
            for col in 0..9 {
                if self.cells[row][col].is_empty() {
                    self.cells[row][col].set_candidates(BitSet::empty());
                }
            }
        }
    }

    /// Validate the entire grid
    pub fn validate(&self) -> ValidationResult {
        let values = self.values();
        let mut invalid_cells = Vec::new();
        let mut violations = Vec::new();

        for row in 0..9 {
            for col in 0..9 {
                if let Some(value) = values[row][col] {
                    let pos = Position::new(row, col);

                    for constraint in self.constraints.iter() {
                        // Temporarily remove the value to check if it conflicts
                        let mut test_values = values;
                        test_values[row][col] = None;

                        if !constraint.validate(&test_values, pos, value) {
                            if !invalid_cells.contains(&pos) {
                                invalid_cells.push(pos);
                            }
                            violations.push(format!(
                                "{} constraint violated at ({}, {})",
                                constraint.name(),
                                row,
                                col
                            ));
                        }
                    }
                }
            }
        }

        if invalid_cells.is_empty() {
            ValidationResult::valid()
        } else {
            ValidationResult::invalid(invalid_cells, violations)
        }
    }

    /// Check if the puzzle is complete (all cells filled and valid)
    pub fn is_complete(&self) -> bool {
        // Check all cells are filled
        for row in 0..9 {
            for col in 0..9 {
                if self.cells[row][col].is_empty() {
                    return false;
                }
            }
        }

        // Check validity
        self.validate().is_valid
    }

    /// Check if the puzzle is solved correctly
    pub fn is_solved(&self) -> bool {
        self.is_complete()
    }

    /// Count the number of empty cells
    pub fn empty_count(&self) -> usize {
        let mut count = 0;
        for row in 0..9 {
            for col in 0..9 {
                if self.cells[row][col].is_empty() {
                    count += 1;
                }
            }
        }
        count
    }

    /// Count the number of given cells
    pub fn given_count(&self) -> usize {
        let mut count = 0;
        for row in 0..9 {
            for col in 0..9 {
                if self.cells[row][col].is_given() {
                    count += 1;
                }
            }
        }
        count
    }

    /// Get all empty positions
    pub fn empty_positions(&self) -> Vec<Position> {
        let mut positions = Vec::new();
        for row in 0..9 {
            for col in 0..9 {
                if self.cells[row][col].is_empty() {
                    positions.push(Position::new(row, col));
                }
            }
        }
        positions
    }

    /// Check if a move would be valid (without making it)
    pub fn is_valid_move(&self, pos: Position, value: u8) -> bool {
        if pos.row >= 9 || pos.col >= 9 || !(1..=9).contains(&value) {
            return false;
        }

        if self.cells[pos.row][pos.col].is_given() {
            return false;
        }

        let values = self.values();
        for constraint in self.constraints.iter() {
            if !constraint.validate(&values, pos, value) {
                return false;
            }
        }

        true
    }

    /// Get the constraints
    pub fn constraints(&self) -> &[ConstraintBox] {
        &self.constraints
    }

    /// Parse a grid from a string (81 characters, 0 or . for empty)
    pub fn from_string(s: &str) -> Option<Self> {
        let chars: Vec<char> = s.chars().filter(|c| !c.is_whitespace()).collect();

        if chars.len() != 81 {
            return None;
        }

        let mut grid = Grid::new_classic();

        for (i, c) in chars.iter().enumerate() {
            let row = i / 9;
            let col = i % 9;
            let pos = Position::new(row, col);

            match c {
                '1'..='9' => {
                    let value = c.to_digit(10)? as u8;
                    grid.set_given(pos, value);
                }
                '0' | '.' => {} // Empty cell
                _ => return None,
            }
        }

        grid.recalculate_candidates();
        Some(grid)
    }

    /// Convert the grid to a string (81 characters)
    pub fn to_string_compact(&self) -> String {
        let mut s = String::with_capacity(81);
        for row in 0..9 {
            for col in 0..9 {
                match self.cells[row][col].value() {
                    Some(v) => s.push(char::from_digit(v as u32, 10).unwrap()),
                    None => s.push('.'),
                }
            }
        }
        s
    }

    /// Create a deep clone with constraints
    pub fn deep_clone(&self) -> Self {
        Self {
            cells: self.cells.clone(),
            // Constraint definitions are immutable and may be custom implementations.
            // Preserve them while keeping every cell and its candidates independent.
            constraints: Arc::clone(&self.constraints),
            standard_sudoku_constraints: self.standard_sudoku_constraints,
            variant: self.variant,
            killer_cages: self.killer_cages.clone(),
        }
    }
}

impl std::fmt::Debug for Grid {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(f, "Grid {{")?;
        for row in 0..9 {
            write!(f, "  ")?;
            for col in 0..9 {
                match self.cells[row][col].value() {
                    Some(v) => write!(f, "{}", v)?,
                    None => write!(f, ".")?,
                }
                if col == 2 || col == 5 {
                    write!(f, "|")?;
                }
            }
            writeln!(f)?;
            if row == 2 || row == 5 {
                writeln!(f, "  ---+---+---")?;
            }
        }
        write!(f, "}}")
    }
}

impl std::fmt::Display for Grid {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for row in 0..9 {
            for col in 0..9 {
                match self.cells[row][col].value() {
                    Some(v) => write!(f, "{}", v)?,
                    None => write!(f, ".")?,
                }
                if col == 2 || col == 5 {
                    write!(f, " ")?;
                }
            }
            writeln!(f)?;
            if row == 2 || row == 5 {
                writeln!(f)?;
            }
        }
        Ok(())
    }
}

#[cfg(test)]
#[path = "grid_tests.rs"]
mod tests;
