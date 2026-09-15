use serde::{Deserialize, Serialize};

/// A position on the Sudoku grid
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct Position {
    pub row: usize,
    pub col: usize,
}

impl Position {
    /// Create a new position
    pub const fn new(row: usize, col: usize) -> Self {
        Self { row, col }
    }

    /// Get the box index for a 9x9 grid (0-8)
    pub fn box_index(&self) -> usize {
        (self.row / 3) * 3 + (self.col / 3)
    }

    /// Get the top-left position of the box containing this position
    pub fn box_origin(&self) -> Position {
        Position::new((self.row / 3) * 3, (self.col / 3) * 3)
    }

    /// Check if this position is on the main diagonal (top-left to bottom-right)
    pub fn is_on_main_diagonal(&self, size: usize) -> bool {
        self.row == self.col && self.row < size
    }

    /// Check if this position is on the anti-diagonal (top-right to bottom-left)
    pub fn is_on_anti_diagonal(&self, size: usize) -> bool {
        self.row < size && self.col < size && self.col == size - 1 - self.row
    }

    /// Iterate over all positions in a 9x9 grid
    pub fn all_9x9() -> impl Iterator<Item = Position> {
        (0..9).flat_map(|row| (0..9).map(move |col| Position::new(row, col)))
    }
}

impl std::fmt::Display for Position {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "({}, {})", self.row, self.col)
    }
}

#[cfg(test)]
#[path = "position_tests.rs"]
mod tests;
