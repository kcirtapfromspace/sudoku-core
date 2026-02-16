use crate::BitSet;
use serde::{Deserialize, Serialize};

/// A single cell in the Sudoku grid
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct Cell {
    /// The current value (None if empty)
    value: Option<u8>,
    /// Candidate values (pencil marks)
    candidates: BitSet,
    /// Whether this cell was part of the original puzzle (given)
    given: bool,
}

impl Cell {
    /// Create a new empty cell with all candidates
    pub fn new_empty() -> Self {
        Self {
            value: None,
            candidates: BitSet::all_9(),
            given: false,
        }
    }

    /// Create a cell with a given value (part of the puzzle)
    pub fn new_given(value: u8) -> Self {
        Self {
            value: Some(value),
            candidates: BitSet::empty(),
            given: true,
        }
    }

    /// Create a cell with a user-entered value
    pub fn new_filled(value: u8) -> Self {
        Self {
            value: Some(value),
            candidates: BitSet::empty(),
            given: false,
        }
    }

    /// Get the cell's value
    pub fn value(&self) -> Option<u8> {
        self.value
    }

    /// Set the cell's value
    pub fn set_value(&mut self, value: Option<u8>) {
        self.value = value;
        if value.is_some() {
            self.candidates = BitSet::empty();
        }
    }

    /// Check if this cell has a value
    pub fn is_filled(&self) -> bool {
        self.value.is_some()
    }

    /// Check if this cell is empty
    pub fn is_empty(&self) -> bool {
        self.value.is_none()
    }

    /// Check if this is a given (puzzle) cell
    pub fn is_given(&self) -> bool {
        self.given
    }

    /// Mark this cell as given
    pub fn set_given(&mut self, given: bool) {
        self.given = given;
    }

    /// Get the candidates for this cell
    pub fn candidates(&self) -> BitSet {
        self.candidates
    }

    /// Set all candidates
    pub fn set_candidates(&mut self, candidates: BitSet) {
        self.candidates = candidates;
    }

    /// Add a candidate
    pub fn add_candidate(&mut self, value: u8) {
        self.candidates.insert(value);
    }

    /// Remove a candidate
    pub fn remove_candidate(&mut self, value: u8) {
        self.candidates.remove(value);
    }

    /// Toggle a candidate
    pub fn toggle_candidate(&mut self, value: u8) {
        self.candidates.toggle(value);
    }

    /// Check if a value is a candidate
    pub fn has_candidate(&self, value: u8) -> bool {
        self.candidates.contains(value)
    }

    /// Get the number of candidates
    pub fn candidate_count(&self) -> u32 {
        self.candidates.count()
    }

    /// Clear the cell (remove value and reset candidates)
    pub fn clear(&mut self) {
        if !self.given {
            self.value = None;
            self.candidates = BitSet::all_9();
        }
    }
}

impl Default for Cell {
    fn default() -> Self {
        Self::new_empty()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_new_empty() {
        let cell = Cell::new_empty();
        assert!(cell.is_empty());
        assert!(!cell.is_given());
        assert_eq!(cell.candidate_count(), 9);
    }

    #[test]
    fn test_new_given() {
        let cell = Cell::new_given(5);
        assert!(cell.is_filled());
        assert!(cell.is_given());
        assert_eq!(cell.value(), Some(5));
        assert_eq!(cell.candidate_count(), 0);
    }

    #[test]
    fn test_clear_given() {
        let mut cell = Cell::new_given(5);
        cell.clear();
        // Given cells should not be clearable
        assert!(cell.is_filled());
        assert_eq!(cell.value(), Some(5));
    }

    #[test]
    fn test_clear_non_given() {
        let mut cell = Cell::new_filled(5);
        cell.clear();
        assert!(cell.is_empty());
        assert_eq!(cell.candidate_count(), 9);
    }

    #[test]
    fn test_candidates() {
        let mut cell = Cell::new_empty();
        cell.remove_candidate(3);
        cell.remove_candidate(7);
        assert!(!cell.has_candidate(3));
        assert!(!cell.has_candidate(7));
        assert!(cell.has_candidate(5));
        assert_eq!(cell.candidate_count(), 7);

        cell.toggle_candidate(3);
        assert!(cell.has_candidate(3));
    }
}
