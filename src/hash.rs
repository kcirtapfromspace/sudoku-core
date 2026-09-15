use crate::Grid;
use sha2::{Digest, Sha256};

/// Compute a canonical SHA-256 hash of a puzzle grid.
/// Uses `Grid::to_string_compact()` (digits for filled cells, `.` for empties)
/// as the canonical form. Returns a 64-character lowercase hex string.
pub fn canonical_puzzle_hash(grid: &Grid) -> String {
    let canonical = grid.to_string_compact();
    canonical_puzzle_hash_str(&canonical)
}

/// Compute SHA-256 hash from an 81-character puzzle string directly.
/// The string should use `.` for empties (canonical form).
pub fn canonical_puzzle_hash_str(puzzle_string: &str) -> String {
    let hash = Sha256::digest(puzzle_string.as_bytes());
    format!("{:x}", hash)
}

#[cfg(test)]
#[path = "hash_tests.rs"]
mod tests;
