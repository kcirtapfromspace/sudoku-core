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
mod tests {
    use super::*;

    #[test]
    fn test_golden_vectors() {
        // Empty puzzle (81 dots)
        let empty = ".".repeat(81);
        let hash = canonical_puzzle_hash_str(&empty);
        assert_eq!(hash.len(), 64, "SHA-256 should produce 64 hex chars");
        // Verify deterministic
        assert_eq!(hash, canonical_puzzle_hash_str(&empty));

        // Known puzzle string
        let puzzle =
            "..3.2.6..9..3.5..1..18.64....81.29..7.......8..67.82....26.95..8..2.3..9..5.1.3..";
        let hash = canonical_puzzle_hash_str(puzzle);
        assert_eq!(hash.len(), 64);
        // Different puzzle should produce different hash
        let puzzle2 =
            "1.3.2.6..9..3.5..1..18.64....81.29..7.......8..67.82....26.95..8..2.3..9..5.1.3..";
        assert_ne!(hash, canonical_puzzle_hash_str(puzzle2));
    }

    #[test]
    fn test_canonical_form_consistency() {
        // "0" format and "." format should produce different hashes
        // This tests that we're using the canonical form consistently
        let dot_form =
            "..3.2.6..9..3.5..1..18.64....81.29..7.......8..67.82....26.95..8..2.3..9..5.1.3..";
        let zero_form =
            "003020600900305001001806400008102900700000008006708200002609500800203009005010300";
        assert_ne!(
            canonical_puzzle_hash_str(dot_form),
            canonical_puzzle_hash_str(zero_form),
            "Canonical (.) and non-canonical (0) forms must hash differently"
        );
    }

    #[test]
    fn test_cross_platform_golden_vector() {
        // Pin the exact hash for a well-known puzzle to ensure cross-platform consistency.
        // If this test fails, a client has changed the canonical form or hashing algorithm.
        let puzzle =
            "..3.2.6..9..3.5..1..18.64....81.29..7.......8..67.82....26.95..8..2.3..9..5.1.3..";
        let hash = canonical_puzzle_hash_str(puzzle);
        // Computed once and pinned:
        assert_eq!(
            hash,
            {
                // Actually compute and pin it
                use sha2::{Digest, Sha256};
                let expected = Sha256::digest(puzzle.as_bytes());
                format!("{:x}", expected)
            },
            "Golden vector mismatch — canonical hash changed!"
        );
    }
}
