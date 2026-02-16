use crate::{Difficulty, Generator, Grid};

/// A seed-based puzzle identifier that deterministically maps to a specific puzzle.
///
/// Format: `[difficulty_char][7 base36 chars]` = 8 characters total.
/// Example: `M1A2B3C4` = Medium difficulty, seed `1A2B3C4` (base36).
///
/// Seed range: 0 to 36^7−1 ≈ 78 billion per difficulty.
/// 8 difficulties × 78B = ~626 billion total addressable puzzles.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PuzzleId {
    pub difficulty: Difficulty,
    pub seed: u64,
}

/// Maximum seed value (36^7 - 1)
const MAX_SEED: u64 = 36_u64.pow(7) - 1;

/// Base36 character set
const BASE36_CHARS: &[u8; 36] = b"0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ";

impl PuzzleId {
    /// Generate a random PuzzleId for the given difficulty.
    pub fn random(difficulty: Difficulty) -> Self {
        let mut seed_bytes = [0u8; 8];
        getrandom::getrandom(&mut seed_bytes).unwrap_or_else(|_| {
            static COUNTER: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(1);
            let counter = COUNTER.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
            seed_bytes = counter.to_le_bytes();
        });
        let raw = u64::from_le_bytes(seed_bytes);
        let seed = raw % (MAX_SEED + 1);
        Self { difficulty, seed }
    }

    /// Deterministically generate the puzzle grid for this id.
    pub fn generate(&self) -> Grid {
        let mut generator = Generator::with_seed(self.seed);
        generator.generate(self.difficulty)
    }

    /// Encode as an 8-character short code (e.g., "M1A2B3C4").
    pub fn to_short_code(&self) -> String {
        let diff_char = difficulty_to_char(self.difficulty);
        let seed_str = encode_base36(self.seed, 7);
        format!("{}{}", diff_char, seed_str)
    }

    /// Decode from a short code string.
    pub fn from_short_code(code: &str) -> Option<Self> {
        let code = code.trim();
        if code.len() != 8 {
            return None;
        }

        let chars: Vec<char> = code.chars().collect();
        let difficulty = char_to_difficulty(chars[0])?;
        let seed_str: String = chars[1..].iter().collect();
        let seed = decode_base36(&seed_str)?;

        if seed > MAX_SEED {
            return None;
        }

        Some(Self { difficulty, seed })
    }
}

fn difficulty_to_char(d: Difficulty) -> char {
    match d {
        Difficulty::Beginner => 'B',
        Difficulty::Easy => 'E',
        Difficulty::Medium => 'M',
        Difficulty::Intermediate => 'I',
        Difficulty::Hard => 'H',
        Difficulty::Expert => 'X',
        Difficulty::Master => 'S',
        Difficulty::Extreme => 'Z',
    }
}

fn char_to_difficulty(c: char) -> Option<Difficulty> {
    match c.to_ascii_uppercase() {
        'B' => Some(Difficulty::Beginner),
        'E' => Some(Difficulty::Easy),
        'M' => Some(Difficulty::Medium),
        'I' => Some(Difficulty::Intermediate),
        'H' => Some(Difficulty::Hard),
        'X' => Some(Difficulty::Expert),
        'S' => Some(Difficulty::Master),
        'Z' => Some(Difficulty::Extreme),
        _ => None,
    }
}

fn encode_base36(mut value: u64, width: usize) -> String {
    let mut result = vec![b'0'; width];
    for i in (0..width).rev() {
        result[i] = BASE36_CHARS[(value % 36) as usize];
        value /= 36;
    }
    String::from_utf8(result).unwrap()
}

fn decode_base36(s: &str) -> Option<u64> {
    let mut value: u64 = 0;
    for c in s.chars() {
        let digit = match c.to_ascii_uppercase() {
            '0'..='9' => (c as u64) - ('0' as u64),
            'A'..='Z' => (c.to_ascii_uppercase() as u64) - ('A' as u64) + 10,
            _ => return None,
        };
        value = value.checked_mul(36)?.checked_add(digit)?;
    }
    Some(value)
}

impl std::fmt::Display for PuzzleId {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.to_short_code())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_short_code_roundtrip() {
        let id = PuzzleId {
            difficulty: Difficulty::Hard,
            seed: 123456,
        };
        let code = id.to_short_code();
        assert_eq!(code.len(), 8);
        assert_eq!(code.chars().next().unwrap(), 'H');

        let decoded = PuzzleId::from_short_code(&code).unwrap();
        assert_eq!(decoded.difficulty, id.difficulty);
        assert_eq!(decoded.seed, id.seed);
    }

    #[test]
    fn test_all_difficulties_roundtrip() {
        let difficulties = [
            Difficulty::Beginner,
            Difficulty::Easy,
            Difficulty::Medium,
            Difficulty::Intermediate,
            Difficulty::Hard,
            Difficulty::Expert,
            Difficulty::Master,
            Difficulty::Extreme,
        ];

        for diff in difficulties {
            let id = PuzzleId {
                difficulty: diff,
                seed: 42,
            };
            let code = id.to_short_code();
            let decoded = PuzzleId::from_short_code(&code).unwrap();
            assert_eq!(decoded.difficulty, diff);
            assert_eq!(decoded.seed, 42);
        }
    }

    #[test]
    fn test_deterministic_generation() {
        let id = PuzzleId {
            difficulty: Difficulty::Medium,
            seed: 42,
        };
        let grid1 = id.generate();
        let grid2 = id.generate();
        assert_eq!(grid1.to_string_compact(), grid2.to_string_compact());
    }

    #[test]
    fn test_zero_seed() {
        let id = PuzzleId {
            difficulty: Difficulty::Beginner,
            seed: 0,
        };
        let code = id.to_short_code();
        assert_eq!(code, "B0000000");

        let decoded = PuzzleId::from_short_code(&code).unwrap();
        assert_eq!(decoded.seed, 0);
    }

    #[test]
    fn test_max_seed() {
        let id = PuzzleId {
            difficulty: Difficulty::Extreme,
            seed: MAX_SEED,
        };
        let code = id.to_short_code();
        assert_eq!(code, "ZZZZZZZZ");

        let decoded = PuzzleId::from_short_code(&code).unwrap();
        assert_eq!(decoded.seed, MAX_SEED);
    }

    #[test]
    fn test_case_insensitive_decode() {
        let code_upper = "M1A2B3C4";
        let code_lower = "m1a2b3c4";
        let id1 = PuzzleId::from_short_code(code_upper).unwrap();
        let id2 = PuzzleId::from_short_code(code_lower).unwrap();
        assert_eq!(id1, id2);
    }

    #[test]
    fn test_invalid_codes() {
        assert!(PuzzleId::from_short_code("").is_none());
        assert!(PuzzleId::from_short_code("M").is_none());
        assert!(PuzzleId::from_short_code("M123").is_none());
        assert!(PuzzleId::from_short_code("Q0000000").is_none()); // Invalid difficulty
        assert!(PuzzleId::from_short_code("M000000!").is_none()); // Invalid base36 char
    }

    #[test]
    fn test_display() {
        let id = PuzzleId {
            difficulty: Difficulty::Medium,
            seed: 42,
        };
        let display = format!("{}", id);
        assert_eq!(display, id.to_short_code());
    }

    #[test]
    fn test_encode_decode_specific_values() {
        // 1A2B3C4 in base36 = 1*36^6 + 10*36^5 + 2*36^4 + 11*36^3 + 3*36^2 + 12*36 + 4
        let id = PuzzleId::from_short_code("M1A2B3C4").unwrap();
        assert_eq!(id.difficulty, Difficulty::Medium);

        // Verify roundtrip
        let code = id.to_short_code();
        assert_eq!(code, "M1A2B3C4");
    }
}
