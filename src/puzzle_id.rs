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
#[path = "puzzle_id_tests.rs"]
mod tests;
