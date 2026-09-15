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
