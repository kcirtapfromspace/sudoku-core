use super::*;

#[test]
fn test_generate_easy() {
    let mut generator = Generator::with_seed(42);
    let grid = generator.generate(Difficulty::Easy);

    assert!(grid.given_count() >= 30);
    assert!(grid.given_count() <= 50);

    let solver = Solver::new();
    assert!(solver.has_unique_solution(&grid));
}

#[test]
fn test_generate_medium() {
    let mut generator = Generator::with_seed(42);
    let grid = generator.generate(Difficulty::Medium);

    let solver = Solver::new();
    assert!(solver.has_unique_solution(&grid));
}

#[test]
fn test_for_se_rating_config() {
    // Low SE → Beginner tier, many givens
    let config = GeneratorConfig::for_se_rating(1.5);
    assert_eq!(config.difficulty, Difficulty::Beginner);
    assert!(config.min_givens >= 42);
    assert!(config.min_se_rating.unwrap() >= 1.2);
    assert!(config.max_se_rating.unwrap() <= 1.8);

    // Mid SE → Hard tier
    let config = GeneratorConfig::for_se_rating(4.0);
    assert_eq!(config.difficulty, Difficulty::Hard);
    assert_eq!(config.symmetry, SymmetryType::Rotational180);

    // High SE → no symmetry
    let config = GeneratorConfig::for_se_rating(7.0);
    assert_eq!(config.difficulty, Difficulty::Extreme);
    assert_eq!(config.symmetry, SymmetryType::None);

    // Extreme SE → clamped to 11.0
    let config = GeneratorConfig::for_se_rating(15.0);
    assert!(config.min_se_rating.unwrap() <= 11.0);

    // Very low → clamped to 1.5
    let config = GeneratorConfig::for_se_rating(0.5);
    assert!(config.min_se_rating.unwrap() >= 1.2);
}

#[test]
fn test_generate_for_se() {
    let mut generator = Generator::with_seed(42);
    let grid = generator.generate_for_se(3.0);

    // Should produce a valid puzzle with unique solution
    let solver = Solver::new();
    assert!(solver.has_unique_solution(&grid));
    assert!(grid.given_count() >= 17);
}

#[test]
fn test_symmetry() {
    let mut generator = Generator::with_seed(42);
    generator.config.symmetry = SymmetryType::Rotational180;
    let grid = generator.generate(Difficulty::Easy);

    // Check rotational symmetry
    for row in 0..9 {
        for col in 0..9 {
            let pos1 = Position::new(row, col);
            let pos2 = Position::new(8 - row, 8 - col);

            let has1 = grid.cell(pos1).is_given();
            let has2 = grid.cell(pos2).is_given();

            assert_eq!(has1, has2, "Symmetry broken at {:?} and {:?}", pos1, pos2);
        }
    }
}

#[test]
fn configured_symmetries_preserve_clue_patterns_and_unique_solutions() {
    for symmetry in [
        SymmetryType::None,
        SymmetryType::Rotational180,
        SymmetryType::Rotational90,
        SymmetryType::Horizontal,
        SymmetryType::Vertical,
        SymmetryType::Diagonal,
    ] {
        let mut generator = Generator::with_seed(81);
        generator.config = GeneratorConfig {
            symmetry,
            max_attempts: 1,
            min_givens: 50,
            max_givens: 60,
            ..GeneratorConfig::beginner()
        };
        let puzzle = generator.generate_with_config();
        assert!(puzzle.validate().is_valid, "{symmetry:?}");
        assert!(Solver::new().has_unique_solution(&puzzle), "{symmetry:?}");
        assert!((50..=60).contains(&puzzle.given_count()), "{symmetry:?}");
        for pos in Position::all_9x9() {
            // Independent transforms assert the public symmetry contract.
            let counterpart = match symmetry {
                SymmetryType::None => pos,
                SymmetryType::Rotational180 => Position::new(8 - pos.row, 8 - pos.col),
                SymmetryType::Rotational90 => Position::new(pos.col, 8 - pos.row),
                SymmetryType::Horizontal => Position::new(8 - pos.row, pos.col),
                SymmetryType::Vertical => Position::new(pos.row, 8 - pos.col),
                SymmetryType::Diagonal => Position::new(pos.col, pos.row),
            };
            assert_eq!(
                puzzle.cell(pos).is_given(),
                puzzle.cell(counterpart).is_given(),
                "{symmetry:?} broke at {pos:?}"
            );
        }
    }
}

#[test]
fn exhausted_generation_budgets_return_reproducible_unique_puzzles() {
    // Impossible SE bounds exercise documented best-candidate fallback; a zero
    // budget exercises the last-resort path without an unbounded search.
    for (attempts, minimum, maximum) in
        [(0, None, None), (2, Some(11.0), None), (2, None, Some(0.0))]
    {
        let generate = || {
            let mut generator = Generator::with_seed(410);
            generator.config = GeneratorConfig {
                max_attempts: attempts,
                min_givens: 50,
                max_givens: 60,
                min_se_rating: minimum,
                max_se_rating: maximum,
                ..GeneratorConfig::beginner()
            };
            generator.generate_with_config()
        };
        let first = generate();
        let second = generate();
        assert_eq!(first.to_string(), second.to_string());
        assert!(first.validate().is_valid);
        assert!(Solver::new().has_unique_solution(&first));
        assert!((50..=60).contains(&first.given_count()));
    }
}

#[test]
fn se_search_configs_bound_work_and_bracket_the_requested_rating() {
    for requested in [1.5_f32, 2.3, 3.0, 3.6, 4.0, 4.7, 5.8, 6.5, 7.5, 9.0, 11.0] {
        let config = GeneratorConfig::for_se_rating(requested);
        assert!(config.min_givens >= 17);
        assert!(config.min_givens <= config.max_givens);
        assert!(config.max_givens <= 55);
        assert!((1..=3000).contains(&config.max_attempts));
        assert!(config.min_se_rating.unwrap() <= requested);
        assert!(config.max_se_rating.unwrap() >= requested);
        assert!(config.max_se_rating.unwrap() - config.min_se_rating.unwrap() <= 1.01);
    }
}

#[test]
fn exact_clue_bounds_allow_small_symmetry_orbits_after_oversized_ones() {
    for (symmetry, givens) in [
        (SymmetryType::Rotational90, 80),
        (SymmetryType::Rotational90, 77),
        (SymmetryType::Rotational180, 79),
        (SymmetryType::None, 80),
    ] {
        let mut generator = Generator::with_seed(81);
        generator.config = GeneratorConfig {
            symmetry,
            min_givens: givens,
            max_givens: givens,
            max_attempts: 1,
            ..GeneratorConfig::beginner()
        };
        let puzzle = generator.generate_with_config();
        assert_eq!(puzzle.given_count(), givens, "{symmetry:?}");
        assert!(Solver::new().has_unique_solution(&puzzle));
    }
}
