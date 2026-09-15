use std::collections::HashSet;

use sudoku_core::{
    Difficulty, Grid, Hint, HintType, Position, ProofCertificate, Solver, Technique,
};

const EASY: &str =
    "530070000600195000098000060800060003400803001700020006060000280000419005000080079";
const SOLUTION: &str =
    "534678912672195348198342567859761423426853791713924856961537284287419635345286179";
const HARD: &str =
    "800000000003600000070090200050007000000045700000100030001000068008500010090000400";
const HARD_SOLUTION: &str =
    "812753649943682175675491283154237896369845721287169534521974368438526917796318452";

fn assert_solution_extends(solution: &Grid, puzzle: &Grid) {
    assert!(solution.is_complete());
    for row in 0..9 {
        let mut row_digits = HashSet::new();
        let mut col_digits = HashSet::new();
        let mut box_digits = HashSet::new();
        for col in 0..9 {
            let pos = Position::new(row, col);
            let digit = solution.get(pos).unwrap();
            assert!((1..=9).contains(&digit));
            assert!(row_digits.insert(digit));
            assert!(col_digits.insert(solution.get(Position::new(col, row)).unwrap()));
            assert!(box_digits.insert(
                solution
                    .get(Position::new(row / 3 * 3 + col / 3, row % 3 * 3 + col % 3,))
                    .unwrap()
            ));
            if let Some(given) = puzzle.get(pos) {
                assert_eq!(digit, given, "solver changed clue at {pos:?}");
            }
            assert_eq!(solution.cell(pos).is_given(), puzzle.cell(pos).is_given());
        }
    }
}

fn ambiguous_rectangle() -> Grid {
    let mut digits = SOLUTION.as_bytes().to_vec();
    // Swapping 6 and 7 in this rectangle produces exactly two completions.
    for idx in [3, 4, 30, 31] {
        digits[idx] = b'.';
    }
    Grid::from_string(std::str::from_utf8(&digits).unwrap()).unwrap()
}

#[test]
fn solves_reference_puzzles_without_mutating_input_or_clues() {
    for (puzzle, expected) in [(EASY, SOLUTION), (HARD, HARD_SOLUTION)] {
        let mut grid = Grid::from_string(puzzle).unwrap();
        // Stale pencil marks must not make an otherwise solvable puzzle fail.
        for pos in grid.empty_positions() {
            for value in 1..=9 {
                grid.cell_mut(pos).remove_candidate(value);
            }
        }
        let before = serde_json::to_value(&grid).unwrap();
        let solver = Solver::default();
        let solved = solver.solve(&grid).expect("reference puzzle is solvable");
        assert_eq!(solved.to_string_compact(), expected);
        assert_solution_extends(&solved, &grid);
        assert_eq!(solver.count_solutions(&grid, 2), 1);
        assert!(solver.has_unique_solution(&grid));
        assert_eq!(serde_json::to_value(&grid).unwrap(), before);
    }
}

#[test]
fn solution_count_obeys_zero_one_and_multiple_solution_limits() {
    let solver = Solver::new();
    let ambiguous = ambiguous_rectangle();
    for (limit, count) in [(0, 0), (1, 1), (2, 2), (3, 2)] {
        assert_eq!(solver.count_solutions(&ambiguous, limit), count);
    }
    assert!(!solver.has_unique_solution(&ambiguous));
    assert_solution_extends(&solver.solve(&ambiguous).unwrap(), &ambiguous);

    let completed = Grid::from_string(SOLUTION).unwrap();
    assert_eq!(solver.count_solutions(&completed, 0), 0);
    assert_eq!(solver.count_solutions(&completed, 10), 1);
}

#[test]
fn ambiguous_rectangle_uses_a_backtracking_hint_and_profile() {
    let solver = Solver::new();
    let grid = ambiguous_rectangle();
    let hint = solver
        .get_hint(&grid)
        .expect("an ambiguous but solvable grid still supports a trial hint");
    assert_eq!(
        hint.technique,
        Technique::Backtracking,
        "unexpected deduction: {hint:?}"
    );
    let HintType::SetValue { pos, value } = hint.hint_type else {
        panic!("backtracking must supply a placement")
    };
    let mut with_hint = grid.deep_clone();
    with_hint.set_cell(pos, value).unwrap();
    assert_eq!(solver.count_solutions(&with_hint, 2), 1);
    let (profile, hardest) = solver.collect_technique_profile(&grid).unwrap();
    assert_eq!(hardest, Technique::Backtracking);
    assert_eq!(profile.len(), 1);
    assert_eq!(profile["Backtracking"], 1);
    assert_eq!(solver.analyze(&grid), (Difficulty::Extreme, 11.0));
}

#[test]
fn inconsistent_full_and_zero_candidate_grids_have_no_solution() {
    let solver = Solver::new();
    let mut duplicate = Grid::from_string(SOLUTION).unwrap();
    duplicate.set_cell_unchecked(Position::new(0, 0), Some(3));
    let mut impossible = Grid::new_classic();
    // r1c1 needs 9 to complete its row, but its column already contains 9.
    for col in 1..9 {
        impossible.set_given(Position::new(0, col), col as u8);
    }
    impossible.set_given(Position::new(1, 0), 9);
    impossible.recalculate_candidates();
    assert!(impossible.validate().is_valid);
    assert!(impossible.get_candidates(Position::new(0, 0)).is_empty());
    for grid in [&duplicate, &impossible] {
        assert!(solver.solve(grid).is_none());
        assert_eq!(solver.count_solutions(grid, 2), 0);
        assert!(!solver.has_unique_solution(grid));
        assert!(solver.get_next_placement(grid).is_none());
    }
    assert!(solver.collect_technique_profile(&duplicate).is_none());
}

#[test]
fn solved_grid_has_no_hint_and_an_empty_technique_profile() {
    let grid = Grid::from_string(SOLUTION).unwrap();
    let solver = Solver::new();
    assert!(solver.get_hint(&grid).is_none());
    assert!(solver.get_next_placement(&grid).is_none());
    let (profile, hardest) = solver.collect_technique_profile(&grid).unwrap();
    assert!(profile.is_empty());
    assert_eq!(hardest, Technique::NakedSingle);
    assert_eq!(solver.analyze(&grid), (Difficulty::Beginner, 2.3));
}

#[test]
fn single_hint_has_coordinates_explanation_proof_and_stable_wire_shape() {
    let mut puzzle = SOLUTION.as_bytes().to_vec();
    puzzle[80] = b'.';
    let grid = Grid::from_string(std::str::from_utf8(&puzzle).unwrap()).unwrap();
    let hint = Solver::new().get_hint(&grid).unwrap();
    assert_eq!(hint.technique, Technique::NakedSingle);
    assert!(matches!(
        hint.hint_type,
        HintType::SetValue {
            pos: Position { row: 8, col: 8 },
            value: 9
        }
    ));
    assert_eq!(hint.involved_cells, vec![Position::new(8, 8)]);
    assert_eq!(
        hint.explanation,
        "Cell (9, 9) can only be 9 - it's the only candidate left."
    );
    assert!(matches!(hint.proof, Some(ProofCertificate::Basic { .. })));
    let wire = serde_json::to_value(&hint).unwrap();
    assert!(
        wire.get("proof").is_none(),
        "proof certificates are not part of the public wire format"
    );
    assert_eq!(
        wire["hint_type"],
        serde_json::json!({"SetValue": {"pos": {"row": 8, "col": 8}, "value": 9}})
    );
    let restored: Hint = serde_json::from_value(wire.clone()).unwrap();
    assert!(restored.proof.is_none());
    assert_eq!(serde_json::to_value(restored).unwrap(), wire);
}

#[test]
fn placement_hint_chain_completes_hard_puzzle_with_reference_solution() {
    let solver = Solver::new();
    let mut working = Grid::from_string(HARD).unwrap();
    let reference = Grid::from_string(HARD_SOLUTION).unwrap();
    let initial_empty = working.empty_count();
    for step in 0..initial_empty {
        let before = serde_json::to_value(&working).unwrap();
        let hint = solver
            .get_next_placement(&working)
            .expect("every unsolved reference state needs a placement");
        assert_eq!(
            serde_json::to_value(&working).unwrap(),
            before,
            "hint mutated caller at step {step}"
        );
        assert!(!hint.explanation.is_empty());
        assert!(!hint.involved_cells.is_empty());
        assert!(hint.proof.is_some());
        let HintType::SetValue { pos, value } = hint.hint_type else {
            panic!("placement API returned an elimination");
        };
        assert!(working.cell(pos).is_empty());
        assert_eq!(
            reference.get(pos),
            Some(value),
            "incorrect {:?} hint at step {step}",
            hint.technique
        );
        working.set_cell(pos, value).unwrap();
        assert_eq!(working.empty_count(), initial_empty - step - 1);
    }
    assert_eq!(working.to_string_compact(), HARD_SOLUTION);
    assert!(solver.get_next_placement(&working).is_none());
}

#[test]
fn profile_and_analysis_agree_for_simple_and_advanced_reference_puzzles() {
    let solver = Solver::new();
    let cases = [
        (EASY, Difficulty::Easy, Technique::NakedSingle),
        (
            "000000010400000000020000000000050407008000300001090000300400200050100000000806000",
            Difficulty::Medium,
            Technique::HiddenSingle,
        ),
        (
            "030008002000190040108040000809060003400000790010920800061030000000000035000006000",
            Difficulty::Expert,
            Technique::EmptyRectangle,
        ),
    ];
    for (puzzle, difficulty, expected_hardest) in cases {
        let grid = Grid::from_string(puzzle).unwrap();
        assert!(solver.has_unique_solution(&grid));
        let before = serde_json::to_value(&grid).unwrap();
        let (profile, hardest) = solver.collect_technique_profile(&grid).unwrap();
        assert_eq!(hardest, expected_hardest);
        assert!(profile.values().all(|&count| count > 0));
        assert!(!profile.contains_key("Backtracking"));
        assert!(profile.contains_key(&hardest.to_string()));
        let placements = profile.get("Naked Single").copied().unwrap_or(0)
            + profile.get("Hidden Single").copied().unwrap_or(0);
        assert_eq!(placements as usize, grid.empty_count());
        assert_eq!(solver.analyze(&grid), (difficulty, hardest.se_rating()));
        assert_eq!(solver.rate_difficulty(&grid), difficulty);
        assert_eq!(solver.rate_se(&grid), hardest.se_rating());
        assert_eq!(serde_json::to_value(&grid).unwrap(), before);
    }
}

#[test]
fn difficulty_metadata_covers_every_public_tier() {
    let expected = [
        (
            Difficulty::Beginner,
            "Beginner",
            Technique::NakedSingle,
            false,
            (1.5, 2.0),
        ),
        (
            Difficulty::Easy,
            "Easy",
            Technique::NakedSingle,
            false,
            (2.0, 2.5),
        ),
        (
            Difficulty::Medium,
            "Medium",
            Technique::HiddenSingle,
            false,
            (2.5, 3.4),
        ),
        (
            Difficulty::Intermediate,
            "Intermediate",
            Technique::HiddenTriple,
            false,
            (3.4, 3.8),
        ),
        (
            Difficulty::Hard,
            "Hard",
            Technique::BoxLineReduction,
            false,
            (3.8, 4.5),
        ),
        (
            Difficulty::Expert,
            "Expert",
            Technique::HiddenRectangle,
            false,
            (4.5, 5.5),
        ),
        (
            Difficulty::Master,
            "Master",
            Technique::BivalueUniversalGrave,
            true,
            (5.5, 7.0),
        ),
        (
            Difficulty::Extreme,
            "Extreme",
            Technique::Backtracking,
            true,
            (7.0, 11.0),
        ),
    ];
    assert_eq!(
        Difficulty::all_levels(),
        expected.iter().map(|case| case.0).collect::<Vec<_>>()
    );
    assert_eq!(
        Difficulty::standard_levels(),
        &Difficulty::all_levels()[..6]
    );
    for (difficulty, label, max, secret, range) in expected {
        assert_eq!(difficulty.to_string(), label);
        assert_eq!(difficulty.max_technique(), max);
        assert_eq!(difficulty.is_secret(), secret);
        assert_eq!(difficulty.se_range(), range);
        assert_eq!(difficulty.se_default(), (range.0 + range.1) / 2.0);
        assert!(!difficulty.technique_hint().is_empty());
        let wire = serde_json::to_value(difficulty).unwrap();
        assert_eq!(wire, serde_json::json!(label));
        assert_eq!(
            serde_json::from_value::<Difficulty>(wire).unwrap(),
            difficulty
        );
    }
}

#[test]
#[allow(deprecated)]
fn technique_names_ratings_and_families_remain_compatible() {
    use Technique::*;
    let techniques = [
        (NakedSingle, "Naked Single", 2.3),
        (HiddenSingle, "Hidden Single", 1.5),
        (NakedPair, "Naked Pair", 3.0),
        (HiddenPair, "Hidden Pair", 3.4),
        (NakedTriple, "Naked Triple", 3.6),
        (HiddenTriple, "Hidden Triple", 3.8),
        (PointingPair, "Pointing Pair", 2.6),
        (BoxLineReduction, "Box/Line Reduction", 2.8),
        (XWing, "X-Wing", 3.2),
        (FinnedXWing, "Finned X-Wing", 3.4),
        (Swordfish, "Swordfish", 3.8),
        (FinnedSwordfish, "Finned Swordfish", 4.0),
        (Jellyfish, "Jellyfish", 5.2),
        (FinnedJellyfish, "Finned Jellyfish", 5.4),
        (NakedQuad, "Naked Quad", 5.0),
        (HiddenQuad, "Hidden Quad", 5.4),
        (EmptyRectangle, "Empty Rectangle", 4.6),
        (AvoidableRectangle, "Avoidable Rectangle", 4.6),
        (UniqueRectangle, "Unique Rectangle", 4.6),
        (HiddenRectangle, "Hidden Rectangle", 4.7),
        (XYWing, "XY-Wing", 4.2),
        (XYZWing, "XYZ-Wing", 4.4),
        (WXYZWing, "WXYZ-Wing", 4.6),
        (WWing, "W-Wing", 4.4),
        (XChain, "X-Chain", 4.5),
        (ThreeDMedusa, "3D Medusa", 5.0),
        (SueDeCoq, "Sue de Coq", 5.0),
        (AIC, "AIC", 6.0),
        (FrankenFish, "Franken Fish", 5.5),
        (SiameseFish, "Siamese Fish", 5.5),
        (AlsXz, "ALS-XZ", 5.5),
        (ExtendedUniqueRectangle, "Extended Unique Rectangle", 5.5),
        (BivalueUniversalGrave, "BUG+1", 5.6),
        (AlsXyWing, "ALS-XY-Wing", 7.0),
        (AlsChain, "ALS Chain", 7.5),
        (MutantFish, "Mutant Fish", 6.5),
        (AlignedPairExclusion, "Aligned Pair Exclusion", 6.2),
        (AlignedTripletExclusion, "Aligned Triplet Exclusion", 7.5),
        (DeathBlossom, "Death Blossom", 8.5),
        (ArithmeticCounting, "Arithmetic Counting", 8.5),
        (NishioForcingChain, "Nishio Forcing Chain", 7.5),
        (KrakenFish, "Kraken Fish", 8.0),
        (RegionForcingChain, "Region Forcing Chain", 8.5),
        (CellForcingChain, "Cell Forcing Chain", 8.3),
        (DynamicForcingChain, "Dynamic Forcing Chain", 9.3),
        (Backtracking, "Backtracking", 11.0),
    ];
    let mut names = HashSet::new();
    for (technique, label, rating) in techniques {
        assert!(names.insert(label));
        assert_eq!(technique.to_string(), label);
        assert_eq!(technique.se_rating(), rating);
        assert_eq!(
            technique.is_legacy(),
            [ThreeDMedusa, AlignedPairExclusion, AlignedTripletExclusion].contains(&technique)
        );
        assert_eq!(
            technique.is_als_xz_family(),
            [XYWing, XYZWing, WXYZWing, AlsXz].contains(&technique)
        );
        assert_eq!(
            serde_json::from_value::<Technique>(serde_json::to_value(technique).unwrap()).unwrap(),
            technique
        );
    }
}

#[test]
fn arithmetic_technique_has_a_stable_wire_name_and_advanced_priority() {
    assert!(Technique::DeathBlossom < Technique::ArithmeticCounting);
    assert!(Technique::ArithmeticCounting < Technique::NishioForcingChain);
    assert_eq!(
        serde_json::to_value(Technique::ArithmeticCounting).unwrap(),
        serde_json::json!("ArithmeticCounting")
    );
    assert!(!Technique::ArithmeticCounting.is_legacy());
}

#[test]
fn arithmetic_hint_rejects_invalid_candidate_bits_before_building_solver_state() {
    for raw in [0, 1, 1 << 10, (1 << 3) | (1 << 10)] {
        let mut grid = Grid::new_classic();
        grid.cell_mut(Position::new(0, 0))
            .set_candidates(sudoku_core::BitSet::from_raw(raw));
        let before = serde_json::to_value(&grid).unwrap();
        assert!(Solver::new().get_arithmetic_hint(&grid).is_none());
        assert_eq!(serde_json::to_value(&grid).unwrap(), before);
    }
}
