use super::*;

// A multi-step technique soundness check must preserve eliminations between
// steps. Public get_hint intentionally recomputes from placed values; its
// one-shot behavior is tested separately. This exercises the same dispatch
// and fallback while retaining the candidate state being proved sound.
fn retained_candidate_hint(solver: &Solver, grid: &Grid) -> Option<Hint> {
    solver
        .find_first_technique(grid)
        .map(|finding| finding.to_hint())
        .or_else(|| backtrack::find_backtracking_hint(grid).map(|finding| finding.to_hint()))
}

#[test]
fn test_solve_easy() {
    let puzzle =
        "530070000600195000098000060800060003400803001700020006060000280000419005000080079";
    let grid = Grid::from_string(puzzle).unwrap();
    let solver = Solver::new();
    let solution = solver.solve(&grid).unwrap();
    assert!(solution.is_complete());
}

#[test]
fn test_unique_solution() {
    let puzzle =
        "530070000600195000098000060800060003400803001700020006060000280000419005000080079";
    let grid = Grid::from_string(puzzle).unwrap();
    let solver = Solver::new();
    assert!(solver.has_unique_solution(&grid));
}

#[test]
fn test_get_hint() {
    let puzzle =
        "530070000600195000098000060800060003400803001700020006060000280000419005000080079";
    let grid = Grid::from_string(puzzle).unwrap();
    let solver = Solver::new();
    let hint = solver.get_hint(&grid);
    assert!(hint.is_some());
}

#[test]
fn test_difficulty_rating() {
    let puzzle =
        "530070000600195000098000060800060003400803001700020006060000280000419005000080079";
    let grid = Grid::from_string(puzzle).unwrap();
    let solver = Solver::new();
    let difficulty = solver.rate_difficulty(&grid);
    assert!(difficulty >= Difficulty::Easy);
}

#[test]
fn test_solve_with_techniques_regression() {
    let puzzle =
        "530070000600195000098000060800060003400803001700020006060000280000419005000080079";
    let grid = Grid::from_string(puzzle).unwrap();
    let solver = Solver::new();
    let mut working = grid.deep_clone();
    let max_tech = solver.solve_with_techniques(&mut working);
    assert!(max_tech < Technique::Backtracking);
    assert!(working.is_complete());
}

/// Soundness test: verify that every elimination/placement returned by hints
/// is consistent with the unique solution.
#[test]
fn test_hint_soundness() {
    let puzzles = [
        // Easy (naked/hidden singles)
        "530070000600195000098000060800060003400803001700020006060000280000419005000080079",
        // Medium
        "000000010400000000020000000000050407008000300001090000300400200050100000000806000",
        // Arto Inkala (requires advanced techniques)
        "800000000003600000070090200050007000000045700000100030001000068008500010090000400",
    ];

    let solver = Solver::new();

    for puzzle_str in &puzzles {
        let grid = Grid::from_string(puzzle_str).unwrap();
        assert!(
            solver.has_unique_solution(&grid),
            "Soundness fixture must be unique"
        );
        let solution = solver
            .solve(&grid)
            .expect("Soundness fixture must be solvable");
        assert!(solution.is_complete());

        let mut working = grid.deep_clone();
        working.recalculate_candidates();

        let mut steps = 0;
        while !working.is_complete() && steps < 300 {
            let hint = match retained_candidate_hint(&solver, &working) {
                Some(h) => h,
                None => break,
            };

            match &hint.hint_type {
                HintType::SetValue { pos, value } => {
                    let sol_val = solution.get(*pos);
                    assert_eq!(
                        sol_val,
                        Some(*value),
                        "Unsound placement by {:?}: ({},{}) = {}, solution has {:?}. Puzzle: {}",
                        hint.technique,
                        pos.row + 1,
                        pos.col + 1,
                        value,
                        sol_val,
                        puzzle_str
                    );
                    working.set_cell_unchecked(*pos, Some(*value));
                    working.recalculate_candidates();
                }
                HintType::EliminateCandidates { pos, values } => {
                    let sol_val = solution.get(*pos).expect("Position should have solution");
                    for &v in values {
                        assert_ne!(
                            v, sol_val,
                            "Unsound elimination by {:?}: removing {} from ({},{}) but solution needs it. Puzzle: {}",
                            hint.technique, v, pos.row + 1, pos.col + 1, puzzle_str
                        );
                    }
                    for &v in values {
                        working.cell_mut(*pos).remove_candidate(v);
                    }
                }
            }
            steps += 1;
        }
        assert!(
            working.is_complete(),
            "Hint chain must finish within its step limit"
        );
        assert_eq!(working.to_string_compact(), solution.to_string_compact());
    }
}

/// Soundness test for generated puzzles at every difficulty tier.
/// Uses deterministic seeds so failures are reproducible.
#[test]
fn test_hint_soundness_all_tiers() {
    use crate::Generator;

    let solver = Solver::new();

    for (seed, difficulty) in [
        (42, Difficulty::Easy),
        (42, Difficulty::Medium),
        (42, Difficulty::Intermediate),
        (42, Difficulty::Hard),
        (42, Difficulty::Expert),
        (42, Difficulty::Master),
        (42, Difficulty::Extreme),
        (99, Difficulty::Expert),
        (99, Difficulty::Master),
        (99, Difficulty::Extreme),
        (123, Difficulty::Extreme),
        (200, Difficulty::Extreme),
        (314, Difficulty::Extreme),
        (500, Difficulty::Extreme),
        (777, Difficulty::Extreme),
        (1000, Difficulty::Extreme),
        (1234, Difficulty::Extreme),
        (2024, Difficulty::Extreme),
        (9999, Difficulty::Extreme),
    ] {
        let mut gen = Generator::with_seed(seed);
        let grid = gen.generate(difficulty);

        assert!(
            solver.has_unique_solution(&grid),
            "Soundness fixture must be unique"
        );
        let solution = solver
            .solve(&grid)
            .expect("Soundness fixture must be solvable");
        assert!(solution.is_complete());

        let mut working = grid.deep_clone();
        working.recalculate_candidates();

        let mut steps = 0;
        while !working.is_complete() && steps < 500 {
            let hint = match retained_candidate_hint(&solver, &working) {
                Some(h) => h,
                None => break,
            };

            match &hint.hint_type {
                HintType::SetValue { pos, value } => {
                    let sol_val = solution.get(*pos);
                    assert_eq!(
                        sol_val,
                        Some(*value),
                        "Unsound placement by {:?} (SE {:.1}): ({},{}) = {}, solution has {:?}. Seed={} Diff={:?} Step={}",
                        hint.technique, hint.technique.se_rating(),
                        pos.row + 1, pos.col + 1, value, sol_val,
                        seed, difficulty, steps
                    );
                    working.set_cell_unchecked(*pos, Some(*value));
                    working.recalculate_candidates();
                }
                HintType::EliminateCandidates { pos, values } => {
                    let sol_val = solution.get(*pos).expect("Position should have solution");
                    for &v in values {
                        assert_ne!(
                            v, sol_val,
                            "Unsound elimination by {:?} (SE {:.1}): removing {} from ({},{}) but solution needs it. Seed={} Diff={:?} Step={}",
                            hint.technique, hint.technique.se_rating(),
                            v, pos.row + 1, pos.col + 1,
                            seed, difficulty, steps
                        );
                    }
                    for &v in values {
                        working.cell_mut(*pos).remove_candidate(v);
                    }
                }
            }
            steps += 1;
        }
        assert!(
            working.is_complete(),
            "Hint chain must finish within its step limit"
        );
        assert_eq!(working.to_string_compact(), solution.to_string_compact());
    }
}

/// Soundness test for get_next_placement() which chains eliminations.
/// This is what the WASM apply_hint() uses.
#[test]
fn test_next_placement_soundness() {
    use crate::Generator;

    let solver = Solver::new();

    for (seed, difficulty) in [
        (42, Difficulty::Expert),
        (42, Difficulty::Master),
        (42, Difficulty::Extreme),
        (99, Difficulty::Extreme),
        (123, Difficulty::Extreme),
        (200, Difficulty::Extreme),
        (314, Difficulty::Extreme),
        (500, Difficulty::Extreme),
    ] {
        let mut gen = Generator::with_seed(seed);
        let grid = gen.generate(difficulty);

        assert!(
            solver.has_unique_solution(&grid),
            "Soundness fixture must be unique"
        );
        let solution = solver
            .solve(&grid)
            .expect("Soundness fixture must be solvable");
        assert!(solution.is_complete());

        let mut working = grid.deep_clone();
        working.recalculate_candidates();

        let mut steps = 0;
        while !working.is_complete() && steps < 200 {
            let hint = match solver.get_next_placement(&working) {
                Some(h) => h,
                None => break,
            };

            match &hint.hint_type {
                HintType::SetValue { pos, value } => {
                    let sol_val = solution.get(*pos);
                    assert_eq!(
                        sol_val,
                        Some(*value),
                        "Unsound placement from get_next_placement by {:?}: ({},{}) = {}, solution has {:?}. Seed={} Diff={:?} Step={}",
                        hint.technique, pos.row + 1, pos.col + 1, value, sol_val,
                        seed, difficulty, steps
                    );
                    working.set_cell_unchecked(*pos, Some(*value));
                    working.recalculate_candidates();
                }
                HintType::EliminateCandidates { .. } => {
                    panic!("get_next_placement returned elimination instead of placement at seed={} diff={:?} step={}", seed, difficulty, steps);
                }
            }
            steps += 1;
        }

        assert!(
            working.is_complete(),
            "Failed to complete puzzle: seed={} diff={:?} stopped at step={}",
            seed,
            difficulty,
            steps
        );
    }
}

/// Verify that hints from Expert+ puzzles carry ProofCertificate data.
#[test]
fn test_hint_carries_proof_certificate() {
    let solver = Solver::new();

    // Expert puzzle that requires Empty Rectangle / fish — should produce proof certificates
    let puzzle =
        "030008002000190040108040000809060003400000790010920800061030000000000035000006000";
    let grid = Grid::from_string(puzzle).unwrap();
    assert!(solver.has_unique_solution(&grid));
    let mut working = grid.deep_clone();
    working.recalculate_candidates();

    let mut found_proof = false;
    let mut steps = 0;
    while !working.is_complete() && steps < 300 {
        let hint = match retained_candidate_hint(&solver, &working) {
            Some(h) => h,
            None => break,
        };

        if hint.proof.is_some() {
            found_proof = true;
        }

        match &hint.hint_type {
            HintType::SetValue { pos, value } => {
                working.set_cell_unchecked(*pos, Some(*value));
                working.recalculate_candidates();
            }
            HintType::EliminateCandidates { pos, values } => {
                for &v in values {
                    working.cell_mut(*pos).remove_candidate(v);
                }
            }
        }
        steps += 1;
    }

    assert!(
        working.is_complete(),
        "Proof-carrying hint chain must complete"
    );
    assert!(
        found_proof,
        "Expected at least one hint with a ProofCertificate on an Expert puzzle"
    );
}

/// Regression test: pin SE ratings for known puzzles so engine changes
/// that alter difficulty classification are caught immediately.
#[test]
fn test_se_rating_regression() {
    let cases: &[(&str, f32, Difficulty)] = &[
        // Naked singles only
        (
            "530070000600195000098000060800060003400803001700020006060000280000419005000080079",
            2.3,
            Difficulty::Easy,
        ),
        // Hidden singles only
        (
            "000000010400000000020000000000050407008000300001090000300400200050100000000806000",
            1.5,
            Difficulty::Medium,
        ),
        // Expert tier (Empty Rectangle)
        (
            "030008002000190040108040000809060003400000790010920800061030000000000035000006000",
            4.6,
            Difficulty::Expert,
        ),
    ];

    let solver = Solver::new();
    for (puzzle_str, expected_se, expected_diff) in cases {
        let grid = Grid::from_string(puzzle_str).unwrap();
        assert!(solver.has_unique_solution(&grid));
        let se = solver.rate_se(&grid);
        let diff = solver.rate_difficulty(&grid);
        assert!(
            (se - expected_se).abs() < 0.01,
            "SE regression: expected {:.1} got {:.1} for puzzle {}",
            expected_se,
            se,
            puzzle_str
        );
        assert_eq!(
            diff, *expected_diff,
            "Difficulty regression: expected {:?} got {:?} for puzzle {}",
            expected_diff, diff, puzzle_str
        );
    }
}

/// Verify that solve_with_techniques produces consistent max-technique
/// for known puzzles (guards against dispatch order regressions).
#[test]
fn test_max_technique_regression() {
    let solver = Solver::new();

    // This puzzle needs only naked singles (SE 2.3)
    let grid = Grid::from_string(
        "530070000600195000098000060800060003400803001700020006060000280000419005000080079",
    )
    .unwrap();
    let mut w = grid.deep_clone();
    let tech = solver.solve_with_techniques(&mut w);
    assert!(w.is_complete(), "Puzzle should be fully solved");
    assert!(
        tech <= Technique::NakedSingle,
        "Expected NakedSingle, got {:?}",
        tech
    );

    // This puzzle needs hidden singles (SE 1.5)
    let grid = Grid::from_string(
        "000000010400000000020000000000050407008000300001090000300400200050100000000806000",
    )
    .unwrap();
    let mut w = grid.deep_clone();
    let tech = solver.solve_with_techniques(&mut w);
    assert!(w.is_complete(), "Puzzle should be fully solved");
    assert!(
        tech <= Technique::HiddenSingle,
        "Expected ≤HiddenSingle, got {:?}",
        tech
    );
}

/// Technique coverage: verify each engine can find its expected technique
/// type on a suitable puzzle state.
#[test]
fn test_technique_coverage_basic() {
    // Verify basic techniques fire on appropriate puzzle states
    let solver = Solver::new();

    // The soundness test already covers that hints are correct.
    // Here we verify the engine can solve all three reference puzzles
    // to completion using human techniques.
    let easy = Grid::from_string(
        "530070000600195000098000060800060003400803001700020006060000280000419005000080079",
    )
    .unwrap();
    let mut w = easy.deep_clone();
    let tech = solver.solve_with_techniques(&mut w);
    assert!(w.is_complete());
    assert!(
        tech < Technique::Backtracking,
        "Easy puzzle should not need backtracking"
    );

    let medium = Grid::from_string(
        "000000010400000000020000000000050407008000300001090000300400200050100000000806000",
    )
    .unwrap();
    let mut w = medium.deep_clone();
    let tech = solver.solve_with_techniques(&mut w);
    assert!(w.is_complete());
    assert!(
        tech < Technique::Backtracking,
        "Medium puzzle should not need backtracking"
    );
}

/// Verify every technique in the dispatch chain is reachable by collecting
/// all techniques used across a battery of generated puzzles at each tier.
#[test]
fn test_technique_tier_coverage() {
    use crate::Generator;
    use std::collections::HashSet;

    let solver = Solver::new();
    let mut seen_techniques = HashSet::new();

    // Generate puzzles at different difficulty tiers using deterministic seeds
    for (seed, difficulty) in [
        (42, Difficulty::Easy),
        (42, Difficulty::Medium),
        (42, Difficulty::Intermediate),
        (42, Difficulty::Hard),
        (42, Difficulty::Expert),
        (99, Difficulty::Easy),
        (99, Difficulty::Medium),
        (99, Difficulty::Hard),
        (99, Difficulty::Expert),
    ] {
        let mut gen = Generator::with_seed(seed);
        let grid = gen.generate(difficulty);
        let mut working = grid.deep_clone();
        working.recalculate_candidates();

        // Solve step by step, collecting every technique used
        let mut steps = 0;
        while !working.is_complete() && steps < 300 {
            let hint = match retained_candidate_hint(&solver, &working) {
                Some(h) => h,
                None => break,
            };
            seen_techniques.insert(hint.technique);

            match &hint.hint_type {
                HintType::SetValue { pos, value } => {
                    working.set_cell_unchecked(*pos, Some(*value));
                    working.recalculate_candidates();
                }
                HintType::EliminateCandidates { pos, values } => {
                    for &v in values {
                        working.cell_mut(*pos).remove_candidate(v);
                    }
                }
            }
            steps += 1;
        }
    }

    // At minimum, the basic techniques should be exercised
    assert!(
        seen_techniques.contains(&Technique::NakedSingle),
        "NakedSingle never fired"
    );
    assert!(
        seen_techniques.contains(&Technique::HiddenSingle),
        "HiddenSingle never fired"
    );

    // Print coverage for diagnostics
    let mut sorted: Vec<_> = seen_techniques.iter().collect();
    sorted.sort();
    println!("Technique coverage ({} observed):", sorted.len());
    for t in &sorted {
        println!("  {:?} (SE {:.1})", t, t.se_rating());
    }
}

#[test]
fn test_collect_technique_profile() {
    let solver = Solver::new();

    // Easy puzzle (naked singles only)
    let grid = Grid::from_string(
        "530070000600195000098000060800060003400803001700020006060000280000419005000080079",
    )
    .unwrap();
    let (techs, max) = solver
        .collect_technique_profile(&grid)
        .expect("should solve easy");
    assert!(max.se_rating() <= 3.0);
    assert!(techs.contains_key("Naked Single"));

    // Expert puzzle (requires elimination techniques like Empty Rectangle / fish)
    let grid = Grid::from_string(
        "030008002000190040108040000809060003400000790010920800061030000000000035000006000",
    )
    .unwrap();
    assert!(solver.has_unique_solution(&grid));
    let (techs, max) = solver
        .collect_technique_profile(&grid)
        .expect("should solve expert");
    assert!(
        max.se_rating() > 3.0,
        "expert puzzle should have SE > 3.0, got {}",
        max.se_rating()
    );
    assert!(
        techs.len() > 2,
        "expert puzzle should use multiple techniques, got: {:?}",
        techs.keys().collect::<Vec<_>>()
    );
}

/// Historical extreme fixtures must still produce a complete, sound profile.
/// Corrected AIC/ALS deductions now solve the first three without guessing;
/// the independently solved Arto fixture keeps the backtracking regression.
#[test]
fn test_collect_technique_profile_extreme() {
    let solver = Solver::new();
    let cases = [
        (
            "020000000000380000095020810000000000840000070003601004002000000000010092600005003",
            "428157369167389425395426817719548236846932571253671984932764158574813692681295743",
            false,
        ),
        (
            "400016005390000000000080070000000200009200700080004003502000060000500000000940100",
            "478316925391725684625489371754693218139258746286174593512837469947561832863942157",
            false,
        ),
        (
            "701600000600000370000000008000100002090007800000405001500000100070300040002009000",
            "781643925625918374943572618357186492194237856268495731539824167876351249412769583",
            false,
        ),
        (
            "800000000003600000070090200050007000000045700000100030001000068008500010090000400",
            "812753649943682175675491283154237896369845721287169534521974368438526917796318452",
            true,
        ),
    ];

    for (puzzle, expected_solution, needs_backtracking) in cases {
        let grid = Grid::from_string(puzzle).unwrap();
        assert!(
            solver.has_unique_solution(&grid),
            "fixture must be unique: {puzzle}"
        );
        let (techniques, max) = solver
            .collect_technique_profile(&grid)
            .expect("profile collection must finish on every extreme fixture");
        assert!(techniques.values().all(|&count| count > 0));
        assert!(techniques.contains_key(&max.to_string()));
        if needs_backtracking {
            assert_eq!(max, Technique::Backtracking);
            assert_eq!(techniques.get("Backtracking"), Some(&1));
        } else {
            assert!(max > Technique::HiddenSingle && max < Technique::Backtracking);
            assert!(!techniques.contains_key("Backtracking"));
        }

        // An improved classification is acceptable only when every deduction
        // remains consistent with the independently enumerated unique solution.
        let mut working = grid.deep_clone();
        working.recalculate_candidates();
        for step in 0..500 {
            if working.is_complete() {
                break;
            }
            let Some(finding) = solver.find_first_technique(&working) else {
                assert!(backtrack::solve_recursive(&mut working));
                break;
            };
            match &finding.inference {
                InferenceResult::Placement { cell, value } => {
                    assert_eq!(
                        expected_solution.as_bytes()[*cell] - b'0',
                        *value,
                        "incorrect placement at step {step}: {finding:?}"
                    );
                }
                InferenceResult::Elimination { cell, values } => {
                    assert!(
                        !values.contains(&(expected_solution.as_bytes()[*cell] - b'0')),
                        "incorrect elimination at step {step}: {finding:?}"
                    );
                }
            }
            apply_finding(&mut working, &finding);
        }
        assert!(
            working.is_complete(),
            "technique trace did not finish: {puzzle}"
        );
        assert_eq!(working.to_string_compact(), expected_solution);
    }
}
