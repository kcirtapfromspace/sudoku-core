use super::*;
use crate::{Grid, Position, Solver};

const SOLUTION: &str =
    "534678912672195348198342567859761423426853791713924856961537284287419635345286179";
const TYPE1: &str =
    "034078912072090348198342567859761423426853791713924856961037284287419635345286179";
const TYPE2: &str =
    "034078912002090348098342567859761423426853791703924856961037284287419635345286179";
const TYPE3: &str =
    "004008902002090348098342567809061423426803791713924856961537284287419635345286179";
const TYPE4: &str =
    "034008902072095348098342567859061423426853791713924856961537284287419635345286179";
const EXTENDED: &str =
    "034070912072090348198342567859761423426853791713924856961037284287419635345280179";

fn transform(cell: usize, transpose: bool) -> usize {
    if transpose {
        cell % 9 * 9 + cell / 9
    } else {
        cell
    }
}

// Every fixture starts from a unique puzzle, with all manually retained candidates
// compatible with its givens and its known solution. Removing other candidates is
// a solution-preserving prior deduction, not an invented impossible pencil state.
fn fixture(puzzle: &str, candidates: &[(usize, &[u8])], transpose: bool) -> (Grid, Grid) {
    let mut text = vec![b'0'; 81];
    let mut answer = vec![b'0'; 81];
    for cell in 0..81 {
        text[transform(cell, transpose)] = puzzle.as_bytes()[cell];
        answer[transform(cell, transpose)] = SOLUTION.as_bytes()[cell];
    }
    let mut grid = Grid::from_string(std::str::from_utf8(&text).unwrap()).unwrap();
    let solution = Grid::from_string(std::str::from_utf8(&answer).unwrap()).unwrap();
    let solver = Solver::new();
    assert!(
        solver.has_unique_solution(&grid),
        "Fixture must really be unique"
    );
    assert_eq!(
        solver.solve(&grid).unwrap().to_string_compact(),
        solution.to_string_compact()
    );
    let initial = grid.deep_clone();
    for pos in grid.empty_positions() {
        grid.cell_mut(pos)
            .set_candidates(BitSet::from_slice(&[solution.get(pos).unwrap()]));
    }
    for &(cell, values) in candidates {
        let pos = idx_to_pos(transform(cell, transpose));
        let selected = BitSet::from_slice(values);
        assert!(initial.cell(pos).is_empty());
        assert!(
            selected.difference(&initial.get_candidates(pos)).is_empty(),
            "Candidates must satisfy givens at {cell}"
        );
        assert!(
            selected.contains(solution.get(pos).unwrap()),
            "Candidates must preserve solution at {cell}"
        );
        grid.cell_mut(pos).set_candidates(selected);
    }
    (grid, solution)
}

fn assert_elimination(finding: &Finding, grid: &Grid, solution: &Grid, pattern: &str) {
    let InferenceResult::Elimination { cell, values } = &finding.inference else {
        panic!("Expected elimination: {finding:?}");
    };
    let pos = idx_to_pos(*cell);
    assert!(!values.is_empty());
    for &value in values {
        assert!(grid.get_candidates(pos).contains(value));
        assert_ne!(
            solution.get(pos),
            Some(value),
            "Unsound finding: {finding:?}"
        );
    }
    let Some(ProofCertificate::Uniqueness {
        pattern: actual,
        floor_cells,
        roof_cells,
    }) = &finding.proof
    else {
        panic!("Missing uniqueness proof: {finding:?}");
    };
    assert_eq!(actual, pattern);
    assert!(!floor_cells.is_empty());
    assert!(!roof_cells.is_empty());
    assert!(floor_cells.iter().all(|cell| !roof_cells.contains(cell)));
}

#[test]
fn unique_rectangle_type_one_preserves_solution_in_both_roof_orders() {
    let candidates: &[(usize, &[u8])] =
        &[(0, &[5, 6]), (3, &[5, 6]), (9, &[5, 6]), (12, &[1, 5, 6])];
    for transpose in [false, true] {
        let (grid, solution) = fixture(TYPE1, candidates, transpose);
        let fab = CandidateFabric::from_grid(&grid);
        for reverse in [false, true] {
            let (third, fourth) = if reverse { (12, 9) } else { (9, 12) };
            let finding = try_ur_hint(
                &fab,
                transform(0, transpose),
                transform(3, transpose),
                transform(third, transpose),
                transform(fourth, transpose),
                5,
                6,
            )
            .unwrap();
            assert_elimination(&finding, &grid, &solution, "UR Type 1");
        }
        let finding = find_unique_rectangle(&fab).unwrap();
        assert_elimination(&finding, &grid, &solution, "UR Type 1");
    }
}

#[test]
fn unique_rectangle_type_two_excludes_shared_extra_from_a_common_peer() {
    for transpose in [false, true] {
        let (grid, solution) = fixture(
            TYPE2,
            &[
                (0, &[5, 6]),
                (3, &[5, 6]),
                (9, &[1, 5, 6]),
                (12, &[1, 5, 6]),
                (10, &[1, 7]),
            ],
            transpose,
        );
        let fab = CandidateFabric::from_grid(&grid);
        let finding = try_ur_hint(
            &fab,
            transform(0, transpose),
            transform(3, transpose),
            transform(9, transpose),
            transform(12, transpose),
            5,
            6,
        )
        .unwrap();
        assert_elimination(&finding, &grid, &solution, "UR Type 2");
        assert!(
            matches!(finding.inference,InferenceResult::Elimination{cell,ref values} if cell==transform(10,transpose) && values==&vec![1])
        );
    }
}

#[test]
fn unique_rectangle_type_three_combines_roof_extras_with_a_naked_pair() {
    for transpose in [false, true] {
        let (grid, solution) = fixture(
            TYPE3,
            &[
                (9, &[1, 6]),
                (12, &[1, 6]),
                (0, &[1, 5, 6]),
                (3, &[1, 6, 7]),
                (4, &[5, 7]),
                (1, &[3, 5, 7]),
            ],
            transpose,
        );
        let fab = CandidateFabric::from_grid(&grid);
        let finding = try_ur_hint(
            &fab,
            transform(9, transpose),
            transform(12, transpose),
            transform(0, transpose),
            transform(3, transpose),
            6,
            1,
        )
        .unwrap();
        assert_elimination(&finding, &grid, &solution, "UR Type 3");
        assert!(
            matches!(finding.inference,InferenceResult::Elimination{cell,ref values} if cell==transform(1,transpose) && values==&vec![5,7])
        );
    }
}

#[test]
fn unique_rectangle_type_four_uses_a_real_conjugate_pair() {
    for transpose in [false, true] {
        let (grid, solution) = fixture(
            TYPE4,
            &[
                (9, &[1, 6]),
                (12, &[1, 6]),
                (0, &[1, 5, 6]),
                (3, &[1, 6, 7]),
            ],
            transpose,
        );
        let fab = CandidateFabric::from_grid(&grid);
        let finding = try_ur_hint(
            &fab,
            transform(9, transpose),
            transform(12, transpose),
            transform(0, transpose),
            transform(3, transpose),
            6,
            1,
        )
        .unwrap();
        assert_elimination(&finding, &grid, &solution, "UR Type 4");
    }
}

#[test]
fn extended_rectangle_checks_both_orientations_with_unique_puzzles() {
    for transpose in [false, true] {
        let (grid, solution) = fixture(
            EXTENDED,
            &[
                (0, &[5, 6]),
                (3, &[5, 6]),
                (5, &[5, 6, 8]),
                (9, &[5, 6]),
                (12, &[1, 5, 6]),
                (14, &[5, 6]),
            ],
            transpose,
        );
        let finding = find_extended_unique_rectangle(&CandidateFabric::from_grid(&grid)).unwrap();
        assert_elimination(&finding, &grid, &solution, "Extended UR");
        assert_eq!(finding.involved_cells.len(), 6);
    }
}

#[test]
fn avoidable_rectangle_requires_all_four_corners_to_be_non_given() {
    for transpose in [false, true] {
        let (mut grid, solution) = fixture(
            TYPE1,
            &[(0, &[5, 6]), (3, &[5, 6]), (9, &[5, 6]), (12, &[1, 5, 6])],
            transpose,
        );
        for cell in [0, 3, 9] {
            let pos = idx_to_pos(transform(cell, transpose));
            grid.set_cell_unchecked(pos, solution.get(pos));
        }
        let finding = find_avoidable_rectangle(&CandidateFabric::from_grid(&grid)).unwrap();
        assert_elimination(&finding, &grid, &solution, "Avoidable Rectangle");
        for cell in [0, 3, 9] {
            let mut fixed = grid.deep_clone();
            fixed
                .cell_mut(idx_to_pos(transform(cell, transpose)))
                .set_given(true);
            assert!(find_avoidable_rectangle(&CandidateFabric::from_grid(&fixed)).is_none());
        }
    }
}

#[test]
fn absent_or_incomplete_uniqueness_patterns_do_not_produce_hints() {
    let solved = Grid::from_string(SOLUTION).unwrap();
    let fab = CandidateFabric::from_grid(&solved);
    assert!(find_avoidable_rectangle(&fab).is_none());
    assert!(find_unique_rectangle(&fab).is_none());
    assert!(find_hidden_rectangle(&fab).is_none());
    assert!(find_extended_unique_rectangle(&fab).is_none());
    assert!(find_bug(&fab).is_none());
    assert!(try_ur_hint(&fab, 0, 3, 9, 12, 5, 6).is_none());
    let (mut grid, _) = fixture(
        TYPE1,
        &[(0, &[5, 6]), (3, &[5, 6]), (9, &[5, 6]), (12, &[1, 5, 6])],
        false,
    );
    grid.cell_mut(Position::new(1, 3)).remove_candidate(5);
    assert!(try_ur_hint(&CandidateFabric::from_grid(&grid), 0, 3, 9, 12, 5, 6).is_none());
    assert!(!sees(0, 0));
    assert!(combinations(&[1, 2], 0).is_empty());
    assert!(combinations(&[1, 2], 3).is_empty());
    assert_eq!(
        combinations(&[1, 2, 3], 2),
        vec![vec![1, 2], vec![1, 3], vec![2, 3]]
    );
}

#[test]
fn arbitrary_bivalue_reductions_cannot_justify_an_unsound_bug_hint() {
    let puzzle =
        "800000000003600000070090200050007000000045700000100030001000068008500010090000400";
    let mut original = Grid::from_string(puzzle).unwrap();
    let solver = Solver::new();
    assert!(solver.has_unique_solution(&original));
    let solution = solver.solve(&original).unwrap();
    loop {
        let singles: Vec<_> = original
            .empty_positions()
            .into_iter()
            .filter(|&p| original.get_candidates(p).count() == 1)
            .collect();
        if singles.is_empty() {
            break;
        }
        for pos in singles {
            original.set_cell_unchecked(pos, solution.get(pos));
        }
        original.recalculate_candidates();
    }
    for offset in 0..12 {
        let mut grid = original.deep_clone();
        let mut triple = false;
        for pos in original.empty_positions() {
            let correct = solution.get(pos).unwrap();
            let others: Vec<_> = original
                .get_candidates(pos)
                .iter()
                .filter(|&v| v != correct)
                .collect();
            assert!(!others.is_empty());
            let mut values = vec![
                correct,
                others[(pos.row * 9 + pos.col + offset) % others.len()],
            ];
            if !triple && others.len() >= 2 {
                values.push(others[(pos.row * 9 + pos.col + offset + 1) % others.len()]);
                triple = true;
            }
            grid.cell_mut(pos)
                .set_candidates(BitSet::from_slice(&values));
        }
        assert!(triple);
        if let Some(finding) = find_bug(&CandidateFabric::from_grid(&grid)) {
            match finding.inference {
                InferenceResult::Placement { cell, value } => assert_eq!(
                    solution.get(idx_to_pos(cell)),
                    Some(value),
                    "offset={offset}, {finding:?}"
                ),
                InferenceResult::Elimination { cell, ref values } => assert!(
                    !values.contains(&solution.get(idx_to_pos(cell)).unwrap()),
                    "offset={offset}, {finding:?}"
                ),
            }
        }
    }
}

#[test]
fn natural_bug_plus_one_has_a_valid_bivalue_remainder() {
    let puzzle =
        "534678912672195348198342567859761423426853791713924856060500284280400600040286100";
    let grid = Grid::from_string(puzzle).unwrap();
    let solver = Solver::new();
    assert!(solver.has_unique_solution(&grid));
    let solution = solver.solve(&grid).unwrap();
    let fab = CandidateFabric::from_grid(&grid);
    let finding = find_bug(&fab).expect("Actual BUG+1 must remain supported");
    assert!(matches!(
        finding.inference,
        InferenceResult::Placement { cell: 65, value: 7 }
    ));
    assert_eq!(solution.get(idx_to_pos(65)), Some(7));
    assert!(is_bug_remainder(&fab, &[(65, vec![7])]));
    assert!(!is_bug_remainder(&fab, &[(65, vec![1])]));
    assert!(!is_bug_remainder(&fab, &[]));
}

#[test]
fn natural_bug_candidates_are_checked_against_the_unique_solution() {
    let puzzles = [
        "030078912072190040198342567859761020006853791713924050901037080087419030300280179",
        "500678912672195348198342000809061400426003091710904000961537284287419635300006179",
        "504078912000105348108040567859700423426853791713924856901530084080400035345080179",
    ];
    for puzzle in puzzles {
        let grid = Grid::from_string(puzzle).unwrap();
        let solver = Solver::new();
        assert!(solver.has_unique_solution(&grid));
        let solution = solver.solve(&grid).unwrap();
        if let Some(finding) = find_bug(&CandidateFabric::from_grid(&grid)) {
            match &finding.inference {
                InferenceResult::Placement { cell, value } => {
                    assert_eq!(solution.get(idx_to_pos(*cell)), Some(*value))
                }
                InferenceResult::Elimination { cell, values } => {
                    assert!(!values.contains(&solution.get(idx_to_pos(*cell)).unwrap()))
                }
            }
        }
    }
}

#[test]
fn hidden_rectangle_requires_both_conjugate_links_and_an_opposite_bivalue() {
    for transpose in [false, true] {
        let (grid, solution) = fixture(
            TYPE4,
            &[
                (9, &[1, 6]),
                (12, &[1, 6]),
                (0, &[1, 5, 6]),
                (3, &[1, 6, 7]),
            ],
            transpose,
        );
        let finding = find_hidden_rectangle(&CandidateFabric::from_grid(&grid))
            .expect("Valid hidden rectangle");
        assert_elimination(&finding, &grid, &solution, "Hidden Rectangle");
    }
}

#[test]
fn bug_plus_two_eliminates_only_from_peers_of_both_extra_candidates() {
    let puzzle =
        "534678912600195348198342567859761423400853791713924856000500284280400600040286100";
    for transpose in [false, true] {
        let mut chars = vec![b'0'; 81];
        for (cell, value) in puzzle.bytes().enumerate() {
            chars[transform(cell, transpose)] = value;
        }
        let grid = Grid::from_string(std::str::from_utf8(&chars).unwrap()).unwrap();
        let solver = Solver::new();
        assert!(solver.has_unique_solution(&grid));
        let solution = solver.solve(&grid).unwrap();
        let fab = CandidateFabric::from_grid(&grid);
        assert!(is_bug_remainder(
            &fab,
            &[
                (transform(56, transpose), vec![7]),
                (transform(65, transpose), vec![7])
            ]
        ));
        let finding = find_bug(&fab).expect("Real BUG+2 remains supported");
        assert_elimination(&finding, &grid, &solution, "BUG+2");
        let InferenceResult::Elimination { cell, values } = finding.inference else {
            unreachable!()
        };
        assert_eq!(values, vec![7]);
        assert!(sees(cell, transform(56, transpose)) && sees(cell, transform(65, transpose)));
        assert!(![transform(56, transpose), transform(65, transpose)].contains(&cell));
    }
}

#[test]
fn bug_extras_without_a_common_unit_do_not_eliminate_other_candidates() {
    let puzzle =
        "534678912672195348198342567859061023026853091013920856060500284280000600040286100";
    let grid = Grid::from_string(puzzle).unwrap();
    assert!(Solver::new().has_unique_solution(&grid));
    let fab = CandidateFabric::from_grid(&grid);
    assert!(is_bug_remainder(
        &fab,
        &[(54, vec![7]), (65, vec![7]), (68, vec![7]), (72, vec![7])]
    ));
    assert!(find_bug(&fab).is_none());
}
