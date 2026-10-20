use super::*;
use crate::grid::GridVariant;
use crate::{Constraint, KillerCageConstraint, Solver};

const SOLUTION: &str =
    "534678912672195348198342567859761423426853791713924856961537284287419635345286179";

#[test]
fn classic_counts_match_general_search_across_clue_densities_and_caps() {
    let solver = Solver::new();
    // Independent general-purpose search remains an oracle for the mask-based
    // implementation, including branch restoration and early solution caps.
    for seed in 0..16_u64 {
        let mut state = seed + 1;
        let mut indices: [usize; 81] = std::array::from_fn(|i| i);
        for i in (1..81).rev() {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            indices.swap(i, (state >> 32) as usize % (i + 1));
        }
        for removed in [0, 4, 20, 40, 50, 60] {
            let mut grid = Grid::from_string(SOLUTION).unwrap();
            for &index in &indices[..removed] {
                grid.set_cell_unchecked(Position::new(index / 9, index % 9), None);
            }
            for limit in [0, 1, 2, 3] {
                let mut reference = grid.deep_clone();
                reference.recalculate_candidates();
                let mut expected = 0;
                count_solutions_recursive(&mut reference, &mut expected, limit);
                assert_eq!(
                    solver.count_solutions(&grid, limit),
                    expected,
                    "seed={seed}, removed={removed}, limit={limit}"
                );
            }
        }
    }
}

#[test]
fn classic_counts_reject_conflicting_and_impossible_values() {
    let solver = Solver::new();
    for peer in [
        Position::new(0, 4),
        Position::new(4, 0),
        Position::new(1, 1),
    ] {
        let mut grid = Grid::new_classic();
        grid.set_given(Position::new(0, 0), 5);
        grid.set_given(peer, 5);
        assert_eq!(solver.count_solutions(&grid, 2), 0);
    }
    for value in [0, 10, 255] {
        let mut grid = Grid::new_classic();
        grid.set_cell_unchecked(Position::new(0, 0), Some(value));
        assert_eq!(solver.count_solutions(&grid, 2), 0);
    }

    let mut grid = Grid::new_classic();
    for col in 1..9 {
        grid.set_given(Position::new(0, col), col as u8);
    }
    grid.set_given(Position::new(1, 0), 9);
    assert!(grid.validate().is_valid);
    assert_eq!(solver.count_solutions(&grid, 2), 0);
}

#[test]
fn counting_uses_all_values_ignores_pencil_marks_and_preserves_cell_state() {
    let mut grid = Grid::from_string(SOLUTION).unwrap();
    for index in [3, 4, 30, 31] {
        grid.set_cell_unchecked(Position::new(index / 9, index % 9), None);
    }
    grid.clear_all_candidates();
    for pos in Position::all_9x9() {
        grid.cell_mut(pos).set_given(false);
    }
    let before = serde_json::to_value(&grid).unwrap();
    let solver = Solver::new();
    assert_eq!(solver.count_solutions(&grid, 1), 1);
    assert_eq!(solver.count_solutions(&grid, 10), 2);
    assert_eq!(serde_json::to_value(&grid).unwrap(), before);
    assert_eq!(solver.count_solutions(&Grid::new_classic(), 2), 2);
}

struct PermissiveRow;

impl Constraint for PermissiveRow {
    fn name(&self) -> &'static str {
        "Row"
    }

    fn validate(&self, _: &[[Option<u8>; 9]; 9], _: Position, _: u8) -> bool {
        true
    }

    fn affected_cells(&self, _: Position) -> Vec<Position> {
        Vec::new()
    }
}

#[test]
fn variant_custom_and_unrestored_rules_never_use_classic_counting() {
    let solver = Solver::new();
    let mut diagonal = Grid::new_x_sudoku();
    let mut killer = Grid::new_killer(vec![KillerCageConstraint::new(
        vec![Position::new(0, 0)],
        1,
    )]);
    for grid in [&mut diagonal, &mut killer] {
        for (index, digit) in SOLUTION.bytes().enumerate() {
            grid.set_given(Position::new(index / 9, index % 9), digit - b'0');
        }
        // A complete classic solution violates these additional rules.
        assert!(!grid.has_classic_sudoku_constraints());
        assert_eq!(solver.count_solutions(grid, 2), 0);
    }

    // A custom constraint can share a built-in name and a Classic variant tag
    // while accepting boards that would be invalid under classic rules.
    let mut custom = Grid::new_with_constraints(vec![Box::new(PermissiveRow)]);
    for pos in Position::all_9x9() {
        custom.set_given(pos, 1);
    }
    assert_eq!(custom.variant(), GridVariant::Classic);
    for grid in [&custom, &custom.clone()] {
        assert!(!grid.has_classic_sudoku_constraints());
        assert_eq!(solver.count_solutions(grid, 2), 1);
    }

    let mut restored: Grid =
        serde_json::from_value(serde_json::to_value(&custom).unwrap()).unwrap();
    assert!(!restored.has_classic_sudoku_constraints());
    assert_eq!(solver.count_solutions(&restored, 2), 1);
    restored.restore_constraints();
    assert!(restored.has_classic_sudoku_constraints());
    assert_eq!(solver.count_solutions(&restored, 2), 0);
}
