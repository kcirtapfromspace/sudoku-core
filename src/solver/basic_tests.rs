use super::*;
use crate::Grid;

#[test]
fn test_naked_single() {
    let puzzle =
        "530070000600195000098000060800060003400803001700020006060000280000419005000080079";
    let grid = Grid::from_string(puzzle).unwrap();
    let fab = CandidateFabric::from_grid(&grid);
    let finding = find_naked_single(&fab);
    assert!(finding.is_some());
    let f = finding.unwrap();
    assert_eq!(f.technique, Technique::NakedSingle);
}

#[test]
fn test_hidden_single() {
    // A uniquely solvable puzzle requiring hidden singles.
    let puzzle =
        "000000010400000000020000000000050407008000300001090000300400200050100000000806000";
    let grid = Grid::from_string(puzzle).unwrap();
    let solver = crate::Solver::new();
    assert!(solver.has_unique_solution(&grid));
    let fab = CandidateFabric::from_grid(&grid);
    let finding = find_hidden_single(&fab).expect("the fixture contains a hidden single");
    assert_eq!(finding.technique, Technique::HiddenSingle);
    let InferenceResult::Placement { cell, value } = finding.inference else {
        panic!("a hidden single must place a value");
    };
    let solution = solver.solve(&grid).unwrap();
    assert_eq!(
        solution.get(crate::Position::new(cell / 9, cell % 9)),
        Some(value)
    );
}

#[test]
fn test_combinations() {
    let items = vec![1usize, 2, 3, 4];
    let combos = combinations_idx(&items, 2);
    assert_eq!(combos.len(), 6);
    assert!(combos.contains(&vec![1, 2]));
    assert!(combos.contains(&vec![3, 4]));
}
