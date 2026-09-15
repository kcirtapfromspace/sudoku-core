use super::*;

#[test]
fn test_new_classic() {
    let grid = Grid::new_classic();
    assert_eq!(grid.empty_count(), 81);
    assert_eq!(grid.given_count(), 0);
}

#[test]
fn test_set_cell() {
    let mut grid = Grid::new_classic();
    let pos = Position::new(0, 0);

    assert!(grid.set_cell(pos, 5).is_ok());
    assert_eq!(grid.get(pos), Some(5));

    // Can't set duplicate in same row
    assert!(grid.set_cell(Position::new(0, 5), 5).is_err());
}

#[test]
fn test_given_cell() {
    let mut grid = Grid::new_classic();
    let pos = Position::new(0, 0);

    grid.set_given(pos, 5);
    assert!(grid.cell(pos).is_given());
    assert!(grid.set_cell(pos, 3).is_err());
}

#[test]
fn test_from_string() {
    let puzzle =
        "530070000600195000098000060800060003400803001700020006060000280000419005000080079";
    let grid = Grid::from_string(puzzle).unwrap();

    assert_eq!(grid.get(Position::new(0, 0)), Some(5));
    assert_eq!(grid.get(Position::new(0, 1)), Some(3));
    assert_eq!(grid.get(Position::new(0, 2)), None);
}

#[test]
fn test_candidates() {
    let mut grid = Grid::new_classic();
    grid.set_given(Position::new(0, 0), 5);

    // Cell in same row should not have 5 as candidate
    let candidates = grid.get_candidates(Position::new(0, 5));
    assert!(!candidates.contains(5));
    assert!(candidates.contains(3));
}

#[test]
fn test_is_complete() {
    let solved =
        "534678912672195348198342567859761423426853791713924856961537284287419635345286179";
    let grid = Grid::from_string(solved).unwrap();
    assert!(grid.is_complete());
}
