use super::*;

fn empty_grid() -> [[Option<u8>; 9]; 9] {
    [[None; 9]; 9]
}

#[test]
fn test_row_constraint() {
    let mut grid = empty_grid();
    grid[0][0] = Some(5);

    let constraint = RowConstraint;

    // Same value in same row should fail
    assert!(!constraint.validate(&grid, Position::new(0, 5), 5));

    // Same value in different row should pass
    assert!(constraint.validate(&grid, Position::new(1, 5), 5));

    // Different value in same row should pass
    assert!(constraint.validate(&grid, Position::new(0, 5), 3));
}

#[test]
fn test_column_constraint() {
    let mut grid = empty_grid();
    grid[0][0] = Some(5);

    let constraint = ColumnConstraint;

    // Same value in same column should fail
    assert!(!constraint.validate(&grid, Position::new(5, 0), 5));

    // Same value in different column should pass
    assert!(constraint.validate(&grid, Position::new(5, 1), 5));
}

#[test]
fn test_box_constraint() {
    let mut grid = empty_grid();
    grid[0][0] = Some(5);

    let constraint = BoxConstraint;

    // Same value in same box should fail
    assert!(!constraint.validate(&grid, Position::new(2, 2), 5));

    // Same value in different box should pass
    assert!(constraint.validate(&grid, Position::new(3, 3), 5));
}

#[test]
fn test_diagonal_constraint() {
    let mut grid = empty_grid();
    grid[0][0] = Some(5);
    grid[0][8] = Some(3);

    let constraint = DiagonalConstraint;

    // Same value on main diagonal should fail
    assert!(!constraint.validate(&grid, Position::new(4, 4), 5));

    // Same value on anti-diagonal should fail
    assert!(!constraint.validate(&grid, Position::new(4, 4), 3));

    // Different value on diagonals should pass
    assert!(constraint.validate(&grid, Position::new(4, 4), 7));

    // Cell not on diagonal shouldn't care about diagonal values
    assert!(constraint.validate(&grid, Position::new(0, 1), 5));
}

#[test]
fn test_killer_cage_constraint() {
    let mut grid = empty_grid();
    let cage = KillerCageConstraint::new(
        vec![
            Position::new(0, 0),
            Position::new(0, 1),
            Position::new(1, 0),
        ],
        15, // Target sum is 15
    );

    grid[0][0] = Some(3);
    grid[0][1] = Some(5);

    // Sum would be 3 + 5 + 7 = 15, and all unique
    assert!(cage.validate(&grid, Position::new(1, 0), 7));

    // Sum would be 3 + 5 + 8 = 16 > 15
    assert!(!cage.validate(&grid, Position::new(1, 0), 8));

    // Duplicate value in cage (3 already exists)
    assert!(!cage.validate(&grid, Position::new(1, 0), 3));
}

#[test]
fn test_thermo_constraint() {
    let mut grid = empty_grid();
    let thermo = ThermoConstraint::new(vec![
        Position::new(0, 0),
        Position::new(0, 1),
        Position::new(0, 2),
    ]);

    grid[0][0] = Some(2);
    grid[0][2] = Some(6);

    // Value must be > 2 and < 6
    assert!(thermo.validate(&grid, Position::new(0, 1), 4));
    assert!(!thermo.validate(&grid, Position::new(0, 1), 2));
    assert!(!thermo.validate(&grid, Position::new(0, 1), 6));
    assert!(!thermo.validate(&grid, Position::new(0, 1), 1));
    assert!(!thermo.validate(&grid, Position::new(0, 1), 7));
}
