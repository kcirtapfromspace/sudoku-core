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

struct NamedNoopConstraint(&'static str);

impl crate::Constraint for NamedNoopConstraint {
    fn validate(&self, _: &[[Option<u8>; 9]; 9], _: Position, _: u8) -> bool {
        true
    }

    fn affected_cells(&self, _: Position) -> Vec<Position> {
        Vec::new()
    }

    fn name(&self) -> &'static str {
        self.0
    }
}

#[test]
fn built_in_grids_and_clones_have_standard_constraints() {
    for grid in [
        Grid::new_classic(),
        Grid::new_x_sudoku(),
        Grid::new_killer(vec![KillerCageConstraint::new(
            vec![Position::new(0, 0), Position::new(0, 1)],
            3,
        )]),
    ] {
        assert!(grid.has_standard_sudoku_constraints());
        assert!(grid.clone().has_standard_sudoku_constraints());
        assert!(grid.deep_clone().has_standard_sudoku_constraints());
    }
}

#[test]
fn custom_grids_cannot_claim_standard_constraints_by_variant_or_rule_name() {
    let cases: Vec<Vec<ConstraintBox>> = vec![
        Vec::new(),
        classic_constraints(),
        vec![
            Box::new(NamedNoopConstraint("Row")),
            Box::new(NamedNoopConstraint("Column")),
            Box::new(NamedNoopConstraint("Box")),
        ],
    ];
    for constraints in cases {
        let mut grid = Grid::new_with_constraints(constraints);
        assert_eq!(grid.variant(), GridVariant::Classic);
        assert!(!grid.has_standard_sudoku_constraints());
        assert!(!grid.clone().has_standard_sudoku_constraints());
        assert!(!grid.deep_clone().has_standard_sudoku_constraints());

        // Restoration actually replaces custom rules with the built-in rules.
        grid.restore_constraints();
        assert!(grid.has_standard_sudoku_constraints());
        grid.set_given(Position::new(0, 0), 5);
        assert!(!grid.is_valid_move(Position::new(0, 1), 5));
    }
}

#[test]
fn deserialized_grids_and_clones_require_constraint_restoration() {
    for mut original in [
        Grid::new_classic(),
        Grid::new_x_sudoku(),
        Grid::new_killer(vec![KillerCageConstraint::new(
            vec![Position::new(0, 0), Position::new(0, 1)],
            3,
        )]),
    ] {
        original.set_given(Position::new(8, 0), 5);
        let wire = serde_json::to_value(&original).unwrap();
        let fields = wire.as_object().unwrap();
        assert_eq!(fields.len(), 3);
        assert!(fields.contains_key("cells"));
        assert!(fields.contains_key("variant"));
        assert!(fields.contains_key("killer_cages"));

        let mut restored: Grid = serde_json::from_value(wire).unwrap();
        assert!(restored.constraints().is_empty());
        assert!(!restored.has_standard_sudoku_constraints());
        assert!(!restored.clone().has_standard_sudoku_constraints());
        assert!(!restored.deep_clone().has_standard_sudoku_constraints());
        assert!(restored.is_valid_move(Position::new(8, 1), 5));

        restored.restore_constraints();
        assert!(restored.has_standard_sudoku_constraints());
        assert!(restored.clone().has_standard_sudoku_constraints());
        assert!(restored.deep_clone().has_standard_sudoku_constraints());
        assert!(!restored.is_valid_move(Position::new(8, 1), 5));
        assert_eq!(restored.constraints().len(), original.constraints().len());
    }
}

#[test]
fn serialized_input_cannot_set_constraint_provenance() {
    let mut wire = serde_json::to_value(Grid::new_classic()).unwrap();
    wire["standard_sudoku_constraints"] = serde_json::json!(true);
    let grid: Grid = serde_json::from_value(wire).unwrap();
    assert!(grid.constraints().is_empty());
    assert!(!grid.has_standard_sudoku_constraints());
}
