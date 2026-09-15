use std::collections::HashSet;
use sudoku_core::{
    BoxConstraint, ColumnConstraint, Constraint, DiagonalConstraint, Grid, KillerCageConstraint,
    Position, RowConstraint, ThermoConstraint,
};
fn p(row: usize, col: usize) -> Position {
    Position::new(row, col)
}

#[test]
fn classic_constraints_report_exact_peers_and_allow_replacing_the_same_cell() {
    let position = p(4, 4);
    let cases: Vec<(Box<dyn Constraint>, &str, Vec<Position>, Position)> = vec![
        (
            Box::new(RowConstraint),
            "Row",
            (0..9)
                .filter(|&col| col != 4)
                .map(|col| p(4, col))
                .collect(),
            p(0, 0),
        ),
        (
            Box::new(ColumnConstraint),
            "Column",
            (0..9)
                .filter(|&row| row != 4)
                .map(|row| p(row, 4))
                .collect(),
            p(0, 0),
        ),
        (
            Box::new(BoxConstraint),
            "Box",
            (3..6)
                .flat_map(|row| (3..6).map(move |col| p(row, col)))
                .filter(|&pos| pos != position)
                .collect(),
            p(0, 0),
        ),
    ];
    for (constraint, name, peers, outside) in cases {
        assert_eq!(constraint.name(), name);
        assert_eq!(constraint.affected_cells(position), peers);
        assert_eq!(peers.iter().collect::<HashSet<_>>().len(), 8);
        let mut cells = [[None; 9]; 9];
        cells[position.row][position.col] = Some(4);
        cells[outside.row][outside.col] = Some(4);
        assert!(constraint.validate(&cells, position, 4));
        for peer in peers {
            cells[peer.row][peer.col] = Some(4);
            assert!(
                !constraint.validate(&cells, position, 4),
                "{name} ignored peer {peer}"
            );
            assert!(constraint.validate(&cells, position, 5));
            cells[peer.row][peer.col] = None;
        }
    }
}

#[test]
fn diagonal_constraint_handles_center_corners_and_non_diagonal_positions() {
    let constraint = DiagonalConstraint;
    assert_eq!(constraint.name(), "Diagonal");
    assert!(constraint.affected_cells(p(0, 1)).is_empty());
    assert_eq!(
        constraint.affected_cells(p(0, 0)),
        (1..9).map(|i| p(i, i)).collect::<Vec<_>>()
    );
    assert_eq!(
        constraint.affected_cells(p(0, 8)),
        (1..9).map(|i| p(i, 8 - i)).collect::<Vec<_>>()
    );
    let peers = constraint.affected_cells(p(4, 4));
    assert_eq!(peers.len(), 16);
    assert_eq!(peers.iter().collect::<HashSet<_>>().len(), 16);
    let mut cells = [[None; 9]; 9];
    cells[0][0] = Some(3);
    cells[0][8] = Some(7);
    assert!(!constraint.validate(&cells, p(4, 4), 3));
    assert!(!constraint.validate(&cells, p(4, 4), 7));
    assert!(constraint.validate(&cells, p(4, 4), 5));
    assert!(constraint.validate(&cells, p(0, 1), 3));
    cells[4][4] = Some(5);
    assert!(constraint.validate(&cells, p(4, 4), 5));
}

#[test]
fn killer_cages_enforce_distinct_values_partial_sums_and_exact_final_sum() {
    let positions = vec![p(0, 0), p(3, 3), p(6, 6)];
    let cage = KillerCageConstraint::new(positions.clone(), 15);
    assert_eq!(cage.name(), "KillerCage");
    assert!(cage.contains(p(3, 3)));
    assert!(!cage.contains(p(0, 1)));
    assert_eq!(cage.affected_cells(p(3, 3)), vec![p(0, 0), p(6, 6)]);
    assert!(cage.affected_cells(p(0, 1)).is_empty());
    let mut cells = [[None; 9]; 9];
    cells[0][0] = Some(9);
    assert!(cage.validate(&cells, p(3, 3), 3));
    assert!(!cage.validate(&cells, p(3, 3), 7));
    assert!(!cage.validate(&cells, p(3, 3), 9));
    assert!(cage.validate(&cells, p(0, 1), 9));
    cells[0][0] = Some(3);
    cells[6][6] = Some(5);
    assert!(cage.validate(&cells, p(3, 3), 7));
    assert!(!cage.validate(&cells, p(3, 3), 6));
    assert!(!cage.validate(&cells, p(3, 3), 8));
    let clone = cage.clone();
    let round_trip: KillerCageConstraint =
        serde_json::from_str(&serde_json::to_string(&cage).unwrap()).unwrap();
    assert_eq!(round_trip.cells, positions);
    assert_eq!(round_trip.sum, 15);
    assert!(clone.validate(&cells, p(3, 3), 7));
    let single = KillerCageConstraint::new(vec![p(0, 0)], 4);
    assert!(single.validate(&[[None; 9]; 9], p(0, 0), 4));
    assert!(!single.validate(&[[None; 9]; 9], p(0, 0), 5));
}

#[test]
fn thermometers_require_strict_order_in_both_directions_and_ignore_other_cells() {
    let positions = vec![p(0, 0), p(3, 3), p(6, 6), p(8, 8)];
    let thermo = ThermoConstraint::new(positions.clone());
    assert_eq!(thermo.name(), "Thermo");
    assert_eq!(
        thermo.affected_cells(p(3, 3)),
        vec![p(0, 0), p(6, 6), p(8, 8)]
    );
    assert!(thermo.affected_cells(p(0, 1)).is_empty());
    let mut cells = [[None; 9]; 9];
    assert!(thermo.validate(&cells, p(3, 3), 4));
    cells[0][0] = Some(2);
    cells[8][8] = Some(8);
    assert!(thermo.validate(&cells, p(3, 3), 4));
    for invalid in [1, 2, 8, 9] {
        assert!(!thermo.validate(&cells, p(3, 3), invalid));
    }
    assert!(thermo.validate(&cells, p(0, 1), 1));
    cells[6][6] = Some(6);
    assert!(!thermo.validate(&cells, p(3, 3), 6));
    assert!(!thermo.validate(&cells, p(8, 8), 6));
    let round_trip: ThermoConstraint =
        serde_json::from_str(&serde_json::to_string(&thermo).unwrap()).unwrap();
    assert_eq!(round_trip.path, positions);
    assert!(round_trip.validate(&cells, p(3, 3), 4));
    assert!(thermo.clone().validate(&cells, p(3, 3), 4));
}

#[test]
fn custom_constraint_grids_enforce_rules_and_clone_them_without_changing_cells() {
    let mut custom = Grid::new_with_constraints(vec![
        Box::new(DiagonalConstraint),
        Box::new(ThermoConstraint::new(vec![p(0, 1), p(4, 4), p(8, 7)])),
    ]);
    custom.set_given(p(0, 0), 5);
    custom.set_given(p(0, 1), 3);
    assert!(!custom.is_valid_move(p(4, 4), 5));
    assert!(!custom.is_valid_move(p(4, 4), 2));
    assert!(custom.is_valid_move(p(4, 4), 4));
    // A custom grid has only its requested rules, so equal values in a row may be valid.
    assert!(custom.is_valid_move(p(0, 2), 5));
    for mut copy in [custom.clone(), custom.deep_clone()] {
        assert_eq!(
            copy.constraints()
                .iter()
                .map(|c| c.name())
                .collect::<Vec<_>>(),
            ["Diagonal", "Thermo"]
        );
        assert!(!copy.is_valid_move(p(4, 4), 5));
        assert!(!copy.is_valid_move(p(4, 4), 2));
        assert!(copy.is_valid_move(p(0, 2), 5));
        copy.set_cell(p(4, 4), 4).unwrap();
    }
    assert_eq!(custom.get(p(4, 4)), None);
}

#[test]
fn solving_a_custom_thermometer_keeps_the_rule_in_every_working_copy() {
    let solution =
        "534678912672195348198342567859761423426853791713924856961537284287419635345286179";
    for (path, satisfiable) in [
        (vec![p(0, 1), p(0, 2), p(0, 0)], true),
        (vec![p(0, 0), p(0, 2), p(0, 1)], false),
    ] {
        let mut grid = Grid::new_with_constraints(vec![
            Box::new(RowConstraint),
            Box::new(ColumnConstraint),
            Box::new(BoxConstraint),
            Box::new(ThermoConstraint::new(path)),
        ]);
        for (index, value) in solution.bytes().enumerate() {
            if index != 2 {
                grid.set_given(p(index / 9, index % 9), value - b'0');
            }
        }
        let result = sudoku_core::Solver::new().solve(&grid);
        if satisfiable {
            let solved = result.expect("the standard solution obeys this thermometer");
            assert_eq!(solved.to_string_compact(), solution);
            assert_eq!(solved.constraints().last().unwrap().name(), "Thermo");
            assert!(solved.validate().is_valid);
        } else {
            assert!(
                result.is_none(),
                "the reversed thermometer must not disappear during cloning"
            );
        }
        assert_eq!(grid.empty_count(), 1);
    }
}
