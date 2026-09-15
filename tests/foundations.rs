use sudoku_core::{
    canonical_puzzle_hash, canonical_puzzle_hash_str, BitSet, Cell, Difficulty, Grid,
    KillerCageConstraint, MoveError, Position, PuzzleId, ValidationResult,
};

const PUZZLE: &str =
    "530070000600195000098000060800060003400803001700020006060000280000419005000080079";
const SOLUTION: &str =
    "534678912672195348198342567859761423426853791713924856961537284287419635345286179";
fn p(row: usize, col: usize) -> Position {
    Position::new(row, col)
}

#[test]
fn bitset_operations_match_set_membership_for_every_standard_mask() {
    for mask in 0..512_u16 {
        let members: Vec<u8> = (1..=9).filter(|v| mask & (1 << (v - 1)) != 0).collect();
        let set = BitSet::from_slice(&members);
        assert_eq!(set.as_raw(), mask << 1);
        assert_eq!(BitSet::from_raw(mask << 1), set);
        assert_eq!(set.to_vec(), members);
        assert_eq!(set.count() as usize, members.len());
        assert_eq!(set.is_empty(), members.is_empty());
        assert_eq!(
            set.single_value(),
            if members.len() == 1 {
                Some(members[0])
            } else {
                None
            }
        );
        let complement = BitSet::all_9().difference(&set);
        assert_eq!(set.union(&complement), BitSet::all_9());
        assert_eq!(set.intersection(&complement), BitSet::empty());
        for value in 1..=9 {
            let mut changed = set;
            changed.toggle(value);
            assert_eq!(changed.contains(value), !set.contains(value));
            changed.toggle(value);
            assert_eq!(changed, set);
            changed.insert(value);
            assert!(changed.contains(value));
            changed.remove(value);
            assert!(!changed.contains(value));
            assert_eq!(
                changed.union(&BitSet::single(value)),
                set.union(&BitSet::single(value))
            );
        }
        let encoded = serde_json::to_string(&set).unwrap();
        assert_eq!(serde_json::from_str::<BitSet>(&encoded).unwrap(), set);
    }
    assert_eq!(
        format!("{}", BitSet::from_slice(&[9, 1, 3, 3])),
        "{1, 3, 9}"
    );
    assert_eq!(format!("{}", BitSet::default()), "{}");
    let mut iter = BitSet::single(15).iter();
    assert_eq!(iter.next(), Some(15));
    assert_eq!(iter.next(), None);
    assert_eq!(iter.next(), None);
}

#[test]
fn bitset_all_supports_the_full_representable_value_range() {
    for limit in 0..=15 {
        assert_eq!(BitSet::all(limit).to_vec(), (1..=limit).collect::<Vec<_>>());
    }
}

#[test]
fn cell_values_givens_and_notes_round_trip_without_aliasing() {
    let mut cell = Cell::default();
    assert!(cell.is_empty());
    assert!(!cell.is_given());
    assert_eq!(cell.candidate_count(), 9);
    cell.set_candidates(BitSet::from_slice(&[2, 4]));
    cell.add_candidate(6);
    cell.remove_candidate(4);
    cell.toggle_candidate(2);
    assert_eq!(cell.candidates(), BitSet::single(6));
    assert!(cell.has_candidate(6));
    cell.set_value(Some(6));
    assert!(cell.is_filled());
    assert_eq!(cell.value(), Some(6));
    assert_eq!(cell.candidate_count(), 0);
    cell.set_given(true);
    cell.clear();
    assert_eq!(cell.value(), Some(6));
    let restored: Cell = serde_json::from_str(&serde_json::to_string(&cell).unwrap()).unwrap();
    assert_eq!(restored, cell);
    let mut copy = cell.clone();
    copy.set_given(false);
    copy.clear();
    assert_eq!(copy, Cell::new_empty());
    assert_eq!(cell.value(), Some(6));
    copy.set_value(None);
    assert!(copy.is_empty());
    assert_eq!(Cell::new_filled(7).value(), Some(7));
    assert!(!Cell::new_filled(7).is_given());
    assert!(Cell::new_given(9).is_given());
}

#[test]
fn positions_cover_the_board_once_and_locate_boxes_and_diagonals() {
    let positions = Position::all_9x9().collect::<Vec<_>>();
    assert_eq!(positions.len(), 81);
    assert_eq!(positions[0], p(0, 0));
    assert_eq!(positions[80], p(8, 8));
    let unique = positions.iter().collect::<std::collections::HashSet<_>>();
    assert_eq!(unique.len(), 81);
    let mut boxes = [0; 9];
    for position in positions {
        boxes[position.box_index()] += 1;
        let origin = position.box_origin();
        assert!(position.row >= origin.row && position.row < origin.row + 3);
        assert!(position.col >= origin.col && position.col < origin.col + 3);
        assert_eq!(
            position.is_on_main_diagonal(9),
            position.row == position.col
        );
        assert_eq!(
            position.is_on_anti_diagonal(9),
            position.row + position.col == 8
        );
        let json = serde_json::to_string(&position).unwrap();
        assert_eq!(serde_json::from_str::<Position>(&json).unwrap(), position);
    }
    assert_eq!(boxes, [9; 9]);
    assert_eq!(p(4, 7).to_string(), "(4, 7)");
    assert!(!p(9, 9).is_on_main_diagonal(9));
    assert!(!p(9, 0).is_on_anti_diagonal(9));
}

#[test]
fn zero_sized_grids_have_no_diagonal_positions() {
    assert!(!p(0, 0).is_on_main_diagonal(0));
    assert!(!p(0, 0).is_on_anti_diagonal(0));
    assert!(!p(usize::MAX, usize::MAX).is_on_anti_diagonal(9));
}

#[test]
fn grid_parsing_canonical_output_and_human_display_preserve_values() {
    let grid = Grid::from_string(PUZZLE).unwrap();
    assert_eq!(grid.to_string_compact(), PUZZLE.replace('0', "."));
    assert_eq!(grid.given_count(), 30);
    assert_eq!(grid.empty_count(), 51);
    assert_eq!(grid.empty_positions().len(), 51);
    assert_eq!(grid.values()[0][0], Some(5));
    assert_eq!(grid.values()[0][2], None);
    assert_eq!(format!("{:?}", grid.variant()), "Classic");
    assert_eq!(
        grid.constraints()
            .iter()
            .map(|c| c.name())
            .collect::<Vec<_>>(),
        ["Row", "Column", "Box"]
    );
    let readable = grid.to_string();
    assert!(readable.starts_with("53. .7. ...\n"));
    assert_eq!(
        Grid::from_string(&readable).unwrap().to_string_compact(),
        grid.to_string_compact()
    );
    let debug = format!("{:?}", grid);
    assert!(debug.starts_with("Grid {\n  53.|.7.|...\n"));
    assert!(debug.contains("---+---+---"));
    for invalid in [
        "".to_string(),
        ".".repeat(80),
        ".".repeat(82),
        format!("{}?", ".".repeat(80)),
        format!("{}é", ".".repeat(80)),
    ] {
        assert!(Grid::from_string(&invalid).is_none());
    }
}

#[test]
fn checked_moves_reject_invalid_input_without_changing_the_grid() {
    let mut grid = Grid::from_string(PUZZLE).unwrap();
    let original = grid.to_string_compact();
    for (position, value, expected) in [
        (p(9, 0), 1, MoveError::PositionOutOfBounds),
        (p(0, 9), 1, MoveError::PositionOutOfBounds),
        (p(0, 2), 0, MoveError::ValueOutOfRange),
        (p(0, 2), 10, MoveError::ValueOutOfRange),
        (p(0, 0), 1, MoveError::CellIsGiven),
        (p(0, 2), 5, MoveError::ConstraintViolation("Row".into())),
    ] {
        assert!(!grid.is_valid_move(position, value));
        assert_eq!(grid.set_cell(position, value), Err(expected.clone()));
        let encoded = serde_json::to_string(&expected).unwrap();
        assert_eq!(
            serde_json::from_str::<MoveError>(&encoded).unwrap(),
            expected
        );
        assert_eq!(grid.to_string_compact(), original);
    }
    assert_eq!(grid.clear_cell(p(0, 0)), Err(MoveError::CellIsGiven));
    assert_eq!(
        MoveError::CellIsGiven.to_string(),
        "Cannot modify a given cell"
    );
    assert_eq!(
        MoveError::ValueOutOfRange.to_string(),
        "Value must be between 1 and 9"
    );
    assert_eq!(
        MoveError::PositionOutOfBounds.to_string(),
        "Position is out of bounds"
    );
    assert_eq!(
        MoveError::ConstraintViolation("Row".into()).to_string(),
        "Constraint violation: Row"
    );
    assert!(grid.is_valid_move(p(0, 2), 4));
    grid.set_cell(p(0, 2), 4).unwrap();
    assert_eq!(grid.get(p(0, 2)), Some(4));
    assert!(grid.get_candidates(p(0, 2)).is_empty());
    grid.clear_cell(p(0, 2)).unwrap();
    grid.clear_cell(p(0, 2)).unwrap();
    assert_eq!(grid.to_string_compact(), original);
}

#[test]
fn clearing_an_out_of_bounds_position_returns_a_move_error() {
    let mut grid = Grid::new_classic();
    for position in [p(9, 0), p(0, 9), p(usize::MAX, usize::MAX)] {
        assert_eq!(
            grid.clear_cell(position),
            Err(MoveError::PositionOutOfBounds)
        );
    }
}

#[test]
fn clearing_a_move_restores_only_candidates_allowed_by_all_constraints() {
    let mut grid = Grid::new_classic();
    grid.set_given(p(0, 8), 5);
    grid.set_cell(p(4, 4), 5).unwrap();
    grid.clear_cell(p(4, 4)).unwrap();
    // Row 0 still contains a 5 even though clearing (4,4) releases column 4.
    assert!(!grid.get_candidates(p(0, 4)).contains(5));
    assert!(grid.get_candidates(p(8, 4)).contains(5));
    for position in Position::all_9x9() {
        assert_eq!(
            grid.get_candidates(position),
            grid.compute_candidates(position),
            "candidates diverged at {position}"
        );
    }
    let mut puzzle = Grid::from_string(PUZZLE).unwrap();
    puzzle.set_cell(p(0, 2), 4).unwrap();
    puzzle.clear_cell(p(0, 2)).unwrap();
    assert_eq!(puzzle.get_candidates(p(0, 2)).to_vec(), [1, 2, 4]);
}

#[test]
fn replacing_an_entered_value_releases_candidates_for_the_previous_value() {
    let mut grid = Grid::new_classic();
    grid.set_cell(p(0, 0), 1).unwrap();
    grid.set_cell(p(0, 0), 2).unwrap();
    assert!(grid.get_candidates(p(0, 1)).contains(1));
    assert!(!grid.get_candidates(p(0, 1)).contains(2));
    for position in Position::all_9x9() {
        assert_eq!(
            grid.get_candidates(position),
            grid.compute_candidates(position)
        );
    }
}

#[test]
fn candidate_recalculation_and_manual_notes_do_not_change_the_puzzle() {
    let mut grid = Grid::from_string(PUZZLE).unwrap();
    assert_eq!(grid.compute_candidates(p(0, 2)).to_vec(), [1, 2, 4]);
    assert!(grid.compute_candidates(p(0, 0)).is_empty());
    grid.clear_all_candidates();
    assert!(Position::all_9x9().all(|pos| grid.get_candidates(pos).is_empty()));
    grid.cell_mut(p(0, 2)).add_candidate(4);
    assert_eq!(grid.get_candidates(p(0, 2)), BitSet::single(4));
    grid.recalculate_candidates();
    for pos in Position::all_9x9() {
        assert_eq!(grid.get_candidates(pos), grid.compute_candidates(pos));
    }
    assert_eq!(grid.to_string_compact(), PUZZLE.replace('0', "."));
}

#[test]
fn validation_identifies_conflicting_cells_and_complete_solutions() {
    let valid = ValidationResult::valid();
    assert!(valid.is_valid && valid.invalid_cells.is_empty() && valid.violations.is_empty());
    let mut grid = Grid::from_string(SOLUTION).unwrap();
    assert!(grid.is_complete());
    assert!(grid.is_solved());
    assert_eq!(grid.validate(), valid);
    grid.set_cell_unchecked(p(0, 0), Some(3));
    let invalid = grid.validate();
    assert!(!invalid.is_valid);
    assert!(!grid.is_complete());
    assert!(invalid.invalid_cells.contains(&p(0, 0)));
    assert!(invalid.invalid_cells.contains(&p(0, 1)));
    assert!(invalid
        .violations
        .contains(&"Row constraint violated at (0, 0)".to_string()));
    assert_eq!(
        invalid
            .invalid_cells
            .iter()
            .collect::<std::collections::HashSet<_>>()
            .len(),
        invalid.invalid_cells.len()
    );
    let round_trip: ValidationResult =
        serde_json::from_str(&serde_json::to_string(&invalid).unwrap()).unwrap();
    assert_eq!(round_trip, invalid);
    grid.set_cell_unchecked(p(0, 0), None);
    assert!(!grid.is_complete());
}

#[test]
fn classic_diagonal_and_killer_rules_survive_clone_and_serialization_restore() {
    let grids = [
        Grid::new_classic(),
        Grid::new_x_sudoku(),
        Grid::new_killer(vec![KillerCageConstraint::new(vec![p(0, 0), p(4, 4)], 10)]),
    ];
    for (index, mut original) in grids.into_iter().enumerate() {
        original.set_given(p(0, 0), 3);
        original
            .cell_mut(p(8, 7))
            .set_candidates(BitSet::from_slice(&[1, 9]));
        let mut cloned = original.clone();
        let mut restored: Grid =
            serde_json::from_str(&serde_json::to_string(&original).unwrap()).unwrap();
        restored.restore_constraints();
        for copy in [&mut cloned, &mut restored] {
            assert_eq!(copy.to_string_compact(), original.to_string_compact());
            assert_eq!(copy.get_candidates(p(8, 7)).to_vec(), [1, 9]);
            assert_eq!(copy.constraints().len(), if index == 0 { 3 } else { 4 });
            match index {
                0 => assert!(copy.is_valid_move(p(4, 4), 3)),
                1 => assert!(!copy.is_valid_move(p(4, 4), 3)),
                _ => {
                    assert!(!copy.is_valid_move(p(4, 4), 6));
                    assert!(copy.is_valid_move(p(4, 4), 7));
                }
            }
            copy.set_cell(p(8, 7), 1).unwrap();
        }
        assert_eq!(original.get(p(8, 7)), None);
        assert_eq!(original.get_candidates(p(8, 7)).to_vec(), [1, 9]);
    }
}

#[test]
fn short_codes_round_trip_boundaries_all_difficulties_and_case() {
    let maximum = 36_u64.pow(7) - 1;
    for (difficulty, prefix) in [
        (Difficulty::Beginner, 'B'),
        (Difficulty::Easy, 'E'),
        (Difficulty::Medium, 'M'),
        (Difficulty::Intermediate, 'I'),
        (Difficulty::Hard, 'H'),
        (Difficulty::Expert, 'X'),
        (Difficulty::Master, 'S'),
        (Difficulty::Extreme, 'Z'),
    ] {
        for seed in [0, 1, 35, 36, 1295, 1296, maximum] {
            let id = PuzzleId { difficulty, seed };
            let code = id.to_short_code();
            assert_eq!(code.len(), 8);
            assert!(code.starts_with(prefix));
            assert_eq!(PuzzleId::from_short_code(&code), Some(id.clone()));
            assert_eq!(
                PuzzleId::from_short_code(&format!(" \n{}\t", code.to_ascii_lowercase())),
                Some(id.clone())
            );
            assert_eq!(id.to_string(), code);
        }
        for _ in 0..8 {
            let random = PuzzleId::random(difficulty);
            assert_eq!(random.difficulty, difficulty);
            assert!(random.seed <= maximum);
            assert_eq!(
                PuzzleId::from_short_code(&random.to_short_code()),
                Some(random)
            );
        }
    }
    assert_eq!(
        PuzzleId {
            difficulty: Difficulty::Extreme,
            seed: maximum
        }
        .to_short_code(),
        "ZZZZZZZZ"
    );
    for invalid in [
        "",
        "H",
        "H000000",
        "H00000000",
        "Q0000000",
        "H000000!",
        "M00000é",
        "é000000",
        "😀0000",
    ] {
        assert!(
            PuzzleId::from_short_code(invalid).is_none(),
            "accepted invalid ID {invalid:?}"
        );
    }
    assert_eq!(PuzzleId::from_short_code("M0000010").unwrap().seed, 36);
}

#[test]
fn canonical_hashes_are_pinned_and_normalize_grid_inputs() {
    // Independent SHA-256 published test vector and pinned Sudoku bytes.
    assert_eq!(
        canonical_puzzle_hash_str("abc"),
        "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"
    );
    let empty = Grid::new_classic();
    assert_eq!(
        canonical_puzzle_hash(&empty),
        "eb7b0297e5b6b33f6680cf9c6b471903832a5caceda3debf8e52cda0a40ab9c3"
    );
    let zeros = Grid::from_string(PUZZLE).unwrap();
    let dots = Grid::from_string(&PUZZLE.replace('0', ".")).unwrap();
    assert_eq!(canonical_puzzle_hash(&zeros), canonical_puzzle_hash(&dots));
    assert_eq!(
        canonical_puzzle_hash(&zeros),
        canonical_puzzle_hash_str(&dots.to_string_compact())
    );
    assert_ne!(
        canonical_puzzle_hash_str(PUZZLE),
        canonical_puzzle_hash(&zeros)
    );
}
