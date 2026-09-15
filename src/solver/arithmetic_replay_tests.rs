use super::*;
use crate::{ArithmeticRequirement, ArithmeticTerm};

fn naked_single() -> (Grid, ArithmeticProof) {
    let mut grid = Grid::new_classic();
    grid.cell_mut(position(0)).set_candidates(BitSet::single(1));
    grid.cell_mut(position(80)).remove_candidate(9);
    let proof = ArithmeticProof::from_terms(
        &grid,
        vec![ArithmeticTerm {
            requirement: ArithmeticRequirement::Cell { cell: 0 },
            weight: 1,
        }],
        0,
        1,
        true,
    )
    .unwrap();
    (grid, proof)
}

#[test]
fn round_trip_reinstalls_native_constraints_and_preserves_all_masks() {
    let (grid, proof) = naked_single();
    let replay = ArithmeticReplay::capture(&grid, &proof).unwrap();
    let json = serde_json::to_string(&replay).unwrap();
    let restored: ArithmeticReplay = serde_json::from_str(&json).unwrap();
    assert!(restored.verify());
    let rebuilt = restored.state.to_grid().unwrap();
    assert!(rebuilt.has_standard_sudoku_constraints());
    assert_eq!(ArithmeticState::capture(&rebuilt).unwrap(), replay.state);
    assert!(!rebuilt.get_candidates(position(80)).contains(9));
}

#[test]
fn application_preserves_prior_eliminations_and_cannot_replay_twice() {
    let (mut grid, proof) = naked_single();
    let replay = ArithmeticReplay::capture(&grid, &proof).unwrap();
    replay.apply(&mut grid).unwrap();
    assert_eq!(grid.get(position(0)), Some(1));
    assert!(!grid.get_candidates(position(1)).contains(1));
    assert!(!grid.get_candidates(position(80)).contains(9));
    let after = ArithmeticState::capture(&grid).unwrap();
    assert!(replay.apply(&mut grid).is_err());
    assert_eq!(ArithmeticState::capture(&grid).unwrap(), after);
}

#[test]
fn changed_premises_or_certificate_fail_without_mutation() {
    let (mut grid, proof) = naked_single();
    let replay = ArithmeticReplay::capture(&grid, &proof).unwrap();
    let mut tampered = replay.clone();
    tampered.state.domains[80] = ALL;
    assert!(!tampered.verify());
    tampered = replay.clone();
    tampered.proof.value = false;
    assert!(!tampered.verify());
    grid.cell_mut(position(79)).remove_candidate(8);
    let before = ArithmeticState::capture(&grid).unwrap();
    assert!(replay.apply(&mut grid).is_err());
    assert_eq!(ArithmeticState::capture(&grid).unwrap(), before);
}

#[test]
fn rejects_malformed_states_and_uninstalled_constraints() {
    let (grid, proof) = naked_single();
    let replay = ArithmeticReplay::capture(&grid, &proof).unwrap();
    for (cell, domain) in [(0, 0), (1, 512)] {
        let mut state = replay.state.clone();
        state.domains[cell] = domain;
        assert!(state.to_grid().is_err());
    }
    let mut state = replay.state.clone();
    state.values[1] = 10;
    assert!(state.to_grid().is_err());
    state = replay.state.clone();
    state.values[1] = 1;
    assert!(state.to_grid().is_err()); // not a singleton domain
    state = replay.state.clone();
    state.values[1] = 1;
    state.domains[1] = 1;
    state.values[2] = 1;
    state.domains[2] = 1;
    assert!(state.to_grid().is_err()); // duplicate placed digit
    state = replay.state.clone();
    state.values.pop();
    assert!(state.to_grid().is_err());
    state = replay.state.clone();
    state.version = 2;
    assert!(state.to_grid().is_err());
    state = replay.state.clone();
    for domain in &mut state.domains[..9] {
        *domain &= !256;
    }
    assert!(state.to_grid().is_err()); // no 9 in the first row
    let raw_grid: Grid = serde_json::from_str(&serde_json::to_string(&grid).unwrap()).unwrap();
    assert!(ArithmeticState::capture(&raw_grid).is_err());
    let custom = Grid::new_with_constraints(vec![]);
    assert!(ArithmeticState::capture(&custom).is_err());
    let mut malformed = Grid::new_classic();
    malformed.set_cell_unchecked(position(0), Some(0));
    malformed
        .cell_mut(position(0))
        .set_candidates(BitSet::all_9());
    assert!(ArithmeticState::capture(&malformed).is_err());
    let mut bad_replay = replay;
    bad_replay.version = 2;
    assert!(!bad_replay.verify());
}

#[test]
fn modular_elimination_round_trips_and_applies_to_exact_state() {
    let mut grid = Grid::new_classic();
    let sectors = [12, 3, 10, 18, 1];
    let retained = [3, 30, 28, 1, 9, 12];
    for sector in sectors {
        for cell in super::super::fabric::sector_cells(sector) {
            if !retained.contains(&cell) {
                grid.cell_mut(position(cell)).remove_candidate(1);
            }
        }
    }
    grid.cell_mut(position(3))
        .set_candidates(BitSet::from_slice(&[1, 2, 3]));
    let mut terms: Vec<_> = sectors
        .into_iter()
        .map(|sector| ArithmeticTerm {
            requirement: ArithmeticRequirement::SectorDigit { sector, digit: 1 },
            weight: 2,
        })
        .collect();
    terms.push(ArithmeticTerm {
        requirement: ArithmeticRequirement::Cell { cell: 3 },
        weight: 1,
    });
    let proof = ArithmeticProof::from_residue_terms(&grid, terms, 3, 2, false, 4).unwrap();
    let replay = ArithmeticReplay::capture(&grid, &proof).unwrap();
    let restored: ArithmeticReplay =
        serde_json::from_str(&serde_json::to_string(&replay).unwrap()).unwrap();
    assert!(matches!(
        restored.check().unwrap(),
        ArithmeticCheck::Residue { modulus: 4, .. }
    ));
    restored.apply(&mut grid).unwrap();
    assert_eq!(grid.get_candidates(position(3)).to_vec(), vec![1, 3]);
    assert!(!grid.get_candidates(position(0)).contains(1));
}

#[test]
fn locally_contradictory_application_rolls_back_atomically() {
    let (mut grid, _) = naked_single();
    grid.cell_mut(position(1)).set_candidates(BitSet::single(1));
    // The two singleton premises are globally inconsistent. A conditional proof
    // can still verify, but apply must not install a state with an empty domain.
    let proof = ArithmeticProof::from_terms(
        &grid,
        vec![ArithmeticTerm {
            requirement: ArithmeticRequirement::Cell { cell: 0 },
            weight: 1,
        }],
        0,
        1,
        true,
    )
    .unwrap();
    let replay = ArithmeticReplay::capture(&grid, &proof).unwrap();
    assert!(replay.verify());
    let before = ArithmeticState::capture(&grid).unwrap();
    assert!(replay.apply(&mut grid).is_err());
    assert_eq!(ArithmeticState::capture(&grid).unwrap(), before);
}

#[test]
fn filled_cells_round_trip_with_empty_and_singleton_distinguished() {
    let mut grid = Grid::new_classic();
    grid.set_given(position(40), 5);
    grid.cell_mut(position(0)).set_candidates(BitSet::single(1));
    let state = ArithmeticState::capture(&grid).unwrap();
    assert_eq!((state.values[40], state.domains[40]), (5, 16));
    assert_eq!((state.values[0], state.domains[0]), (0, 1));
    assert_eq!(
        ArithmeticState::capture(&state.to_grid().unwrap()).unwrap(),
        state
    );
}
