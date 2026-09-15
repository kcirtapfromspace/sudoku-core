use super::{
    explain::{Finding, InferenceResult, ProofCertificate},
    fabric::CandidateFabric,
    fish_engine, Technique,
};
use crate::{Grid, Position};

const SOLUTION: &str =
    "534678912672195348198342567859761423426853791713924856961537284287419635345286179";

fn digit_fish(size: usize, transpose: bool, fin: bool) -> Grid {
    let mut grid = Grid::new_classic();
    let covers = [0, 5, 6, 1]; // The known solution's 5 positions in rows 1–4.
    for row in 0..size {
        for col in 0..9 {
            if !covers[..size].contains(&col) && !(fin && row == 0 && col == 3) {
                let pos = if transpose {
                    Position::new(col, row)
                } else {
                    Position::new(row, col)
                };
                grid.cell_mut(pos).remove_candidate(5);
            }
        }
    }
    grid
}

fn assert_sound_elimination(finding: &Finding, grid: &Grid, transpose: bool) {
    let InferenceResult::Elimination { cell, values } = &finding.inference else {
        panic!("Expected an elimination, got {finding:?}");
    };
    let pos = Position::new(cell / 9, cell % 9);
    let solution_idx = if transpose {
        pos.col * 9 + pos.row
    } else {
        *cell
    };
    assert!(!values.is_empty());
    for &digit in values {
        assert!(grid.cell(pos).has_candidate(digit));
        assert_ne!(
            SOLUTION.as_bytes()[solution_idx] - b'0',
            digit,
            "{finding:?}"
        );
    }
}

#[test]
fn basic_fish_proofs_eliminate_only_outside_the_bases_in_both_orientations() {
    for (size, technique) in [
        (2, Technique::XWing),
        (3, Technique::Swordfish),
        (4, Technique::Jellyfish),
    ] {
        for transpose in [false, true] {
            let grid = digit_fish(size, transpose, false);
            let finding =
                fish_engine::find_basic_fish(&CandidateFabric::from_grid(&grid), size).unwrap();
            assert_eq!(finding.technique, technique);
            assert_sound_elimination(&finding, &grid, transpose);
            let Some(ProofCertificate::Fish {
                digit,
                base_sectors,
                cover_sectors,
                fins,
            }) = finding.proof
            else {
                panic!("Fish proof missing");
            };
            assert_eq!(digit, 5);
            assert_eq!(base_sectors.len(), size);
            assert_eq!(cover_sectors.len(), size);
            assert!(fins.is_empty());
            assert!(base_sectors
                .iter()
                .all(|base| !cover_sectors.contains(base)));
        }
    }
}

#[test]
fn finned_fish_proof_keeps_eliminations_in_the_fin_box() {
    for transpose in [false, true] {
        let grid = digit_fish(2, transpose, true);
        let finding = fish_engine::find_finned_fish(&CandidateFabric::from_grid(&grid), 2).unwrap();
        assert_eq!(finding.technique, Technique::FinnedXWing);
        assert_sound_elimination(&finding, &grid, transpose);
        let InferenceResult::Elimination { cell, .. } = finding.inference else {
            unreachable!()
        };
        let Some(ProofCertificate::Fish { digit, fins, .. }) = finding.proof else {
            panic!("Fish proof missing")
        };
        assert_eq!(digit, 5);
        assert!(!fins.is_empty());
        let target_box = Position::new(cell / 9, cell % 9).box_index();
        assert!(fins
            .iter()
            .all(|fin| Position::new(fin / 9, fin % 9).box_index() == target_box));
    }
}

#[test]
fn generalized_fish_cannot_remove_either_completion_of_an_ambiguous_rectangle() {
    let mut grid = Grid::new_classic();
    for (idx, byte) in SOLUTION.bytes().enumerate() {
        if ![3, 4, 30, 31].contains(&idx) {
            grid.set_given(Position::new(idx / 9, idx % 9), byte - b'0');
        }
    }
    grid.recalculate_candidates();
    assert_eq!(super::Solver::new().count_solutions(&grid, 3), 2);
    let fab = CandidateFabric::from_grid(&grid);
    assert!(fish_engine::find_franken_fish(&fab).is_none());
    assert!(fish_engine::find_mutant_fish(&fab).is_none());
    assert!(fish_engine::find_siamese_fish(&fab).is_none());
}

#[test]
fn generalized_fish_with_independent_groups_still_produce_sound_proofs() {
    let mut grid = Grid::new_classic();
    for row in 0..2 {
        for col in 0..9 {
            if ![0, 3, 4, 5].contains(&col) {
                grid.cell_mut(Position::new(row, col)).remove_candidate(5);
            }
        }
    }
    let fabric = CandidateFabric::from_grid(&grid);
    for (finding, expected) in [
        (
            fish_engine::find_franken_fish(&fabric),
            Technique::FrankenFish,
        ),
        (
            fish_engine::find_mutant_fish(&fabric),
            Technique::MutantFish,
        ),
    ] {
        let finding = finding.expect("Independent line/box groups form a valid fish");
        assert_eq!(finding.technique, expected);
        assert_sound_elimination(&finding, &grid, false);
        let Some(ProofCertificate::Fish {
            base_sectors,
            cover_sectors,
            ..
        }) = finding.proof
        else {
            panic!("Missing proof")
        };
        assert!(base_sectors
            .iter()
            .chain(&cover_sectors)
            .any(|&sector| sector >= 18));
        for sectors in [&base_sectors, &cover_sectors] {
            let mut occurrences = [0u8; 81];
            for &sector in sectors {
                for cell in fabric.sector_cells_with_candidate(sector, 5) {
                    occurrences[cell] += 1;
                    assert_eq!(
                        occurrences[cell], 1,
                        "A placement cannot satisfy two groups"
                    );
                }
            }
        }
    }
}

#[test]
fn larger_finned_fish_keep_the_solution_in_each_orientation() {
    for (size, expected) in [
        (3, Technique::FinnedSwordfish),
        (4, Technique::FinnedJellyfish),
    ] {
        for transpose in [false, true] {
            let mut grid = Grid::new_classic();
            let bases = [0, 3, 6, 7];
            let covers = [0, 1, 3, 8];
            for &row in &bases[..size] {
                for col in 0..9 {
                    if !covers[..size].contains(&col) && !(row == 0 && col == 2) {
                        let pos = if transpose {
                            Position::new(col, row)
                        } else {
                            Position::new(row, col)
                        };
                        grid.cell_mut(pos).remove_candidate(5);
                    }
                }
            }
            let finding =
                fish_engine::find_finned_fish(&CandidateFabric::from_grid(&grid), size).unwrap();
            assert_eq!(finding.technique, expected);
            assert_sound_elimination(&finding, &grid, transpose);
        }
    }
}

fn pencil_grid_with_five(transpose: bool) -> Grid {
    let mut grid = Grid::new_classic();
    for (idx, byte) in SOLUTION.bytes().enumerate() {
        let pos = if transpose {
            Position::new(idx % 9, idx / 9)
        } else {
            Position::new(idx / 9, idx % 9)
        };
        grid.cell_mut(pos)
            .set_candidates(crate::BitSet::from_slice(&[5, byte - b'0']));
    }
    grid
}

#[test]
fn overlapping_finned_fish_preserve_the_solution_in_both_orientations() {
    for transpose in [false, true] {
        let mut grid = pencil_grid_with_five(transpose);
        for row in [0, 3, 8] {
            for col in 0..9 {
                if ![0, 1].contains(&col) && !(row == 8 && col == 2) {
                    let pos = if transpose {
                        Position::new(col, row)
                    } else {
                        Position::new(row, col)
                    };
                    grid.cell_mut(pos).remove_candidate(5);
                }
            }
        }
        let finding = fish_engine::find_siamese_fish(&CandidateFabric::from_grid(&grid)).unwrap();
        assert_eq!(finding.technique, Technique::SiameseFish);
        assert_sound_elimination(&finding, &grid, transpose);
        let Some(ProofCertificate::Fish { digit, fins, .. }) = finding.proof else {
            panic!("Siamese proof missing")
        };
        assert_eq!(digit, 5);
        assert!(!fins.is_empty());
    }
}

#[test]
fn box_line_intersections_eliminate_only_outside_the_intersection() {
    for transpose in [false, true] {
        let mut pointing = pencil_grid_with_five(transpose);
        for row in 1..3 {
            for col in 0..3 {
                let pos = if transpose {
                    Position::new(col, row)
                } else {
                    Position::new(row, col)
                };
                pointing.cell_mut(pos).remove_candidate(5);
            }
        }
        let finding =
            fish_engine::find_pointing_pair(&CandidateFabric::from_grid(&pointing)).unwrap();
        assert_eq!(finding.technique, Technique::PointingPair);
        assert_sound_elimination(&finding, &pointing, transpose);
        let mut claiming = pencil_grid_with_five(transpose);
        for col in 3..9 {
            let pos = if transpose {
                Position::new(col, 0)
            } else {
                Position::new(0, col)
            };
            claiming.cell_mut(pos).remove_candidate(5);
        }
        let finding =
            fish_engine::find_box_line_reduction(&CandidateFabric::from_grid(&claiming)).unwrap();
        assert_eq!(finding.technique, Technique::BoxLineReduction);
        assert_sound_elimination(&finding, &claiming, transpose);
    }
}
