use super::super::tests::{assert_forced, pos, restrict_one, sector_terms, xwing};
use super::*;
use crate::{ArithmeticReplay, ArithmeticState, BitSet, Solver};

const PARTITIONS: [ArithmeticPartition; 4] = [
    ArithmeticPartition::Cells,
    ArithmeticPartition::RowDigits,
    ArithmeticPartition::ColumnDigits,
    ArithmeticPartition::BoxDigits,
];

fn cell_pair() -> (Grid, Vec<ArithmeticTerm>) {
    let mut grid = Grid::new_classic();
    for digit in [1, 2] {
        for cell in 0..9 {
            if cell != 0 && cell != 3 {
                grid.cell_mut(pos(cell)).remove_candidate(digit);
            }
        }
    }
    grid.cell_mut(pos(0))
        .set_candidates(BitSet::from_slice(&[1, 2, 3]));
    grid.cell_mut(pos(3))
        .set_candidates(BitSet::from_slice(&[1, 2]));
    let terms = [1, 2]
        .into_iter()
        .map(|digit| ArithmeticTerm {
            requirement: ArithmeticRequirement::SectorDigit { sector: 0, digit },
            weight: 1,
        })
        .collect();
    (grid, terms)
}

fn proof(
    grid: &Grid,
    terms: Vec<ArithmeticTerm>,
    cell: usize,
    digit: u8,
    value: bool,
    partition: ArithmeticPartition,
    modulus: u8,
) -> ArithmeticProof {
    ArithmeticProof::from_terms_with_terminal(
        grid,
        terms,
        cell,
        digit,
        value,
        ArithmeticTerminal::PartitionResidue { modulus, partition },
    )
    .unwrap()
}

#[test]
fn partitions_are_disjoint_complete_native_requirements() {
    let snapshot = Snapshot::from_grid(&Grid::new_classic()).unwrap();
    for partition in PARTITIONS {
        let groups = groups(&snapshot, partition);
        let mut visits = [0; VARIABLE_COUNT];
        for (index, group) in groups.iter().enumerate() {
            assert_eq!(group.len(), 9);
            let requirement = match partition {
                ArithmeticPartition::Cells => ArithmeticRequirement::Cell { cell: index },
                _ => ArithmeticRequirement::SectorDigit {
                    sector: index / 9
                        + match partition {
                            ArithmeticPartition::RowDigits => 0,
                            ArithmeticPartition::ColumnDigits => 9,
                            ArithmeticPartition::BoxDigits => 18,
                            _ => unreachable!(),
                        },
                    digit: (index % 9 + 1) as u8,
                },
            };
            assert_eq!(
                *group,
                snapshot.requirement_variables(requirement.id().unwrap())
            );
            for &variable in group {
                visits[variable] += 1;
            }
        }
        assert_eq!(visits, [1; VARIABLE_COUNT]);
    }
}

#[test]
fn grouped_remainders_match_exhaustive_integer_choices() {
    let choices: &[&[i16]] = &[
        &[],
        &[-2],
        &[-1],
        &[0],
        &[1],
        &[2],
        &[-1, 1],
        &[0, 0],
        &[1, 1],
        &[-2, 0, 2],
    ];
    for length in 0..=4u32 {
        for mut encoding in 0..choices.len().pow(length) {
            let groups: Vec<_> = (0..length)
                .map(|_| {
                    let group = choices[encoding % choices.len()];
                    encoding /= choices.len();
                    group
                })
                .collect();
            // Enumerate integer sums before taking any remainders. No bitsets
            // or production recurrence are used by this reference calculation.
            let mut sums = vec![0i32];
            for group in &groups {
                sums = sums
                    .iter()
                    .flat_map(|sum| group.iter().map(move |&v| sum + i32::from(v)))
                    .collect();
            }
            for modulus in 2..=16 {
                let expected = sums.iter().fold(0u16, |set, &s| {
                    set | (1 << s.rem_euclid(i32::from(modulus)))
                });
                let actual = groups.iter().fold(1, |set, group| {
                    let choices = group
                        .iter()
                        .fold(0, |bits, &c| bits | (1 << residue(c, modulus)));
                    sum_residues(set, choices, modulus)
                });
                assert_eq!(actual, expected, "groups={groups:?}, modulus={modulus}");
            }
        }
    }
    // Exactly one of these two equal coefficients contributes: never zero or two.
    assert_eq!(sum_residues(1, 1 << 1, 4), 1 << 1);
    assert_eq!(sum_residues(0b1111, 0, 4), 0);
}

#[test]
fn zero_coefficient_target_is_strictly_beyond_both_old_terminals() {
    let (grid, terms) = cell_pair();
    assert!(ArithmeticProof::from_terms(&grid, terms.clone(), 0, 3, false).is_none());
    for modulus in 2..=16 {
        assert!(
            ArithmeticProof::from_residue_terms(&grid, terms.clone(), 0, 3, false, modulus)
                .is_none()
        );
        let proof = proof(
            &grid,
            terms.clone(),
            0,
            3,
            false,
            ArithmeticPartition::Cells,
            modulus,
        );
        assert_eq!(proof.version, 3);
    }
    let proof = proof(
        &grid,
        terms.clone(),
        0,
        3,
        false,
        ArithmeticPartition::Cells,
        2,
    );
    assert_eq!(
        proof.check(&grid),
        Some(ArithmeticCheck::PartitionResidue {
            rhs: 2,
            modulus: 2,
            partition: ArithmeticPartition::Cells,
            required_residue: 0,
            reachable_residues: vec![1],
        })
    );
    // The separate exact-cover oracle verifies satisfiability and rules out the
    // opposite step using all original domains and all 324 native requirements.
    assert_forced(&grid, 0, 3, false);
    assert!(ArithmeticProof::from_terms_with_terminal(
        &grid,
        terms,
        0,
        3,
        false,
        ArithmeticTerminal::PartitionResidue {
            modulus: 2,
            partition: ArithmeticPartition::RowDigits
        }
    )
    .is_none());
}

#[test]
fn all_partition_deductions_have_independent_full_sudoku_witnesses() {
    let (grid, terms) = cell_pair();
    let mut fixtures = vec![(grid, terms, 0, 3, ArithmeticPartition::Cells)];
    fixtures.push((
        xwing(),
        sector_terms(&[0, 1], &[1, 1]),
        27,
        1,
        ArithmeticPartition::ColumnDigits,
    ));
    let mut grid = Grid::new_classic();
    restrict_one(&mut grid, 9, &[0, 54]);
    restrict_one(&mut grid, 10, &[1, 55]);
    fixtures.push((
        grid,
        sector_terms(&[9, 10], &[1, 1]),
        3,
        1,
        ArithmeticPartition::RowDigits,
    ));
    let mut grid = Grid::new_classic();
    restrict_one(&mut grid, 0, &[0, 1]);
    fixtures.push((
        grid,
        sector_terms(&[0], &[1]),
        9,
        1,
        ArithmeticPartition::BoxDigits,
    ));
    for (grid, terms, cell, digit, partition) in fixtures {
        if matches!(
            partition,
            ArithmeticPartition::RowDigits | ArithmeticPartition::ColumnDigits
        ) {
            // The required sum 2 aliases the attainable sum 0 modulo 2.
            // Modulus 3 separates them; a weaker modulus must not invent a step.
            assert!(ArithmeticProof::from_terms_with_terminal(
                &grid,
                terms.clone(),
                cell,
                digit,
                false,
                ArithmeticTerminal::PartitionResidue {
                    modulus: 2,
                    partition
                },
            )
            .is_none());
        }
        let certificate = proof(&grid, terms, cell, digit, false, partition, 3);
        assert_forced(&grid, cell, digit, false);
        let snapshot = Snapshot::from_grid(&grid).unwrap();
        let hint = make_finding(
            &snapshot,
            certificate.clone(),
            certificate.check(&grid).unwrap(),
        )
        .to_hint();
        assert!(hint.explanation.contains("exactly one candidate per"));
        assert!(hint.involved_cells.contains(&pos(cell)));
        let options = ArithmeticSearchOptions {
            max_sources: 2,
            ..Default::default()
        };
        let result = Solver::new().search_partition_arithmetic(&grid, &options, partition, 3);
        let hint = result
            .hint
            .expect("bounded search discovers a partition deduction");
        let Some(ProofCertificate::Arithmetic(found)) = hint.proof else {
            panic!("missing proof")
        };
        assert!(matches!(
            found.terminal,
            ArithmeticTerminal::PartitionResidue { .. }
        ));
        assert!(found.verify(&grid));
        assert_forced(&grid, found.cell, found.digit, found.value);
    }
}

#[test]
fn partition_proof_replay_round_trips_and_preserves_prior_eliminations() {
    let (mut grid, terms) = cell_pair();
    let certificate = proof(&grid, terms, 0, 3, false, ArithmeticPartition::Cells, 2);
    let replay = ArithmeticReplay::capture(&grid, &certificate).unwrap();
    let json = serde_json::to_string(&replay).unwrap();
    let decoded: ArithmeticReplay = serde_json::from_str(&json).unwrap();
    assert_eq!(decoded, replay);
    assert!(decoded.verify());
    let before = ArithmeticState::capture(&grid).unwrap();
    decoded.apply(&mut grid).unwrap();
    let mut expected = before;
    expected.domains[0] &= !(1 << 2);
    assert_eq!(ArithmeticState::capture(&grid).unwrap(), expected);
    assert!(decoded.apply(&mut grid).is_err());
    assert!(!certificate.verify(&grid));
}

#[test]
fn malformed_and_both_endpoint_partition_proofs_are_rejected() {
    let (grid, terms) = cell_pair();
    let certificate = proof(&grid, terms, 0, 3, false, ArithmeticPartition::Cells, 2);
    for version in [0, 1, 2, 4, u8::MAX] {
        let mut bad = certificate.clone();
        bad.version = version;
        assert!(!bad.verify(&grid));
    }
    for modulus in [0, 1, 17, u8::MAX] {
        let mut bad = certificate.clone();
        bad.terminal = ArithmeticTerminal::PartitionResidue {
            modulus,
            partition: ArithmeticPartition::Cells,
        };
        assert!(!bad.verify(&grid));
    }
    for weight in [0, 9, -9, i8::MIN, i8::MAX] {
        let mut bad = certificate.clone();
        bad.terms[0].weight = weight;
        assert!(!bad.verify(&grid));
    }
    let mut bad = certificate.clone();
    bad.value = true;
    assert!(!bad.verify(&grid));
    bad = certificate.clone();
    bad.terms[1] = bad.terms[0].clone();
    assert!(!bad.verify(&grid));
    let mut stale = grid.clone();
    stale.cell_mut(pos(80)).remove_candidate(9);
    assert!(!certificate.verify(&stale));

    let mut inconsistent = Grid::new_classic();
    for cell in [0, 1] {
        inconsistent
            .cell_mut(pos(cell))
            .set_candidates(BitSet::from_slice(&[1]));
    }
    restrict_one(&mut inconsistent, 0, &[0, 1]);
    let terms = vec![
        ArithmeticTerm {
            requirement: ArithmeticRequirement::Cell { cell: 0 },
            weight: 1,
        },
        ArithmeticTerm {
            requirement: ArithmeticRequirement::Cell { cell: 1 },
            weight: 1,
        },
        ArithmeticTerm {
            requirement: ArithmeticRequirement::SectorDigit {
                sector: 0,
                digit: 1,
            },
            weight: -1,
        },
    ];
    for partition in PARTITIONS {
        assert!(ArithmeticProof::from_terms_with_terminal(
            &inconsistent,
            terms.clone(),
            2,
            2,
            false,
            ArithmeticTerminal::PartitionResidue {
                modulus: 2,
                partition
            }
        )
        .is_none());
    }
}

#[test]
fn negative_and_repeated_coefficients_preserve_exclusive_choices() {
    let (grid, terms) = cell_pair();
    for scale in [-8, -2, -1, 1, 2, 8] {
        let terms = terms
            .iter()
            .map(|t| ArithmeticTerm {
                requirement: t.requirement.clone(),
                weight: scale,
            })
            .collect();
        let certificate = proof(&grid, terms, 0, 3, false, ArithmeticPartition::Cells, 3);
        assert!(certificate.verify(&grid));
    }
    // Saturation in earlier groups cannot hide an empty last group under x_t=0.
    let mut grid = Grid::new_classic();
    grid.cell_mut(pos(80))
        .set_candidates(BitSet::from_slice(&[9]));
    let certificate = proof(
        &grid,
        sector_terms(&[0], &[1]),
        80,
        9,
        true,
        ArithmeticPartition::Cells,
        2,
    );
    assert_forced(&grid, 80, 9, true);
    assert!(
        matches!(certificate.check(&grid), Some(ArithmeticCheck::PartitionResidue { reachable_residues, .. }) if reachable_residues.is_empty())
    );
}

#[test]
fn partition_search_is_deterministic_bounded_and_preserves_masks() {
    let (grid, _) = cell_pair();
    let before = ArithmeticState::capture(&grid).unwrap();
    let options = ArithmeticSearchOptions {
        max_sources: 2,
        ..Default::default()
    };
    let solver = Solver::new();
    let first = solver.search_partition_arithmetic(&grid, &options, ArithmeticPartition::Cells, 2);
    let second = solver.search_partition_arithmetic(&grid, &options, ArithmeticPartition::Cells, 2);
    assert_eq!(first.tested_combinations, second.tested_combinations);
    assert_eq!(
        serde_json::to_value(&first.hint).unwrap(),
        serde_json::to_value(&second.hint).unwrap()
    );
    assert!(first.hint.is_some());
    assert!(first.tested_combinations <= options.max_combinations);
    assert_eq!(ArithmeticState::capture(&grid).unwrap(), before);
    for modulus in [0, 1, 17, u8::MAX] {
        let result = solver.search_partition_arithmetic(
            &grid,
            &options,
            ArithmeticPartition::Cells,
            modulus,
        );
        assert!(result.error.unwrap().contains("modulus"));
        assert_eq!(result.tested_combinations, 0);
    }
    for budget in [0, 1] {
        let options = ArithmeticSearchOptions {
            max_combinations: budget,
            ..options.clone()
        };
        let result =
            solver.search_partition_arithmetic(&grid, &options, ArithmeticPartition::Cells, 2);
        assert!(result.hint.is_none());
        assert!(result.budget_exhausted);
        assert_eq!(result.tested_combinations, budget);
    }
    let result = solver.search_partition_arithmetic(
        &Grid::new_classic(),
        &options,
        ArithmeticPartition::Cells,
        2,
    );
    assert!(result.hint.is_none());
    let options = ArithmeticSearchOptions {
        max_sources: 7,
        ..options
    };
    assert!(solver
        .search_partition_arithmetic(&grid, &options, ArithmeticPartition::Cells, 2)
        .error
        .is_some());
}
