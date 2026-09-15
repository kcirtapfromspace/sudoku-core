use super::*;
use crate::{Grid, Position};

fn fabric(pattern: &[(usize, &[u8])]) -> CandidateFabric {
    fabric_with_values(pattern, &[])
}

fn fabric_with_values(pattern: &[(usize, &[u8])], placed: &[(usize, u8)]) -> CandidateFabric {
    let mut grid = Grid::new_classic();
    for &(cell, value) in placed {
        grid.set_cell_unchecked(Position::new(cell / 9, cell % 9), Some(value));
    }
    for &(cell, candidates) in pattern {
        grid.cell_mut(Position::new(cell / 9, cell % 9))
            .set_candidates(BitSet::from_slice(candidates));
    }
    CandidateFabric::from_grid(&grid)
}

fn share_unit(a: usize, b: usize) -> bool {
    a / 9 == b / 9 || a % 9 == b % 9 || (a / 27, a % 9 / 3) == (b / 27, b % 9 / 3)
}

// Independent finite assignment oracle for the candidate pattern. Unrestricted
// cells add no assumptions; every constrained cell and the elimination target
// must admit an assignment respecting standard Sudoku peers.
fn assignment_exists(fab: &CandidateFabric, fixed: Option<(usize, u8)>) -> bool {
    let mut cells: Vec<usize> = (0..81)
        .filter(|&c| fab.values[c].is_none() && fab.cell_cands[c].count() < 9)
        .collect();
    if let Some((target, _)) = fixed {
        if !cells.contains(&target) {
            cells.push(target);
        }
    }
    cells.sort_by_key(|&c| {
        if fixed.is_some_and(|(target, _)| target == c) {
            0
        } else {
            fab.cell_cands[c].count()
        }
    });
    fn search(
        fab: &CandidateFabric,
        cells: &[usize],
        assigned: &mut Vec<(usize, u8)>,
        fixed: Option<(usize, u8)>,
    ) -> bool {
        let Some((&cell, remaining)) = cells.split_first() else {
            return true;
        };
        let candidates: Vec<u8> = match fixed {
            Some((target, value)) if target == cell => vec![value],
            _ => fab.cell_cands[cell].iter().collect(),
        };
        for value in candidates {
            if assigned
                .iter()
                .any(|&(other, existing)| existing == value && share_unit(cell, other))
            {
                continue;
            }
            assigned.push((cell, value));
            if search(fab, remaining, assigned, fixed) {
                return true;
            }
            assigned.pop();
        }
        false
    }
    search(fab, &cells, &mut Vec::new(), fixed)
}

fn assert_sound(fab: &CandidateFabric, finding: &Finding) {
    assert!(
        assignment_exists(fab, None),
        "fixture has no valid local assignment"
    );
    let InferenceResult::Elimination { cell, values } = &finding.inference else {
        panic!("ALS should eliminate candidates");
    };
    assert!(!values.is_empty());
    for &value in values {
        assert!(fab.cell_cands[*cell].contains(value));
        assert!(
            !assignment_exists(fab, Some((*cell, value))),
            "{:?} eliminated feasible {value} from cell {cell}; proof {:?}",
            finding.technique,
            finding.proof
        );
    }
    assert!(matches!(finding.proof, Some(ProofCertificate::Als { .. })));
}

fn assert_descriptor(fab: &CandidateFabric, descriptor: &AlsProofDescriptor) {
    assert!(!descriptor.cells.is_empty());
    let mut union = BitSet::empty();
    for &cell in &descriptor.cells {
        assert!(
            sector_cells(descriptor.sector).contains(&cell),
            "wrong proof sector for {descriptor:?}"
        );
        union = union.union(&fab.cell_cands[cell]);
    }
    assert_eq!(union.to_vec(), descriptor.candidates);
}

#[test]
fn catalog_enumerates_each_almost_locked_subset_once() {
    let fab = fabric(&[(0, &[1, 2]), (1, &[2, 3]), (2, &[1, 3]), (9, &[4, 5])]);
    let catalog = enumerate_als(&fab);
    for als in &catalog {
        assert_eq!(als.candidates.count() as usize, als.cells.len() + 1);
        assert_descriptor(&fab, &als_to_descriptor(als));
    }
    let cell_sets: std::collections::HashSet<_> = catalog.iter().map(|a| a.cells.clone()).collect();
    assert_eq!(cell_sets.len(), catalog.len());
    assert!(cell_sets.contains(&vec![0]));
    assert!(cell_sets.contains(&vec![0, 1]));
    assert!(cell_sets.contains(&vec![0, 2]));
    assert!(cell_sets.contains(&vec![1, 2]));
    assert!(!cell_sets.contains(&vec![0, 1, 2])); // three cells, only three digits: locked, not ALS
    assert!(enumerate_als(&fabric(&[])).is_empty());
}

#[test]
fn restricted_candidates_require_disjoint_sets_and_every_occurrence_to_see_the_other_set() {
    let fab = fabric(&[(0, &[1, 2]), (3, &[1, 3]), (30, &[1, 4])]);
    let a = Als {
        cells: vec![0],
        candidates: BitSet::from_slice(&[1, 2]),
        sector: 0,
    };
    let b = Als {
        cells: vec![3],
        candidates: BitSet::from_slice(&[1, 3]),
        sector: 0,
    };
    let distant = Als {
        cells: vec![30],
        candidates: BitSet::from_slice(&[1, 4]),
        sector: 3,
    };
    assert_eq!(find_rccs(&fab, &a, &b), vec![1]);
    assert!(find_rccs(&fab, &a, &a).is_empty());
    assert!(!is_rcc(&fab, &a, &distant, 1));
    assert!(!is_rcc(&fab, &a, &b, 9));
    let multiple = Als {
        cells: vec![3, 30],
        candidates: BitSet::from_slice(&[1, 3, 4]),
        sector: 12,
    };
    assert!(!is_rcc(&fab, &a, &multiple, 1));
}

#[test]
fn combination_enumeration_is_complete_without_duplicates_or_invalid_subsets() {
    assert!(combinations_usize(&[2, 4, 6], 0).is_empty());
    assert!(combinations_usize(&[2, 4, 6], 4).is_empty());
    assert_eq!(
        combinations_usize(&[2, 4, 6], 2),
        vec![vec![2, 4], vec![2, 6], vec![4, 6]]
    );
    assert_eq!(combinations_usize(&[2, 4, 6], 3), vec![vec![2, 4, 6]]);
}

#[test]
#[allow(deprecated)]
fn pair_and_triplet_exclusion_follow_all_distinct_assignment_logic() {
    let pair = fabric(&[(0, &[1, 2]), (1, &[1, 2]), (2, &[1, 2, 3])]);
    let finding = find_aligned_pair_exclusion(&pair).expect("two cells lock digits 1 and 2");
    assert_eq!(finding.technique, Technique::AlignedPairExclusion);
    assert_sound(&pair, &finding);
    let triple = fabric(&[(0, &[1, 2]), (1, &[2, 3]), (2, &[1, 3]), (3, &[1, 2, 3, 4])]);
    let finding = find_aligned_triplet_exclusion(&triple).expect("three cells lock digits 1,2,3");
    assert_eq!(finding.technique, Technique::AlignedTripletExclusion);
    assert_sound(&triple, &finding);
    let unlocked = fabric(&[(0, &[1, 2]), (1, &[3, 4]), (30, &[1, 2, 3, 4])]);
    assert!(find_aligned_pair_exclusion(&unlocked).is_none());
    assert!(find_aligned_triplet_exclusion(&unlocked).is_none());
    let impossible_triple = fabric(&[(0, &[1, 2]), (1, &[1, 2]), (2, &[1, 2])]);
    assert!(find_aligned_triplet_exclusion(&impossible_triple).is_none());
}

#[test]
fn wing_and_general_pair_patterns_emit_sound_proofs() {
    for (pattern, expected) in [
        (vec![(0, vec![1, 2]), (1, vec![1, 2])], Technique::XYWing),
        (
            vec![(0, vec![1, 2]), (1, vec![1, 3]), (2, vec![2, 3])],
            Technique::XYZWing,
        ),
        (
            vec![
                (0, vec![1, 2]),
                (1, vec![1, 3]),
                (2, vec![3, 4]),
                (3, vec![2, 4]),
            ],
            Technique::WXYZWing,
        ),
        (
            vec![
                (0, vec![1, 2]),
                (1, vec![1, 3]),
                (2, vec![3, 4]),
                (3, vec![4, 5]),
                (4, vec![2, 5]),
            ],
            Technique::AlsXz,
        ),
    ] {
        let input: Vec<_> = pattern.iter().map(|(c, v)| (*c, v.as_slice())).collect();
        let fab = fabric(&input);
        let finding = match expected {
            Technique::XYWing => find_xy_wing(&fab),
            Technique::XYZWing => find_xyz_wing(&fab),
            Technique::WXYZWing => find_wxyz_wing(&fab),
            _ => find_als_xz(&fab),
        }
        .expect("the ring forms a locked subset via two disjoint ALS");
        assert_eq!(finding.technique, expected);
        assert_sound(&fab, &finding);
        if let Some(ProofCertificate::Als {
            als_chain,
            rcc_values,
            z_value,
        }) = &finding.proof
        {
            assert_eq!(als_chain.len(), 2);
            assert_eq!(rcc_values.len(), 1);
            assert_ne!(*z_value, Some(rcc_values[0]));
            for descriptor in als_chain {
                assert_descriptor(&fab, descriptor);
            }
        }
    }
    assert!(find_als_xz_filtered(&fabric(&[]), None).is_none());
    assert_eq!(technique_name(Technique::HiddenSingle), "ALS");
    assert_eq!(technique_name(Technique::AlsXyWing), "ALS-XY-Wing");
    assert_eq!(technique_name(Technique::AlsChain), "ALS Chain");
    assert_eq!(technique_name(Technique::SueDeCoq), "Sue de Coq");
    assert_eq!(technique_name(Technique::DeathBlossom), "Death Blossom");
}

#[test]
fn three_and_four_als_chains_produce_sound_eliminations() {
    let patterns: Vec<Vec<(usize, Vec<u8>)>> = vec![
        vec![
            (0, vec![1, 3]),
            (3, vec![1, 2]),
            (30, vec![2, 3]),
            (27, vec![3, 4]),
        ],
        vec![
            (0, vec![1, 4]),
            (3, vec![1, 2]),
            (30, vec![2, 3]),
            (27, vec![3, 4]),
            (9, vec![4, 5]),
        ],
    ];
    for (index, pattern) in patterns.iter().enumerate() {
        let input: Vec<_> = pattern.iter().map(|(c, v)| (*c, v.as_slice())).collect();
        let fab = fabric(&input);
        let finding = if index == 0 {
            find_als_xy_wing(&fab)
        } else {
            find_als_chain(&fab)
        }
        .expect("closed endpoint chain");
        assert_sound(&fab, &finding);
        if let Some(ProofCertificate::Als {
            als_chain,
            rcc_values,
            ..
        }) = &finding.proof
        {
            assert_eq!(als_chain.len(), index + 3);
            assert_eq!(rcc_values.len(), index + 2);
            for descriptor in als_chain {
                assert_descriptor(&fab, descriptor);
            }
        }
    }
    let no_links = fabric(&[(0, &[1, 2]), (40, &[3, 4]), (80, &[5, 6])]);
    assert!(find_als_xy_wing(&no_links).is_none());
    assert!(find_als_chain(&no_links).is_none());
}

#[test]
fn sue_de_coq_locks_disjoint_digit_groups_in_box_and_line() {
    // Two intersection cells choose one digit each from {1,2} and {3,4};
    // one outside ALS in each unit consumes the remaining digit of its group.
    for (box_target, line_target) in [(true, false), (false, true)] {
        let mut pattern: Vec<(usize, &[u8])> = vec![
            (0, &[1, 2, 3, 4]),
            (1, &[1, 2, 3, 4]),
            (9, &[1, 2]),
            (3, &[3, 4]),
        ];
        // The third overlap cell is filled in this local pattern, so the DDS
        // intersection contains exactly the two cells that supply the freedom.
        if box_target {
            pattern.push((10, &[1, 2, 5]));
        }
        if line_target {
            pattern.push((4, &[3, 4, 5]));
        }
        let fab = fabric_with_values(&pattern, &[(2, 9)]);
        let finding = find_sue_de_coq(&fab).expect("distributed locked subsets");
        assert_eq!(finding.technique, Technique::SueDeCoq);
        assert_sound(&fab, &finding);
        if let Some(ProofCertificate::Als { als_chain, .. }) = &finding.proof {
            for descriptor in als_chain {
                assert_descriptor(&fab, descriptor);
            }
        }
    }
    assert!(find_sue_de_coq(&fabric(&[])).is_none());
}

#[test]
fn death_blossom_petals_force_a_shared_digit_away_from_common_peers() {
    let fab = fabric(&[(0, &[1, 2]), (3, &[1, 3]), (27, &[2, 3]), (30, &[3, 4])]);
    let finding = find_death_blossom(&fab).expect("either stem value forces 3 into a petal");
    assert_eq!(finding.technique, Technique::DeathBlossom);
    assert_sound(&fab, &finding);
    if let Some(ProofCertificate::Als {
        als_chain,
        rcc_values,
        z_value,
    }) = &finding.proof
    {
        assert_eq!(als_chain.len(), 3);
        assert_eq!(rcc_values.len(), 2);
        assert!(z_value.is_some());
    }
    let no_petals = fabric(&[(0, &[1, 2]), (40, &[3, 4])]);
    assert!(find_death_blossom(&no_petals).is_none());
    assert!(find_death_blossom(&fabric(&[])).is_none());
}

#[test]
#[allow(deprecated)]
fn als_engines_preserve_both_completions_of_an_ambiguous_rectangle() {
    let solution =
        "534678912672195348198342567859761423426853791713924856961537284287419635345286179";
    let mut grid = Grid::from_string(solution).unwrap();
    for cell in [3, 4, 30, 31] {
        let pos = Position::new(cell / 9, cell % 9);
        grid.set_cell_unchecked(pos, None);
        grid.cell_mut(pos).set_given(false);
    }
    grid.recalculate_candidates();
    let fab = CandidateFabric::from_grid(&grid);
    for cell in [3, 4, 30, 31] {
        assert_eq!(fab.cell_cands[cell].to_vec(), vec![6, 7]);
    }
    for (name, find) in [
        (
            "xy_wing",
            find_xy_wing as fn(&CandidateFabric) -> Option<Finding>,
        ),
        ("xyz_wing", find_xyz_wing),
        ("wxyz_wing", find_wxyz_wing),
        ("als_xz", find_als_xz),
        ("als_xy_wing", find_als_xy_wing),
        ("als_chain", find_als_chain),
        ("sue_de_coq", find_sue_de_coq),
        ("death_blossom", find_death_blossom),
        ("aligned_pair", find_aligned_pair_exclusion),
        ("aligned_triplet", find_aligned_triplet_exclusion),
    ] {
        assert!(
            find(&fab).is_none(),
            "{name} must preserve both valid 6/7 completions"
        );
    }
}

#[test]
fn death_blossom_requires_the_stem_to_see_every_linking_digit_in_a_petal() {
    // The 3-cell column ALS has candidates {1,3,4,5}, but the stem sees
    // only its first occurrence of 1. Stem=1 therefore does not lock this
    // petal: (stem,A1,A2,A3,B,target)=(1,4,1,5,2,3) is a valid assignment.
    let fab = fabric(&[
        (0, &[1, 2]),
        (3, &[1, 4, 5]),
        (30, &[1, 3, 5]),
        (57, &[3, 4, 5]),
        (9, &[2, 3]),
        (12, &[3, 6]),
    ]);
    assert!(assignment_exists(&fab, Some((12, 3))));
    if let Some(finding) = find_death_blossom(&fab) {
        assert_sound(&fab, &finding);
    }
}
