use super::*;
use crate::{BitSet, Position, Solver};

fn fixture() -> Grid {
    let json: serde_json::Value = serde_json::from_str(include_str!(
        "../../research/nonmodular/three-house-847.json"
    ))
    .unwrap();
    let mut grid = Grid::new_classic();
    for (cell, mask) in json["domains"].as_array().unwrap().iter().enumerate() {
        grid.cell_mut(idx_to_pos(cell))
            .set_candidates(BitSet::from_raw((mask.as_u64().unwrap() as u16) << 1));
    }
    grid
}
fn certificate() -> (Grid, ThreeHouseProof) {
    let grid = fixture();
    let proof = ThreeHouseProof::for_elimination(&grid, 0, 3, 1).unwrap();
    (grid, proof)
}
fn hint_proof(hint: &Hint) -> &ThreeHouseProof {
    match hint.proof.as_ref().unwrap() {
        ProofCertificate::ThreeHouse(proof) => proof,
        other => panic!("unexpected proof: {other:?}"),
    }
}

#[test]
fn all_crossing_geometries_and_signature_counts() {
    for center in 0..81 {
        let layout = Layout::new(center).unwrap();
        assert_eq!(
            layout.atoms.each_ref().map(|cells| cells.len()),
            [1, 2, 2, 6, 6, 4]
        );
        assert_eq!(layout.cells.len(), 21);
        let mut union = layout.atoms.concat();
        union.sort();
        assert_eq!(union, layout.cells);
        for sector in layout.houses {
            assert!(sector_cells(sector).iter().all(|cell| union.contains(cell)));
        }
    }
    assert!(Layout::new(81).is_none());
    let mut count = 0;
    assert!(signatures(ALL, |sets| {
        count += 1;
        assert_eq!(sets.map(u16::count_ones), [1, 2, 2, 6, 6, 4]);
        assert_eq!(sets[0] | sets[1] | sets[3], ALL);
        assert_eq!(sets[0] | sets[2] | sets[4], ALL);
        assert_eq!(sets[0] | sets[1] | sets[2] | sets[5], ALL);
        true
    }));
    assert_eq!(count, 3780);
    assert!(!signatures(1, |_| false));
}

// Independent enumerator visits actual injective assignments, without DP.
fn enumerate(domains: &[u16], digits: u16) -> Option<Vec<u16>> {
    fn walk(
        domains: &[u16],
        digits: u16,
        path: &mut Vec<u16>,
        support: &mut [u16],
        count: &mut usize,
    ) {
        if path.len() == domains.len() {
            *count += 1;
            for (i, &bit) in path.iter().enumerate() {
                support[i] |= bit;
            }
            return;
        }
        for digit in 0..9 {
            let bit = 1 << digit;
            if digits & domains[path.len()] & bit != 0 {
                path.push(bit);
                walk(domains, digits ^ bit, path, support, count);
                path.pop();
            }
        }
    }
    let mut count = 0;
    let mut support = vec![0; domains.len()];
    walk(domains, digits, &mut Vec::new(), &mut support, &mut count);
    (count != 0).then_some(support)
}

#[test]
fn subset_matching_agrees_with_exhaustive_assignments_and_hall() {
    let mut seed = 0x19a384u32;
    for n in 1..=6 {
        let cells: Vec<_> = (0..n).collect();
        for _ in 0..150 {
            let domains: Vec<_> = (0..n)
                .map(|_| {
                    seed = seed.wrapping_mul(1664525).wrapping_add(1013904223);
                    (seed >> 11) as u16 & ALL
                })
                .collect();
            let digits = (1 << n) - 1;
            let brute = enumerate(&domains, digits);
            let actual = matching(&domains, digits);
            assert_eq!(actual.as_ref().map(|m| &m.support), brute.as_ref());
            if let Some(found) = actual {
                let mut used = 0;
                for (&digit, domain) in found.witness.iter().zip(&domains) {
                    let bit = 1 << (digit - 1);
                    assert_ne!(domain & bit, 0);
                    assert_eq!(used & bit, 0);
                    used |= bit;
                }
                assert_eq!(used, digits);
                assert!(hall_subset(&cells, &domains, digits).is_none());
            } else {
                assert!(hall_subset(&cells, &domains, digits).is_some());
            }
        }
    }
}

#[test]
fn synthetic_three_house_proof_and_search_preserve_candidates() {
    let (grid, proof) = certificate();
    let state = CandidateState::capture(&grid).unwrap();
    let checked = proof.check(&grid).unwrap();
    assert_eq!(checked.signatures, 2100);
    assert!(checked.hall_obstructions > 0);
    assert!(proof.verify(&grid));
    let hint = Solver::new().get_three_house_hint(&grid).unwrap();
    assert_eq!(hint.technique, Technique::ThreeHouseMatching);
    assert!(hint.explanation.contains("Hall-deficient"));
    assert!(hint_proof(&hint).verify(&grid));
    assert_eq!(CandidateState::capture(&grid).unwrap(), state);
    assert!(ThreeHouseProof::for_elimination(&grid, 0, 3, 4).is_none()); // a known completion uses 4
    for (center, cell, digit) in [
        (81, 3, 1),
        (0, 81, 1),
        (0, 80, 1),
        (0, 3, 0),
        (0, 3, 10),
        (0, 3, 2),
    ] {
        assert!(ThreeHouseProof::for_elimination(&grid, center, cell, digit).is_none());
    }
}

#[test]
fn checker_rejects_tampered_evidence_and_wrong_state() {
    let (grid, proof) = certificate();
    let mut mutations = Vec::new();
    macro_rules! bad {
        ($field:ident, $value:expr) => {{
            let mut p = proof.clone();
            p.$field = $value;
            mutations.push(p);
        }};
    }
    bad!(version, 2);
    bad!(state_hash, "incorrect".into());
    bad!(center, 81);
    bad!(cell, 81);
    bad!(cell, 80);
    bad!(digit, 0);
    bad!(digit, 10);
    bad!(digit, 2);
    bad!(witness, vec![]);
    bad!(witness, vec![1; 21]);
    bad!(blockers, vec![]);
    bad!(blockers, vec![proof.blockers[0].clone(); MAX_BLOCKERS + 1]);
    let mut p = proof.clone();
    p.witness[0] = 0;
    mutations.push(p);
    let mut p = proof.clone();
    p.witness[0] = 10;
    mutations.push(p);
    let mut p = proof.clone();
    p.witness[3] = 1;
    mutations.push(p);
    for (atom, digits, cells) in [
        (0, 3, 1),
        (6, 3, 1),
        (1, 0x200, 1),
        (1, 1, 1),
        (1, 3, 0),
        (1, 3, 4),
        (1, 3, 1),
    ] {
        let mut p = proof.clone();
        p.blockers[0] = HallObstruction {
            atom,
            digits,
            cells,
        };
        mutations.push(p);
    }
    let mut p = proof.clone();
    p.blockers.push(p.blockers[0].clone());
    mutations.push(p);
    for p in mutations {
        assert!(!p.verify(&grid), "tampered proof accepted: {p:?}");
    }
    let mut changed = grid.deep_clone();
    changed.cell_mut(Position::new(8, 8)).remove_candidate(9);
    assert!(!proof.verify(&changed));
    let mut placed = grid.deep_clone();
    placed.set_cell_unchecked(idx_to_pos(3), Some(1));
    assert!(!proof.verify(&placed));
    let mut json = serde_json::to_value(&proof).unwrap();
    json["unknown"] = 1.into();
    assert!(serde_json::from_value::<ThreeHouseProof>(json).is_err());
}

#[test]
fn replay_is_portable_and_atomic() {
    let (mut grid, proof) = certificate();
    let replay = ThreeHouseReplay::capture(&grid, &proof).unwrap();
    assert!(replay.verify());
    let decoded: ThreeHouseReplay =
        serde_json::from_str(&serde_json::to_string(&replay).unwrap()).unwrap();
    assert_eq!(decoded, replay);
    assert_eq!(decoded.check().unwrap(), proof.check(&grid).unwrap());
    let original = CandidateState::capture(&grid).unwrap();
    let mut wrong = grid.deep_clone();
    wrong.cell_mut(idx_to_pos(80)).remove_candidate(9);
    let before = CandidateState::capture(&wrong).unwrap();
    assert!(replay.apply(&mut wrong).is_err());
    assert_eq!(CandidateState::capture(&wrong).unwrap(), before);
    replay.apply(&mut grid).unwrap();
    let after = CandidateState::capture(&grid).unwrap();
    for cell in 0..81 {
        assert_eq!(
            after.domains[cell],
            if cell == 3 {
                original.domains[cell] & !1
            } else {
                original.domains[cell]
            }
        );
    }
    assert!(replay.apply(&mut grid).is_err());
    assert_eq!(CandidateState::capture(&grid).unwrap(), after);
    let mut malformed = replay.clone();
    malformed.version = 99;
    assert!(!malformed.verify());
    assert!(ThreeHouseReplay::capture(&grid, &proof).is_err());
    for field in ["values", "domains"] {
        let mut json = serde_json::to_value(&replay).unwrap();
        json["state"][field] = serde_json::json!([]);
        assert!(!serde_json::from_value::<ThreeHouseReplay>(json)
            .unwrap()
            .verify());
    }
}

#[test]
fn budgets_and_infeasible_relaxations_never_create_vacuous_proofs() {
    let solver = Solver::new();
    let blank = Grid::new_classic();
    let result = solver.search_three_house(&blank, &ThreeHouseSearchOptions { max_centers: 0 });
    assert!(result.hint.is_none() && !result.exhausted && result.error.is_none());
    assert_eq!(result.centers_checked, 0);
    let result = solver.search_three_house(&blank, &ThreeHouseSearchOptions { max_centers: 82 });
    assert!(result.error.is_some());
    let start = std::time::Instant::now();
    let result = solver.search_three_house(&blank, &ThreeHouseSearchOptions::default());
    assert!(result.hint.is_none() && result.exhausted && result.error.is_none());
    assert_eq!(result.centers_checked, 81);
    println!("blank full search: {:?}", start.elapsed());
    let mut bad = blank.deep_clone();
    // Every sector still covers all digits, but these two cells cannot differ.
    for cell in [0, 1] {
        bad.cell_mut(idx_to_pos(cell))
            .set_candidates(BitSet::single(1));
    }
    assert!(CandidateState::capture(&bad).is_ok());
    assert!(ThreeHouseProof::for_elimination(&bad, 0, 3, 1).is_none());
    let result = solver.search_three_house(&bad, &ThreeHouseSearchOptions { max_centers: 1 });
    assert!(result.hint.is_none() && result.error.is_none());
    bad.cell_mut(idx_to_pos(0)).set_candidates(BitSet::empty());
    assert!(solver
        .search_three_house(&bad, &ThreeHouseSearchOptions::default())
        .error
        .is_some());
    assert!(ThreeHouseProof::for_elimination(&bad, 0, 3, 1).is_none());
    assert!(solver
        .search_three_house(
            &Grid::new_with_constraints(vec![]),
            &ThreeHouseSearchOptions::default()
        )
        .error
        .is_some());
}

#[test]
fn native_rule_provenance_variants_and_atomic_variant_rejection() {
    let (grid, proof) = certificate();
    let state = CandidateState::capture(&grid).unwrap();
    let mut variant = Grid::new_x_sudoku();
    for (cell, &mask) in state.domains.iter().enumerate() {
        variant
            .cell_mut(idx_to_pos(cell))
            .set_candidates(BitSet::from_raw(mask << 1));
    }
    // Added restrictions preserve the native theorem; native premises hash alike.
    assert!(proof.verify(&variant));
    ThreeHouseReplay::capture(&variant, &proof)
        .unwrap()
        .apply(&mut variant)
        .unwrap();
    let deserialized: Grid = serde_json::from_str(&serde_json::to_string(&grid).unwrap()).unwrap();
    assert!(!proof.verify(&deserialized)); // constraints were not reinstalled
    assert!(Solver::new().get_three_house_hint(&deserialized).is_none());
    // A pre-existing variant conflict is outside the three selected houses.
    // The native conditional proof still checks; application must reject it.
    let mut variant = Grid::new_x_sudoku();
    for (cell, &mask) in state.domains.iter().enumerate() {
        variant
            .cell_mut(idx_to_pos(cell))
            .set_candidates(BitSet::from_raw(mask << 1));
    }
    variant.set_cell_unchecked(idx_to_pos(40), Some(1));
    variant.set_cell_unchecked(idx_to_pos(80), Some(1));
    let p = ThreeHouseProof::for_elimination(&variant, 0, 3, 1).unwrap();
    let replay = ThreeHouseReplay::capture(&variant, &p).unwrap();
    let before = CandidateState::capture(&variant).unwrap();
    assert!(replay.apply(&mut variant).is_err());
    assert_eq!(CandidateState::capture(&variant).unwrap(), before);
    // Native state outside the cluster can lack a digit after deletion, even
    // when the three-house relaxation itself has a valid opposite witness.
    let mut inconsistent = grid.deep_clone();
    for cell in 4..9 {
        inconsistent.cell_mut(idx_to_pos(cell)).remove_candidate(1);
    }
    for cell in [12, 13, 14, 21, 22, 23] {
        inconsistent.cell_mut(idx_to_pos(cell)).remove_candidate(1);
    }
    // Box 2 has only the target's 1. Native relaxed support still refutes it.
    let p = ThreeHouseProof::for_elimination(&inconsistent, 0, 3, 1).unwrap();
    let replay = ThreeHouseReplay::capture(&inconsistent, &p).unwrap();
    let before = CandidateState::capture(&inconsistent).unwrap();
    assert!(replay.apply(&mut inconsistent).is_err());
    assert_eq!(CandidateState::capture(&inconsistent).unwrap(), before);
}
