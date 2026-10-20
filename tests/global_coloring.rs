//! Givens-derived published and generated controls, with portable exact-state proofs.
use sudoku_core::{
    CandidateState, GlobalColoringProof, GlobalColoringReplay, GlobalColoringSearchOptions, Grid,
    HintType, Position, ProofCertificate, Solver, Technique,
};

const CASES: [&str; 7] = [
    include_str!("../research/global_obstructions/real_4662_after_existing.json"),
    include_str!("../research/global_obstructions/real_4401.json"),
    include_str!("../research/global_obstructions/real_Loki.json"),
    include_str!("../research/global_obstructions/real_four_digit_candidate.json"),
    include_str!("../research/global_obstructions/real_four_digit_unique.json"),
    include_str!("../research/global_obstructions/expanded_unique_nonforcing.json"),
    include_str!("../research/global_obstructions/nonrank_full_old_stall.json"),
];

#[test]
fn global_obstructions_certify_both_published_and_generated_game_states() {
    for text in CASES {
        let case: serde_json::Value = serde_json::from_str(text).unwrap();
        let state: CandidateState = serde_json::from_value(case["state"].clone()).unwrap();
        let grid = state.to_grid().unwrap();
        let givens = Grid::from_string(case["puzzle"].as_str().unwrap()).unwrap();
        let expected_count = case["solutions_capped_at_2"].as_u64().unwrap() as usize;
        assert_eq!(Solver::new().count_solutions(&givens, 2), expected_count);
        let completion = case["solution"]
            .as_str()
            .or_else(|| case["completion"].as_str())
            .unwrap();
        let solved = Grid::from_string(completion).unwrap();
        assert!(solved.is_complete() && solved.validate().is_valid);
        let digits: Vec<u8> = completion.bytes().map(|b| b - b'0').collect();
        for (c, &d) in digits.iter().enumerate() {
            assert_ne!(state.domains[c] & (1 << (d - 1)), 0);
            if let Some(given) = givens.get(Position::new(c / 9, c % 9)) {
                assert_eq!(given, d);
            }
        }
        let cells: Vec<usize> = serde_json::from_value(case["graph_cells"].clone()).unwrap();
        let cell = case["conclusion"]["cell"].as_u64().unwrap() as usize;
        let digit = case["conclusion"]["digit"].as_u64().unwrap() as u8;
        let proof = GlobalColoringProof::from_cells(&grid, cells, cell, digit, true).unwrap();
        assert!(proof.verify(&grid));
        assert_eq!(digits[cell], digit);
        let replay = GlobalColoringReplay::capture(&grid, &proof).unwrap();
        let decoded: GlobalColoringReplay =
            serde_json::from_str(&serde_json::to_string(&replay).unwrap()).unwrap();
        assert!(decoded.verify());
        let mut next = grid.deep_clone();
        decoded.apply(&mut next).unwrap();
        assert_eq!(next.get(Position::new(cell / 9, cell % 9)), Some(digit));
        assert_eq!(CandidateState::capture(&grid).unwrap(), state);
        for c in 0..81 {
            if next.get(Position::new(c / 9, c % 9)).is_none() {
                assert_ne!(
                    next.get_candidates(Position::new(c / 9, c % 9)).as_raw() & (1 << digits[c]),
                    0
                );
                assert_eq!(
                    (next.get_candidates(Position::new(c / 9, c % 9)).as_raw() >> 1)
                        & !state.domains[c],
                    0
                );
            }
        }
        let after = CandidateState::capture(&next).unwrap();
        assert!(decoded.apply(&mut next).is_err());
        assert_eq!(CandidateState::capture(&next).unwrap(), after);
        let result =
            Solver::new().search_global_coloring(&grid, &GlobalColoringSearchOptions::default());
        assert!(result.error.is_none());
        // The unique generated puzzle is an explicit old-method-solvable control;
        // a bounded dedicated scan may still find another valid local obstruction.
        {
            let hint = result.hint.expect("bounded scan must certify this fixture");
            assert_eq!(hint.technique, Technique::GlobalColoring);
            let Some(ProofCertificate::GlobalColoring(p)) = &hint.proof else {
                panic!("wrong proof family")
            };
            assert!(p.verify(&grid));
            match hint.hint_type {
                HintType::SetValue { pos, value } => {
                    assert_eq!(digits[pos.row * 9 + pos.col], value);
                    assert!(p.value);
                }
                HintType::EliminateCandidates { pos, values } => {
                    assert!(!values.contains(&digits[pos.row * 9 + pos.col]));
                    assert!(!p.value);
                }
            }
        }
    }
}

#[test]
fn committed_exports_verify_and_cannot_replay_after_application() {
    for text in [
        include_str!("../research/global_obstructions/replay_4662.json"),
        include_str!("../research/global_obstructions/replay_4401.json"),
        include_str!("../research/global_obstructions/replay_Loki.json"),
        include_str!("../research/global_obstructions/replay_four_digit.json"),
        include_str!("../research/global_obstructions/replay_four_digit_generated.json"),
        include_str!("../research/global_obstructions/replay_four_digit_unique.json"),
        include_str!("../research/global_obstructions/replay_expanded_unique_nonforcing.json"),
        include_str!("../research/global_obstructions/replay_nonrank_full_old_stall.json"),
        include_str!("../research/global_obstructions/replay_native27.json"),
    ] {
        let replay: GlobalColoringReplay = serde_json::from_str(text).unwrap();
        assert!(replay.verify());
        let mut grid = replay.state.to_grid().unwrap();
        replay.apply(&mut grid).unwrap();
        assert!(replay.apply(&mut grid).is_err());
    }
}

#[test]
fn normal_hint_ladder_reaches_the_global_certificate_on_published_givens() {
    let case: serde_json::Value = serde_json::from_str(CASES[1]).unwrap();
    let grid = Grid::from_string(case["puzzle"].as_str().unwrap()).unwrap();
    let hint = Solver::new().get_hint(&grid).unwrap();
    assert_eq!(hint.technique, Technique::GlobalColoring);
    assert!(matches!(hint.hint_type,HintType::SetValue{pos,value:9} if pos==Position::new(3,3)));
    let Some(ProofCertificate::GlobalColoring(proof)) = hint.proof else {
        panic!("missing checked global certificate")
    };
    assert!(proof.verify(&grid));
}
