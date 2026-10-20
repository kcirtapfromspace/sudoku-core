//! Public, givens-derived regressions: exact baseline traces and all recorded
//! center/candidate deductions, rather than only synthetic pencil-mark states.
use sudoku_core::{
    CandidateState, Grid, HintType, Position, ProofCertificate, Solver, Technique, ThreeHouseProof,
    ThreeHouseReplay,
};

#[test]
fn published_puzzles_replay_to_the_certified_candidate_states() {
    let fixtures: serde_json::Value =
        serde_json::from_str(include_str!("../research/nonmodular/real_game_cases.json")).unwrap();
    let exported: ThreeHouseReplay = serde_json::from_str(include_str!(
        "../research/nonmodular/three_house_example_replay.json"
    ))
    .unwrap();
    assert!(exported.verify());
    assert_eq!(
        serde_json::to_value(&exported.state).unwrap(),
        fixtures["cases"][0]["state"]
    );
    let mut count = 0;
    for case in fixtures["cases"].as_array().unwrap() {
        let id = case["id"].as_str().unwrap();
        let puzzle = Grid::from_string(case["puzzle"].as_str().unwrap()).unwrap();
        assert_eq!(Solver::new().count_solutions(&puzzle, 2), 1, "{id}");
        let solution: Vec<u8> = case["solution"]
            .as_str()
            .unwrap()
            .bytes()
            .map(|b| b - b'0')
            .collect();
        let solved = Grid::from_string(case["solution"].as_str().unwrap()).unwrap();
        assert!(solved.is_complete() && solved.validate().is_valid);
        let pos = |cell: usize| Position::new(cell / 9, cell % 9);
        let mut reached = puzzle.deep_clone();
        reached.recalculate_candidates();
        for step in case["baseline_trace"].as_array().unwrap() {
            let cell = step[1].as_u64().unwrap() as usize;
            let digits: Vec<u8> = step[2]
                .as_array()
                .unwrap()
                .iter()
                .map(|v| v.as_u64().unwrap() as u8)
                .collect();
            assert!(reached.get(pos(cell)).is_none());
            assert!(digits
                .iter()
                .all(|&d| reached.get_candidates(pos(cell)).contains(d)));
            if step[3].as_bool().unwrap() {
                assert_eq!(digits, vec![solution[cell]], "{id}");
                reached.set_cell_unchecked(pos(cell), Some(digits[0]));
                reached.update_candidates_after_move(pos(cell), digits[0]);
            } else {
                for digit in digits {
                    assert_ne!(solution[cell], digit, "{id}");
                    reached.cell_mut(pos(cell)).remove_candidate(digit);
                }
            }
        }
        let state: CandidateState = serde_json::from_value(case["state"].clone()).unwrap();
        assert_eq!(CandidateState::capture(&reached).unwrap(), state, "{id}");
        for cell in 0..81 {
            assert_ne!(state.domains[cell] & (1 << (solution[cell] - 1)), 0);
            if let Some(given) = puzzle.get(pos(cell)) {
                assert_eq!(given, solution[cell]);
            }
        }
        let grid = state.to_grid().unwrap();
        let hint = Solver::new().get_three_house_hint(&grid).expect(id);
        assert_eq!(hint.technique, Technique::ThreeHouseMatching);
        let proof = match hint.proof.as_ref().unwrap() {
            ProofCertificate::ThreeHouse(proof) => proof,
            _ => panic!("wrong proof family"),
        };
        match hint.hint_type {
            HintType::EliminateCandidates { pos, values } => {
                assert_eq!(pos.row * 9 + pos.col, proof.cell);
                assert_eq!(values, vec![proof.digit]);
            }
            _ => panic!("wrong hint kind"),
        }
        assert_ne!(proof.digit, solution[proof.cell]);
        assert!(proof.verify(&grid));
        let replay = ThreeHouseReplay::capture(&grid, proof).unwrap();
        let decoded: ThreeHouseReplay =
            serde_json::from_str(&serde_json::to_string(&replay).unwrap()).unwrap();
        assert!(decoded.verify());
        let mut next = grid.deep_clone();
        decoded.apply(&mut next).unwrap();
        let after = CandidateState::capture(&next).unwrap();
        for cell in 0..81 {
            let mask = if cell == proof.cell {
                state.domains[cell] & !(1 << (proof.digit - 1))
            } else {
                state.domains[cell]
            };
            assert_eq!(after.domains[cell], mask);
        }
        assert_eq!(CandidateState::capture(&grid).unwrap(), state);
        assert!(decoded.apply(&mut next).is_err());
        for cluster in case["three_house"].as_array().unwrap() {
            let center = cluster["center"].as_u64().unwrap() as usize;
            for elimination in cluster["eliminations"].as_array().unwrap() {
                let cell = elimination["cell"].as_u64().unwrap() as usize;
                let mask = elimination["mask"].as_u64().unwrap() as u16;
                for digit in 1..=9 {
                    if mask & (1 << (digit - 1)) == 0 {
                        continue;
                    }
                    assert_ne!(solution[cell], digit, "{id}");
                    let proof =
                        ThreeHouseProof::for_elimination(&grid, center, cell, digit).expect(id);
                    assert!(proof.verify(&grid));
                    count += 1;
                }
            }
        }
    }
    assert_eq!(count, 25);
}
