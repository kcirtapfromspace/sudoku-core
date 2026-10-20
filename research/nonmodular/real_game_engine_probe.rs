//! Public API checks for givens-derived positive fixtures. Run through
//! `python3 research/nonmodular/real_game_engine_check.py`.
use serde_json::{json, Value};
use std::{fs, time::Instant};
use sudoku_core::{ArithmeticSearchOptions, CandidateState, HintType, Position, ProofCertificate, Solver,
    ThreeHouseProof, ThreeHouseReplay, ThreeHouseSearchOptions};

fn main() {
    let file = std::env::args().nth(1).expect("fixture path");
    let fixtures: Value = serde_json::from_str(&fs::read_to_string(file).unwrap()).unwrap();
    let mut checked = Vec::new();
    for case in fixtures["cases"].as_array().unwrap() {
        let state: CandidateState = serde_json::from_value(case["state"].clone()).unwrap();
        let mut grid = state.to_grid().unwrap();
        // CandidateState intentionally omits clue flags; restore them for the
        // existing uniqueness detectors used only by this comparison probe.
        for (cell, byte) in case["puzzle"].as_str().unwrap().bytes().enumerate() {
            grid.cell_mut(Position::new(cell/9,cell%9)).set_given((b'1'..=b'9').contains(&byte));
        }
        let solution: Vec<u8> = case["solution"].as_str().unwrap().bytes().map(|b|b-b'0').collect();
        assert_eq!(state, CandidateState::capture(&grid).unwrap());
        let baseline_finding=Solver::new().research_baseline_finding(&grid,true);
        assert!(baseline_finding.is_none(), "stored snapshot is no longer a restricted baseline fixed point: {:?}",baseline_finding);
        let arithmetic_started=Instant::now();
        let arithmetic=Solver::new().search_arithmetic(&grid,&ArithmeticSearchOptions::default());
        let arithmetic_seconds=arithmetic_started.elapsed().as_secs_f64();
        let arithmetic_hint=arithmetic.hint.as_ref().map(|h|serde_json::to_value(h).unwrap());
        let started=Instant::now();
        let result=Solver::new().search_three_house(&grid, &ThreeHouseSearchOptions{max_centers:81});
        let search_seconds=started.elapsed().as_secs_f64();
        assert!(result.error.is_none());
        assert!(!result.exhausted);
        assert_eq!(state, CandidateState::capture(&grid).unwrap(), "search mutated premises");
        let hint=result.hint.expect("positive Python fixture needs production deduction");
        let proof=match hint.proof.as_ref().unwrap() {
            ProofCertificate::ThreeHouse(p) => p,
            _ => panic!("wrong certificate family"),
        };
        let check=proof.check(&grid).expect("independent Hall checker rejected proof");
        assert_ne!(solution[proof.cell],proof.digit);
        match hint.hint_type {
            HintType::EliminateCandidates{pos,values} => {
                assert_eq!(pos.row*9+pos.col,proof.cell);
                assert_eq!(values,vec![proof.digit]);
            }
            _ => panic!("expected elimination"),
        }
        let replay=ThreeHouseReplay::capture(&grid,proof).unwrap();
        if checked.is_empty() {
            if let Some(path) = std::env::args().nth(2) {
                fs::write(path, serde_json::to_string_pretty(&replay).unwrap() + "\n").unwrap();
            }
        }
        let encoded=serde_json::to_string(&replay).unwrap();
        let decoded: ThreeHouseReplay=serde_json::from_str(&encoded).unwrap();
        assert!(decoded.verify());
        let mut next=grid.deep_clone();
        decoded.apply(&mut next).unwrap();
        let after=CandidateState::capture(&next).unwrap();
        for cell in 0..81 {
            let expected=if cell==proof.cell {state.domains[cell]&!(1<<(proof.digit-1))}
                else {state.domains[cell]};
            assert_eq!(after.domains[cell],expected);
            assert_eq!(after.values[cell],state.values[cell]);
        }
        assert!(decoded.apply(&mut next).is_err(), "stale replay accepted");
        let mut individual_proofs=0;
        for cluster in case["three_house"].as_array().unwrap() {
            let center=cluster["center"].as_u64().unwrap() as usize;
            for elimination in cluster["eliminations"].as_array().unwrap() {
                let cell=elimination["cell"].as_u64().unwrap() as usize;
                let mask=elimination["mask"].as_u64().unwrap() as u16;
                for digit in 1u8..=9 {
                    if mask & (1 << (digit-1)) == 0 {continue;}
                    assert_ne!(solution[cell],digit);
                    let p=ThreeHouseProof::for_elimination(&grid,center,cell,digit)
                        .expect("production could not certify Python-supported deduction");
                    assert!(p.verify(&grid));
                    individual_proofs+=1;
                }
            }
        }
        checked.push(json!({"id":case["id"],"centers_checked":result.centers_checked,
            "search_seconds":search_seconds,"first_proof":{"center":proof.center,
            "cell":proof.cell,"digit":proof.digit,"signatures":check.signatures,
            "hall_obstructions":check.hall_obstructions},
            "individually_verified_center_candidate_proofs":individual_proofs,
            "restricted_baseline_fixed_point_rechecked":true,
            "default_arithmetic":{"hint":arithmetic_hint,"seconds":arithmetic_seconds,
                "tested_combinations":arithmetic.tested_combinations,"budget_exhausted":arithmetic.budget_exhausted,
                "beam_pruned":arithmetic.beam_pruned,"error":arithmetic.error},
            "input_unchanged":true,"json_replay_applied":true,"stale_replay_rejected":true}));
    }
    println!("{}",serde_json::to_string_pretty(&json!({"cases":checked})).unwrap());
}
