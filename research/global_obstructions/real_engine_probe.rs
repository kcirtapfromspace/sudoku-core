//! Compile as an example in a temporary copy; uses only the public engine API.
use serde_json::{json, Value};
use std::io::{self, BufRead};
use std::time::Instant;
use sudoku_core::{
    CandidateState, GlobalColoringProof, GlobalColoringReplay,
    GlobalColoringSearchOptions, Position, ProofCertificate, Solver,
};

fn main() {
    for line in io::stdin().lock().lines() {
        let case: Value = serde_json::from_str(&line.unwrap()).unwrap();
        let state: CandidateState = serde_json::from_value(case["state"].clone()).unwrap();
        let grid = state.to_grid().unwrap();
        let cells: Vec<usize> = serde_json::from_value(case["graph_cells"].clone()).unwrap();
        let cell = case["conclusion"]["cell"].as_u64().unwrap() as usize;
        let digit = case["conclusion"]["digit"].as_u64().unwrap() as u8;
        assert!(case["conclusion"]["value"].as_bool().unwrap());
        assert_eq!(case["solution"].as_str().unwrap().as_bytes()[cell] - b'0', digit);
        let proof = GlobalColoringProof::from_cells(&grid, cells, cell, digit, true)
            .expect("explicit fixture graph must produce a placement certificate");
        let check = proof.check(&grid).unwrap();
        let replay = GlobalColoringReplay::capture(&grid, &proof).unwrap();
        let encoded = serde_json::to_string(&replay).unwrap();
        let decoded: GlobalColoringReplay = serde_json::from_str(&encoded).unwrap();
        assert!(decoded.verify());
        let mut applied = grid.clone();
        decoded.apply(&mut applied).unwrap();
        assert_eq!(applied.get(Position::new(cell / 9, cell % 9)), Some(digit));
        // Applying a stale certificate must fail without changing the board.
        let after = CandidateState::capture(&applied).unwrap();
        assert!(decoded.apply(&mut applied).is_err());
        assert_eq!(CandidateState::capture(&applied).unwrap(), after);
        let started = Instant::now();
        let found = Solver::new().search_global_coloring(
            &grid, &GlobalColoringSearchOptions::default());
        let seconds = started.elapsed().as_secs_f64();
        assert!(found.error.is_none());
        let automatic = found.hint.as_ref().expect("default search must find a fixture proof");
        let Some(ProofCertificate::GlobalColoring(automatic_proof)) = &automatic.proof else {
            panic!("automatic hint omitted its certificate");
        };
        assert!(automatic_proof.verify(&grid));
        let solution_digit = case["solution"].as_str().unwrap().as_bytes()[automatic_proof.cell] - b'0';
        assert_eq!(solution_digit == automatic_proof.digit, automatic_proof.value);
        assert_eq!(CandidateState::capture(&grid).unwrap(), state);
        println!("{}", json!({"name":case["name"],"explicit_nodes":check.nodes,
            "explicit_splits":check.splits,"explicit_conflicts":check.conflicts,
            "replay_bytes":encoded.len(),"replay_roundtrip":true,"stale_replay_rejected":true,
            "input_unchanged":true,"automatic_patterns":found.patterns_tested,
            "automatic_nodes":found.nodes_visited,"automatic_seconds":seconds,
            "automatic_cell":automatic_proof.cell,"automatic_digit":automatic_proof.digit,
            "automatic_value":automatic_proof.value,"budget_exhausted":found.budget_exhausted}));
    }
}
