//! Read a portable CandidateState JSON file and export one checked deduction.
//! cargo run --release --example global_coloring -- candidate-state.json
use std::{env, fs, io};
use sudoku_core::{
    CandidateState, GlobalColoringReplay, GlobalColoringSearchOptions, ProofCertificate, Solver,
};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let path = env::args()
        .nth(1)
        .ok_or_else(|| io::Error::other("Pass a CandidateState JSON path"))?;
    let state: CandidateState = serde_json::from_str(&fs::read_to_string(path)?)?;
    let grid = state.to_grid()?;
    let result =
        Solver::new().search_global_coloring(&grid, &GlobalColoringSearchOptions::default());
    let hint = result.hint.ok_or_else(|| io::Error::other(format!(
        "No certified hint: {} patterns, {} visited nodes, budget exhausted={}, oversized={}, error={:?}",
        result.patterns_tested, result.nodes_visited, result.budget_exhausted, result.oversized_components, result.error)))?;
    let proof = match hint.proof.as_ref() {
        Some(ProofCertificate::GlobalColoring(proof)) => proof,
        _ => return Err(io::Error::other("Unexpected certificate family").into()),
    };
    let replay = GlobalColoringReplay::capture(&grid, proof)?;
    println!("{}", serde_json::to_string_pretty(&replay)?);
    Ok(())
}
