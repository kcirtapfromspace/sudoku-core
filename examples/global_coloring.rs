//! Read a portable CandidateState JSON file and export one checked deduction.
//! cargo run --release --example global_coloring -- candidate-state.json
//! cargo run --release --example global_coloring -- --explicit proof-request.json
use std::{env, fs, io};
use sudoku_core::{
    CandidateState, GlobalColoringProof, GlobalColoringReplay, GlobalColoringSearchOptions,
    ProofCertificate, Solver,
};

#[derive(serde::Deserialize)]
#[serde(deny_unknown_fields)]
struct ExplicitRequest {
    state: CandidateState,
    cells: Vec<usize>,
    cell: usize,
    digit: u8,
    value: bool,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut args = env::args().skip(1);
    let path = args
        .next()
        .ok_or_else(|| io::Error::other("Pass a CandidateState JSON path"))?;
    if path == "--explicit" {
        let path = args
            .next()
            .ok_or_else(|| io::Error::other("Pass an explicit proof request JSON path"))?;
        let request: ExplicitRequest = serde_json::from_str(&fs::read_to_string(path)?)?;
        let grid = request.state.to_grid()?;
        let proof = GlobalColoringProof::from_cells(
            &grid,
            request.cells,
            request.cell,
            request.digit,
            request.value,
        )
        .ok_or_else(|| io::Error::other("No complete explicit certificate within the bounds"))?;
        let replay = GlobalColoringReplay::capture(&grid, &proof)?;
        println!("{}", serde_json::to_string_pretty(&replay)?);
        return Ok(());
    }
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
