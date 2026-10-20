//! Offline named-detector coverage search. Successful states still need certificates.
use serde_json::{json, Value};
use std::{collections::HashSet, env, fs};
use sudoku_core::{
    CandidateState, ForcingSearchOptions, Grid, HintType, Position, Solver, Technique,
};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<_> = env::args().collect();
    let input = args
        .get(1)
        .ok_or("Pass input state records and an output directory")?;
    let out = args.get(2).ok_or("Pass an output directory")?;
    fs::create_dir_all(out)?;
    let records: Vec<Value> = serde_json::from_str(&fs::read_to_string(input)?)?;
    let mut found = HashSet::new();
    for t in Technique::all() {
        if std::path::Path::new(&format!("{out}/{t:?}.json")).exists() {
            found.insert(*t);
        }
    }
    let solver = Solver::new();
    let options = ForcingSearchOptions {
        max_depth: 2,
        max_nodes: 4096,
        max_steps: 32768,
        ..Default::default()
    };
    for (i, record) in records.iter().enumerate() {
        if found.len() == Technique::all().len() {
            break;
        }
        let puzzle = record["puzzle"].as_str().ok_or("No original puzzle")?;
        let Some(givens) = Grid::from_string(puzzle) else {
            continue;
        };
        if solver.count_solutions(&givens, 2) != 1 {
            continue;
        }
        let solution = solver.solve(&givens).unwrap();
        let state: CandidateState = serde_json::from_value(record["state"].clone())?;
        if (0..81).any(|c| {
            let p = Position::new(c / 9, c % 9);
            let d = solution.get(p).unwrap();
            state.values[c] != 0 && state.values[c] != d || state.domains[c] & (1 << (d - 1)) == 0
        }) {
            continue;
        }
        let Ok(mut grid) = state.to_grid() else {
            continue;
        };
        for c in 0..81 {
            let p = Position::new(c / 9, c % 9);
            grid.cell_mut(p).set_given(givens.get(p).is_some());
        }
        for &technique in Technique::all() {
            if found.contains(&technique) {
                continue;
            }
            let Some(hint) = solver.get_hint_for_technique(&grid, technique, &options)? else {
                continue;
            };
            let sound = match &hint.hint_type {
                HintType::SetValue { pos, value } => solution.get(*pos) == Some(*value),
                HintType::EliminateCandidates { pos, values } => {
                    !values.contains(&solution.get(*pos).unwrap())
                }
            };
            if !sound {
                eprintln!("Rejected false endpoint: {technique:?} / {}", record["id"]);
                continue;
            }
            let pattern = hint.proof.clone();
            fs::write(
                format!("{out}/{technique:?}.json"),
                serde_json::to_string_pretty(
                    &json!({"source":record,"hint":hint,"pattern":pattern}),
                )?,
            )?;
            found.insert(technique);
            eprintln!("{}/48 {technique:?} in {}", found.len(), record["id"]);
        }
        if i % 25 == 0 {
            eprintln!("Scanned {i}: {} found", found.len());
        }
    }
    eprintln!(
        "Missing: {:?}",
        Technique::all()
            .iter()
            .filter(|t| !found.contains(t))
            .collect::<Vec<_>>()
    );
    Ok(())
}
