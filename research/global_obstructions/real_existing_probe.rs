//! Runs only in an instrumented temporary engine copy.
use serde_json::{json, Value};
use std::io::{self, BufRead};
use std::time::Instant;
use sudoku_core::{ArithmeticSearchOptions, BitSet, CandidateState, Position, Solver,
    ThreeHouseSearchOptions};

fn main() {
    for line in io::stdin().lock().lines() {
        let request: Value = serde_json::from_str(&line.unwrap()).unwrap();
        let state: CandidateState = serde_json::from_value(request["state"].clone()).unwrap();
        let mut grid = state.to_grid().unwrap();
        for (cell, byte) in request["puzzle"].as_str().unwrap().bytes().enumerate() {
            grid.cell_mut(Position::new(cell / 9, cell % 9))
                .set_given((b'1'..=b'9').contains(&byte));
        }
        let solver = Solver::new();
        let unique = request["unique"].as_bool().unwrap_or(true);
        let mut trace = Vec::new();
        let started = Instant::now();
        for _ in 0..1000 {
            let Some(step) = solver.research_baseline_finding(&grid, unique) else { break };
            let (technique, cell, digits, placement) = &step;
            let pos = Position::new(cell / 9, cell % 9);
            if *placement {
                let old = CandidateState::capture(&grid).unwrap();
                grid.set_cell_unchecked(pos, Some(digits[0]));
                grid.recalculate_candidates();
                for i in 0..81 {
                    let p = Position::new(i / 9, i % 9);
                    if grid.get(p).is_none() {
                        let keep = grid.get_candidates(p).as_raw() & (old.domains[i] << 1);
                        grid.cell_mut(p).set_candidates(BitSet::from_raw(keep));
                    }
                }
            } else {
                for &digit in digits { grid.cell_mut(pos).remove_candidate(digit); }
            }
            trace.push(json!([technique, cell, digits, placement]));
        }
        assert!(solver.research_baseline_finding(&grid, unique).is_none());
        let baseline_seconds = started.elapsed().as_secs_f64();
        let arithmetic_started = Instant::now();
        let arithmetic = solver.search_arithmetic(&grid, &ArithmeticSearchOptions::default());
        let arithmetic_seconds = arithmetic_started.elapsed().as_secs_f64();
        let triple = solver.search_three_house(&grid, &ThreeHouseSearchOptions { max_centers: 81 });
        let full_started = Instant::now();
        let full_executed = unique && request["full"].as_bool().unwrap_or(false);
        let full = if full_executed {
            solver.research_full_finding(&grid)
        } else { None };
        println!("{}", json!({"id":request["id"],
            "state":CandidateState::capture(&grid).unwrap(), "trace":trace,
            "baseline_seconds":baseline_seconds,"uniqueness_warranted":unique,
            "arithmetic":{"hint":arithmetic.hint,"tested_combinations":arithmetic.tested_combinations,
                "beam_pruned":arithmetic.beam_pruned,"budget_exhausted":arithmetic.budget_exhausted,
                "error":arithmetic.error,"seconds":arithmetic_seconds},
            "three_house":{"hint":triple.hint,"centers_checked":triple.centers_checked,
                "exhausted":triple.exhausted,"error":triple.error},
            "full_dispatcher_executed":full_executed,
            "full_current_finding":full,"full_current_seconds":full_started.elapsed().as_secs_f64()}));
    }
}
