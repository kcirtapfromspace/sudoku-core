//! Run with: cargo run --release --example arithmetic_counting
use sudoku_core::{
    ArithmeticProof, ArithmeticRequirement, ArithmeticSearchOptions, ArithmeticTerm, Grid,
    Position, ProofCertificate, Solver,
};

fn main() {
    // A candidate-state fixture for the familiar five-equation Guardian.
    let mut grid = Grid::new_classic();
    for (sector, retained) in [
        (0, vec![0, 3]),
        (12, vec![3, 30]),
        (3, vec![30, 28]),
        (10, vec![28, 10]),
        (18, vec![0, 10, 20]),
    ] {
        for cell in 0..81 {
            let belongs = match sector {
                0..=8 => cell / 9 == sector,
                9..=17 => cell % 9 == sector - 9,
                _ => cell / 27 * 3 + cell % 9 / 3 == sector - 18,
            };
            if belongs && !retained.contains(&cell) {
                grid.cell_mut(Position::new(cell / 9, cell % 9))
                    .remove_candidate(1);
            }
        }
    }

    let terms = [0, 12, 3, 10, 18]
        .into_iter()
        .map(|sector| ArithmeticTerm {
            requirement: ArithmeticRequirement::SectorDigit { sector, digit: 1 },
            weight: 1,
        })
        .collect();
    let guardian = ArithmeticProof::from_terms(&grid, terms, 20, 1, true)
        .expect("the five equations force r3c3=1");
    println!("Guardian certificate: {:?}", guardian.check(&grid).unwrap());

    let result = Solver::new().search_arithmetic(&grid, &ArithmeticSearchOptions::default());
    println!(
        "Search: {} combinations, budget exhausted={}, beam pruned={}",
        result.tested_combinations, result.budget_exhausted, result.beam_pruned
    );
    let hint = result
        .hint
        .expect("this fixture has a discoverable deduction");
    let Some(ProofCertificate::Arithmetic(proof)) = hint.proof else {
        panic!("the arithmetic engine must supply its certificate");
    };
    assert!(proof.verify(&grid));
    println!("{}: {}", hint.technique, hint.explanation);
}
