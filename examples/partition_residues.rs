//! Run: cargo run --release --example partition_residues > partition-proof.json
use sudoku_core::{
    ArithmeticPartition, ArithmeticProof, ArithmeticReplay, ArithmeticRequirement,
    ArithmeticSearchOptions, ArithmeticTerm, ArithmeticTerminal, BitSet, Grid, Position,
    ProofCertificate, Solver,
};

fn main() {
    let mut grid = Grid::new_classic();
    // Row 1 places both digits 1 and 2 in A=r1c1 and B=r1c4.
    for col in 0..9 {
        if col != 0 && col != 3 {
            for digit in [1, 2] {
                grid.cell_mut(Position::new(0, col)).remove_candidate(digit);
            }
        }
    }
    grid.cell_mut(Position::new(0, 0))
        .set_candidates(BitSet::from_slice(&[1, 2, 3]));
    grid.cell_mut(Position::new(0, 3))
        .set_candidates(BitSet::from_slice(&[1, 2]));
    let terms: Vec<_> = [1, 2]
        .into_iter()
        .map(|digit| ArithmeticTerm {
            requirement: ArithmeticRequirement::SectorDigit { sector: 0, digit },
            weight: 1,
        })
        .collect();
    // The target x_A3 has coefficient zero: earlier terminals cannot distinguish
    // its endpoints from this equation, regardless of the Boolean modulus.
    assert!(ArithmeticProof::from_terms(&grid, terms.clone(), 0, 3, false).is_none());
    for modulus in 2..=16 {
        assert!(
            ArithmeticProof::from_residue_terms(&grid, terms.clone(), 0, 3, false, modulus)
                .is_none()
        );
    }
    let explicit = ArithmeticProof::from_terms_with_terminal(
        &grid,
        terms,
        0,
        3,
        false,
        ArithmeticTerminal::PartitionResidue {
            modulus: 2,
            partition: ArithmeticPartition::Cells,
        },
    )
    .expect("exactly one choice per cell excludes A=3");
    eprintln!("Explicit certificate: {:?}", explicit.check(&grid).unwrap());

    let options = ArithmeticSearchOptions {
        max_sources: 2,
        ..Default::default()
    };
    let result =
        Solver::new().search_partition_arithmetic(&grid, &options, ArithmeticPartition::Cells, 2);
    let hint = result
        .hint
        .expect("partition-aware search discovers a deduction");
    eprintln!(
        "{} combinations; {}",
        result.tested_combinations, hint.explanation
    );
    let Some(ProofCertificate::Arithmetic(proof)) = hint.proof else {
        panic!("missing certificate")
    };
    let replay = ArithmeticReplay::capture(&grid, &proof).unwrap();
    let json = serde_json::to_string_pretty(&replay).unwrap();
    let decoded: ArithmeticReplay = serde_json::from_str(&json).unwrap();
    assert!(decoded.verify());
    decoded.apply(&mut grid).unwrap();
    println!("{json}");
}
