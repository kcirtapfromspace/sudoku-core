//! Research comparison only; compiled in a temporary copy by compare_engine.py.
//! Uniqueness, forcing/backtracking, and optional partition-residue APIs are
//! outside this precisely bounded baseline. No claim of solver completeness.
use super::*;
use crate::{BitSet, Grid, Position};

#[test]
#[allow(deprecated)]
fn compare_stored_candidate_fixture() {
    let state: serde_json::Value =
        serde_json::from_str(include_str!("three-house-847.json")).unwrap();
    let mut grid = Grid::new_classic();
    for (cell, mask) in state["domains"].as_array().unwrap().iter().enumerate() {
        grid.cell_mut(Position::new(cell / 9, cell % 9))
            .set_candidates(BitSet::from_raw((mask.as_u64().unwrap() as u16) << 1));
    }
    let before = serde_json::to_value(&grid).unwrap();
    let fab = CandidateFabric::from_grid(&grid);
    let graph = aic_engine::build_link_graph(&fab);
    let checks = [
        basic::find_naked_single(&fab),
        basic::find_hidden_single(&fab),
        basic::find_naked_subset(&fab, 2),
        basic::find_hidden_subset(&fab, 2),
        basic::find_naked_subset(&fab, 3),
        basic::find_hidden_subset(&fab, 3),
        basic::find_naked_subset(&fab, 4),
        basic::find_hidden_subset(&fab, 4),
        fish_engine::find_pointing_pair(&fab),
        fish_engine::find_box_line_reduction(&fab),
        fish_engine::find_basic_fish(&fab, 2),
        fish_engine::find_finned_fish(&fab, 2),
        fish_engine::find_basic_fish(&fab, 3),
        fish_engine::find_finned_fish(&fab, 3),
        fish_engine::find_basic_fish(&fab, 4),
        fish_engine::find_finned_fish(&fab, 4),
        aic_engine::find_empty_rectangle(&fab),
        als_engine::find_xy_wing(&fab),
        als_engine::find_xyz_wing(&fab),
        als_engine::find_wxyz_wing(&fab),
        aic_engine::find_w_wing(&fab),
        aic_engine::find_x_chain(&fab, &graph),
        aic_engine::find_medusa(&fab, &graph),
        als_engine::find_sue_de_coq(&fab),
        aic_engine::find_aic(&fab, &graph),
        fish_engine::find_franken_fish(&fab),
        fish_engine::find_siamese_fish(&fab),
        als_engine::find_als_xz(&fab),
        als_engine::find_als_xy_wing(&fab),
        als_engine::find_als_chain(&fab),
        fish_engine::find_mutant_fish(&fab),
        als_engine::find_aligned_pair_exclusion(&fab),
        als_engine::find_aligned_triplet_exclusion(&fab),
        als_engine::find_death_blossom(&fab),
    ];
    let checked = checks.len();
    let findings: Vec<_> = checks
        .into_iter()
        .flatten()
        .map(|finding| format!("{:?}", finding.technique))
        .collect();
    assert!(findings.is_empty(), "Baseline now finds: {findings:?}");
    let result = Solver::new().search_arithmetic(&grid, &ArithmeticSearchOptions::default());
    assert!(result.error.is_none());
    assert!(
        result.hint.is_none(),
        "Arithmetic baseline now finds a step"
    );
    let added = Solver::new().get_three_house_hint(&grid).expect("integrated detector closes the gap");
    let proof = match added.proof.as_ref().unwrap() {
        ProofCertificate::ThreeHouse(proof) => proof,
        other => panic!("unexpected proof: {other:?}"),
    };
    assert!(proof.verify(&grid));
    println!("Integrated Three-House Matching verifies: cell {} != {}", proof.cell, proof.digit);
    assert_eq!(serde_json::to_value(&grid).unwrap(), before);
    println!("No finding from {checked} tested non-forcing, uniqueness-independent detectors.");
    println!(
        "Default arithmetic: {} combinations, exhausted={}, beam_pruned={}",
        result.tested_combinations, result.budget_exhausted, result.beam_pruned
    );
}
