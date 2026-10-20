//! Read-only baseline probe, compiled only in a temporary engine copy.
use super::*;
use crate::{BitSet, Grid, Position};
#[test]
#[allow(deprecated)]
fn structural_old_detector_comparison() {
 let fixture:serde_json::Value=serde_json::from_str(include_str!("structural_four_digit_13.json")).unwrap();
 for state_name in ["base_domains","assumed_domains"] {
  let mut grid=Grid::new_classic();
  for (cell,mask) in fixture[state_name].as_array().unwrap().iter().enumerate(){
   grid.cell_mut(Position::new(cell/9,cell%9)).set_candidates(BitSet::from_raw((mask.as_u64().unwrap() as u16)<<1));
  }
  let before=serde_json::to_value(&grid).unwrap();
  let fab=CandidateFabric::from_grid(&grid);
  let graph=aic_engine::build_link_graph(&fab);
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

  let findings:Vec<_>=checks.into_iter().flatten().map(|f|format!("{:?}",f.technique)).collect();
  let arithmetic=Solver::new().search_arithmetic(&grid,&ArithmeticSearchOptions::default());
  let three_house=Solver::new().search_three_house(&grid,&ThreeHouseSearchOptions::default());
  assert!(arithmetic.error.is_none(),"{:?}",arithmetic.error);
  assert!(three_house.error.is_none(),"{:?}",three_house.error);
  assert_eq!(serde_json::to_value(&grid).unwrap(),before);
  println!("STRUCTURAL_BASELINE {}",serde_json::json!({"state":state_name,"explicit_old_detectors":34,"findings":findings,"arithmetic_hint":arithmetic.hint.map(|h|format!("{:?}",h.technique)),"arithmetic_combinations":arithmetic.tested_combinations,"arithmetic_budget_exhausted":arithmetic.budget_exhausted,"arithmetic_beam_pruned":arithmetic.beam_pruned,"three_house_hint":three_house.hint.map(|h|format!("{:?}",h.technique)),"centers_checked":three_house.centers_checked}));
 }
}
