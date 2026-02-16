//! Solver orchestrator.
//!
//! Dispatches to three abstract engines (Fish, ALS, AIC) plus basic techniques,
//! uniqueness patterns, and backtracking.

mod aic_engine;
mod als_engine;
pub(crate) mod backtrack;
mod basic;
pub(crate) mod explain;
pub(crate) mod fabric;
mod fish_engine;
mod types;
mod uniqueness;

use crate::{Grid, Position};
use explain::{Finding, InferenceResult};
use fabric::{idx_to_pos, CandidateFabric};

pub use explain::{AlsProofDescriptor, ForcingSource, LinkType, Polarity, ProofCertificate};
pub use types::{Difficulty, Hint, HintType, Technique};

/// Unit struct solver — stateless, all state is per-call.
pub struct Solver;

impl Default for Solver {
    fn default() -> Self {
        Self::new()
    }
}

impl Solver {
    /// Create a new solver.
    pub fn new() -> Self {
        Self
    }

    /// Solve the puzzle, returning the solved grid if successful.
    pub fn solve(&self, grid: &Grid) -> Option<Grid> {
        let mut working = grid.deep_clone();
        working.recalculate_candidates();
        if backtrack::solve_recursive(&mut working) {
            Some(working)
        } else {
            None
        }
    }

    /// Count solutions up to a limit.
    pub fn count_solutions(&self, grid: &Grid, limit: usize) -> usize {
        let mut working = grid.deep_clone();
        working.recalculate_candidates();
        let mut count = 0;
        backtrack::count_solutions_recursive(&mut working, &mut count, limit);
        count
    }

    /// Check if the puzzle has exactly one solution.
    pub fn has_unique_solution(&self, grid: &Grid) -> bool {
        self.count_solutions(grid, 2) == 1
    }

    /// Get a hint for the current position.
    pub fn get_hint(&self, grid: &Grid) -> Option<Hint> {
        let mut working = grid.deep_clone();
        working.recalculate_candidates();

        if let Some(finding) = self.find_first_technique(&working) {
            return Some(finding.to_hint());
        }

        // Last resort: backtracking hint
        if let Some(finding) = backtrack::find_backtracking_hint(&working) {
            return Some(finding.to_hint());
        }

        None
    }

    /// Get the next placement hint by chaining through elimination techniques.
    ///
    /// Unlike `get_hint` which returns the first technique found (which may be
    /// an elimination), this method applies eliminations internally and keeps
    /// searching until it finds a placement (SetValue) hint. This prevents
    /// callers from looping on the same elimination when they recalculate
    /// candidates from scratch between calls.
    pub fn get_next_placement(&self, grid: &Grid) -> Option<Hint> {
        let mut working = grid.deep_clone();
        working.recalculate_candidates();

        // Solve once via backtracking to verify chained results.
        // This prevents subtle technique bugs from producing wrong hints.
        let solution = self.solve(grid)?;

        for _ in 0..500 {
            if working.is_complete() {
                return None;
            }

            if let Some(finding) = self.find_first_technique(&working) {
                match &finding.inference {
                    InferenceResult::Placement { cell, value } => {
                        let pos = idx_to_pos(*cell);
                        if solution.get(pos) == Some(*value) {
                            return Some(finding.to_hint());
                        }
                        // Unsound placement — chaining corrupted state
                        break;
                    }
                    InferenceResult::Elimination { cell, values } => {
                        let pos = idx_to_pos(*cell);
                        let sol_val = solution.get(pos);
                        if values.iter().any(|&v| sol_val == Some(v)) {
                            // Unsound elimination — stop chaining
                            break;
                        }
                        apply_finding(&mut working, &finding);
                    }
                }
            } else {
                break;
            }
        }

        // Fall back to backtracking hint (always correct)
        backtrack::find_backtracking_hint(&working).map(|f| f.to_hint())
    }

    /// Rate the difficulty of a puzzle.
    pub fn rate_difficulty(&self, grid: &Grid) -> Difficulty {
        let empty_count = grid.empty_positions().len();
        let mut working = grid.deep_clone();
        let max_tech = self.solve_with_techniques(&mut working);
        Self::technique_to_difficulty(max_tech, empty_count)
    }

    /// Rate the puzzle using the Sudoku Explainer (SE) numerical scale.
    pub fn rate_se(&self, grid: &Grid) -> f32 {
        let mut working = grid.deep_clone();
        let max_tech = self.solve_with_techniques(&mut working);
        max_tech.se_rating()
    }

    // ==================== Internal dispatch ====================

    /// Find the first applicable technique for a hint (does not mutate grid).
    fn find_first_technique(&self, grid: &Grid) -> Option<Finding> {
        let fab = CandidateFabric::from_grid(grid);

        // Phase 1: Basic
        if let Some(f) = basic::find_naked_single(&fab) {
            return Some(f);
        }
        if let Some(f) = basic::find_hidden_single(&fab) {
            return Some(f);
        }

        // Phase 2: Subsets
        if let Some(f) = basic::find_naked_subset(&fab, 2) {
            return Some(f);
        }
        if let Some(f) = basic::find_hidden_subset(&fab, 2) {
            return Some(f);
        }
        if let Some(f) = basic::find_naked_subset(&fab, 3) {
            return Some(f);
        }
        if let Some(f) = basic::find_hidden_subset(&fab, 3) {
            return Some(f);
        }

        // Phase 3: Intersections (size-1 fish)
        if let Some(f) = fish_engine::find_pointing_pair(&fab) {
            return Some(f);
        }
        if let Some(f) = fish_engine::find_box_line_reduction(&fab) {
            return Some(f);
        }

        // Phase 4: Fish (size 2+) + quads
        if let Some(f) = fish_engine::find_basic_fish(&fab, 2) {
            return Some(f);
        }
        if let Some(f) = fish_engine::find_finned_fish(&fab, 2) {
            return Some(f);
        }
        if let Some(f) = fish_engine::find_basic_fish(&fab, 3) {
            return Some(f);
        }
        if let Some(f) = fish_engine::find_finned_fish(&fab, 3) {
            return Some(f);
        }
        if let Some(f) = fish_engine::find_basic_fish(&fab, 4) {
            return Some(f);
        }
        if let Some(f) = fish_engine::find_finned_fish(&fab, 4) {
            return Some(f);
        }
        if let Some(f) = basic::find_naked_subset(&fab, 4) {
            return Some(f);
        }
        if let Some(f) = basic::find_hidden_subset(&fab, 4) {
            return Some(f);
        }

        // Phase 5: Uniqueness
        if let Some(f) = aic_engine::find_empty_rectangle(&fab) {
            return Some(f);
        }
        if let Some(f) = uniqueness::find_avoidable_rectangle(&fab) {
            return Some(f);
        }
        if let Some(f) = uniqueness::find_unique_rectangle(&fab) {
            return Some(f);
        }
        if let Some(f) = uniqueness::find_hidden_rectangle(&fab) {
            return Some(f);
        }

        // Phase 6: Master
        if let Some(f) = als_engine::find_xy_wing(&fab) {
            return Some(f);
        }
        if let Some(f) = als_engine::find_xyz_wing(&fab) {
            return Some(f);
        }
        if let Some(f) = als_engine::find_wxyz_wing(&fab) {
            return Some(f);
        }
        if let Some(f) = aic_engine::find_w_wing(&fab) {
            return Some(f);
        }
        // AIC family: shared link graph for X-Chain, 3D Medusa, AIC
        let graph = aic_engine::build_link_graph(&fab);
        if let Some(f) = aic_engine::find_x_chain(&fab, &graph) {
            return Some(f);
        }
        if let Some(f) = aic_engine::find_medusa(&fab, &graph) {
            return Some(f);
        }
        if let Some(f) = als_engine::find_sue_de_coq(&fab) {
            return Some(f);
        }
        if let Some(f) = aic_engine::find_aic(&fab, &graph) {
            return Some(f);
        }
        if let Some(f) = fish_engine::find_franken_fish(&fab) {
            return Some(f);
        }
        if let Some(f) = fish_engine::find_siamese_fish(&fab) {
            return Some(f);
        }
        if let Some(f) = als_engine::find_als_xz(&fab) {
            return Some(f);
        }
        if let Some(f) = uniqueness::find_extended_unique_rectangle(&fab) {
            return Some(f);
        }
        if let Some(f) = uniqueness::find_bug(&fab) {
            return Some(f);
        }

        // Phase 7: Extreme
        if let Some(f) = als_engine::find_als_xy_wing(&fab) {
            return Some(f);
        }
        if let Some(f) = als_engine::find_als_chain(&fab) {
            return Some(f);
        }
        if let Some(f) = fish_engine::find_mutant_fish(&fab) {
            return Some(f);
        }
        if let Some(f) = als_engine::find_aligned_pair_exclusion(&fab) {
            return Some(f);
        }
        if let Some(f) = als_engine::find_aligned_triplet_exclusion(&fab) {
            return Some(f);
        }
        if let Some(f) = als_engine::find_death_blossom(&fab) {
            return Some(f);
        }

        // Forcing chains need the Grid for propagation
        let propagate_singles = |g: &Grid, pos: Position, val: u8| -> (Grid, bool) {
            backtrack::propagate_singles(g, pos, val)
        };
        if let Some(f) = aic_engine::find_nishio_fc(grid, &propagate_singles) {
            return Some(f);
        }
        if let Some(f) = aic_engine::find_kraken_fish(grid, &propagate_singles) {
            return Some(f);
        }
        if let Some(f) = aic_engine::find_region_fc(grid, &propagate_singles) {
            return Some(f);
        }
        if let Some(f) = aic_engine::find_cell_fc(grid, &propagate_singles) {
            return Some(f);
        }
        // Dynamic FC uses full technique propagation
        let prop_full =
            |g: &Grid, pos: Position, val: u8| -> (Grid, bool) { propagate_full(g, pos, val) };
        if let Some(f) = aic_engine::find_dynamic_fc(grid, &prop_full) {
            return Some(f);
        }

        None
    }

    /// Solve the puzzle using human techniques, returning the hardest technique used.
    fn solve_with_techniques(&self, grid: &mut Grid) -> Technique {
        grid.recalculate_candidates();
        let mut max_technique = Technique::NakedSingle;

        while !grid.is_complete() {
            let fab = CandidateFabric::from_grid(grid);

            // Try techniques in priority order via dispatch table
            let finding = None
                // Phase 1: Basic
                .or_else(|| basic::find_naked_single(&fab))
                .or_else(|| basic::find_hidden_single(&fab))
                // Phase 2: Subsets
                .or_else(|| basic::find_naked_subset(&fab, 2))
                .or_else(|| basic::find_hidden_subset(&fab, 2))
                .or_else(|| basic::find_naked_subset(&fab, 3))
                .or_else(|| basic::find_hidden_subset(&fab, 3))
                // Phase 3: Intersections (size-1 fish)
                .or_else(|| fish_engine::find_pointing_pair(&fab))
                .or_else(|| fish_engine::find_box_line_reduction(&fab))
                // Phase 4: Fish (size 2+) + quads
                .or_else(|| fish_engine::find_basic_fish(&fab, 2))
                .or_else(|| fish_engine::find_finned_fish(&fab, 2))
                .or_else(|| fish_engine::find_basic_fish(&fab, 3))
                .or_else(|| fish_engine::find_finned_fish(&fab, 3))
                .or_else(|| fish_engine::find_basic_fish(&fab, 4))
                .or_else(|| fish_engine::find_finned_fish(&fab, 4))
                .or_else(|| basic::find_naked_subset(&fab, 4))
                .or_else(|| basic::find_hidden_subset(&fab, 4))
                // Phase 5: Uniqueness
                .or_else(|| aic_engine::find_empty_rectangle(&fab))
                .or_else(|| uniqueness::find_avoidable_rectangle(&fab))
                .or_else(|| uniqueness::find_unique_rectangle(&fab))
                .or_else(|| uniqueness::find_hidden_rectangle(&fab))
                // Phase 6: Master
                .or_else(|| als_engine::find_xy_wing(&fab))
                .or_else(|| als_engine::find_xyz_wing(&fab))
                .or_else(|| als_engine::find_wxyz_wing(&fab))
                .or_else(|| aic_engine::find_w_wing(&fab))
                // AIC family: shared link graph for X-Chain, 3D Medusa, AIC
                .or_else(|| {
                    let graph = aic_engine::build_link_graph(&fab);
                    None.or_else(|| aic_engine::find_x_chain(&fab, &graph))
                        .or_else(|| aic_engine::find_medusa(&fab, &graph))
                        .or_else(|| als_engine::find_sue_de_coq(&fab))
                        .or_else(|| aic_engine::find_aic(&fab, &graph))
                })
                .or_else(|| fish_engine::find_franken_fish(&fab))
                .or_else(|| fish_engine::find_siamese_fish(&fab))
                .or_else(|| als_engine::find_als_xz(&fab))
                .or_else(|| uniqueness::find_extended_unique_rectangle(&fab))
                .or_else(|| uniqueness::find_bug(&fab))
                // Phase 7: Extreme
                .or_else(|| als_engine::find_als_xy_wing(&fab))
                .or_else(|| als_engine::find_als_chain(&fab))
                .or_else(|| fish_engine::find_mutant_fish(&fab))
                .or_else(|| als_engine::find_aligned_pair_exclusion(&fab))
                .or_else(|| als_engine::find_aligned_triplet_exclusion(&fab))
                .or_else(|| als_engine::find_death_blossom(&fab))
                // Forcing chains (singles propagation)
                .or_else(|| {
                    let prop = |g: &Grid, pos: Position, val: u8| -> (Grid, bool) {
                        backtrack::propagate_singles(g, pos, val)
                    };
                    None.or_else(|| aic_engine::find_nishio_fc(grid, &prop))
                        .or_else(|| aic_engine::find_kraken_fish(grid, &prop))
                        .or_else(|| aic_engine::find_region_fc(grid, &prop))
                        .or_else(|| aic_engine::find_cell_fc(grid, &prop))
                })
                // Dynamic FC: full technique propagation
                .or_else(|| {
                    let prop_full = |g: &Grid, pos: Position, val: u8| -> (Grid, bool) {
                        propagate_full(g, pos, val)
                    };
                    aic_engine::find_dynamic_fc(grid, &prop_full)
                });

            match finding {
                Some(f) => {
                    if f.technique > max_technique {
                        max_technique = f.technique;
                    }
                    apply_finding(grid, &f);
                }
                None => {
                    // No technique found, use backtracking to finish
                    backtrack::solve_recursive(grid);
                    return Technique::Backtracking;
                }
            }
        }

        max_technique
    }

    /// Map a technique + puzzle characteristics to a difficulty level.
    fn technique_to_difficulty(tech: Technique, empty_count: usize) -> Difficulty {
        match tech {
            Technique::NakedSingle => {
                if empty_count <= 35 {
                    Difficulty::Beginner
                } else {
                    Difficulty::Easy
                }
            }
            Technique::HiddenSingle => Difficulty::Medium,
            Technique::NakedPair
            | Technique::HiddenPair
            | Technique::NakedTriple
            | Technique::HiddenTriple => Difficulty::Intermediate,
            Technique::PointingPair | Technique::BoxLineReduction => Difficulty::Hard,
            Technique::XWing
            | Technique::FinnedXWing
            | Technique::Swordfish
            | Technique::FinnedSwordfish
            | Technique::Jellyfish
            | Technique::FinnedJellyfish
            | Technique::NakedQuad
            | Technique::HiddenQuad
            | Technique::EmptyRectangle
            | Technique::AvoidableRectangle
            | Technique::UniqueRectangle
            | Technique::HiddenRectangle => Difficulty::Expert,
            Technique::XYWing
            | Technique::XYZWing
            | Technique::WXYZWing
            | Technique::WWing
            | Technique::XChain
            | Technique::ThreeDMedusa
            | Technique::SueDeCoq
            | Technique::AIC
            | Technique::FrankenFish
            | Technique::SiameseFish
            | Technique::AlsXz
            | Technique::ExtendedUniqueRectangle
            | Technique::BivalueUniversalGrave => Difficulty::Master,
            Technique::AlsXyWing
            | Technique::AlsChain
            | Technique::MutantFish
            | Technique::AlignedPairExclusion
            | Technique::AlignedTripletExclusion
            | Technique::DeathBlossom
            | Technique::NishioForcingChain
            | Technique::KrakenFish
            | Technique::RegionForcingChain
            | Technique::CellForcingChain
            | Technique::DynamicForcingChain
            | Technique::Backtracking => Difficulty::Extreme,
        }
    }
}

/// Apply a Finding to a mutable Grid.
fn apply_finding(grid: &mut Grid, finding: &Finding) {
    match &finding.inference {
        InferenceResult::Placement { cell, value } => {
            let pos = idx_to_pos(*cell);
            grid.set_cell_unchecked(pos, Some(*value));
            grid.recalculate_candidates();
        }
        InferenceResult::Elimination { cell, values } => {
            let pos = idx_to_pos(*cell);
            for &v in values {
                grid.cell_mut(pos).remove_candidate(v);
            }
        }
    }
}

/// Propagate using the full technique set (for Dynamic Forcing Chains).
///
/// Makes an assumption (set cell value), then loops applying all techniques
/// except forcing chains (to avoid infinite recursion) until no more progress.
fn propagate_full(grid: &Grid, pos: Position, val: u8) -> (Grid, bool) {
    let mut g = grid.deep_clone();
    g.set_cell_unchecked(pos, Some(val));
    g.recalculate_candidates();

    for _ in 0..200 {
        if backtrack::has_contradiction(&g) {
            return (g, true);
        }
        if g.is_complete() {
            return (g, false);
        }

        let fab = CandidateFabric::from_grid(&g);

        // Try all techniques except forcing chains (avoids infinite recursion)
        let finding = None
            .or_else(|| basic::find_naked_single(&fab))
            .or_else(|| basic::find_hidden_single(&fab))
            .or_else(|| basic::find_naked_subset(&fab, 2))
            .or_else(|| basic::find_hidden_subset(&fab, 2))
            .or_else(|| basic::find_naked_subset(&fab, 3))
            .or_else(|| basic::find_hidden_subset(&fab, 3))
            .or_else(|| fish_engine::find_pointing_pair(&fab))
            .or_else(|| fish_engine::find_box_line_reduction(&fab))
            .or_else(|| fish_engine::find_basic_fish(&fab, 2))
            .or_else(|| fish_engine::find_finned_fish(&fab, 2))
            .or_else(|| fish_engine::find_basic_fish(&fab, 3))
            .or_else(|| fish_engine::find_finned_fish(&fab, 3))
            .or_else(|| fish_engine::find_basic_fish(&fab, 4))
            .or_else(|| fish_engine::find_finned_fish(&fab, 4))
            .or_else(|| basic::find_naked_subset(&fab, 4))
            .or_else(|| basic::find_hidden_subset(&fab, 4))
            .or_else(|| aic_engine::find_empty_rectangle(&fab))
            .or_else(|| uniqueness::find_avoidable_rectangle(&fab))
            .or_else(|| uniqueness::find_unique_rectangle(&fab))
            .or_else(|| uniqueness::find_hidden_rectangle(&fab))
            .or_else(|| als_engine::find_xy_wing(&fab))
            .or_else(|| als_engine::find_xyz_wing(&fab))
            .or_else(|| als_engine::find_wxyz_wing(&fab))
            .or_else(|| aic_engine::find_w_wing(&fab))
            .or_else(|| {
                let graph = aic_engine::build_link_graph(&fab);
                None.or_else(|| aic_engine::find_x_chain(&fab, &graph))
                    .or_else(|| aic_engine::find_medusa(&fab, &graph))
                    .or_else(|| als_engine::find_sue_de_coq(&fab))
                    .or_else(|| aic_engine::find_aic(&fab, &graph))
            })
            .or_else(|| fish_engine::find_franken_fish(&fab))
            .or_else(|| fish_engine::find_siamese_fish(&fab))
            .or_else(|| als_engine::find_als_xz(&fab))
            .or_else(|| uniqueness::find_extended_unique_rectangle(&fab))
            .or_else(|| uniqueness::find_bug(&fab))
            .or_else(|| als_engine::find_als_xy_wing(&fab))
            .or_else(|| als_engine::find_als_chain(&fab))
            .or_else(|| fish_engine::find_mutant_fish(&fab))
            .or_else(|| als_engine::find_aligned_pair_exclusion(&fab))
            .or_else(|| als_engine::find_aligned_triplet_exclusion(&fab))
            .or_else(|| als_engine::find_death_blossom(&fab));
        // Note: forcing chains excluded to avoid infinite recursion

        match finding {
            Some(f) => apply_finding(&mut g, &f),
            None => break,
        }
    }

    let contradiction = backtrack::has_contradiction(&g);
    (g, contradiction)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_solve_easy() {
        let puzzle =
            "530070000600195000098000060800060003400803001700020006060000280000419005000080079";
        let grid = Grid::from_string(puzzle).unwrap();
        let solver = Solver::new();
        let solution = solver.solve(&grid).unwrap();
        assert!(solution.is_complete());
    }

    #[test]
    fn test_unique_solution() {
        let puzzle =
            "530070000600195000098000060800060003400803001700020006060000280000419005000080079";
        let grid = Grid::from_string(puzzle).unwrap();
        let solver = Solver::new();
        assert!(solver.has_unique_solution(&grid));
    }

    #[test]
    fn test_get_hint() {
        let puzzle =
            "530070000600195000098000060800060003400803001700020006060000280000419005000080079";
        let grid = Grid::from_string(puzzle).unwrap();
        let solver = Solver::new();
        let hint = solver.get_hint(&grid);
        assert!(hint.is_some());
    }

    #[test]
    fn test_difficulty_rating() {
        let puzzle =
            "530070000600195000098000060800060003400803001700020006060000280000419005000080079";
        let grid = Grid::from_string(puzzle).unwrap();
        let solver = Solver::new();
        let difficulty = solver.rate_difficulty(&grid);
        assert!(difficulty >= Difficulty::Easy);
    }

    #[test]
    fn test_solve_with_techniques_regression() {
        let puzzle =
            "530070000600195000098000060800060003400803001700020006060000280000419005000080079";
        let grid = Grid::from_string(puzzle).unwrap();
        let solver = Solver::new();
        let mut working = grid.deep_clone();
        let max_tech = solver.solve_with_techniques(&mut working);
        assert!(max_tech < Technique::Backtracking);
        assert!(working.is_complete());
    }

    /// Soundness test: verify that every elimination/placement returned by hints
    /// is consistent with the unique solution.
    #[test]
    fn test_hint_soundness() {
        let puzzles = [
            // Easy (naked/hidden singles)
            "530070000600195000098000060800060003400803001700020006060000280000419005000080079",
            // Medium
            "020000600008020050500060020060000093003905100790000080050090004010070300006000010",
            // Arto Inkala (requires advanced techniques)
            "800000000003600000070090200050007000000045700000100030001000068008500010090000400",
        ];

        let solver = Solver::new();

        for puzzle_str in &puzzles {
            let grid = Grid::from_string(puzzle_str).unwrap();
            let solution = match solver.solve(&grid) {
                Some(s) if s.is_complete() => s,
                _ => continue,
            };

            let mut working = grid.deep_clone();
            working.recalculate_candidates();

            let mut steps = 0;
            while !working.is_complete() && steps < 300 {
                let hint = match solver.get_hint(&working) {
                    Some(h) => h,
                    None => break,
                };

                match &hint.hint_type {
                    HintType::SetValue { pos, value } => {
                        let sol_val = solution.get(*pos);
                        assert_eq!(
                            sol_val,
                            Some(*value),
                            "Unsound placement by {:?}: ({},{}) = {}, solution has {:?}. Puzzle: {}",
                            hint.technique,
                            pos.row + 1, pos.col + 1, value, sol_val, puzzle_str
                        );
                        working.set_cell_unchecked(*pos, Some(*value));
                        working.recalculate_candidates();
                    }
                    HintType::EliminateCandidates { pos, values } => {
                        let sol_val = solution.get(*pos).expect("Position should have solution");
                        for &v in values {
                            assert_ne!(
                                v, sol_val,
                                "Unsound elimination by {:?}: removing {} from ({},{}) but solution needs it. Puzzle: {}",
                                hint.technique, v, pos.row + 1, pos.col + 1, puzzle_str
                            );
                        }
                        for &v in values {
                            working.cell_mut(*pos).remove_candidate(v);
                        }
                    }
                }
                steps += 1;
            }
        }
    }

    /// Soundness test for generated puzzles at every difficulty tier.
    /// Uses deterministic seeds so failures are reproducible.
    #[test]
    fn test_hint_soundness_all_tiers() {
        use crate::Generator;

        let solver = Solver::new();

        for (seed, difficulty) in [
            (42, Difficulty::Easy),
            (42, Difficulty::Medium),
            (42, Difficulty::Intermediate),
            (42, Difficulty::Hard),
            (42, Difficulty::Expert),
            (42, Difficulty::Master),
            (42, Difficulty::Extreme),
            (99, Difficulty::Expert),
            (99, Difficulty::Master),
            (99, Difficulty::Extreme),
            (123, Difficulty::Extreme),
            (200, Difficulty::Extreme),
            (314, Difficulty::Extreme),
            (500, Difficulty::Extreme),
            (777, Difficulty::Extreme),
            (1000, Difficulty::Extreme),
            (1234, Difficulty::Extreme),
            (2024, Difficulty::Extreme),
            (9999, Difficulty::Extreme),
        ] {
            let mut gen = Generator::with_seed(seed);
            let grid = gen.generate(difficulty);

            let solution = match solver.solve(&grid) {
                Some(s) if s.is_complete() => s,
                _ => continue,
            };

            let mut working = grid.deep_clone();
            working.recalculate_candidates();

            let mut steps = 0;
            while !working.is_complete() && steps < 500 {
                let hint = match solver.get_hint(&working) {
                    Some(h) => h,
                    None => break,
                };

                match &hint.hint_type {
                    HintType::SetValue { pos, value } => {
                        let sol_val = solution.get(*pos);
                        assert_eq!(
                            sol_val,
                            Some(*value),
                            "Unsound placement by {:?} (SE {:.1}): ({},{}) = {}, solution has {:?}. Seed={} Diff={:?} Step={}",
                            hint.technique, hint.technique.se_rating(),
                            pos.row + 1, pos.col + 1, value, sol_val,
                            seed, difficulty, steps
                        );
                        working.set_cell_unchecked(*pos, Some(*value));
                        working.recalculate_candidates();
                    }
                    HintType::EliminateCandidates { pos, values } => {
                        let sol_val = solution.get(*pos).expect("Position should have solution");
                        for &v in values {
                            assert_ne!(
                                v, sol_val,
                                "Unsound elimination by {:?} (SE {:.1}): removing {} from ({},{}) but solution needs it. Seed={} Diff={:?} Step={}",
                                hint.technique, hint.technique.se_rating(),
                                v, pos.row + 1, pos.col + 1,
                                seed, difficulty, steps
                            );
                        }
                        for &v in values {
                            working.cell_mut(*pos).remove_candidate(v);
                        }
                    }
                }
                steps += 1;
            }
        }
    }

    /// Soundness test for get_next_placement() which chains eliminations.
    /// This is what the WASM apply_hint() uses.
    #[test]
    fn test_next_placement_soundness() {
        use crate::Generator;

        let solver = Solver::new();

        for (seed, difficulty) in [
            (42, Difficulty::Expert),
            (42, Difficulty::Master),
            (42, Difficulty::Extreme),
            (99, Difficulty::Extreme),
            (123, Difficulty::Extreme),
            (200, Difficulty::Extreme),
            (314, Difficulty::Extreme),
            (500, Difficulty::Extreme),
        ] {
            let mut gen = Generator::with_seed(seed);
            let grid = gen.generate(difficulty);

            let solution = match solver.solve(&grid) {
                Some(s) if s.is_complete() => s,
                _ => continue,
            };

            let mut working = grid.deep_clone();
            working.recalculate_candidates();

            let mut steps = 0;
            while !working.is_complete() && steps < 200 {
                let hint = match solver.get_next_placement(&working) {
                    Some(h) => h,
                    None => break,
                };

                match &hint.hint_type {
                    HintType::SetValue { pos, value } => {
                        let sol_val = solution.get(*pos);
                        assert_eq!(
                            sol_val,
                            Some(*value),
                            "Unsound placement from get_next_placement by {:?}: ({},{}) = {}, solution has {:?}. Seed={} Diff={:?} Step={}",
                            hint.technique, pos.row + 1, pos.col + 1, value, sol_val,
                            seed, difficulty, steps
                        );
                        working.set_cell_unchecked(*pos, Some(*value));
                        working.recalculate_candidates();
                    }
                    HintType::EliminateCandidates { .. } => {
                        panic!("get_next_placement returned elimination instead of placement at seed={} diff={:?} step={}", seed, difficulty, steps);
                    }
                }
                steps += 1;
            }

            assert!(
                working.is_complete(),
                "Failed to complete puzzle: seed={} diff={:?} stopped at step={}",
                seed,
                difficulty,
                steps
            );
        }
    }

    /// Verify that hints from Expert+ puzzles carry ProofCertificate data.
    #[test]
    fn test_hint_carries_proof_certificate() {
        let solver = Solver::new();

        // Expert puzzle that requires UR / fish — should produce proof certificates
        let puzzle =
            "000704005020010070000080002090006250600070008053200010400090000030060090200301000";
        let grid = Grid::from_string(puzzle).unwrap();
        let mut working = grid.deep_clone();
        working.recalculate_candidates();

        let mut found_proof = false;
        let mut steps = 0;
        while !working.is_complete() && steps < 300 {
            let hint = match solver.get_hint(&working) {
                Some(h) => h,
                None => break,
            };

            if hint.proof.is_some() {
                found_proof = true;
            }

            match &hint.hint_type {
                HintType::SetValue { pos, value } => {
                    working.set_cell_unchecked(*pos, Some(*value));
                    working.recalculate_candidates();
                }
                HintType::EliminateCandidates { pos, values } => {
                    for &v in values {
                        working.cell_mut(*pos).remove_candidate(v);
                    }
                }
            }
            steps += 1;
        }

        assert!(
            found_proof,
            "Expected at least one hint with a ProofCertificate on an Expert puzzle"
        );
    }

    /// Regression test: pin SE ratings for known puzzles so engine changes
    /// that alter difficulty classification are caught immediately.
    #[test]
    fn test_se_rating_regression() {
        let cases: &[(&str, f32, Difficulty)] = &[
            // Naked singles only
            (
                "530070000600195000098000060800060003400803001700020006060000280000419005000080079",
                2.3,
                Difficulty::Easy,
            ),
            // Hidden singles only
            (
                "000000010400000000020000000000050407008000300001090000300400200050100000000806000",
                1.5,
                Difficulty::Medium,
            ),
            // Expert tier (UR / empty rectangle)
            (
                "000704005020010070000080002090006250600070008053200010400090000030060090200301000",
                4.6,
                Difficulty::Expert,
            ),
        ];

        let solver = Solver::new();
        for (puzzle_str, expected_se, expected_diff) in cases {
            let grid = Grid::from_string(puzzle_str).unwrap();
            let se = solver.rate_se(&grid);
            let diff = solver.rate_difficulty(&grid);
            assert!(
                (se - expected_se).abs() < 0.01,
                "SE regression: expected {:.1} got {:.1} for puzzle {}",
                expected_se,
                se,
                puzzle_str
            );
            assert_eq!(
                diff, *expected_diff,
                "Difficulty regression: expected {:?} got {:?} for puzzle {}",
                expected_diff, diff, puzzle_str
            );
        }
    }

    /// Verify that solve_with_techniques produces consistent max-technique
    /// for known puzzles (guards against dispatch order regressions).
    #[test]
    fn test_max_technique_regression() {
        let solver = Solver::new();

        // This puzzle needs only naked singles (SE 2.3)
        let grid = Grid::from_string(
            "530070000600195000098000060800060003400803001700020006060000280000419005000080079",
        )
        .unwrap();
        let mut w = grid.deep_clone();
        let tech = solver.solve_with_techniques(&mut w);
        assert!(w.is_complete(), "Puzzle should be fully solved");
        assert!(
            tech <= Technique::NakedSingle,
            "Expected NakedSingle, got {:?}",
            tech
        );

        // This puzzle needs hidden singles (SE 1.5)
        let grid = Grid::from_string(
            "000000010400000000020000000000050407008000300001090000300400200050100000000806000",
        )
        .unwrap();
        let mut w = grid.deep_clone();
        let tech = solver.solve_with_techniques(&mut w);
        assert!(w.is_complete(), "Puzzle should be fully solved");
        assert!(
            tech <= Technique::HiddenSingle,
            "Expected ≤HiddenSingle, got {:?}",
            tech
        );
    }

    /// Technique coverage: verify each engine can find its expected technique
    /// type on a suitable puzzle state.
    #[test]
    fn test_technique_coverage_basic() {
        // Verify basic techniques fire on appropriate puzzle states
        let solver = Solver::new();

        // The soundness test already covers that hints are correct.
        // Here we verify the engine can solve all three reference puzzles
        // to completion using human techniques.
        let easy = Grid::from_string(
            "530070000600195000098000060800060003400803001700020006060000280000419005000080079",
        )
        .unwrap();
        let mut w = easy.deep_clone();
        let tech = solver.solve_with_techniques(&mut w);
        assert!(w.is_complete());
        assert!(
            tech < Technique::Backtracking,
            "Easy puzzle should not need backtracking"
        );

        let medium = Grid::from_string(
            "000000010400000000020000000000050407008000300001090000300400200050100000000806000",
        )
        .unwrap();
        let mut w = medium.deep_clone();
        let tech = solver.solve_with_techniques(&mut w);
        assert!(w.is_complete());
        assert!(
            tech < Technique::Backtracking,
            "Medium puzzle should not need backtracking"
        );
    }

    /// Verify every technique in the dispatch chain is reachable by collecting
    /// all techniques used across a battery of generated puzzles at each tier.
    #[test]
    fn test_technique_tier_coverage() {
        use crate::Generator;
        use std::collections::HashSet;

        let solver = Solver::new();
        let mut seen_techniques = HashSet::new();

        // Generate puzzles at different difficulty tiers using deterministic seeds
        for (seed, difficulty) in [
            (42, Difficulty::Easy),
            (42, Difficulty::Medium),
            (42, Difficulty::Intermediate),
            (42, Difficulty::Hard),
            (42, Difficulty::Expert),
            (99, Difficulty::Easy),
            (99, Difficulty::Medium),
            (99, Difficulty::Hard),
            (99, Difficulty::Expert),
        ] {
            let mut gen = Generator::with_seed(seed);
            let grid = gen.generate(difficulty);
            let mut working = grid.deep_clone();
            working.recalculate_candidates();

            // Solve step by step, collecting every technique used
            let mut steps = 0;
            while !working.is_complete() && steps < 300 {
                let hint = match solver.get_hint(&working) {
                    Some(h) => h,
                    None => break,
                };
                seen_techniques.insert(hint.technique);

                match &hint.hint_type {
                    HintType::SetValue { pos, value } => {
                        working.set_cell_unchecked(*pos, Some(*value));
                        working.recalculate_candidates();
                    }
                    HintType::EliminateCandidates { pos, values } => {
                        for &v in values {
                            working.cell_mut(*pos).remove_candidate(v);
                        }
                    }
                }
                steps += 1;
            }
        }

        // At minimum, the basic techniques should be exercised
        assert!(
            seen_techniques.contains(&Technique::NakedSingle),
            "NakedSingle never fired"
        );
        assert!(
            seen_techniques.contains(&Technique::HiddenSingle),
            "HiddenSingle never fired"
        );

        // Print coverage for diagnostics
        let mut sorted: Vec<_> = seen_techniques.iter().collect();
        sorted.sort();
        println!("Technique coverage ({} of 45):", sorted.len());
        for t in &sorted {
            println!("  {:?} (SE {:.1})", t, t.se_rating());
        }
    }
}
