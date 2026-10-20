//! Direct technique selection for tutorials and audits. The normal ladder is unchanged.
use super::*;

impl Solver {
    /// Request one named detector without first applying easier techniques.
    ///
    /// Retains supplied candidate masks. An absent result includes detector limits.
    /// Metadata returned by older detectors is not a checked proof. Uniqueness
    /// methods retain their unique-solution premise. Use independent validation
    /// before treating an arbitrary detector endpoint as established.
    #[allow(deprecated)]
    pub fn get_hint_for_technique(
        &self,
        grid: &Grid,
        technique: Technique,
        options: &ForcingSearchOptions,
    ) -> Result<Option<Hint>, String> {
        solver_state::validate(grid)?;
        options.validate()?;
        let f = CandidateFabric::from_grid(grid);
        use Technique::*;
        let found = match technique {
            NakedSingle => basic::find_naked_single(&f),
            HiddenSingle => basic::find_hidden_single(&f),
            NakedPair => basic::find_naked_subset(&f, 2),
            NakedTriple => basic::find_naked_subset(&f, 3),
            NakedQuad => basic::find_naked_subset(&f, 4),
            HiddenPair => basic::find_hidden_subset(&f, 2),
            HiddenTriple => basic::find_hidden_subset(&f, 3),
            HiddenQuad => basic::find_hidden_subset(&f, 4),
            PointingPair => fish_engine::find_pointing_pair(&f),
            BoxLineReduction => fish_engine::find_box_line_reduction(&f),
            XWing => fish_engine::find_basic_fish(&f, 2),
            Swordfish => fish_engine::find_basic_fish(&f, 3),
            Jellyfish => fish_engine::find_basic_fish(&f, 4),
            FinnedXWing => fish_engine::find_finned_fish(&f, 2),
            FinnedSwordfish => fish_engine::find_finned_fish(&f, 3),
            FinnedJellyfish => fish_engine::find_finned_fish(&f, 4),
            FrankenFish => fish_engine::find_franken_fish(&f),
            SiameseFish => fish_engine::find_siamese_fish(&f),
            MutantFish => fish_engine::find_mutant_fish(&f),
            EmptyRectangle => aic_engine::find_empty_rectangle(&f),
            WWing => aic_engine::find_w_wing(&f),
            XChain => aic_engine::find_x_chain(&f, &aic_engine::build_link_graph(&f)),
            AIC => aic_engine::find_aic(&f, &aic_engine::build_link_graph(&f)),
            ThreeDMedusa => aic_engine::find_medusa(&f, &aic_engine::build_link_graph(&f)),
            AvoidableRectangle => uniqueness::find_avoidable_rectangle(&f),
            UniqueRectangle => uniqueness::find_unique_rectangle(&f),
            HiddenRectangle => uniqueness::find_hidden_rectangle(&f),
            ExtendedUniqueRectangle => uniqueness::find_extended_unique_rectangle(&f),
            BivalueUniversalGrave => uniqueness::find_bug(&f),
            XYWing => als_engine::find_xy_wing(&f),
            XYZWing => als_engine::find_xyz_wing(&f),
            WXYZWing => als_engine::find_wxyz_wing(&f),
            AlsXz => als_engine::find_als_xz(&f),
            AlsXyWing => als_engine::find_als_xy_wing(&f),
            AlsChain => als_engine::find_als_chain(&f),
            SueDeCoq => als_engine::find_sue_de_coq(&f),
            DeathBlossom => als_engine::find_death_blossom(&f),
            AlignedPairExclusion => als_engine::find_aligned_pair_exclusion(&f),
            AlignedTripletExclusion => als_engine::find_aligned_triplet_exclusion(&f),
            ThreeHouseMatching => three_house::find(grid),
            ArithmeticCounting => arithmetic::find(grid, &f),
            GlobalColoring => global_coloring::find(grid),
            Backtracking => backtrack::find_backtracking_hint(grid),
            NishioForcingChain | KrakenFish | RegionForcingChain | CellForcingChain
            | DynamicForcingChain => {
                let budget = std::cell::RefCell::new(forcing::SearchBudget::new(options));
                let policy = if technique == DynamicForcingChain {
                    PropagationPolicy::Logical
                } else {
                    PropagationPolicy::Singles
                };
                let propagate = |g: &Grid, p: Position, d: u8| {
                    forcing::propagate_budgeted(g, p, d, policy, &mut budget.borrow_mut())
                };
                let detector = match technique {
                    NishioForcingChain => aic_engine::find_nishio_fc,
                    KrakenFish => aic_engine::find_kraken_fish,
                    RegionForcingChain => aic_engine::find_region_fc,
                    CellForcingChain => aic_engine::find_cell_fc,
                    _ => aic_engine::find_dynamic_fc,
                };
                detector(grid, &propagate)
            }
        };
        // Some general fish searches also return a finned classification.
        Ok(found
            .filter(|f| f.technique == technique)
            .map(|f| f.to_hint()))
    }
}
