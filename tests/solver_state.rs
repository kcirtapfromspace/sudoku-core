use sudoku_core::{
    BitSet, CandidateState, ForcingSearchOptions, GlobalColoringSearchOptions, Grid,
    GuardianLiteral, GuardianProof, HintType, Position, ProofCertificate, PropagationPolicy,
    Solver, SolverState, Technique,
};
fn pos(c: usize) -> Position {
    Position::new(c / 9, c % 9)
}
fn set_domain(grid: &mut Grid, c: usize, mask: u16) {
    grid.cell_mut(pos(c))
        .set_candidates(BitSet::from_raw(mask << 1));
}

#[test]
fn player_notes_and_retained_deductions_are_explicitly_separate() {
    let mut player = Grid::new_classic();
    set_domain(&mut player, 0, 1);
    let fresh = SolverState::new(&player).unwrap();
    assert_eq!(fresh.grid().get_candidates(pos(0)), BitSet::all_9());
    let retained = SolverState::from_candidates(&player).unwrap();
    let solver = Solver::new();
    let hint = retained.get_hint(&solver).unwrap().unwrap();
    assert!(matches!(hint.hint_type, HintType::SetValue { pos: p, value: 1 } if p == pos(0)));
    assert_eq!(player.get(pos(0)), None);
    assert_eq!(player.get_candidates(pos(0)), BitSet::single(1));
    // The ordinary API rebuilds from values, so it cannot call this a single.
    assert_ne!(
        solver.get_hint(&player).unwrap().technique,
        Technique::NakedSingle
    );
    set_domain(&mut player, 0, 0);
    assert!(solver.get_hint_with_candidates(&player).is_err());
    assert!(SolverState::from_candidates(&player).is_err());
    assert!(SolverState::new(&player).is_ok());
}

#[test]
fn retained_full_ladder_exposes_the_previous_private_stall_regression() {
    let fixture: serde_json::Value = serde_json::from_str(include_str!(
        "../research/global_obstructions/different_lead.json"
    ))
    .unwrap();
    let state: CandidateState = serde_json::from_value(fixture["state"].clone()).unwrap();
    let mut grid = state.to_grid().unwrap();
    // CandidateState is a native proof projection; restore game given flags for
    // the full ladder's existing uniqueness techniques.
    for (c, byte) in fixture["puzzle"].as_str().unwrap().bytes().enumerate() {
        grid.cell_mut(pos(c))
            .set_given((b'1'..=b'9').contains(&byte));
    }
    let hint = Solver::new()
        .get_hint_with_candidates(&grid)
        .unwrap()
        .unwrap();
    assert_eq!(hint.technique, Technique::GlobalColoring);
    assert!(matches!(hint.hint_type, HintType::SetValue { pos: p, value: 7 } if p == pos(68)));
    assert!(matches!(hint.proof, Some(ProofCertificate::GlobalColoring(p)) if p.verify(&grid)));
    assert_eq!(CandidateState::capture(&grid).unwrap(), state);
}

#[test]
fn solver_placements_preserve_nonpeer_deletions_and_edits_rebuild_them() {
    let mut grid = Grid::new_classic();
    set_domain(&mut grid, 0, 1);
    grid.cell_mut(pos(80)).remove_candidate(9);
    let mut state = SolverState::from_candidates(&grid).unwrap();
    state.step(&Solver::new()).unwrap().unwrap();
    assert_eq!(state.grid().get(pos(0)), Some(1));
    assert!(!state.grid().get_candidates(pos(80)).contains(9));
    assert!(!state.grid().get_candidates(pos(1)).contains(1));
    state.edit(pos(0), Some(2)).unwrap();
    assert!(state.grid().get_candidates(pos(80)).contains(9));
    assert!(state.grid().get_candidates(pos(1)).contains(1));
    assert!(!state.grid().get_candidates(pos(1)).contains(2));
    state.edit(pos(0), None).unwrap();
    assert_eq!(state.grid().get_candidates(pos(0)), BitSet::all_9());
}

fn guardian_state() -> (Grid, GuardianProof) {
    let mut grid = Grid::new_classic();
    set_domain(&mut grid, 0, 3);
    set_domain(&mut grid, 1, 7);
    set_domain(&mut grid, 2, 7);
    // An independent later single eliminates just the first guardian.
    set_domain(&mut grid, 28, 4);
    let proof = GuardianProof::from_cells(
        &grid,
        vec![0, 1, 2],
        vec![
            GuardianLiteral { cell: 1, digit: 3 },
            GuardianLiteral { cell: 2, digit: 3 },
        ],
        1000,
    )
    .unwrap();
    (grid, proof)
}

#[test]
fn guardian_or_survives_deductions_and_propagates_common_peers_then_a_unit() {
    let (grid, proof) = guardian_state();
    assert!(proof.verify());
    let decoded: GuardianProof =
        serde_json::from_str(&serde_json::to_string(&proof).unwrap()).unwrap();
    assert_eq!(decoded, proof);
    let mut state = SolverState::from_candidates(&grid).unwrap();
    state.retain_guardian(proof.clone()).unwrap();
    state.retain_guardian(proof).unwrap();
    assert_eq!(state.guardians().len(), 1);
    let solver = Solver::new();
    let mut placed = false;
    for _ in 0..81 {
        let before = state.grid().deep_clone();
        let hint = state.step(&solver).unwrap().unwrap();
        if let Some(ProofCertificate::Guardian(p)) = &hint.proof {
            assert!(p.verify(&before));
            assert!(!p.verify(state.grid()));
            let mut forged = p.clone();
            forged.digit = 0;
            assert!(!forged.verify(&before));
        }
        assert_eq!(state.guardians().len(), 1);
        if let HintType::SetValue { pos: p, value: 3 } = hint.hint_type {
            if p == pos(2) {
                assert_eq!(hint.technique, Technique::GlobalColoring);
                placed = true;
                break;
            }
        }
    }
    assert!(placed, "remaining guardian must become a placement");
    assert_eq!(state.grid().get(pos(28)), Some(3));
    assert!(state.undo());
    assert!(state.guardians().is_empty());
}

#[test]
fn guardian_proofs_and_lifecycle_reject_forged_and_stale_relations() {
    let (grid, proof) = guardian_state();
    let mut altered = proof.clone();
    altered.guardians.pop();
    assert!(!altered.verify());
    altered = proof.clone();
    altered.guardians.push(altered.guardians[0]);
    assert!(!altered.verify());
    altered = proof.clone();
    altered.guardians[0].cell = 81;
    assert!(!altered.verify());
    altered = proof.clone();
    altered.guardians[0].digit = 0;
    assert!(!altered.verify());
    altered = proof.clone();
    altered.nodes.clear();
    assert!(!altered.verify());
    altered = proof.clone();
    altered.witness[0] = 0;
    assert!(!altered.verify());
    altered = proof.clone();
    altered.premise.domains.pop();
    assert!(!altered.verify());
    assert!(GuardianProof::from_cells(&grid, vec![0, 1, 2], proof.guardians.clone(), 0).is_none());
    let mut state = SolverState::from_candidates(&grid).unwrap();
    assert!(state.retain_guardian(altered).is_err());
    assert!(state.guardians().is_empty());
    state.retain_guardian(proof.clone()).unwrap();
    let before = CandidateState::capture(state.grid()).unwrap();
    assert!(state.edit(pos(81), Some(3)).is_err());
    assert_eq!(CandidateState::capture(state.grid()).unwrap(), before);
    assert_eq!(state.guardians().len(), 1);
    state.edit(pos(80), Some(9)).unwrap();
    assert!(state.guardians().is_empty());
    assert!(state.retain_guardian(proof).is_err());
    assert!(state.undo());
    assert_eq!(CandidateState::capture(state.grid()).unwrap(), before);
    assert!(state.guardians().is_empty());
    assert!(!state.undo());
    assert!(state
        .learn_guardian(&GlobalColoringSearchOptions {
            max_patterns: 16385,
            max_nodes: 0
        })
        .is_err());
    assert!(state
        .learn_guardian(&GlobalColoringSearchOptions {
            max_patterns: 0,
            max_nodes: 0
        })
        .unwrap()
        .is_none());
}

#[test]
fn policies_and_invalid_masks_are_checked_at_retained_api_boundary() {
    let mut grid = Grid::new_classic();
    let solver = Solver::new();
    assert!(solver
        .get_hint_with_options(
            &grid,
            PropagationPolicy::Singles,
            &ForcingSearchOptions::default()
        )
        .unwrap()
        .is_none());
    grid.set_given(pos(0), 1);
    grid.cell_mut(pos(1)).add_candidate(1);
    assert!(solver.get_hint_with_candidates(&grid).is_err());
    assert!(solver
        .get_hint_with_candidates(&Grid::new_with_constraints(vec![]))
        .is_err());
    let mut variant = Grid::new_x_sudoku();
    variant.set_given(pos(0), 1);
    variant.set_cell_unchecked(pos(40), Some(1));
    variant.recalculate_candidates();
    assert!(solver.get_hint_with_candidates(&variant).is_err());
}

#[test]
fn automatic_guardian_scan_retains_a_verified_multi_literal_consequence() {
    let mut grid = Grid::new_classic();
    for c in 0..2 {
        set_domain(&mut grid, c, 7);
    }
    for c in 2..4 {
        set_domain(&mut grid, c, 15);
    }
    let mut state = SolverState::from_candidates(&grid).unwrap();
    let proof = state
        .learn_guardian(&GlobalColoringSearchOptions {
            max_patterns: 1,
            max_nodes: 256,
        })
        .unwrap()
        .expect("four row cells cannot all use a three-digit palette");
    assert!(proof.verify());
    assert_eq!(
        proof.guardians,
        vec![
            GuardianLiteral { cell: 2, digit: 4 },
            GuardianLiteral { cell: 3, digit: 4 }
        ]
    );
    assert_eq!(state.guardians().len(), 1);
    let hint = state.get_hint(&Solver::new()).unwrap().unwrap();
    assert!(matches!(hint.proof, Some(ProofCertificate::Guardian(p)) if p.verify(&grid)));
    assert!(state
        .get_hint_with_options(
            &Solver::new(),
            PropagationPolicy::Full,
            &ForcingSearchOptions {
                max_depth: 9,
                ..Default::default()
            }
        )
        .is_err());
}

#[test]
fn state_queries_and_steps_share_policies_and_invalid_step_limits_are_atomic() {
    let solver = Solver::new();
    let mut state = SolverState::new(&Grid::new_classic()).unwrap();
    assert!(state
        .step_with_options(
            &solver,
            PropagationPolicy::Singles,
            &ForcingSearchOptions::default()
        )
        .unwrap()
        .is_none());
    let before = CandidateState::capture(state.grid()).unwrap();
    assert!(state
        .step_with_options(
            &solver,
            PropagationPolicy::Full,
            &ForcingSearchOptions {
                max_depth: 9,
                ..Default::default()
            }
        )
        .is_err());
    assert_eq!(CandidateState::capture(state.grid()).unwrap(), before);
    assert!(!state.undo());
    let (grid, proof) = guardian_state();
    let mut state = SolverState::from_candidates(&grid).unwrap();
    state.retain_guardian(proof).unwrap();
    let query = state.get_hint(&solver).unwrap().unwrap();
    let applied = state.step(&solver).unwrap().unwrap();
    assert_eq!(
        serde_json::to_value(query).unwrap(),
        serde_json::to_value(applied).unwrap()
    );
    let mut impossible = Grid::new_classic();
    set_domain(&mut impossible, 0, 1);
    set_domain(&mut impossible, 1, 1);
    assert!(solver.get_hint_with_candidates(&impossible).is_err());
}
