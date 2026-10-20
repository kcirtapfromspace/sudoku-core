use serde::Deserialize;
use sudoku_core::{
    CandidateState, ForcingLeaf, ForcingProof, ForcingPropagation, ForcingSearchOptions,
    ForcingStep, Grid, Position, Solver, SolverState,
};

#[derive(Deserialize)]
struct Replay {
    state: CandidateState,
    proof: ForcingProof,
}
fn replay() -> Replay {
    serde_json::from_str(include_str!("fixtures/forcing_mith515.json")).unwrap()
}

#[test]
fn seven_node_research_tree_replays_with_checked_hall_and_three_house_leaves() {
    let Replay { state, proof } = replay();
    let mut grid = state.to_grid().unwrap();
    let check = proof
        .check(&grid)
        .expect("independent Hall replay must verify");
    assert_eq!(
        (
            check.nodes,
            check.splits,
            check.hall_leaves,
            check.three_house_leaves,
            check.steps,
            check.depth
        ),
        (7, 3, 1, 3, 84, 2)
    );
    assert_eq!((proof.cell, proof.digit, proof.value), (7, 3, false));
    let decoded: ForcingProof =
        serde_json::from_str(&serde_json::to_string(&proof).unwrap()).unwrap();
    assert_eq!(decoded, proof);
    decoded.apply(&mut grid).unwrap();
    let after = CandidateState::capture(&grid).unwrap();
    assert_eq!(after.domains[7], state.domains[7] & !4);
    for c in 0..81 {
        assert_eq!(after.domains[c] & !state.domains[c], 0);
    }
    assert!(decoded.apply(&mut grid).is_err());
    assert_eq!(CandidateState::capture(&grid).unwrap(), after);
    let mut session = SolverState::from_candidates(&state.to_grid().unwrap()).unwrap();
    session.apply_forcing(&proof).unwrap();
    assert_eq!(CandidateState::capture(session.grid()).unwrap(), after);
    assert!(session.undo());
    assert_eq!(CandidateState::capture(session.grid()).unwrap(), state);
}

#[test]
fn forcing_checker_rejects_incomplete_branches_altered_steps_leaves_and_stale_state() {
    let Replay { state, proof } = replay();
    let grid = state.to_grid().unwrap();
    let mut corruptions = Vec::new();
    macro_rules! alter {
        ($p:ident, $action:expr) => {{
            let mut $p = proof.clone();
            $action;
            corruptions.push($p);
        }};
    }
    alter!(p, p.version = 2);
    alter!(p, p.state_hash = "0".repeat(64));
    alter!(p, p.cell = 81);
    alter!(p, p.digit = 0);
    alter!(p, p.digit = 10);
    alter!(p, p.value = true);
    alter!(p, p.nodes.clear());
    alter!(p, {
        p.nodes.last_mut().unwrap().branches.pop();
    });
    alter!(p, {
        let b = p.nodes.last().unwrap().branches[0].clone();
        p.nodes.last_mut().unwrap().branches.push(b);
    });
    alter!(p, p.nodes.last_mut().unwrap().branches[0].digit = 0);
    alter!(p, p.nodes.last_mut().unwrap().branches[0].digit = 1);
    alter!(p, p.nodes.last_mut().unwrap().branches[0].child = 6);
    alter!(
        p,
        p.nodes.last_mut().unwrap().branches[0].child = usize::MAX
    );
    alter!(p, p.nodes.last_mut().unwrap().cell = Some(81));
    alter!(p, {
        let node = p.nodes[0].clone();
        p.nodes.insert(0, node);
    });
    alter!(p, {
        p.nodes.last_mut().unwrap().steps.pop();
    });
    alter!(p, {
        match &mut p.nodes.last_mut().unwrap().steps[0] {
            ForcingStep::House { edits, .. } | ForcingStep::ThreeHouse { edits, .. } => {
                edits[0].before = 511
            }
        }
    });
    alter!(p, {
        match &mut p.nodes.last_mut().unwrap().steps[0] {
            ForcingStep::House { edits, .. } | ForcingStep::ThreeHouse { edits, .. } => {
                edits[0].after = 511
            }
        }
    });
    alter!(p, {
        match &mut p.nodes.last_mut().unwrap().steps[0] {
            ForcingStep::House { sector, .. } => *sector = 27,
            ForcingStep::ThreeHouse { center, .. } => *center = 81,
        }
    });
    alter!(p, {
        match &mut p.nodes.last_mut().unwrap().steps[0] {
            ForcingStep::House { edits, .. } | ForcingStep::ThreeHouse { edits, .. } => {
                edits.clear()
            }
        }
    });
    let hall = proof
        .nodes
        .iter()
        .position(|n| matches!(n.leaf, Some(ForcingLeaf::Hall { .. })))
        .unwrap();
    let three = proof
        .nodes
        .iter()
        .position(|n| matches!(n.leaf, Some(ForcingLeaf::ThreeHouse { .. })))
        .unwrap();
    alter!(p, p.nodes[hall].cell = Some(0));
    alter!(p, {
        if let Some(ForcingLeaf::Hall { sector, .. }) = &mut p.nodes[hall].leaf {
            *sector = 27;
        }
    });
    alter!(p, {
        if let Some(ForcingLeaf::Hall { cells, .. }) = &mut p.nodes[hall].leaf {
            cells.clear();
        }
    });
    alter!(p, {
        if let Some(ForcingLeaf::Hall { cells, .. }) = &mut p.nodes[hall].leaf {
            cells.push(cells[0]);
        }
    });
    alter!(p, {
        if let Some(ForcingLeaf::Hall { cells, .. }) = &mut p.nodes[hall].leaf {
            cells[0] = 81;
        }
    });
    alter!(p, {
        if let Some(ForcingLeaf::ThreeHouse { blockers, .. }) = &mut p.nodes[three].leaf {
            blockers.clear();
        }
    });
    alter!(p, {
        if let Some(ForcingLeaf::ThreeHouse { blockers, .. }) = &mut p.nodes[three].leaf {
            blockers[0].cells = 0;
        }
    });
    alter!(p, {
        if let Some(ForcingLeaf::ThreeHouse { blockers, .. }) = &mut p.nodes[three].leaf {
            blockers[0].digits = 512;
        }
    });
    alter!(p, {
        if let Some(ForcingLeaf::ThreeHouse { blockers, .. }) = &mut p.nodes[three].leaf {
            blockers[0].atom = 6;
        }
    });
    alter!(p, {
        if let Some(ForcingLeaf::ThreeHouse { center, .. }) = &mut p.nodes[three].leaf {
            *center = 81;
        }
    });
    for (i, p) in corruptions.iter().enumerate() {
        assert!(!p.verify(&grid), "accepted corruption {i}");
    }
    let mut stale = grid.deep_clone();
    stale.cell_mut(Position::new(0, 7)).remove_candidate(3);
    let before = CandidateState::capture(&stale).unwrap();
    assert!(!proof.verify(&stale));
    assert!(proof.apply(&mut stale).is_err());
    assert_eq!(CandidateState::capture(&stale).unwrap(), before);
    assert!(!proof.verify(&Grid::new_with_constraints(vec![])));
    let mut raw = serde_json::to_value(proof).unwrap();
    raw["unexpected"] = true.into();
    assert!(serde_json::from_value::<ForcingProof>(raw).is_err());
}

#[test]
fn bounded_search_returns_only_checked_proofs_and_reports_cutoffs() {
    let grid = Grid::from_string(
        "530070000600195000098000060800060003400803001700020006060000280000419005000080079",
    )
    .unwrap();
    let solver = Solver::new();
    let options = ForcingSearchOptions {
        max_depth: 0,
        ..Default::default()
    };
    let result = solver.search_forcing(&grid, 2, 1, false, &options);
    assert!(result.error.is_none() && !result.budget_exhausted);
    let proof = result
        .proof
        .expect("wrong candidate fails local propagation");
    assert!(proof.verify(&grid));
    for options in [
        ForcingSearchOptions {
            max_nodes: 0,
            ..options.clone()
        },
        ForcingSearchOptions {
            max_steps: 0,
            ..options.clone()
        },
    ] {
        let r = solver.search_forcing(&grid, 2, 1, false, &options);
        assert!(r.proof.is_none() && r.budget_exhausted && r.error.is_none());
        assert!(r.nodes_visited <= options.max_nodes && r.steps_evaluated <= options.max_steps);
    }
    for options in [
        ForcingSearchOptions {
            max_depth: 9,
            ..options.clone()
        },
        ForcingSearchOptions {
            max_nodes: 4097,
            ..options.clone()
        },
        ForcingSearchOptions {
            max_steps: 32769,
            ..options.clone()
        },
    ] {
        assert!(solver
            .search_forcing(&grid, 2, 1, false, &options)
            .error
            .is_some());
    }
    for (cell, digit) in [(81, 1), (0, 1), (2, 0), (2, 10)] {
        assert!(solver
            .search_forcing(&grid, cell, digit, false, &options)
            .error
            .is_some());
    }
    let r = solver.search_forcing(
        &Grid::new_classic(),
        0,
        1,
        false,
        &ForcingSearchOptions {
            max_nodes: 1,
            max_depth: 1,
            ..Default::default()
        },
    );
    assert!(r.proof.is_none() && r.budget_exhausted);
    let r = solver.search_forcing(
        &grid,
        2,
        1,
        false,
        &ForcingSearchOptions {
            propagation: ForcingPropagation::ThreeHouse,
            ..options
        },
    );
    assert!(r.proof.unwrap().verify(&grid));
}

#[test]
fn native_forcing_search_rederives_the_seven_node_endpoint_with_two_split_levels() {
    let Replay { state, .. } = replay();
    let grid = state.to_grid().unwrap();
    let options = ForcingSearchOptions {
        max_depth: 2,
        max_nodes: 4096,
        max_steps: 32768,
        propagation: ForcingPropagation::ThreeHouse,
    };
    let result = Solver::new().search_forcing(&grid, 7, 3, false, &options);
    let proof = result
        .proof
        .expect("bounded native search must refute the known endpoint");
    assert!(proof.verify(&grid));
    assert!(proof.check(&grid).unwrap().depth <= 2);
    assert!(
        result.nodes_visited <= options.max_nodes && result.steps_evaluated <= options.max_steps
    );
}

#[test]
fn forcing_application_rejects_added_variant_conflicts_atomically() {
    use sudoku_core::KillerCageConstraint;
    let grid = Grid::from_string(
        "530070000600195000098000060800060003400803001700020006060000280000419005000080079",
    )
    .unwrap();
    let proof = Solver::new()
        .search_forcing(&grid, 40, 5, true, &ForcingSearchOptions::default())
        .proof
        .unwrap();
    let state = CandidateState::capture(&grid).unwrap();
    let mut variant = Grid::new_killer(vec![KillerCageConstraint::new(
        vec![Position::new(4, 4)],
        6,
    )]);
    for c in 0..81 {
        let p = Position::new(c / 9, c % 9);
        if state.values[c] != 0 {
            variant.set_cell_unchecked(p, Some(state.values[c]));
        } else {
            variant
                .cell_mut(p)
                .set_candidates(sudoku_core::BitSet::from_raw(state.domains[c] << 1));
        }
    }
    assert!(proof.verify(&variant));
    let before = CandidateState::capture(&variant).unwrap();
    assert!(proof.apply(&mut variant).is_err());
    assert_eq!(CandidateState::capture(&variant).unwrap(), before);
}
