use super::*;
use crate::{BitSet, Position, Solver};

fn real_fixture() -> (Grid, Vec<usize>, usize, u8) {
    let json: serde_json::Value = serde_json::from_str(include_str!(
        "../../research/global_obstructions/real_4401.json"
    ))
    .unwrap();
    let state: CandidateState = serde_json::from_value(json["state"].clone()).unwrap();
    let cells = json["graph_cells"]
        .as_array()
        .unwrap()
        .iter()
        .map(|c| c.as_u64().unwrap() as usize)
        .collect();
    (state.to_grid().unwrap(), cells, 30, 9)
}
fn four_fixture() -> (Grid, Vec<usize>, usize, u8) {
    let json: serde_json::Value = serde_json::from_str(include_str!(
        "../../research/global_obstructions/structural_four_digit_13.json"
    ))
    .unwrap();
    let domains: Vec<u16> = serde_json::from_value(json["base_domains"].clone()).unwrap();
    let mut grid = Grid::new_classic();
    for (cell, domain) in domains.into_iter().enumerate() {
        grid.cell_mut(idx_to_pos(cell))
            .set_candidates(BitSet::from_raw(domain << 1));
    }
    let cells: Vec<usize> = serde_json::from_value(json["support"].clone()).unwrap();
    (grid, cells, 16, 1)
}
fn certificate() -> (Grid, GlobalColoringProof) {
    let (grid, cells, cell, digit) = real_fixture();
    let proof = GlobalColoringProof::from_cells(&grid, cells, cell, digit, true).unwrap();
    (grid, proof)
}

#[test]
fn full_candidate_ladder_advances_the_frozen_old_engine_stall() {
    let case: serde_json::Value = serde_json::from_str(include_str!(
        "../../research/global_obstructions/nonrank_full_old_stall.json"
    ))
    .unwrap();
    let state: CandidateState = serde_json::from_value(case["state"].clone()).unwrap();
    let mut grid = state.to_grid().unwrap();
    for (cell, byte) in case["puzzle"].as_str().unwrap().bytes().enumerate() {
        grid.cell_mut(idx_to_pos(cell))
            .set_given((b'1'..=b'9').contains(&byte));
    }
    // Use the actual logical dispatcher on retained candidates. Public get_hint
    // rebuilds masks; the direct coloring API and replay preserve these premises.
    let finding = Solver::new().find_first_technique(&grid).unwrap();
    assert_eq!(finding.technique, Technique::GlobalColoring);
    assert!(matches!(
        finding.inference,
        InferenceResult::Placement { cell: 35, value: 9 }
    ));
    let Some(ProofCertificate::GlobalColoring(proof)) = finding.proof else {
        panic!("full ladder omitted the global certificate")
    };
    assert!(proof.verify(&grid));
    assert_eq!(CandidateState::capture(&grid).unwrap(), state);
}

#[test]
fn real_and_four_digit_obstructions_generate_checked_trees() {
    for (grid, cells, cell, digit) in [real_fixture(), four_fixture()] {
        let before = CandidateState::capture(&grid).unwrap();
        let p = GlobalColoringProof::from_cells(&grid, cells, cell, digit, true).unwrap();
        let check = p.check(&grid).unwrap();
        assert_eq!(check.nodes, p.nodes.len());
        assert!(check.splits > 0 && check.conflicts > 0);
        assert!(p.verify(&grid));
        let start = std::time::Instant::now();
        let result =
            Solver::new().search_global_coloring(&grid, &GlobalColoringSearchOptions::default());
        println!(
            "global search {:?}, {} patterns, {} nodes",
            start.elapsed(),
            result.patterns_tested,
            result.nodes_visited
        );
        assert!(result.error.is_none());
        let hint = result
            .hint
            .expect("bounded guardian search must find the fixture");
        assert_eq!(hint.technique, Technique::GlobalColoring);
        let Some(ProofCertificate::GlobalColoring(proof)) = &hint.proof else {
            panic!("missing global proof")
        };
        assert!(proof.verify(&grid));
        assert_eq!(CandidateState::capture(&grid).unwrap(), before);
    }
}

fn brute_color(cells: &[usize], domains: &[u16]) -> bool {
    fn walk(cells: &[usize], domains: &[u16], assigned: &mut Vec<u8>) -> bool {
        let i = assigned.len();
        if i == cells.len() {
            return true;
        }
        for digit in 1..=9 {
            if domains[i] & (1 << (digit - 1)) == 0
                || (0..i).any(|j| peers(cells[i], cells[j]) && assigned[j] == digit)
            {
                continue;
            }
            assigned.push(digit);
            if walk(cells, domains, assigned) {
                return true;
            }
            assigned.pop();
        }
        false
    }
    walk(cells, domains, &mut Vec::new())
}

#[test]
fn exhaustive_small_assignments_agree_with_search_and_tree_checker() {
    let fabric = super::super::fabric::CandidateFabric::from_grid(&Grid::new_classic());
    for a in 0..81 {
        for b in 0..81 {
            assert_eq!(peers(a, b), a != b && fabric.sees(a, b));
        }
    }
    let mut seed = 0x1ab342u32;
    for n in 1..=7 {
        for _ in 0..120 {
            let mut cells = Vec::new();
            while cells.len() < n {
                seed = seed.wrapping_mul(1664525).wrapping_add(1013904223);
                let cell = (seed >> 9) as usize % 81;
                if !cells.contains(&cell) {
                    cells.push(cell);
                }
            }
            cells.sort_unstable();
            let domains: Vec<u16> = cells
                .iter()
                .map(|_| {
                    seed = seed.wrapping_mul(1664525).wrapping_add(1013904223);
                    ((seed >> 12) as u16 % 7) + 1
                })
                .collect();
            let brute = brute_color(&cells, &domains);
            let mut nodes = Vec::new();
            let mut budget = Budget::new(50_000);
            match solve(
                &Graph::new(&cells),
                &domains,
                0,
                &mut nodes,
                &mut budget,
                &mut 0,
            ) {
                Outcome::Supported(witness) => {
                    assert!(brute);
                    for i in 0..n {
                        assert_ne!(domains[i] & (1 << (witness[i] - 1)), 0);
                        for j in 0..i {
                            if peers(cells[i], cells[j]) {
                                assert_ne!(witness[i], witness[j]);
                            }
                        }
                    }
                }
                Outcome::Refuted(root) => {
                    assert!(!brute);
                    assert_eq!(root + 1, nodes.len());
                    let proof = GlobalColoringProof {
                        version: 1,
                        state_hash: String::new(),
                        cells: cells.clone(),
                        cell: cells[0],
                        digit: 1,
                        value: true,
                        witness: vec![],
                        nodes,
                    };
                    let mut check = GlobalColoringCheck {
                        cells: n,
                        nodes: 0,
                        splits: 0,
                        conflicts: 0,
                    };
                    let mut visited = vec![false; proof.nodes.len()];
                    assert!(proof
                        .check_node(root, &domains, 0, &mut visited, &mut check)
                        .is_some());
                    assert!(visited.into_iter().all(|v| v));
                }
                Outcome::Cutoff => panic!("small exhaustive test exceeded budget"),
            }
        }
    }
}

#[test]
fn malformed_certificates_and_incomplete_branch_coverage_are_rejected() {
    let (grid, proof) = certificate();
    let mut mutations = Vec::new();
    macro_rules! bad {
        ($field:ident, $value:expr) => {{
            let mut p = proof.clone();
            p.$field = $value;
            mutations.push(p);
        }};
    }
    bad!(version, 2);
    bad!(state_hash, "stale".into());
    bad!(cells, vec![]);
    bad!(cells, (0..25).collect());
    bad!(cell, 81);
    bad!(digit, 0);
    bad!(digit, 10);
    bad!(digit, 1);
    bad!(value, false);
    bad!(witness, vec![]);
    bad!(witness, vec![1; proof.cells.len()]);
    bad!(nodes, vec![]);
    bad!(
        nodes,
        vec![ColoringNode::Conflict { cell: 30 }; MAX_PROOF_NODES + 1]
    );
    let mut p = proof.clone();
    p.cells[0] = 81;
    mutations.push(p);
    let mut p = proof.clone();
    p.cells[0] = p.cells[1];
    mutations.push(p);
    let mut p = proof.clone();
    p.cells.reverse();
    mutations.push(p);
    let mut p = proof.clone();
    p.witness[0] = 0;
    mutations.push(p);
    let mut p = proof.clone();
    p.witness[0] = 10;
    mutations.push(p);
    let mut p = proof.clone();
    p.nodes[0] = ColoringNode::Conflict { cell: 0 };
    mutations.push(p);
    // Replace a genuine conflict with a nonempty domain (root cannot be a conflict).
    bad!(nodes, vec![ColoringNode::Conflict { cell: proof.cell }]);
    let split = proof
        .nodes
        .iter()
        .position(|n| matches!(n,ColoringNode::Split{branches,..} if branches.len()>1))
        .unwrap();
    for kind in 0..9 {
        let mut p = proof.clone();
        let ColoringNode::Split { cell, branches } = &mut p.nodes[split] else {
            unreachable!()
        };
        match kind {
            0 => {
                branches.pop();
            }
            1 => {
                branches.push(branches[0].clone());
            }
            2 => {
                branches[0].digit = 0;
            }
            3 => {
                branches[0].digit = 10;
            }
            4 => {
                branches[0].child = split;
            }
            5 => {
                branches[0].child = usize::MAX;
            }
            6 => {
                branches[1].child = branches[0].child;
            }
            7 => {
                *cell = 81;
            }
            _ => {
                branches.clear();
            }
        }
        mutations.push(p);
    }
    let mut p = proof.clone();
    p.nodes.insert(0, ColoringNode::Conflict { cell: 0 });
    for node in &mut p.nodes[1..] {
        if let ColoringNode::Split { branches, .. } = node {
            for b in branches {
                b.child += 1;
            }
        }
    }
    mutations.push(p); // unused node
    for p in mutations {
        assert!(!p.verify(&grid), "malformed certificate accepted");
    }
    let mut changed = grid.deep_clone();
    changed.cell_mut(idx_to_pos(0)).remove_candidate(1);
    assert!(!proof.verify(&changed));
    let mut json = serde_json::to_value(&proof).unwrap();
    json["unknown"] = true.into();
    assert!(serde_json::from_value::<GlobalColoringProof>(json).is_err());
    assert!(GlobalColoringProof::from_cells(&grid, vec![30, 30], 30, 9, true).is_none());
    assert!(GlobalColoringProof::from_cells(&grid, vec![30], 30, 9, false).is_none());
    assert!(GlobalColoringProof::from_cells(&grid, vec![30], 30, 1, true).is_none());
}

#[test]
fn portable_replay_preserves_masks_for_both_result_polarities() {
    let (grid, proof) = certificate();
    let mut cells = proof.cells.clone();
    cells.push(3);
    let elimination = GlobalColoringProof::from_cells(&grid, cells, 3, 9, false).unwrap();
    for p in [proof, elimination] {
        assert!(p.verify(&grid));
        let replay = GlobalColoringReplay::capture(&grid, &p).unwrap();
        let decoded: GlobalColoringReplay =
            serde_json::from_str(&serde_json::to_string(&replay).unwrap()).unwrap();
        assert_eq!(decoded, replay);
        assert!(decoded.verify());
        let before = CandidateState::capture(&grid).unwrap();
        let mut next = grid.deep_clone();
        decoded.apply(&mut next).unwrap();
        let after = CandidateState::capture(&next).unwrap();
        for cell in 0..81 {
            let expected = if cell == p.cell {
                if p.value {
                    1 << (p.digit - 1)
                } else {
                    before.domains[cell] & !(1 << (p.digit - 1))
                }
            } else if p.value && before.values[cell] == 0 && peers(cell, p.cell) {
                before.domains[cell] & !(1 << (p.digit - 1))
            } else {
                before.domains[cell]
            };
            assert_eq!(after.domains[cell], expected);
        }
        assert_eq!(after.values[p.cell], if p.value { p.digit } else { 0 });
        assert!(decoded.apply(&mut next).is_err());
        assert_eq!(CandidateState::capture(&next).unwrap(), after);
        let mut stale = grid.deep_clone();
        stale.cell_mut(idx_to_pos(0)).remove_candidate(1);
        let stale_before = CandidateState::capture(&stale).unwrap();
        assert!(decoded.apply(&mut stale).is_err());
        assert_eq!(CandidateState::capture(&stale).unwrap(), stale_before);
        let mut malformed = replay.clone();
        malformed.version = 99;
        assert!(!malformed.verify());
        let mut malformed = replay.clone();
        malformed.state.values.clear();
        assert!(!malformed.verify());
        let mut malformed = replay;
        malformed.state.domains.clear();
        assert!(!malformed.verify());
    }
}

#[test]
fn bounds_and_inconsistent_unforced_graphs_never_yield_proofs() {
    let (grid, cells, cell, digit) = real_fixture();
    for options in [
        GlobalColoringSearchOptions {
            max_patterns: 0,
            max_nodes: 50_000,
        },
        GlobalColoringSearchOptions {
            max_patterns: 512,
            max_nodes: 0,
        },
        GlobalColoringSearchOptions {
            max_patterns: 512,
            max_nodes: 2,
        },
    ] {
        let result = Solver::new().search_global_coloring(&grid, &options);
        assert!(result.hint.is_none() && result.budget_exhausted && result.error.is_none());
    }
    for options in [
        GlobalColoringSearchOptions {
            max_patterns: MAX_PATTERNS + 1,
            max_nodes: 0,
        },
        GlobalColoringSearchOptions {
            max_patterns: 0,
            max_nodes: MAX_SEARCH_NODES + 1,
        },
    ] {
        assert!(Solver::new()
            .search_global_coloring(&grid, &options)
            .error
            .is_some());
    }
    let mut budget = Budget::new(50_000);
    let mut calls = MAX_PROOF_NODES;
    assert!(matches!(
        solve(
            &Graph::new(&cells),
            &vec![0x1ff; cells.len()],
            0,
            &mut Vec::new(),
            &mut budget,
            &mut calls
        ),
        Outcome::Cutoff
    ));
    assert!(budget.exhausted);
    let blank = Grid::new_classic();
    assert!(Solver::new().get_global_coloring_hint(&blank).is_none());
    let mut large = Grid::new_classic();
    for c in 0..81 {
        let row = c / 9;
        let col = c % 9;
        // Leave one full-domain cell per row, column and box.
        if col != (row * 3 + row / 3) % 9 {
            large
                .cell_mut(idx_to_pos(c))
                .set_candidates(BitSet::from_slice(&[1, 2, 3]));
        }
    }
    large
        .cell_mut(idx_to_pos(1))
        .set_candidates(BitSet::from_slice(&[1, 2, 3, 4]));
    let result =
        Solver::new().search_global_coloring(&large, &GlobalColoringSearchOptions::default());
    assert!(result.error.is_none() && result.hint.is_none() && result.oversized_components > 0);
    let custom = Grid::new_with_constraints(vec![]);
    assert!(Solver::new()
        .search_global_coloring(&custom, &GlobalColoringSearchOptions::default())
        .error
        .is_some());
    assert!(GlobalColoringProof::from_cells(&custom, cells.clone(), cell, digit, true).is_none());
    let mut bad = blank.deep_clone();
    bad.cell_mut(idx_to_pos(0))
        .set_candidates(BitSet::from_slice(&[1, 2]));
    for c in [1, 2] {
        bad.cell_mut(idx_to_pos(c))
            .set_candidates(BitSet::single(1));
    }
    assert!(CandidateState::capture(&bad).is_ok());
    assert!(GlobalColoringProof::from_cells(&bad, vec![0, 1, 2], 0, 2, true).is_none());
    bad.cell_mut(idx_to_pos(0)).set_candidates(BitSet::empty());
    assert!(Solver::new()
        .search_global_coloring(&bad, &GlobalColoringSearchOptions::default())
        .error
        .is_some());
}

#[test]
fn native_provenance_variant_and_local_contradiction_rejections_are_atomic() {
    let (grid, proof) = certificate();
    let raw: Grid = serde_json::from_str(&serde_json::to_string(&grid).unwrap()).unwrap();
    assert!(!proof.verify(&raw));
    let (grid, cells, cell, digit) = four_fixture();
    let state = CandidateState::capture(&grid).unwrap();
    let mut variant = Grid::new_x_sudoku();
    for c in 0..81 {
        variant
            .cell_mut(idx_to_pos(c))
            .set_candidates(BitSet::from_raw(state.domains[c] << 1));
    }
    let p = GlobalColoringProof::from_cells(&variant, cells.clone(), cell, digit, true).unwrap();
    assert!(p.verify(&variant));
    GlobalColoringReplay::capture(&variant, &p)
        .unwrap()
        .apply(&mut variant)
        .unwrap();
    let mut variant = Grid::new_x_sudoku();
    for c in 0..81 {
        variant
            .cell_mut(idx_to_pos(c))
            .set_candidates(BitSet::from_raw(state.domains[c] << 1));
    }
    variant.set_cell_unchecked(Position::new(3, 5), Some(1)); // same anti-diagonal, no native house shared
    let p = GlobalColoringProof::from_cells(&variant, cells, cell, digit, true).unwrap();
    let replay = GlobalColoringReplay::capture(&variant, &p).unwrap();
    let before = CandidateState::capture(&variant).unwrap();
    assert!(replay.apply(&mut variant).is_err());
    assert_eq!(CandidateState::capture(&variant).unwrap(), before);
    // A conflict outside the selected graph must not allow a partial mutation.
    let (mut grid, proof) = certificate();
    grid.cell_mut(idx_to_pos(3))
        .set_candidates(BitSet::single(9));
    let p = GlobalColoringProof::from_cells(&grid, proof.cells, 30, 9, true).unwrap();
    let replay = GlobalColoringReplay::capture(&grid, &p).unwrap();
    let before = CandidateState::capture(&grid).unwrap();
    assert!(replay.apply(&mut grid).is_err());
    assert_eq!(CandidateState::capture(&grid).unwrap(), before);
}
