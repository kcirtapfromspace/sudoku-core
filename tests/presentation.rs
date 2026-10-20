use sudoku_core::{
    AlsProofDescriptor, BitSet, CandidateState, EvidenceStatus, ForcingProof, ForcingSearchOptions,
    FrameKind, GlobalColoringReplay, Grid, GridScene, HighlightRole, Hint, HintType, LinkType,
    Polarity, Position, ProofCertificate, Solver, Technique, ThreeHouseReplay,
};

const EASY: &str =
    "530070000600195000098000060800060003400803001700020006060000280000419005000080079";

fn hint(cell: usize, digit: u8, value: bool, proof: Option<ProofCertificate>) -> Hint {
    let pos = Position::new(cell / 9, cell % 9);
    Hint {
        technique: Technique::DynamicForcingChain,
        hint_type: if value {
            HintType::SetValue { pos, value: digit }
        } else {
            HintType::EliminateCandidates {
                pos,
                values: vec![digit],
            }
        },
        explanation: "Inspect the supplied evidence.".into(),
        involved_cells: vec![pos],
        proof,
    }
}

fn forcing() -> (Grid, Hint) {
    #[derive(serde::Deserialize)]
    struct Replay {
        state: CandidateState,
        proof: ForcingProof,
    }
    let r: Replay = serde_json::from_str(include_str!("fixtures/forcing_mith515.json")).unwrap();
    (
        r.state.to_grid().unwrap(),
        hint(
            r.proof.cell,
            r.proof.digit,
            r.proof.value,
            Some(ProofCertificate::ForcingProof(r.proof)),
        ),
    )
}

#[test]
fn every_registered_technique_renders_without_claiming_metadata_is_a_proof() {
    // Exercise the public wire names, including retired methods and fallback.
    let names = "NakedSingle HiddenSingle NakedPair HiddenPair NakedTriple HiddenTriple PointingPair BoxLineReduction XWing FinnedXWing Swordfish FinnedSwordfish Jellyfish FinnedJellyfish NakedQuad HiddenQuad EmptyRectangle AvoidableRectangle UniqueRectangle HiddenRectangle XYWing XYZWing WXYZWing WWing XChain ThreeDMedusa SueDeCoq AIC FrankenFish SiameseFish AlsXz ExtendedUniqueRectangle BivalueUniversalGrave AlsXyWing AlsChain MutantFish AlignedPairExclusion AlignedTripletExclusion DeathBlossom ThreeHouseMatching ArithmeticCounting GlobalColoring NishioForcingChain KrakenFish RegionForcingChain CellForcingChain DynamicForcingChain Backtracking";
    let mut grid = Grid::from_string(EASY).unwrap();
    // Retain an earlier deletion and ensure presentation does not restore it.
    grid.cell_mut(Position::new(0, 2)).remove_candidate(4);
    let before = CandidateState::capture(&grid).unwrap();
    assert_eq!(names.split_whitespace().count(), 48);
    for name in names.split_whitespace() {
        let mut h = hint(2, 1, false, None);
        h.technique = serde_json::from_value(name.into()).unwrap();
        let scene = GridScene::from_hint(&grid, &h).unwrap();
        assert_eq!(
            scene.evidence.status,
            EvidenceStatus::MetadataOnly,
            "{name}"
        );
        assert_eq!(scene.frames.len(), 2);
        assert_eq!(scene.frames[0].cells[2].candidates, vec![1, 2]);
        assert_eq!(scene.frames[1].cells[2].candidates, vec![2]);
        let json = serde_json::to_value(&scene).unwrap();
        assert_eq!(json["technique"], name);
        assert_eq!(json["frames"][0]["cells"].as_array().unwrap().len(), 81);
        assert!(scene.frames[0]
            .highlights
            .iter()
            .any(|h| h.cell == 2 && h.digit == Some(1) && h.role == HighlightRole::Elimination));
        assert_eq!(CandidateState::capture(&grid).unwrap(), before);
    }
}

#[test]
fn forcing_trace_contains_all_steps_and_branches_without_mutating_premises() {
    let (grid, h) = forcing();
    let before = CandidateState::capture(&grid).unwrap();
    let scene = GridScene::from_hint(&grid, &h).unwrap();
    assert_eq!(scene.evidence.status, EvidenceStatus::CheckedCertificate);
    assert!(!scene.trace_truncated);
    assert_eq!(
        scene
            .frames
            .iter()
            .filter(|f| f.kind == FrameKind::Propagation)
            .count(),
        84
    );
    assert_eq!(
        scene
            .frames
            .iter()
            .filter(|f| f.kind == FrameKind::Contradiction)
            .count(),
        4
    );
    assert_eq!(
        scene
            .frames
            .iter()
            .filter(|f| f.kind == FrameKind::Assumption)
            .count(),
        7
    );
    assert_eq!(scene.frames.last().unwrap().kind, FrameKind::Conclusion);
    assert!(!scene.frames.last().unwrap().cells[7]
        .candidates
        .contains(&3));
    assert_eq!(CandidateState::capture(&grid).unwrap(), before);
    // A proof of a different endpoint must never lend this hint a green badge.
    let mut mismatch = h.clone();
    mismatch.hint_type = HintType::EliminateCandidates {
        pos: Position::new(0, 7),
        values: vec![3, 8],
    };
    assert_eq!(
        GridScene::from_hint(&grid, &mismatch)
            .unwrap()
            .evidence
            .status,
        EvidenceStatus::Rejected
    );
    let mut corrupt = h;
    if let Some(ProofCertificate::ForcingProof(p)) = &mut corrupt.proof {
        p.nodes.last_mut().unwrap().branches.pop();
    }
    let rejected = GridScene::from_hint(&grid, &corrupt).unwrap();
    assert_eq!(rejected.evidence.status, EvidenceStatus::Rejected);
    assert_eq!(rejected.frames.last().unwrap().kind, FrameKind::Validation);
    assert_eq!(
        rejected.frames.last().unwrap().cells,
        rejected.frames[0].cells
    );
}

#[test]
fn native_coloring_and_three_house_scenes_are_checked_and_large_trees_are_bounded() {
    let replay: ThreeHouseReplay = serde_json::from_str(include_str!(
        "../research/nonmodular/three_house_example_replay.json"
    ))
    .unwrap();
    let p = replay.proof;
    let scene = GridScene::from_hint(
        &replay.state.to_grid().unwrap(),
        &hint(
            p.cell,
            p.digit,
            false,
            Some(ProofCertificate::ThreeHouse(p)),
        ),
    )
    .unwrap();
    assert_eq!(scene.evidence.status, EvidenceStatus::CheckedCertificate);
    assert!(scene
        .frames
        .iter()
        .any(|f| f.kind == FrameKind::Contradiction));
    for raw in [
        include_str!("../research/global_obstructions/replay_different_lead.json"),
        include_str!("../research/global_obstructions/replay_native27.json"),
    ] {
        let replay: GlobalColoringReplay = serde_json::from_str(raw).unwrap();
        let p = replay.proof;
        let large = p.nodes.len() > 1000;
        let grid = replay.state.to_grid().unwrap();
        let scene = GridScene::from_hint(
            &grid,
            &hint(
                p.cell,
                p.digit,
                p.value,
                Some(ProofCertificate::GlobalColoring(p)),
            ),
        )
        .unwrap();
        assert_eq!(scene.evidence.status, EvidenceStatus::CheckedCertificate);
        assert_eq!(scene.trace_truncated, large);
        assert!(scene.frames.len() <= 256);
        assert_eq!(scene.frames.last().unwrap().kind, FrameKind::Conclusion);
        assert!(scene
            .frames
            .iter()
            .any(|f| f.kind == FrameKind::Contradiction));
        if large {
            assert_eq!(scene.frames.len(), 256);
        }
    }
}

#[test]
fn structure_annotations_preserve_candidate_polarity_links_and_pattern_roles() {
    let grid = Grid::new_classic();
    let structures = [
        ProofCertificate::Aic {
            chain: vec![
                (0, 1, Polarity::On),
                (1, 1, Polarity::Off),
                (1, 2, Polarity::On),
            ],
            link_types: vec![LinkType::WeakInference, LinkType::Strong],
        },
        ProofCertificate::Fish {
            digit: 1,
            base_sectors: vec![0, 1],
            cover_sectors: vec![9, 12],
            fins: vec![2],
        },
        ProofCertificate::Als {
            als_chain: vec![AlsProofDescriptor {
                cells: vec![0, 1],
                candidates: vec![1, 2, 3],
                sector: 0,
            }],
            rcc_values: vec![1],
            z_value: Some(3),
        },
        ProofCertificate::Uniqueness {
            pattern: "rectangle".into(),
            floor_cells: vec![0, 1],
            roof_cells: vec![9, 10],
        },
    ];
    let mut roles = vec![];
    for proof in structures {
        let scene = GridScene::from_hint(&grid, &hint(20, 1, false, Some(proof))).unwrap();
        assert_eq!(scene.evidence.status, EvidenceStatus::MetadataOnly);
        roles.extend(scene.frames[0].highlights.iter().map(|h| h.role));
        if !scene.frames[0].links.is_empty() {
            assert_eq!(scene.frames[0].links.len(), 2);
            assert_eq!(
                scene.frames[0].links[0].kind,
                sudoku_core::SceneLinkKind::Weak
            );
            assert_eq!(
                scene.frames[0].links[1].kind,
                sudoku_core::SceneLinkKind::Strong
            );
        }
    }
    for role in [
        HighlightRole::On,
        HighlightRole::Off,
        HighlightRole::Fin,
        HighlightRole::Base,
        HighlightRole::Cover,
        HighlightRole::Floor,
        HighlightRole::Roof,
    ] {
        assert!(roles.contains(&role));
    }
    // Malformed structural metadata remains unverified and cannot export bad SVG coordinates.
    let bad = ProofCertificate::Aic {
        chain: vec![(usize::MAX, 255, Polarity::On), (0, 1, Polarity::Off)],
        link_types: vec![LinkType::Strong],
    };
    let scene = GridScene::from_hint(&grid, &hint(20, 1, false, Some(bad))).unwrap();
    assert!(scene.frames[0].links.is_empty());
    assert!(scene.frames[0].highlights.iter().all(|h| h.cell < 81));
    let mut restricted = grid.clone();
    restricted
        .cell_mut(Position::new(0, 0))
        .set_candidates(BitSet::from_slice(&[1, 2]));
    let als = ProofCertificate::Als {
        als_chain: vec![AlsProofDescriptor {
            cells: vec![0, 1],
            candidates: vec![1, 2, 3],
            sector: 0,
        }],
        rcc_values: vec![1],
        z_value: Some(3),
    };
    let scene = GridScene::from_hint(&restricted, &hint(20, 1, false, Some(als))).unwrap();
    assert!(!scene.frames[0]
        .highlights
        .iter()
        .any(|h| h.cell == 0 && h.digit == Some(3)));
}

#[test]
fn board_validation_marks_the_actual_conflicts_and_rejects_malformed_endpoints() {
    let mut grid = Grid::from_string(EASY).unwrap();
    let valid = GridScene::from_validation(&grid).unwrap();
    assert_eq!(valid.evidence.status, EvidenceStatus::RulesSatisfied);
    assert_eq!(valid.frames[0].cells[0].candidates, vec![5]);
    grid.set_cell_unchecked(Position::new(0, 2), Some(5));
    let invalid = GridScene::from_validation(&grid).unwrap();
    assert_eq!(invalid.evidence.status, EvidenceStatus::RuleViolation);
    for c in [0, 2] {
        assert!(invalid.frames[0]
            .highlights
            .iter()
            .any(|h| h.cell == c && h.role == HighlightRole::Conflict));
    }
    let grid = Grid::new_classic();
    for (cell, digit) in [(81, 1), (0, 0), (0, 10)] {
        assert!(GridScene::from_hint(&grid, &hint(cell, digit, false, None)).is_err());
    }
    for values in [vec![], vec![1, 1]] {
        let mut h = hint(0, 1, false, None);
        h.hint_type = HintType::EliminateCandidates {
            pos: Position::new(0, 0),
            values,
        };
        assert!(GridScene::from_hint(&grid, &h).is_err());
    }
    let mut malformed = grid.clone();
    malformed
        .cell_mut(Position::new(0, 0))
        .set_candidates(BitSet::from_raw(1));
    assert!(GridScene::from_validation(&malformed).is_err());
}

#[test]
fn a_native_certificate_does_not_hide_an_added_variant_rule_violation() {
    let grid = Grid::from_string(EASY).unwrap();
    let p = Solver::new()
        .search_forcing(&grid, 40, 5, true, &ForcingSearchOptions::default())
        .proof
        .unwrap();
    let state = CandidateState::capture(&grid).unwrap();
    let mut variant = Grid::new_killer(vec![sudoku_core::KillerCageConstraint::new(
        vec![Position::new(4, 4)],
        6,
    )]);
    for c in 0..81 {
        let pos = Position::new(c / 9, c % 9);
        if state.values[c] != 0 {
            variant.set_cell_unchecked(pos, Some(state.values[c]));
        }
    }
    for c in 0..81 {
        if state.values[c] == 0 {
            variant
                .cell_mut(Position::new(c / 9, c % 9))
                .set_candidates(BitSet::from_raw(state.domains[c] << 1));
        }
    }
    assert!(p.verify(&variant));
    let scene = GridScene::from_hint(
        &variant,
        &hint(40, 5, true, Some(ProofCertificate::ForcingProof(p))),
    )
    .unwrap();
    assert_eq!(scene.evidence.status, EvidenceStatus::Rejected);
    assert!(scene.frames.iter().all(|f| f.kind != FrameKind::Conclusion));
}

#[test]
fn arithmetic_terminals_and_retained_guardian_propagation_have_checked_grid_evidence() {
    use sudoku_core::{
        ArithmeticProof, ArithmeticRequirement, ArithmeticTerm, GlobalColoringSearchOptions,
        SolverState,
    };
    let terms = |sectors: &[usize], weights: &[i8]| {
        sectors
            .iter()
            .zip(weights)
            .map(|(&sector, &weight)| ArithmeticTerm {
                requirement: ArithmeticRequirement::SectorDigit { sector, digit: 1 },
                weight,
            })
            .collect::<Vec<_>>()
    };
    let restrict_one = |grid: &mut Grid, sector: usize, retained: &[usize]| {
        for c in 0..81 {
            let belongs = if sector < 9 {
                c / 9 == sector
            } else if sector < 18 {
                c % 9 == sector - 9
            } else {
                Position::new(c / 9, c % 9).box_index() == sector - 18
            };
            if belongs && !retained.contains(&c) {
                grid.cell_mut(Position::new(c / 9, c % 9))
                    .remove_candidate(1);
            }
        }
    };
    let mut guardian = Grid::new_classic();
    for (s, cells) in [
        (0, vec![0, 3]),
        (12, vec![3, 30]),
        (3, vec![30, 28]),
        (10, vec![28, 10]),
        (18, vec![0, 10, 20]),
    ] {
        restrict_one(&mut guardian, s, &cells);
    }
    let divisibility =
        ArithmeticProof::from_terms(&guardian, terms(&[0, 12, 3, 10, 18], &[1; 5]), 20, 1, true)
            .unwrap();
    let mut xwing = Grid::new_classic();
    restrict_one(&mut xwing, 0, &[0, 6]);
    restrict_one(&mut xwing, 1, &[9, 15]);
    let interval =
        ArithmeticProof::from_terms(&xwing, terms(&[9, 15, 0, 1], &[1, 1, -1, -1]), 27, 1, false)
            .unwrap();
    let mut residue_grid = Grid::new_classic();
    for s in [12, 3, 10, 18, 1] {
        restrict_one(&mut residue_grid, s, &[3, 30, 28, 1, 9, 12]);
    }
    residue_grid
        .cell_mut(Position::new(0, 3))
        .set_candidates(BitSet::from_slice(&[1, 2, 3]));
    let mut residue_terms = terms(&[12, 3, 10, 18, 1], &[2; 5]);
    residue_terms.push(ArithmeticTerm {
        requirement: ArithmeticRequirement::Cell { cell: 3 },
        weight: 1,
    });
    let residue =
        ArithmeticProof::from_residue_terms(&residue_grid, residue_terms, 3, 2, false, 4).unwrap();
    for (grid, proof, expected) in [
        (guardian, divisibility, "divisible"),
        (xwing, interval, "interval"),
        (residue_grid, residue, "remainder"),
    ] {
        let h = hint(
            proof.cell,
            proof.digit,
            proof.value,
            Some(ProofCertificate::Arithmetic(proof)),
        );
        let scene = GridScene::from_hint(&grid, &h).unwrap();
        assert_eq!(scene.evidence.status, EvidenceStatus::CheckedCertificate);
        assert!(scene
            .frames
            .iter()
            .any(|f| f.kind == FrameKind::Contradiction && f.description.contains(expected)));
    }
    let mut grid = Grid::new_classic();
    for c in 0..4 {
        grid.cell_mut(Position::new(0, c))
            .set_candidates(BitSet::from_raw(if c < 2 { 7 << 1 } else { 15 << 1 }));
    }
    let mut state = SolverState::from_candidates(&grid).unwrap();
    state
        .learn_guardian(&GlobalColoringSearchOptions {
            max_patterns: 1,
            max_nodes: 256,
        })
        .unwrap()
        .unwrap();
    let h = state.get_hint(&Solver::new()).unwrap().unwrap();
    assert!(matches!(h.proof, Some(ProofCertificate::Guardian(_))));
    let scene = GridScene::from_hint(&grid, &h).unwrap();
    assert_eq!(scene.evidence.status, EvidenceStatus::CheckedCertificate);
    assert!(scene.frames.iter().any(|f| f.title.contains("OR relation")));
}
