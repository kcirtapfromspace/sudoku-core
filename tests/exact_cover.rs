use sudoku_core::{
    CandidateState, EvidenceStatus, ExactCoverProof, ExactCoverSearchOptions, ExactNode,
    ExactPurpose, ForcingSearchOptions, FrameKind, Grid, GridScene, Hint, HintType, Position,
    ProofCertificate, Solver, Technique, TechniqueExampleCorpus,
};

const EASY: &str =
    "530070000600195000098000060800060003400803001700020006060000280000419005000080079";
fn corpus() -> TechniqueExampleCorpus {
    serde_json::from_str(include_str!("fixtures/technique_examples.json")).unwrap()
}
fn hint(values: Vec<u8>, placement: bool) -> Hint {
    Hint {
        technique: Technique::Backtracking,
        hint_type: if placement {
            HintType::SetValue {
                pos: Position::new(0, 2),
                value: values[0],
            }
        } else {
            HintType::EliminateCandidates {
                pos: Position::new(0, 2),
                values,
            }
        },
        explanation: "Independent endpoint test".into(),
        involved_cells: vec![Position::new(0, 2)],
        proof: None,
    }
}

#[test]
fn all_48_real_examples_have_replayable_deductions_and_authenticated_setup() {
    let corpus = corpus();
    corpus.verify().unwrap();
    assert_eq!(corpus.examples.len(), 48);
    assert_eq!(corpus.puzzles.len(), 11);
    let solver = Solver::new();
    for case in &corpus.examples {
        let grid = case.grid(&corpus.puzzles[case.puzzle]).unwrap();
        let state = CandidateState::capture(&grid).unwrap();
        let hint = case.checked_hint();
        let scene = GridScene::from_hint(&grid, &hint).unwrap();
        assert_eq!(
            scene.evidence.status,
            EvidenceStatus::CheckedPatternAndDeduction
        );
        assert!(scene
            .frames
            .iter()
            .any(|f| f.kind == FrameKind::Contradiction));
        assert_eq!(scene.frames.last().unwrap().kind, FrameKind::Conclusion);
        for frame in scene
            .frames
            .iter()
            .filter(|f| f.kind == FrameKind::Contradiction)
        {
            for (cell, snapshot) in frame.cells.iter().enumerate() {
                if state.values[cell] == 0 && snapshot.value.is_some() {
                    assert!(frame.highlights.iter().any(|h| h.cell == cell
                        && h.digit == snapshot.value
                        && h.role == sudoku_core::HighlightRole::Assumption));
                }
            }
        }
        // The recorded named method must still be produced by its actual detector.
        // This is a regression test, separate from replay-only proof verification.
        let current = solver
            .get_hint_for_technique(
                &grid,
                hint.technique,
                &ForcingSearchOptions {
                    max_depth: 2,
                    max_nodes: 4096,
                    max_steps: 32768,
                    ..Default::default()
                },
            )
            .unwrap()
            .unwrap_or_else(|| panic!("Frozen state must trigger {:?}", hint.technique));
        assert_eq!(current.technique, hint.technique);
        ExactCoverProof::certify_hint(&grid, &current, &ExactCoverSearchOptions::default())
            .unwrap();
        assert_eq!(CandidateState::capture(&grid).unwrap(), state);
    }
    let roundtrip: TechniqueExampleCorpus =
        serde_json::from_str(&serde_json::to_string(&corpus).unwrap()).unwrap();
    roundtrip.verify().unwrap();
}

#[test]
fn independent_replay_rejects_missing_cases_wrong_requirements_and_forged_witnesses() {
    let corpus = corpus();
    let puzzle = &corpus.puzzles[0];
    let grid = Grid::from_string(&puzzle.givens).unwrap();
    let proof = &puzzle.uniqueness;
    let mut bad = vec![];
    macro_rules! alter {
        ($p:ident, $body:expr) => {{
            let mut $p = proof.clone();
            $body;
            bad.push($p);
        }};
    }
    alter!(p, p.version = 2);
    alter!(p, p.state_hash = "0".repeat(64));
    alter!(p, p.nodes.clear());
    alter!(p, p.witness[0] = 10);
    alter!(p, p.witness.pop());
    alter!(p, p.witness[0] = p.witness[1]);
    alter!(p, *p.nodes.last_mut().unwrap() = ExactNode::Solution);
    alter!(
        p,
        *p.nodes.last_mut().unwrap() = ExactNode::Conflict { requirement: 324 }
    );
    alter!(
        p,
        *p.nodes.last_mut().unwrap() = ExactNode::Conflict { requirement: 0 }
    );
    alter!(p, {
        if let ExactNode::Split { branches, .. } = p.nodes.last_mut().unwrap() {
            branches.pop();
        }
    });
    alter!(p, {
        if let ExactNode::Split { branches, .. } = p.nodes.last_mut().unwrap() {
            branches[0].candidate = 729;
        }
    });
    alter!(p, {
        let root = p.nodes.len() - 1;
        if let ExactNode::Split { branches, .. } = p.nodes.last_mut().unwrap() {
            branches[0].child = root;
        }
    });
    alter!(p, {
        p.nodes.push(p.nodes.last().unwrap().clone());
    });
    for (i, p) in bad.iter().enumerate() {
        assert!(!p.verify(&grid), "accepted corruption {i}");
    }
    let mut changed = grid.clone();
    let cell = changed.empty_positions()[0];
    changed
        .cell_mut(cell)
        .remove_candidate(proof.witness[cell.row * 9 + cell.col]);
    assert!(!proof.verify(&changed));
}

#[test]
fn cutoffs_counterexamples_and_inconsistent_premises_never_become_proofs() {
    let grid = Grid::from_string(EASY).unwrap();
    let options = ExactCoverSearchOptions::default();
    let wrong = ExactCoverProof::search_hint(&grid, &hint(vec![1], true), &options);
    assert!(wrong.proof.is_none() && wrong.counterexample.is_some() && !wrong.budget_exhausted);
    let ambiguous = ExactCoverProof::search_unique(&Grid::new_classic(), &options);
    assert!(ambiguous.proof.is_none() && ambiguous.counterexample.is_some());
    let cutoff = ExactCoverProof::search_unique(&grid, &ExactCoverSearchOptions { max_nodes: 0 });
    assert!(cutoff.proof.is_none() && cutoff.budget_exhausted);
    assert!(ExactCoverProof::search_unique(
        &grid,
        &ExactCoverSearchOptions {
            max_nodes: 2_000_001
        }
    )
    .error
    .is_some());
    assert!(
        ExactCoverProof::search_unique(&Grid::new_x_sudoku(), &options)
            .error
            .is_some()
    );
    let mut impossible = Grid::new_classic();
    for c in 0..2 {
        impossible
            .cell_mut(Position::new(0, c))
            .set_candidates(sudoku_core::BitSet::from_slice(&[1]));
    }
    assert!(
        ExactCoverProof::search_hint(&impossible, &hint(vec![1], false), &options)
            .proof
            .is_none()
    );
    let empty = Grid::new_classic();
    let fallback = Solver::new()
        .get_hint_for_technique(&empty, Technique::Backtracking, &Default::default())
        .unwrap()
        .unwrap();
    assert!(ExactCoverProof::certify_hint(&empty, &fallback, &options).is_err());
    assert!(
        ExactCoverProof::search_hint(&grid, &hint(vec![0], true), &options)
            .error
            .is_some()
    );
    assert!(
        ExactCoverProof::search_hint(&grid, &hint(vec![1, 1], false), &options)
            .error
            .is_some()
    );
    let mut out_of_bounds = hint(vec![1], true);
    out_of_bounds.hint_type = HintType::SetValue {
        pos: Position::new(9, 0),
        value: 1,
    };
    assert!(
        ExactCoverProof::search_hint(&grid, &out_of_bounds, &options)
            .error
            .is_some()
    );
}

#[test]
fn multi_candidate_proofs_cover_the_entire_opposite_and_bind_the_exact_endpoint() {
    let grid = Grid::from_string(EASY).unwrap();
    let checked =
        ExactCoverProof::certify_hint(&grid, &hint(vec![1, 2], false), &Default::default())
            .unwrap();
    let scene = GridScene::from_hint(&grid, &checked).unwrap();
    assert_eq!(scene.evidence.status, EvidenceStatus::CheckedDeduction);
    assert_eq!(scene.frames[1].cells[2].candidates, vec![1, 2]);
    let mut mismatch = checked.clone();
    mismatch.hint_type = HintType::EliminateCandidates {
        pos: Position::new(0, 2),
        values: vec![1],
    };
    assert_eq!(
        GridScene::from_hint(&grid, &mismatch)
            .unwrap()
            .evidence
            .status,
        EvidenceStatus::Rejected
    );
    let ProofCertificate::Certified { deduction, .. } = checked.proof.as_ref().unwrap() else {
        panic!()
    };
    for purpose in [
        ExactPurpose::Deduction {
            cell: 81,
            digits: vec![1],
            placement: false,
        },
        ExactPurpose::Deduction {
            cell: 2,
            digits: vec![0],
            placement: false,
        },
        ExactPurpose::Deduction {
            cell: 2,
            digits: vec![1, 1],
            placement: false,
        },
        ExactPurpose::Unique,
    ] {
        let mut p = deduction.clone();
        p.purpose = purpose;
        assert!(!p.verify(&grid));
    }
    let mut corrupted = checked;
    if let Some(ProofCertificate::Certified { deduction, .. }) = &mut corrupted.proof {
        deduction.nodes.pop();
    }
    assert_eq!(
        GridScene::from_hint(&grid, &corrupted)
            .unwrap()
            .evidence
            .status,
        EvidenceStatus::Rejected
    );
}

#[test]
fn corpus_verification_rejects_coverage_gaps_and_forged_setup_even_with_a_valid_endpoint() {
    let original = corpus();
    let mut partial = original.clone();
    partial.examples.pop();
    assert!(partial.verify().is_err());
    let mut duplicate = original.clone();
    duplicate.examples[1] = duplicate.examples[0].clone();
    assert!(duplicate.verify().is_err());
    let mut setup = original.clone();
    let example = &mut setup.examples[0];
    let witness = &setup.puzzles[example.puzzle].uniqueness.witness;
    let cell = example.state.values.iter().position(|&v| v == 0).unwrap();
    example.state.domains[cell] &= !(1 << (witness[cell] - 1));
    assert!(setup.verify().is_err());
    let mut nested = original.clone();
    if let ProofCertificate::Certified {
        pattern, deduction, ..
    } = &mut nested.examples[0].certificate
    {
        *pattern = Some(Box::new(ProofCertificate::Certified {
            named: None,
            pattern: None,
            deduction: deduction.clone(),
        }));
    }
    assert!(nested.verify().is_err());
    let example = &nested.examples[0];
    let grid = example.grid(&nested.puzzles[example.puzzle]).unwrap();
    assert_eq!(
        GridScene::from_hint(&grid, &example.checked_hint())
            .unwrap()
            .evidence
            .status,
        EvidenceStatus::Rejected
    );
    let mut single = original.examples[0].clone();
    let puzzle = original.puzzles[single.puzzle].clone();
    single.puzzle = 0;
    let bundle = TechniqueExampleCorpus {
        version: 1,
        puzzles: vec![puzzle],
        examples: vec![single],
    };
    bundle.verify_subset().unwrap();
    assert!(bundle.verify().is_err());
    let mut raw = serde_json::to_value(bundle).unwrap();
    raw["unexpected"] = true.into();
    assert!(serde_json::from_value::<TechniqueExampleCorpus>(raw).is_err());
}

#[test]
fn a_missing_branch_in_a_real_bug_refutation_cannot_pass_replay() {
    let corpus = corpus();
    let e = corpus
        .examples
        .iter()
        .find(|e| e.hint.technique == Technique::BivalueUniversalGrave)
        .unwrap();
    let grid = e.grid(&corpus.puzzles[e.puzzle]).unwrap();
    let ProofCertificate::Certified { deduction, .. } = &e.certificate else {
        panic!()
    };
    let mut missing = deduction.clone();
    let ExactNode::Split { branches, .. } = missing.nodes.last_mut().unwrap() else {
        panic!()
    };
    assert!(branches.len() > 1);
    branches.pop();
    assert!(!missing.verify(&grid));
    let mut repeated = deduction.clone();
    let ExactNode::Split { branches, .. } = repeated.nodes.last_mut().unwrap() else {
        panic!()
    };
    branches[1] = branches[0].clone();
    assert!(!repeated.verify(&grid));
}
