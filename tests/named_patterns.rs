use sudoku_core::{
    EvidenceStatus, GridScene, HintType, NamedPatternProof, PatternWitness, ProofCertificate,
    Technique, TechniqueExampleCorpus,
};
fn corpus() -> TechniqueExampleCorpus {
    serde_json::from_str(include_str!("fixtures/technique_examples.json")).unwrap()
}
#[test]
fn all_named_proofs_bind_the_state_label_action_and_pattern() {
    let corpus = corpus();
    corpus.verify().unwrap();
    for e in &corpus.examples {
        let grid = e.grid(&corpus.puzzles[e.puzzle]).unwrap();
        let ProofCertificate::Certified {
            named: Some(proof), ..
        } = &e.certificate
        else {
            panic!("Missing named proof")
        };
        proof.check(&grid, &e.hint).unwrap();
        let mut bad = proof.as_ref().clone();
        bad.version = 2;
        assert!(!bad.verify(&grid, &e.hint));
        bad = proof.as_ref().clone();
        bad.state_hash.clear();
        assert!(!bad.verify(&grid, &e.hint));
        bad = proof.as_ref().clone();
        bad.technique = if e.hint.technique == Technique::Backtracking {
            Technique::NakedSingle
        } else {
            Technique::Backtracking
        };
        assert!(!bad.verify(&grid, &e.hint));
        let mut relabeled = e.hint.clone();
        relabeled.technique = bad.technique;
        assert!(
            !bad.verify(&grid, &relabeled),
            "Wrong family accepted for {:?}",
            e.hint.technique
        );
        let mut other = e.hint.clone();
        match &mut other.hint_type {
            HintType::SetValue { value, .. } => *value = if *value == 9 { 1 } else { *value + 1 },
            HintType::EliminateCandidates { values, .. } => values.clear(),
        };
        assert!(!proof.verify(&grid, &other));
        bad = proof.as_ref().clone();
        bad.uniqueness = if bad.uniqueness.is_some() {
            None
        } else {
            corpus.examples.iter().find_map(|e| {
                if let ProofCertificate::Certified { named: Some(p), .. } = &e.certificate {
                    p.uniqueness.clone()
                } else {
                    None
                }
            })
        };
        assert!(!bad.verify(&grid, &e.hint));
        bad = proof.as_ref().clone();
        match &mut bad.witness {
            PatternWitness::Subset { house, .. } => *house = 27,
            PatternWitness::Fish { parts } => parts[0].covers[0] = 27,
            PatternWitness::Wing { cells, .. } => cells[0] = cells[1],
            PatternWitness::Als { sets, .. } => sets[0].digits.clear(),
            PatternWitness::Chain { nodes, .. } => nodes[0] = 729,
            PatternWitness::Trade { cells, .. } => cells[0] = 81,
            PatternWitness::Bug => bad.uniqueness.as_mut().unwrap().proof.nodes.clear(),
            PatternWitness::Aligned { exclusions, .. } => exclusions.clear(),
            PatternWitness::Forcing { branches, .. } => branches[0].steps.clear(),
            PatternWitness::Native { proof } => match proof.as_mut() {
                ProofCertificate::Arithmetic(p) => p.state_hash.clear(),
                ProofCertificate::ThreeHouse(p) => p.blockers.clear(),
                ProofCertificate::GlobalColoring(p) => p.nodes.clear(),
                _ => panic!(),
            },
            PatternWitness::Search { proof } => proof.nodes.clear(),
        }
        assert!(
            !bad.verify(&grid, &e.hint),
            "Malformed witness accepted for {:?}",
            e.hint.technique
        );
        let mut corrupt = e.checked_hint();
        if let Some(ProofCertificate::Certified { named, .. }) = &mut corrupt.proof {
            *named = Some(Box::new(bad));
        }
        assert_eq!(
            GridScene::from_hint(&grid, &corrupt)
                .unwrap()
                .evidence
                .status,
            EvidenceStatus::Rejected
        );
    }
}
#[test]
fn a_valid_endpoint_does_not_establish_a_named_pattern() {
    let mut corpus = corpus();
    let ProofCertificate::Certified { named, .. } = &mut corpus.examples[0].certificate else {
        panic!()
    };
    *named = None;
    assert!(corpus.verify().is_err());
    assert!(corpus.verify_subset().is_err());
    let e = &corpus.examples[0];
    let grid = e.grid(&corpus.puzzles[e.puzzle]).unwrap();
    assert_eq!(
        GridScene::from_hint(&grid, &e.checked_hint())
            .unwrap()
            .evidence
            .status,
        EvidenceStatus::CheckedDeduction
    );
    corpus.version = 1;
    corpus.verify_subset().unwrap();
    assert!(corpus.verify().is_err());
}
#[test]
fn pattern_deserialization_rejects_unknown_fields() {
    let c = corpus();
    let ProofCertificate::Certified { named: Some(p), .. } = &c.examples[0].certificate else {
        panic!()
    };
    let mut value = serde_json::to_value(p).unwrap();
    value["unverified_theorem"] = true.into();
    assert!(serde_json::from_value::<NamedPatternProof>(value).is_err());
    let mut value = serde_json::to_value(p).unwrap();
    value["witness"]["Subset"]["unchecked"] = true.into();
    assert!(serde_json::from_value::<NamedPatternProof>(value).is_err());
}
