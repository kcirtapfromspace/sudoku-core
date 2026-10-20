//! Freeze mined detector examples with exact uniqueness and deduction certificates.
use std::{collections::HashMap, env, fs};
use sudoku_core::{
    CandidateState, ExactCoverProof, ExactCoverSearchOptions, Grid, Hint, NamedPatternProof,
    PatternUniqueness, PatternWitness, ProofCertificate, Technique, TechniqueExample,
    TechniqueExampleCorpus, TechniqueExamplePuzzle,
};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<_> = env::args().collect();
    let input = args
        .get(1)
        .ok_or("Pass the mined-example directory and output JSON path")?;
    let output = args.get(2).ok_or("Pass an output JSON path")?;
    let mut corpus = TechniqueExampleCorpus {
        version: 2,
        puzzles: vec![],
        examples: vec![],
    };
    let mut puzzles = HashMap::new();
    let options = ExactCoverSearchOptions {
        max_nodes: 2_000_000,
    };
    for &technique in Technique::all() {
        let record: serde_json::Value =
            serde_json::from_str(&fs::read_to_string(format!("{input}/{technique:?}.json"))?)?;
        let source = &record["source"];
        let givens = source["puzzle"].as_str().ok_or("No original givens")?;
        let index = if let Some(&i) = puzzles.get(givens) {
            i
        } else {
            let grid = Grid::from_string(givens).ok_or("Invalid puzzle")?;
            let mut result = ExactCoverProof::search_unique(&grid, &options);
            let proof = result.proof.take().ok_or_else(|| {
                format!(
                    "Uniqueness not certified for {:?}: {result:?}",
                    source["id"]
                )
            })?;
            eprintln!("Unique {}: {} exact nodes", source["id"], proof.nodes.len());
            let index = corpus.puzzles.len();
            corpus.puzzles.push(TechniqueExamplePuzzle {
                id: source["id"].as_str().ok_or("No source id")?.into(),
                givens: givens.into(),
                source_url: source["source_url"].as_str().unwrap_or("").into(),
                attribution: source["attribution"].as_str().unwrap_or("").into(),
                source_sha256: source["source_sha256"].as_str().map(str::to_owned),
                source_index_1based: source["source_index_1based"].as_u64().map(|n| n as usize),
                uniqueness: proof,
            });
            puzzles.insert(givens.to_string(), index);
            index
        };
        let state: CandidateState = serde_json::from_value(source["state"].clone())?;
        let mut hint: Hint = serde_json::from_value(record["hint"].clone())?;
        hint.proof = serde_json::from_value::<Option<ProofCertificate>>(record["pattern"].clone())?;
        let stub = TechniqueExample {
            puzzle: index,
            state,
            hint: hint.clone(),
            certificate: ProofCertificate::Backtracking,
        };
        let grid = stub.grid(&corpus.puzzles[index])?;
        let checked = ExactCoverProof::certify_hint(&grid, &hint, &options)?;
        let mut certificate = checked.proof.unwrap();
        if let Some(raw) = record.get("named_pattern").filter(|v| !v.is_null()) {
            let witness: PatternWitness = if raw.as_str() == Some("Search") {
                let ProofCertificate::Certified { deduction, .. } = &certificate else {
                    unreachable!()
                };
                PatternWitness::Search {
                    proof: deduction.clone(),
                }
            } else {
                serde_json::from_value(raw.clone())?
            };
            if let ProofCertificate::Certified {
                named, deduction, ..
            } = &mut certificate
            {
                let unique = matches!(
                    technique,
                    Technique::UniqueRectangle
                        | Technique::HiddenRectangle
                        | Technique::AvoidableRectangle
                        | Technique::ExtendedUniqueRectangle
                        | Technique::BivalueUniversalGrave
                );
                let proof = NamedPatternProof {
                    version: 1,
                    state_hash: deduction.state_hash.clone(),
                    technique,
                    conclusion: deduction.purpose.clone(),
                    witness,
                    uniqueness: unique.then(|| PatternUniqueness {
                        givens: givens.into(),
                        proof: corpus.puzzles[index].uniqueness.clone(),
                    }),
                };
                proof.check(&grid, &hint)?;
                *named = Some(Box::new(proof));
            }
        }
        if let ProofCertificate::Certified { deduction, .. } = &certificate {
            eprintln!("{technique:?}: {} refutation nodes", deduction.nodes.len());
        }
        corpus.examples.push(TechniqueExample {
            certificate,
            ..stub
        });
    }
    corpus.verify()?;
    fs::write(output, serde_json::to_string_pretty(&corpus)?)?;
    eprintln!(
        "Certified all {} techniques in {} original puzzles",
        corpus.examples.len(),
        corpus.puzzles.len()
    );
    Ok(())
}
