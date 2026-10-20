//! Replay proof JSON without running any named detector or proof producer.
use std::{env, fs};
use sudoku_core::TechniqueExampleCorpus;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<_> = env::args().skip(1).collect();
    let source = if let Some(path) = args.iter().find(|a| !a.starts_with("--")) {
        fs::read_to_string(path)?
    } else {
        include_str!("../tests/fixtures/technique_examples.json").to_owned()
    };
    let corpus: TechniqueExampleCorpus = serde_json::from_str(&source)?;
    if args.iter().any(|a| a == "--subset") {
        corpus.verify_subset()?;
    } else {
        corpus.verify()?;
    }
    println!(
        "Verified {} technique examples and {} original-givens uniqueness proofs",
        corpus.examples.len(),
        corpus.puzzles.len()
    );
    Ok(())
}
