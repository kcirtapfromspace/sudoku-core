//! Export real engine scenes into a self-contained, offline browser viewer.
use serde::{Deserialize, Serialize};
use std::{collections::HashSet, env, fs, path::PathBuf};
use sudoku_core::{
    CandidateState, ForcingProof, ForcingSearchOptions, GlobalColoringReplay, Grid, GridScene,
    Hint, HintType, Position, ProofCertificate, PropagationPolicy, Solver, SolverState, Technique,
    ThreeHouseReplay,
};

const EASY: &str =
    "530070000600195000098000060800060003400803001700020006060000280000419005000080079";

#[derive(Serialize)]
struct Example {
    title: String,
    origin: String,
    source_url: String,
    original_puzzle: String,
    scene: GridScene,
    proof_bundle: Option<sudoku_core::TechniqueExampleCorpus>,
}

fn example(title: &str, origin: &str, source_url: &str, puzzle: &str, scene: GridScene) -> Example {
    Example {
        title: title.into(),
        origin: origin.into(),
        source_url: source_url.into(),
        original_puzzle: puzzle.into(),
        scene,
        proof_bundle: None,
    }
}

fn hint(
    cell: usize,
    digit: u8,
    value: bool,
    technique: Technique,
    proof: ProofCertificate,
) -> Hint {
    let pos = Position::new(cell / 9, cell % 9);
    Hint {
        technique,
        hint_type: if value {
            HintType::SetValue { pos, value: digit }
        } else {
            HintType::EliminateCandidates {
                pos,
                values: vec![digit],
            }
        },
        explanation: format!(
            "{} {} at r{}c{}.",
            if value { "Place" } else { "Remove" },
            digit,
            pos.row + 1,
            pos.col + 1
        ),
        involved_cells: vec![pos],
        proof: Some(proof),
    }
}

fn set_given_flags(grid: &mut Grid, puzzle: &str) {
    for (cell, c) in puzzle.bytes().enumerate() {
        grid.cell_mut(Position::new(cell / 9, cell % 9))
            .set_given((b'1'..=b'9').contains(&c));
    }
}

// The completion is an additional example check, not a technique certificate.
fn check_endpoint(hint: &Hint, solution: &Grid) {
    match &hint.hint_type {
        HintType::SetValue { pos, value } => assert_eq!(solution.get(*pos), Some(*value)),
        HintType::EliminateCandidates { pos, values } => {
            assert!(!values.contains(&solution.get(*pos).unwrap()))
        }
    }
}

fn puzzle_examples(
    puzzle: &str,
    name: &str,
    source: &str,
    all_steps: bool,
    examples: &mut Vec<Example>,
    seen: &mut HashSet<Technique>,
) -> Result<(), String> {
    let grid = Grid::from_string(puzzle).ok_or("Invalid puzzle string")?;
    let solver = Solver::new();
    if solver.count_solutions(&grid, 2) != 1 {
        return Err(format!("{name}: example requires a unique puzzle"));
    }
    let solution = solver.solve(&grid).ok_or("Example has no completion")?;
    let uniqueness = sudoku_core::ExactCoverProof::search_unique(
        &grid,
        &sudoku_core::ExactCoverSearchOptions::default(),
    )
    .proof
    .ok_or("Original-puzzle uniqueness certificate exceeded its budget")?;
    let source = sudoku_core::TechniqueExamplePuzzle {
        id: name.into(),
        givens: puzzle.into(),
        source_url: source.into(),
        attribution: "User-supplied unique puzzle".into(),
        source_sha256: None,
        source_index_1based: None,
        uniqueness,
    };
    let mut session = SolverState::new(&grid)?;
    for n in 0..256 {
        let before = session.grid().clone();
        let Some(hint) = session.step_with_options(
            &solver,
            PropagationPolicy::Logical,
            &ForcingSearchOptions::default(),
        )?
        else {
            break;
        };
        check_endpoint(&hint, &solution);
        if all_steps || seen.insert(hint.technique) {
            let hint = sudoku_core::ExactCoverProof::certify_hint(
                &before,
                &hint,
                &sudoku_core::ExactCoverSearchOptions::default(),
            )?;
            let case = sudoku_core::TechniqueExample {
                puzzle: 0,
                state: CandidateState::capture(&before)?,
                hint: hint.clone(),
                certificate: hint
                    .proof
                    .clone()
                    .ok_or("Missing independent certificate")?,
            };
            let mut item = example(
                &format!("{name} · step {} · {}", n + 1, hint.technique),
                "Unique givens-only puzzle; independently checked setup and deduction.",
                &source.source_url,
                puzzle,
                GridScene::from_hint(&before, &hint)?,
            );
            let bundle = sudoku_core::TechniqueExampleCorpus {
                version: 1,
                puzzles: vec![source.clone()],
                examples: vec![case],
            };
            bundle.verify_subset()?;
            item.proof_bundle = Some(bundle);
            examples.push(item);
        }
    }
    Ok(())
}

fn gallery() -> Result<Vec<Example>, Box<dyn std::error::Error>> {
    let corpus: sudoku_core::TechniqueExampleCorpus =
        serde_json::from_str(include_str!("../tests/fixtures/technique_examples.json"))?;
    corpus.verify()?;
    let mut examples = vec![];
    for entry in &corpus.examples {
        let source = &corpus.puzzles[entry.puzzle];
        let grid = entry.grid(source)?;
        let h = entry.checked_hint();
        let start = CandidateState::capture(
            &Grid::from_string(&source.givens).ok_or("Invalid original givens")?,
        )?;
        let placements = (0..81)
            .filter(|&c| start.values[c] == 0 && entry.state.values[c] != 0)
            .count();
        let removals: u32 = (0..81)
            .filter(|&c| entry.state.values[c] == 0)
            .map(|c| (start.domains[c] & !entry.state.domains[c]).count_ones())
            .sum();
        let origin = format!("{}. Original-givens uniqueness proof authenticates {placements} setup placements and {removals} candidate removals. This selected detector can coexist with easier steps.", source.attribution.trim_end_matches('.'));
        let mut item = example(
            &format!("{} · {}", h.technique, source.id),
            &origin,
            &source.source_url,
            &source.givens,
            GridScene::from_hint(&grid, &h)?,
        );
        let mut single = entry.clone();
        single.puzzle = 0;
        item.proof_bundle = Some(sudoku_core::TechniqueExampleCorpus {
            version: 2,
            puzzles: vec![source.clone()],
            examples: vec![single],
        });
        examples.push(item);
    }
    let cases: serde_json::Value =
        serde_json::from_str(include_str!("../research/nonmodular/real_game_cases.json"))?;
    let replay: ThreeHouseReplay = serde_json::from_str(include_str!(
        "../research/nonmodular/three_house_example_replay.json"
    ))?;
    let case = &cases["cases"][0];
    let puzzle = case["puzzle"].as_str().unwrap();
    let mut grid = replay.state.to_grid()?;
    set_given_flags(&mut grid, puzzle);
    let p = replay.proof;
    let h = hint(
        p.cell,
        p.digit,
        false,
        Technique::ThreeHouseMatching,
        ProofCertificate::ThreeHouse(p),
    );
    check_endpoint(
        &h,
        &Solver::new()
            .solve(&Grid::from_string(puzzle).ok_or("Invalid puzzle string")?)
            .unwrap(),
    );
    examples.push(example("Three-house obstruction · Norvig top95 #20", "Unique published puzzle. Exact retained baseline trace: research/nonmodular/real_game_cases.json.", case["source_url"].as_str().unwrap(), puzzle, GridScene::from_hint(&grid, &h)?));

    let lead: serde_json::Value = serde_json::from_str(include_str!(
        "../research/global_obstructions/different_lead.json"
    ))?;
    let puzzle = lead["puzzle"].as_str().unwrap();
    let source = lead["source"]["source_url"].as_str().unwrap();
    let givens = Grid::from_string(puzzle).ok_or("Invalid puzzle string")?;
    assert_eq!(Solver::new().count_solutions(&givens, 2), 1);
    let solution = Solver::new().solve(&givens).unwrap();
    let replay: GlobalColoringReplay = serde_json::from_str(include_str!(
        "../research/global_obstructions/replay_different_lead.json"
    ))?;
    let mut grid = replay.state.to_grid()?;
    set_given_flags(&mut grid, puzzle);
    let p = replay.proof;
    let h = hint(
        p.cell,
        p.digit,
        p.value,
        Technique::GlobalColoring,
        ProofCertificate::GlobalColoring(p),
    );
    check_endpoint(&h, &solution);
    examples.push(example("20-cell coloring · Mith515", "Unique published puzzle at a retained research state. Known-method benchmark; no originality claim.", source, puzzle, GridScene::from_hint(&grid, &h)?));

    #[derive(Deserialize)]
    struct ForcingReplay {
        state: CandidateState,
        proof: ForcingProof,
    }
    let replay: ForcingReplay =
        serde_json::from_str(include_str!("../tests/fixtures/forcing_mith515.json"))?;
    let mut grid = replay.state.to_grid()?;
    set_given_flags(&mut grid, puzzle);
    let p = replay.proof;
    let h = hint(
        p.cell,
        p.digit,
        p.value,
        Technique::DynamicForcingChain,
        ProofCertificate::ForcingProof(p),
    );
    check_endpoint(&h, &solution);
    examples.push(example("Seven-node forcing proof · Mith515", "Same unique puzzle and candidate premises. Inspect every branch and all 84 recorded propagation events.", source, puzzle, GridScene::from_hint(&grid, &h)?));

    let mut broken = h.clone();
    if let Some(ProofCertificate::ForcingProof(p)) = &mut broken.proof {
        p.nodes.last_mut().unwrap().branches.pop();
    }
    examples.push(example("Rejected certificate · missing branch", "Deliberately altered proof on the real Mith515 state. This example demonstrates checker rejection.", source, puzzle, GridScene::from_hint(&grid, &broken)?));

    let easy = Grid::from_string(EASY).ok_or("Invalid example puzzle")?;
    let proof = Solver::new()
        .search_forcing(
            &easy,
            2,
            1,
            false,
            &ForcingSearchOptions {
                max_depth: 0,
                ..Default::default()
            },
        )
        .proof
        .ok_or("Missing easy forcing proof")?;
    let h = hint(
        proof.cell,
        proof.digit,
        proof.value,
        Technique::NishioForcingChain,
        ProofCertificate::ForcingProof(proof),
    );
    examples.push(example(
        "Candidate invalidation · r1c3≠1",
        "Unique givens-only puzzle. The opposite assumption fails native Hall propagation.",
        "",
        EASY,
        GridScene::from_hint(&easy, &h)?,
    ));

    let mut invalid = easy.clone();
    invalid.set_cell_unchecked(Position::new(0, 2), Some(5));
    examples.push(example("Invalid move · duplicate 5", "Deliberately incorrect move in the classic example. Red cells show installed-rule conflicts.", "", EASY, GridScene::from_validation(&invalid)?));
    examples.push(example("Board validation · placed values", "Valid partial board. This check concerns installed rules, not uniqueness or candidate notes.", "", EASY, GridScene::from_validation(&easy)?));

    // A concrete synthetic fixture supplements the unique-puzzle examples.
    let mut grid = Grid::new_classic();
    for (sector, keep) in [
        (0, vec![0, 3]),
        (12, vec![3, 30]),
        (3, vec![30, 28]),
        (10, vec![28, 10]),
        (18, vec![0, 10, 20]),
    ] {
        for cell in 0..81 {
            let belongs = if sector < 9 {
                cell / 9 == sector
            } else if sector < 18 {
                cell % 9 == sector - 9
            } else {
                Position::new(cell / 9, cell % 9).box_index() == sector - 18
            };
            if belongs && !keep.contains(&cell) {
                grid.cell_mut(Position::new(cell / 9, cell % 9))
                    .remove_candidate(1);
            }
        }
    }
    let p = sudoku_core::ArithmeticProof::from_terms(
        &grid,
        [0, 12, 3, 10, 18]
            .into_iter()
            .map(|sector| sudoku_core::ArithmeticTerm {
                requirement: sudoku_core::ArithmeticRequirement::SectorDigit { sector, digit: 1 },
                weight: 1,
            })
            .collect(),
        20,
        1,
        true,
    )
    .ok_or("Missing arithmetic fixture proof")?;
    let h = hint(
        p.cell,
        p.digit,
        p.value,
        Technique::ArithmeticCounting,
        ProofCertificate::Arithmetic(p),
    );
    examples.push(example("Arithmetic · five-equation guardian", "Synthetic candidate-state fixture. These pencil marks are supplied premises, not a unique givens-only puzzle.", "", "", GridScene::from_hint(&grid, &h)?));
    Ok(examples)
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut output = PathBuf::from("target/grid-viewer.html");
    let mut puzzle = None;
    let mut args = env::args().skip(1);
    while let Some(arg) = args.next() {
        match arg.as_str() {
            "--output" => output = args.next().ok_or("--output needs a path")?.into(),
            "--puzzle" => puzzle = Some(args.next().ok_or("--puzzle needs 81 cells")?),
            _ => {
                return Err(
                    format!("Unknown option: {arg}. Use --output PATH or --puzzle STRING.").into(),
                )
            }
        }
    }
    let examples = if let Some(puzzle) = puzzle {
        let mut examples = vec![];
        puzzle_examples(
            &puzzle,
            "Your puzzle",
            "",
            true,
            &mut examples,
            &mut HashSet::new(),
        )?;
        if examples.is_empty() {
            return Err("No logical steps found for this puzzle".into());
        }
        examples
    } else {
        gallery()?
    };
    // Script-safe JSON. The viewer uses textContent for all scene text.
    let json = serde_json::to_string(&examples)?
        .replace('<', "\\u003c")
        .replace('>', "\\u003e")
        .replace('&', "\\u0026");
    let html = include_str!("grid_viewer/viewer.html")
        .replace("/* VIEWER_CSS */", include_str!("grid_viewer/viewer.css"))
        .replace("/* VIEWER_JS */", include_str!("grid_viewer/viewer.js"))
        .replace("<!-- SCENE_DATA -->", &json);
    if let Some(parent) = output.parent().filter(|p| !p.as_os_str().is_empty()) {
        fs::create_dir_all(parent)?;
    }
    fs::write(&output, html)?;
    println!("Wrote {} examples to {}", examples.len(), output.display());
    Ok(())
}
