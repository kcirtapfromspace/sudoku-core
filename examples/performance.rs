//! Repeatable, dependency-free timings for the generation validation path.
use std::hint::black_box;
use std::time::Instant;
use sudoku_core::{Difficulty, Generator, Grid, Position, Solver};

const EASY: &str =
    "530070000600195000098000060800060003400803001700020006060000280000419005000080079";
const HARD: &str =
    "800000000003600000070090200050007000000045700000100030001000068008500010090000400";
const SOLUTION: &str =
    "534678912672195348198342567859761423426853791713924856961537284287419635345286179";

fn main() {
    let mode = std::env::args().nth(1).unwrap_or_else(|| "all".into());
    assert!(matches!(mode.as_str(), "all" | "count" | "generate"));
    let solver = Solver::new();

    if mode != "generate" {
        let mut ambiguous = Grid::from_string(SOLUTION).unwrap();
        for index in [3, 4, 30, 31] {
            ambiguous.set_cell_unchecked(Position::new(index / 9, index % 9), None);
        }
        for (name, grid, iterations, expected) in [
            ("easy", Grid::from_string(EASY).unwrap(), 1000, 1),
            ("hard", Grid::from_string(HARD).unwrap(), 100, 1),
            ("ambiguous", ambiguous, 1000, 2),
            ("empty", Grid::new_classic(), 100, 2),
        ] {
            let start = Instant::now();
            for _ in 0..iterations {
                assert_eq!(
                    black_box(solver.count_solutions(black_box(&grid), 2)),
                    expected
                );
            }
            println!(
                "count/{name}: {:.2} us/puzzle ({iterations} iterations)",
                start.elapsed().as_secs_f64() * 1e6 / iterations as f64
            );
        }
    }

    if mode != "count" {
        for difficulty in [
            Difficulty::Beginner,
            Difficulty::Easy,
            Difficulty::Medium,
            Difficulty::Hard,
            Difficulty::Expert,
        ] {
            // Time generation alone. Check correctness and report puzzle identity
            // outside the timed interval so revisions can be compared fairly.
            let mut elapsed = 0.0;
            let mut fingerprint = 0_u64;
            for seed in [42, 99, 314] {
                let mut generator = Generator::with_seed(seed);
                let start = Instant::now();
                let puzzle = black_box(generator.generate(difficulty));
                elapsed += start.elapsed().as_secs_f64();
                assert!(puzzle.validate().is_valid);
                assert!(solver.has_unique_solution(&puzzle));
                for byte in puzzle.to_string_compact().bytes() {
                    fingerprint = fingerprint.wrapping_mul(31).wrapping_add(byte as u64);
                }
            }
            println!(
                "generate/{difficulty:?}: {:.2} ms/puzzle (3 seeds, fingerprint {fingerprint:016x})",
                elapsed * 1e3 / 3.0
            );
        }
    }
}
