//! Basic example of using the Sudoku engine

use sudoku_core::{Difficulty, Generator, Grid, Solver};

fn main() {
    // Generate a puzzle
    println!("Generating a Medium difficulty puzzle...\n");
    let mut generator = Generator::new();
    let puzzle = generator.generate(Difficulty::Medium);

    println!("Generated puzzle:");
    println!("{}", puzzle);

    // Show some stats
    println!("Given cells: {}", puzzle.given_count());
    println!("Empty cells: {}", puzzle.empty_count());

    // Rate the difficulty
    let solver = Solver::new();
    let actual_difficulty = solver.rate_difficulty(&puzzle);
    println!("Rated difficulty: {}\n", actual_difficulty);

    // Solve it
    println!("Solving...\n");
    if let Some(solution) = solver.solve(&puzzle) {
        println!("Solution:");
        println!("{}", solution);
    } else {
        println!("No solution found (this shouldn't happen for a generated puzzle!)");
    }

    // Get a hint for the puzzle
    println!("\nGetting a hint for the original puzzle:");
    if let Some(hint) = solver.get_hint(&puzzle) {
        println!("Technique: {}", hint.technique);
        println!("Explanation: {}", hint.explanation);
    }

    // Parse a puzzle from a string
    println!("\n--- Parsing a puzzle from string ---\n");
    let puzzle_string =
        "530070000600195000098000060800060003400803001700020006060000280000419005000080079";
    if let Some(grid) = Grid::from_string(puzzle_string) {
        println!("Parsed puzzle:");
        println!("{}", grid);

        // Check uniqueness
        let solutions = solver.count_solutions(&grid, 2);
        println!("Number of solutions (up to 2): {}", solutions);
    }
}
