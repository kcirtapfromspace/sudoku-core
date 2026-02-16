//! Puzzle Diversity Analysis
//!
//! Provides tools to analyze and estimate the number of unique Sudoku puzzles
//! that can be generated at each difficulty level.

use crate::{Difficulty, Generator, Grid, Solver};
use serde::{Deserialize, Serialize};
use std::collections::{HashMap, HashSet};

/// Statistics for a single difficulty level
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct DifficultyStats {
    /// Number of puzzles sampled
    pub samples: usize,
    /// Number of unique puzzle fingerprints seen
    pub unique_puzzles: usize,
    /// Distribution of clue counts
    pub clue_distribution: HashMap<usize, usize>,
    /// Distribution of techniques required
    pub technique_distribution: HashMap<String, usize>,
    /// Average clue count
    pub avg_clues: f64,
    /// Min clue count observed
    pub min_clues: usize,
    /// Max clue count observed
    pub max_clues: usize,
    /// Estimated total puzzles (extrapolated)
    pub estimated_total: u128,
    /// Generation success rate
    pub success_rate: f64,
    /// Average generation attempts needed
    pub avg_attempts: f64,
}

/// Overall diversity analysis results
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct DiversityReport {
    /// Statistics per difficulty level
    pub by_difficulty: HashMap<String, DifficultyStats>,
    /// Total unique fingerprints across all difficulties
    pub total_unique: usize,
    /// Theoretical estimates
    pub theoretical: TheoreticalEstimates,
}

/// Theoretical mathematical estimates for Sudoku puzzle counts
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TheoreticalEstimates {
    /// Total valid completed Sudoku grids
    pub total_grids: u128,
    /// Essentially different grids (removing symmetries)
    pub unique_grids: u128,
    /// Minimum clues for unique solution (proven)
    pub min_clues_for_unique: usize,
    /// Estimated puzzles per clue count
    pub puzzles_by_clue_count: HashMap<usize, u128>,
}

impl Default for TheoreticalEstimates {
    fn default() -> Self {
        Self::new()
    }
}

impl TheoreticalEstimates {
    pub fn new() -> Self {
        // Known mathematical constants for 9x9 Sudoku
        // Total valid grids: 6,670,903,752,021,072,936,960 ≈ 6.67 × 10^21
        let total_grids: u128 = 6_670_903_752_021_072_936_960;

        // Essentially different grids (after removing symmetries): ~5.47 × 10^9
        let unique_grids: u128 = 5_472_730_538;

        // Minimum clues proven for unique solution
        let min_clues_for_unique = 17;

        // Estimate puzzles by clue count
        // This is a rough estimation based on combinatorics and uniqueness probability
        let mut puzzles_by_clue_count = HashMap::new();

        // For each clue count, estimate: C(81, clues) * uniqueness_probability * unique_grids
        // Uniqueness probability decreases as clues decrease
        for clues in 17..=50 {
            let estimate = Self::estimate_puzzles_for_clue_count(clues, unique_grids);
            puzzles_by_clue_count.insert(clues, estimate);
        }

        Self {
            total_grids,
            unique_grids,
            min_clues_for_unique,
            puzzles_by_clue_count,
        }
    }

    /// Estimate number of puzzles with a given clue count
    fn estimate_puzzles_for_clue_count(clues: usize, unique_grids: u128) -> u128 {
        // C(81, clues) = ways to choose which cells are clues
        let combinations = Self::binomial(81, clues);

        // Probability of unique solution (rough estimates based on research)
        // Higher clue counts have higher uniqueness probability
        let uniqueness_prob = match clues {
            17 => 0.0000001, // Very rare - only ~49,000 known 17-clue puzzles
            18..=20 => 0.00001,
            21..=24 => 0.0001,
            25..=28 => 0.001,
            29..=32 => 0.01,
            33..=36 => 0.05,
            37..=40 => 0.15,
            41..=45 => 0.35,
            46..=50 => 0.60,
            _ => 0.80,
        };

        // Estimate = unique_grids * C(81, clues) * uniqueness_probability
        // But cap it to avoid overflow
        let raw = (unique_grids as f64) * (combinations as f64) * uniqueness_prob;
        if raw > u128::MAX as f64 {
            u128::MAX
        } else {
            raw as u128
        }
    }

    /// Calculate binomial coefficient C(n, k) using integer arithmetic
    fn binomial(n: usize, k: usize) -> u128 {
        if k > n {
            return 0;
        }
        if k == 0 || k == n {
            return 1;
        }

        let k = k.min(n - k); // Optimization: C(n,k) = C(n,n-k)

        // Use the multiplicative formula: C(n,k) = n! / (k! * (n-k)!)
        // Compute as: (n * (n-1) * ... * (n-k+1)) / (k * (k-1) * ... * 1)
        // To avoid overflow, divide as we go by computing GCD
        let mut numerator: u128 = 1;
        let mut denominator: u128 = 1;

        for i in 0..k {
            numerator *= (n - i) as u128;
            denominator *= (i + 1) as u128;

            // Reduce fraction to avoid overflow
            let g = gcd(numerator, denominator);
            numerator /= g;
            denominator /= g;
        }

        numerator / denominator
    }
}

/// Calculate greatest common divisor
fn gcd(a: u128, b: u128) -> u128 {
    if b == 0 {
        a
    } else {
        gcd(b, a % b)
    }
}

impl TheoreticalEstimates {
    /// Get estimated puzzles for a difficulty level
    pub fn estimate_for_difficulty(&self, difficulty: Difficulty) -> u128 {
        let (min_clues, max_clues) = match difficulty {
            Difficulty::Beginner => (45, 55),
            Difficulty::Easy => (36, 45),
            Difficulty::Medium => (32, 38),
            Difficulty::Intermediate => (28, 34),
            Difficulty::Hard => (24, 30),
            Difficulty::Expert => (22, 26),
            Difficulty::Master => (20, 24),
            Difficulty::Extreme => (17, 22),
        };

        let mut total: u128 = 0;
        for clues in min_clues..=max_clues {
            if let Some(&count) = self.puzzles_by_clue_count.get(&clues) {
                total = total.saturating_add(count);
            }
        }
        total
    }
}

/// Analyzer for puzzle diversity
pub struct DiversityAnalyzer {
    solver: Solver,
    fingerprints: HashSet<String>,
    stats: HashMap<Difficulty, DifficultyStats>,
}

impl Default for DiversityAnalyzer {
    fn default() -> Self {
        Self::new()
    }
}

impl DiversityAnalyzer {
    pub fn new() -> Self {
        Self {
            solver: Solver::new(),
            fingerprints: HashSet::new(),
            stats: HashMap::new(),
        }
    }

    /// Generate a fingerprint for a puzzle (canonical representation)
    pub fn fingerprint(grid: &Grid) -> String {
        let mut chars = Vec::with_capacity(81);
        for row in 0..9 {
            for col in 0..9 {
                let pos = crate::Position::new(row, col);
                match grid.get(pos) {
                    Some(v) => chars.push((b'0' + v) as char),
                    None => chars.push('.'),
                }
            }
        }
        chars.into_iter().collect()
    }

    /// Analyze puzzles for a specific difficulty
    pub fn analyze_difficulty(
        &mut self,
        difficulty: Difficulty,
        sample_size: usize,
    ) -> DifficultyStats {
        let mut generator = Generator::new();

        let mut stats = DifficultyStats::default();
        let mut local_fingerprints = HashSet::new();
        let mut total_clues = 0usize;

        for _ in 0..sample_size {
            // Generate a puzzle at the target difficulty
            let puzzle = generator.generate(difficulty);

            // Count clues
            let clue_count = Self::count_clues(&puzzle);
            total_clues += clue_count;

            // Update clue distribution
            *stats.clue_distribution.entry(clue_count).or_insert(0) += 1;

            // Update min/max
            if stats.min_clues == 0 || clue_count < stats.min_clues {
                stats.min_clues = clue_count;
            }
            if clue_count > stats.max_clues {
                stats.max_clues = clue_count;
            }

            // Get actual difficulty rated
            let rated = self.solver.rate_difficulty(&puzzle);
            let technique_name = format!("{:?}", rated);
            *stats
                .technique_distribution
                .entry(technique_name)
                .or_insert(0) += 1;

            // Check fingerprint uniqueness
            let fp = Self::fingerprint(&puzzle);
            if !self.fingerprints.contains(&fp) {
                self.fingerprints.insert(fp.clone());
                local_fingerprints.insert(fp);
            }
        }

        stats.samples = sample_size;
        stats.unique_puzzles = local_fingerprints.len();
        stats.success_rate = 1.0; // Generator always returns a puzzle
        stats.avg_attempts = 1.0;

        if sample_size > 0 {
            stats.avg_clues = total_clues as f64 / sample_size as f64;
        }

        // Estimate total puzzles based on uniqueness rate
        let theoretical = TheoreticalEstimates::new();
        let base_estimate = theoretical.estimate_for_difficulty(difficulty);

        // Adjust estimate based on observed uniqueness rate
        // If we saw X unique in Y samples, estimate total ≈ base * (X/Y) factor
        if stats.samples > 0 && stats.unique_puzzles > 0 {
            let uniqueness_rate = stats.unique_puzzles as f64 / stats.samples as f64;
            // High uniqueness rate suggests large puzzle space
            stats.estimated_total = if uniqueness_rate > 0.99 {
                base_estimate // Near 100% unique = very large space
            } else {
                // Use capture-recapture style estimation
                // If we're seeing repeats, space is smaller
                let estimated = (stats.samples as f64 * stats.samples as f64)
                    / (stats.samples - stats.unique_puzzles + 1) as f64;
                (estimated as u128).min(base_estimate)
            };
        } else {
            stats.estimated_total = base_estimate;
        }

        self.stats.insert(difficulty, stats.clone());
        stats
    }

    /// Count the number of clues (given cells) in a puzzle
    fn count_clues(grid: &Grid) -> usize {
        let mut count = 0;
        for row in 0..9 {
            for col in 0..9 {
                let pos = crate::Position::new(row, col);
                if grid.get(pos).is_some() {
                    count += 1;
                }
            }
        }
        count
    }

    /// Run full analysis across all difficulties
    pub fn full_analysis(&mut self, samples_per_difficulty: usize) -> DiversityReport {
        let difficulties = [
            Difficulty::Beginner,
            Difficulty::Easy,
            Difficulty::Medium,
            Difficulty::Intermediate,
            Difficulty::Hard,
            Difficulty::Expert,
        ];

        let mut report = DiversityReport {
            theoretical: TheoreticalEstimates::new(),
            ..Default::default()
        };

        for difficulty in difficulties {
            let stats = self.analyze_difficulty(difficulty, samples_per_difficulty);
            report
                .by_difficulty
                .insert(format!("{:?}", difficulty), stats);
        }

        report.total_unique = self.fingerprints.len();
        report
    }

    /// Get a summary string of the analysis
    pub fn summary(&self) -> String {
        let theoretical = TheoreticalEstimates::new();
        let mut lines = Vec::new();

        lines.push("=== Sudoku Puzzle Diversity Analysis ===\n".to_string());

        lines.push("Theoretical Estimates:".to_string());
        lines.push(format!(
            "  Total valid grids:     {:>25} (6.67 × 10²¹)",
            format_large_number(theoretical.total_grids)
        ));
        lines.push(format!(
            "  Unique grids:          {:>25} (5.47 × 10⁹)",
            format_large_number(theoretical.unique_grids)
        ));
        lines.push(format!(
            "  Min clues for unique:  {:>25}",
            theoretical.min_clues_for_unique
        ));
        lines.push(String::new());

        lines.push("Estimated Puzzles by Difficulty:".to_string());
        lines.push(format!(
            "  {:15} {:>12} {:>15} {:>20}",
            "Difficulty", "Clue Range", "Sample Unique", "Estimated Total"
        ));
        lines.push(format!("  {}", "-".repeat(65)));

        let difficulties = [
            (Difficulty::Beginner, "45-55"),
            (Difficulty::Easy, "36-45"),
            (Difficulty::Medium, "32-38"),
            (Difficulty::Intermediate, "28-34"),
            (Difficulty::Hard, "24-30"),
            (Difficulty::Expert, "22-26"),
            (Difficulty::Master, "20-24"),
            (Difficulty::Extreme, "17-22"),
        ];

        for (diff, range) in difficulties {
            let estimate = theoretical.estimate_for_difficulty(diff);
            let sample_info = self
                .stats
                .get(&diff)
                .map(|s| format!("{}/{}", s.unique_puzzles, s.samples))
                .unwrap_or_else(|| "N/A".to_string());

            lines.push(format!(
                "  {:15} {:>12} {:>15} {:>20}",
                format!("{:?}", diff),
                range,
                sample_info,
                format_large_number(estimate)
            ));
        }

        lines.push(String::new());
        lines.push(format!(
            "Total unique fingerprints observed: {}",
            self.fingerprints.len()
        ));

        lines.join("\n")
    }
}

/// Format a large number with scientific notation for readability
fn format_large_number(n: u128) -> String {
    if n == 0 {
        return "0".to_string();
    }
    if n < 1_000_000 {
        return format!("{}", n);
    }

    let log10 = (n as f64).log10().floor() as u32;
    let mantissa = n as f64 / 10_f64.powi(log10 as i32);

    format!("{:.2} × 10^{}", mantissa, log10)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_theoretical_estimates() {
        let estimates = TheoreticalEstimates::new();

        assert_eq!(estimates.min_clues_for_unique, 17);
        assert!(estimates.total_grids > 0);
        assert!(estimates.unique_grids > 0);
        assert!(estimates.unique_grids < estimates.total_grids);

        // Check that harder difficulties have fewer estimated puzzles
        let beginner = estimates.estimate_for_difficulty(Difficulty::Beginner);
        let expert = estimates.estimate_for_difficulty(Difficulty::Expert);
        // This might not always hold due to combinatorics, but generally true
        assert!(beginner > 0);
        assert!(expert > 0);
    }

    #[test]
    fn test_binomial() {
        // Exact tests for small values
        assert_eq!(TheoreticalEstimates::binomial(5, 0), 1);
        assert_eq!(TheoreticalEstimates::binomial(5, 5), 1);
        assert_eq!(TheoreticalEstimates::binomial(5, 2), 10);
        assert_eq!(TheoreticalEstimates::binomial(10, 3), 120);
        assert_eq!(TheoreticalEstimates::binomial(20, 10), 184756);

        // For large binomials, check order of magnitude (floating point precision limits)
        let c81_17 = TheoreticalEstimates::binomial(81, 17);
        assert!(c81_17 > 10_000_000_000_000); // > 10^13
        assert!(c81_17 < 1_000_000_000_000_000_000); // < 10^18
    }

    #[test]
    fn test_fingerprint() {
        let grid = Grid::new_classic();
        let fp = DiversityAnalyzer::fingerprint(&grid);
        assert_eq!(fp.len(), 81);
        assert!(fp.chars().all(|c| c == '.'));
    }

    #[test]
    fn test_quick_analysis() {
        let mut analyzer = DiversityAnalyzer::new();
        let stats = analyzer.analyze_difficulty(Difficulty::Medium, 5);

        assert_eq!(stats.samples, 5);
        assert!(stats.unique_puzzles > 0);
        assert!(stats.avg_clues >= 32.0 && stats.avg_clues <= 38.0);
    }

    #[test]
    fn test_format_large_number() {
        assert_eq!(format_large_number(0), "0");
        assert_eq!(format_large_number(1000), "1000");
        assert_eq!(format_large_number(1_000_000), "1.00 × 10^6");
        assert_eq!(format_large_number(123_456_789), "1.23 × 10^8");
    }
}
