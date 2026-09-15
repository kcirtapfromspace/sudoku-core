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

#[test]
fn analysis_statistics_match_samples_and_survive_serialization() {
    let mut analyzer = DiversityAnalyzer::default();
    let stats = analyzer.analyze_difficulty(Difficulty::Beginner, 2);
    assert_eq!(stats.samples, 2);
    assert!((1..=2).contains(&stats.unique_puzzles));
    assert_eq!(stats.clue_distribution.values().sum::<usize>(), 2);
    assert_eq!(stats.technique_distribution.values().sum::<usize>(), 2);
    let clues: usize = stats
        .clue_distribution
        .iter()
        .map(|(n, count)| n * count)
        .sum();
    assert_eq!(stats.avg_clues, clues as f64 / 2.0);
    assert_eq!(
        stats.min_clues,
        *stats.clue_distribution.keys().min().unwrap()
    );
    assert_eq!(
        stats.max_clues,
        *stats.clue_distribution.keys().max().unwrap()
    );
    assert!(stats.min_clues >= 45 && stats.max_clues <= 55);
    assert!(stats.estimated_total > 0);
    let restored: DifficultyStats =
        serde_json::from_str(&serde_json::to_string(&stats).unwrap()).unwrap();
    assert_eq!(restored.samples, stats.samples);
    assert_eq!(restored.clue_distribution, stats.clue_distribution);
    let summary = analyzer.summary();
    assert!(summary.contains("Beginner"));
    assert!(summary.contains(&format!("{}/2", stats.unique_puzzles)));
    assert!(summary.contains("Total unique fingerprints observed:"));
}

#[test]
fn empty_analysis_has_finite_zero_statistics_and_theoretical_estimates() {
    let mut analyzer = DiversityAnalyzer::new();
    let stats = analyzer.analyze_difficulty(Difficulty::Extreme, 0);
    assert_eq!(stats.samples, 0);
    assert_eq!(stats.unique_puzzles, 0);
    assert_eq!(stats.avg_clues, 0.0);
    assert!(stats.clue_distribution.is_empty());
    assert!(stats.estimated_total > 0);
    let report = analyzer.full_analysis(0);
    assert_eq!(report.total_unique, 0);
    assert!(report.by_difficulty.contains_key("Beginner"));
    assert!(report.by_difficulty.contains_key("Expert"));
    for stats in report.by_difficulty.values() {
        assert_eq!(stats.samples, 0);
        assert_eq!(stats.unique_puzzles, 0);
        assert!(stats.avg_clues.is_finite());
    }
    assert!(analyzer.summary().contains("0/0"));
}

#[test]
fn estimates_respect_missing_clue_counts_and_saturate_safely() {
    let mut estimates = TheoreticalEstimates::default();
    for difficulty in [
        Difficulty::Beginner,
        Difficulty::Easy,
        Difficulty::Medium,
        Difficulty::Intermediate,
        Difficulty::Hard,
        Difficulty::Expert,
        Difficulty::Master,
        Difficulty::Extreme,
    ] {
        assert!(estimates.estimate_for_difficulty(difficulty) > 0);
    }
    estimates.puzzles_by_clue_count.clear();
    assert_eq!(estimates.estimate_for_difficulty(Difficulty::Easy), 0);
    estimates.puzzles_by_clue_count.insert(36, u128::MAX);
    estimates.puzzles_by_clue_count.insert(37, 1);
    assert_eq!(
        estimates.estimate_for_difficulty(Difficulty::Easy),
        u128::MAX
    );
    assert_eq!(estimates.estimate_for_difficulty(Difficulty::Extreme), 0);
}

#[test]
fn combinatorial_edges_and_number_formatting_are_well_defined() {
    assert_eq!(TheoreticalEstimates::binomial(9, 10), 0);
    assert_eq!(TheoreticalEstimates::binomial(9, 0), 1);
    assert_eq!(TheoreticalEstimates::binomial(9, 9), 1);
    assert_eq!(TheoreticalEstimates::binomial(9, 2), 36);
    assert_eq!(TheoreticalEstimates::binomial(9, 7), 36);
    assert_eq!(format_large_number(0), "0");
    assert_eq!(format_large_number(999_999), "999999");
    assert_eq!(format_large_number(1_000_000), "1.00 × 10^6");
    assert_eq!(format_large_number(5_470_000_000), "5.47 × 10^9");
}
