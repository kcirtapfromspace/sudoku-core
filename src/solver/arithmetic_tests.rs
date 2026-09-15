use super::*;
use crate::{BitSet, Grid, Hint, HintType, Position, ProofCertificate, Solver};

const ALL: u16 = 0x1ff;

fn pos(cell: usize) -> Position {
    Position::new(cell / 9, cell % 9)
}

fn sector_contains(sector: usize, cell: usize) -> bool {
    match sector {
        0..=8 => cell / 9 == sector,
        9..=17 => cell % 9 == sector - 9,
        18..=26 => (cell / 27) * 3 + cell % 9 / 3 == sector - 18,
        _ => false,
    }
}

fn restrict_one(grid: &mut Grid, sector: usize, retained: &[usize]) {
    for cell in 0..81 {
        if sector_contains(sector, cell) && !retained.contains(&cell) {
            grid.cell_mut(pos(cell)).remove_candidate(1);
        }
    }
}

fn sector_terms(sectors: &[usize], weights: &[i8]) -> Vec<ArithmeticTerm> {
    assert_eq!(sectors.len(), weights.len());
    sectors
        .iter()
        .zip(weights)
        .map(|(&sector, &weight)| ArithmeticTerm {
            requirement: ArithmeticRequirement::SectorDigit { sector, digit: 1 },
            weight,
        })
        .collect()
}

fn guardian() -> Grid {
    let mut grid = Grid::new_classic();
    for (sector, cells) in [
        (0, vec![0, 3]),
        (12, vec![3, 30]),
        (3, vec![30, 28]),
        (10, vec![28, 10]),
        (18, vec![0, 10, 20]),
    ] {
        restrict_one(&mut grid, sector, &cells);
    }
    grid
}

fn xwing() -> Grid {
    let mut grid = Grid::new_classic();
    restrict_one(&mut grid, 0, &[0, 6]);
    restrict_one(&mut grid, 1, &[9, 15]);
    grid
}

fn weight_two_fixture() -> Grid {
    let mut grid = Grid::new_classic();
    let retained = [4, 3, 30, 28, 1, 12, 9];
    // Native source IDs [81,189,108,171,243,90] are respectively
    // digit 1 in R1,C4,R4,C2,B1,R2 (sector indices below).
    for sector in [0, 12, 3, 10, 18, 1] {
        restrict_one(&mut grid, sector, &retained);
    }
    grid
}

fn cell_obstruction() -> Grid {
    let mut grid = Grid::new_classic();
    let retained = [3, 30, 28, 1, 9, 12];
    for sector in [12, 3, 10, 18, 1] {
        restrict_one(&mut grid, sector, &retained);
    }
    grid.cell_mut(pos(3))
        .set_candidates(BitSet::from_slice(&[1, 2, 3]));
    grid
}

fn canonical_completion() -> [u8; 81] {
    std::array::from_fn(|cell| (1 + (cell / 9 * 3 + cell / 27 + cell % 9) % 9) as u8)
}

fn parse_completion(text: &str) -> [u8; 81] {
    assert_eq!(text.len(), 81);
    std::array::from_fn(|i| text.as_bytes()[i] - b'0')
}

fn domains(grid: &Grid) -> [u16; 81] {
    std::array::from_fn(|cell| match grid.get(pos(cell)) {
        Some(value) => 1 << (value - 1),
        None => grid
            .get_candidates(pos(cell))
            .iter()
            .fold(0, |mask, digit| mask | (1 << (digit - 1))),
    })
}

fn assert_completion(grid: &Grid, values: &[u8; 81]) {
    let allowed = domains(grid);
    for cell in 0..81 {
        assert!((1..=9).contains(&values[cell]));
        assert_ne!(allowed[cell] & (1 << (values[cell] - 1)), 0, "cell {cell}");
    }
    for sector in 0..27 {
        let mask = (0..81)
            .filter(|&cell| sector_contains(sector, cell))
            .fold(0, |mask, cell| mask | (1 << (values[cell] - 1)));
        assert_eq!(mask, ALL, "sector {sector}");
    }
}

// Independent exact-cover oracle. Each recursive decision chooses the smallest
// unsatisfied cell or sector/digit requirement. It never calls Solver::solve,
// recalculates candidates, or loses the caller's candidate deletions.
#[derive(Debug)]
enum OracleResult {
    Sat([u8; 81]),
    Unsat,
    Budget,
}

fn oracle(grid: &Grid, assumption: Option<(usize, u8, bool)>) -> OracleResult {
    let mut allowed = domains(grid);
    if let Some((cell, digit, value)) = assumption {
        if value {
            allowed[cell] &= 1 << (digit - 1);
        } else {
            allowed[cell] &= !(1 << (digit - 1));
        }
    }
    let mut remaining = 100_000;
    oracle_search(&allowed, &mut [0; 81], &mut remaining)
}

fn oracle_search(
    allowed: &[u16; 81],
    values: &mut [u8; 81],
    remaining: &mut usize,
) -> OracleResult {
    if *remaining == 0 {
        return OracleResult::Budget;
    }
    *remaining -= 1;
    let mut used = [0u16; 27];
    for (cell, &digit) in values.iter().enumerate() {
        if digit == 0 {
            continue;
        }
        let bit = 1 << (digit - 1);
        for sector in [cell / 9, 9 + cell % 9, 18 + (cell / 27) * 3 + cell % 9 / 3] {
            assert_eq!(used[sector] & bit, 0);
            used[sector] |= bit;
        }
    }
    let mut candidates = [0u16; 81];
    let mut best: Vec<(usize, u8)> = Vec::new();
    let mut best_len = 10;
    for cell in 0..81 {
        if values[cell] != 0 {
            continue;
        }
        candidates[cell] = allowed[cell]
            & !(used[cell / 9] | used[9 + cell % 9] | used[18 + (cell / 27) * 3 + cell % 9 / 3]);
        let count = candidates[cell].count_ones() as usize;
        if count == 0 {
            return OracleResult::Unsat;
        }
        if count < best_len {
            best = (1..=9)
                .filter(|digit| candidates[cell] & (1 << (digit - 1)) != 0)
                .map(|digit| (cell, digit))
                .collect();
            best_len = count;
        }
    }
    if best_len == 10 {
        return OracleResult::Sat(*values);
    }
    for (sector, &occupied) in used.iter().enumerate() {
        for digit in 1..=9 {
            let bit = 1 << (digit - 1);
            if occupied & bit != 0 {
                continue;
            }
            let choices: Vec<_> = (0..81)
                .filter(|&cell| sector_contains(sector, cell) && candidates[cell] & bit != 0)
                .map(|cell| (cell, digit))
                .collect();
            if choices.is_empty() {
                return OracleResult::Unsat;
            }
            if choices.len() < best_len {
                best_len = choices.len();
                best = choices;
            }
        }
    }
    for (cell, digit) in best {
        values[cell] = digit;
        let result = oracle_search(allowed, values, remaining);
        values[cell] = 0;
        match result {
            OracleResult::Unsat => {}
            other => return other,
        }
    }
    OracleResult::Unsat
}

fn assert_forced(grid: &Grid, cell: usize, digit: u8, value: bool) {
    match oracle(grid, Some((cell, digit, value))) {
        OracleResult::Sat(completion) => assert_completion(grid, &completion),
        other => panic!("expected a compatible full completion, got {other:?}"),
    }
    assert!(matches!(
        oracle(grid, Some((cell, digit, !value))),
        OracleResult::Unsat
    ));
}

fn hint_proof(hint: &Hint) -> &ArithmeticProof {
    match hint.proof.as_ref() {
        Some(ProofCertificate::Arithmetic(proof)) => proof,
        other => panic!("expected an arithmetic proof, got {other:?}"),
    }
}

#[test]
fn guardian_parity_certificate_and_opposite_are_independently_checked() {
    let grid = guardian();
    let proof = ArithmeticProof::from_terms(
        &grid,
        sector_terms(&[0, 12, 3, 10, 18], &[1; 5]),
        20,
        1,
        true,
    )
    .expect("five source Guardian placement");
    assert!(proof.verify(&grid));
    assert!(matches!(
        proof.check(&grid),
        Some(ArithmeticCheck::Divisibility {
            residual: 5,
            divisor: 2
        })
    ));
    assert_completion(
        &grid,
        &parse_completion(
            "234156789567489123891237456315742698678913245429568317153694872982371564746825931",
        ),
    );
    assert_forced(&grid, 20, 1, true);
}

#[test]
fn xwing_interval_certificate_and_opposite_are_independently_checked() {
    let grid = xwing();
    let proof = ArithmeticProof::from_terms(
        &grid,
        sector_terms(&[9, 15, 0, 1], &[1, 1, -1, -1]),
        27,
        1,
        false,
    )
    .expect("ordinary X-wing elimination");
    assert!(matches!(
        proof.check(&grid),
        Some(ArithmeticCheck::Interval { .. })
    ));
    assert_completion(&grid, &canonical_completion());
    assert_forced(&grid, 27, 1, false);
}

#[test]
fn six_source_weight_two_certificate_is_verified_against_full_masks() {
    let grid = weight_two_fixture();
    let proof = ArithmeticProof::from_terms(
        &grid,
        sector_terms(&[0, 12, 3, 10, 18, 1], &[2, -1, 1, -1, -1, 1]),
        4,
        1,
        false,
    )
    .expect("2t+a=1 must exclude Boolean t=1");
    assert!(proof.verify(&grid));
    assert_completion(
        &grid,
        &parse_completion(
            "234156789156798234789234156412567398593812467678349512921485673367921845845673921",
        ),
    );
    assert_forced(&grid, 4, 1, false);
}

#[test]
fn overlapping_fish_trap_rejects_an_elimination_with_a_real_completion() {
    let mut grid = Grid::new_classic();
    let solution = canonical_completion();
    for (cell, &digit) in solution.iter().enumerate() {
        if digit != 1 && cell != 1 && cell != 9 {
            grid.cell_mut(pos(cell)).remove_candidate(1);
        }
    }
    assert_eq!(solution[15], 1);
    assert_completion(&grid, &solution);
    assert!(ArithmeticProof::from_terms(
        &grid,
        sector_terms(&[1, 18, 0, 9], &[1, 1, -1, -1]),
        15,
        1,
        false,
    )
    .is_none());
}

#[test]
fn every_unit_combination_misses_the_fixed_six_source_cell_obstruction() {
    let grid = cell_obstruction();
    assert_forced(&grid, 3, 2, false);
    let mut sources: Vec<_> = [12, 3, 10, 18, 1]
        .into_iter()
        .map(|sector| ArithmeticRequirement::SectorDigit { sector, digit: 1 })
        .collect();
    sources.push(ArithmeticRequirement::Cell { cell: 3 });
    let mut checked = 0;
    for mut code in 1..729 {
        let mut weights = [0i8; 6];
        for weight in &mut weights {
            *weight = match code % 3 {
                0 => 0,
                1 => 1,
                _ => -1,
            };
            code /= 3;
        }
        if weights.iter().find(|&&weight| weight != 0) != Some(&1) {
            continue;
        }
        let terms = sources
            .iter()
            .zip(weights)
            .filter(|(_, weight)| *weight != 0)
            .map(|(requirement, weight)| ArithmeticTerm {
                requirement: requirement.clone(),
                weight,
            })
            .collect();
        assert!(ArithmeticProof::from_terms(&grid, terms, 3, 2, false).is_none());
        checked += 1;
    }
    assert_eq!(checked, 364);
}

#[test]
fn verifier_rejects_tampering_and_candidate_additions() {
    let grid = guardian();
    let proof = ArithmeticProof::from_terms(
        &grid,
        sector_terms(&[0, 12, 3, 10, 18], &[1; 5]),
        20,
        1,
        true,
    )
    .unwrap();
    let mut tampered = proof.clone();
    tampered.value = false;
    assert!(!tampered.verify(&grid));
    tampered = proof.clone();
    tampered.cell = 0;
    assert!(!tampered.verify(&grid));
    tampered = proof.clone();
    tampered.version = u8::MAX;
    assert!(!tampered.verify(&grid));
    tampered = proof.clone();
    tampered.state_hash.push('0');
    assert!(!tampered.verify(&grid));

    let mut changed = grid.clone();
    // This candidate is outside every proof source. The arithmetic remains
    // valid, but the proof must still be bound to the exact original state.
    changed.cell_mut(pos(80)).remove_candidate(9);
    let fresh = ArithmeticProof::from_terms(
        &changed,
        proof.terms.clone(),
        proof.cell,
        proof.digit,
        proof.value,
    )
    .unwrap();
    assert!(!proof.verify(&changed));
    changed.cell_mut(pos(80)).add_candidate(9);
    assert!(!fresh.verify(&changed));
    assert!(proof.verify(&changed));
}

#[test]
fn malformed_identifiers_and_weights_are_rejected_without_panicking() {
    let grid = guardian();
    let proof = ArithmeticProof::from_terms(
        &grid,
        sector_terms(&[0, 12, 3, 10, 18], &[1; 5]),
        20,
        1,
        true,
    )
    .unwrap();
    let invalid_requirements = [
        ArithmeticRequirement::Cell { cell: 81 },
        ArithmeticRequirement::Cell { cell: usize::MAX },
        ArithmeticRequirement::SectorDigit {
            sector: 27,
            digit: 1,
        },
        ArithmeticRequirement::SectorDigit {
            sector: usize::MAX,
            digit: 1,
        },
        ArithmeticRequirement::SectorDigit {
            sector: 0,
            digit: 0,
        },
        ArithmeticRequirement::SectorDigit {
            sector: 0,
            digit: 10,
        },
        ArithmeticRequirement::SectorDigit {
            sector: 0,
            digit: u8::MAX,
        },
    ];
    for requirement in invalid_requirements {
        let mut bad = proof.clone();
        bad.terms[0].requirement = requirement;
        assert!(!bad.verify(&grid));
    }
    for weight in [i8::MIN, i8::MAX, 0] {
        let mut bad = proof.clone();
        bad.terms[0].weight = weight;
        assert!(!bad.verify(&grid));
    }
    for (cell, digit) in [(81, 1), (usize::MAX, 1), (20, 0), (20, 10), (20, u8::MAX)] {
        let mut bad = proof.clone();
        bad.cell = cell;
        bad.digit = digit;
        assert!(!bad.verify(&grid));
    }
    let mut duplicate = proof.clone();
    duplicate.terms.push(duplicate.terms[0].clone());
    assert!(!duplicate.verify(&grid));
    assert!(ArithmeticProof::from_terms(&grid, Vec::new(), 20, 1, true).is_none());
}

#[test]
fn placed_values_are_normalized_to_singletons_in_source_equations() {
    let mut grid = Grid::new_classic();
    grid.set_cell_unchecked(pos(0), Some(1));
    assert!(grid.cell(pos(0)).candidates().is_empty());
    // Including the placed candidate gives x_(0,1)+sum(other row ones)=1.
    // Its cell equation then proves another row candidate false. Treating the
    // placed cell's stored empty candidates as an empty equation would fail.
    let terms = vec![
        ArithmeticTerm {
            requirement: ArithmeticRequirement::SectorDigit {
                sector: 0,
                digit: 1,
            },
            weight: 1,
        },
        ArithmeticTerm {
            requirement: ArithmeticRequirement::Cell { cell: 0 },
            weight: -1,
        },
    ];
    let proof = ArithmeticProof::from_terms(&grid, terms, 1, 1, false)
        .expect("placed values must contribute their singleton to native equations");
    assert!(proof.verify(&grid));
    assert_forced(&grid, 1, 1, false);
}

#[test]
fn bounded_search_is_deterministic_preserves_masks_and_returns_sound_proofs() {
    let grid = guardian();
    let before = serde_json::to_value(&grid).unwrap();
    let options = ArithmeticSearchOptions::default();
    let solver = Solver::new();
    let first = solver.search_arithmetic(&grid, &options);
    let second = solver.search_arithmetic(&grid, &options);
    assert!(first.error.is_none(), "{:?}", first.error);
    assert!(second.error.is_none(), "{:?}", second.error);
    assert!(first.tested_combinations <= options.max_combinations);
    assert_eq!(first.tested_combinations, second.tested_combinations);
    assert_eq!(first.budget_exhausted, second.budget_exhausted);
    assert_eq!(first.beam_pruned, second.beam_pruned);
    let first_hint = first
        .hint
        .expect("default search must find a deduction in the Guardian fixture");
    let second_hint = second.hint.expect("repeat search must agree");
    let proof = hint_proof(&first_hint);
    assert!(proof.verify(&grid));
    assert_eq!(
        serde_json::to_value(proof).unwrap(),
        serde_json::to_value(hint_proof(&second_hint)).unwrap()
    );
    assert_forced(&grid, proof.cell, proof.digit, proof.value);
    match &first_hint.hint_type {
        HintType::SetValue {
            pos: position,
            value,
        } => {
            assert!(proof.value);
            assert_eq!(*position, pos(proof.cell));
            assert_eq!(*value, proof.digit);
        }
        HintType::EliminateCandidates {
            pos: position,
            values,
        } => {
            assert!(!proof.value);
            assert_eq!(*position, pos(proof.cell));
            assert_eq!(values, &vec![proof.digit]);
        }
    }
    let convenience = solver
        .get_arithmetic_hint(&grid)
        .expect("same default search");
    assert!(hint_proof(&convenience).verify(&grid));
    assert_eq!(serde_json::to_value(&grid).unwrap(), before);
}

#[test]
fn zero_search_budget_reports_no_deduction_and_no_completion_claim() {
    let grid = guardian();
    let options = ArithmeticSearchOptions {
        max_combinations: 0,
        ..Default::default()
    };
    let result = Solver::new().search_arithmetic(&grid, &options);
    assert!(result.hint.is_none());
    assert_eq!(result.tested_combinations, 0);
    assert!(result.budget_exhausted);
}

fn singleton_grid() -> Grid {
    let mut grid = Grid::new_classic();
    grid.cell_mut(pos(80))
        .set_candidates(BitSet::from_slice(&[9]));
    grid
}

fn singleton_terms() -> Vec<ArithmeticTerm> {
    vec![ArithmeticTerm {
        requirement: ArithmeticRequirement::Cell { cell: 80 },
        weight: 1,
    }]
}

fn assert_arithmetic_rejects_grid(grid: &Grid, expected_error: &str) {
    let before = serde_json::to_value(grid).unwrap();
    let solver = Solver::new();
    let result = solver.search_arithmetic(grid, &ArithmeticSearchOptions::default());
    assert!(result.hint.is_none());
    assert_eq!(result.tested_combinations, 0);
    assert!(!result.budget_exhausted);
    assert!(!result.beam_pruned);
    assert!(
        result
            .error
            .as_deref()
            .is_some_and(|error| error.contains(expected_error)),
        "expected {expected_error:?}, got {:?}",
        result.error
    );
    // The convenience API must validate before constructing other solver data,
    // including when a public BitSet contains digits outside the Sudoku range.
    assert!(solver.get_arithmetic_hint(grid).is_none());
    assert!(ArithmeticProof::from_terms(grid, singleton_terms(), 80, 9, true).is_none());
    assert_eq!(serde_json::to_value(grid).unwrap(), before);
}

#[test]
fn invalid_search_options_are_rejected_before_any_evaluation() {
    let grid = singleton_grid();
    let mut invalid = Vec::new();
    for max_sources in [0, 7, usize::MAX] {
        invalid.push(ArithmeticSearchOptions {
            max_sources,
            ..Default::default()
        });
    }
    for max_weight in [i8::MIN, -1, 0, 3, i8::MAX] {
        invalid.push(ArithmeticSearchOptions {
            max_weight,
            ..Default::default()
        });
    }
    for beam_width in [0, 513, usize::MAX] {
        invalid.push(ArithmeticSearchOptions {
            beam_width,
            ..Default::default()
        });
    }
    for max_combinations in [1_000_001, usize::MAX] {
        invalid.push(ArithmeticSearchOptions {
            max_combinations,
            ..Default::default()
        });
    }
    for options in invalid {
        let result = Solver::new().search_arithmetic(&grid, &options);
        assert!(result.error.is_some(), "options {options:?}");
        assert!(result.hint.is_none());
        assert_eq!(result.tested_combinations, 0);
        assert!(!result.budget_exhausted);
        assert!(!result.beam_pruned);
    }
}

#[test]
fn valid_search_option_endpoints_accept_an_immediate_verified_deduction() {
    let grid = singleton_grid();
    for options in [
        ArithmeticSearchOptions {
            max_sources: 1,
            max_weight: 1,
            beam_width: 1,
            max_combinations: 1,
        },
        ArithmeticSearchOptions {
            max_sources: 6,
            max_weight: 2,
            beam_width: 512,
            max_combinations: 1_000_000,
        },
    ] {
        let result = Solver::new().search_arithmetic(&grid, &options);
        assert!(result.error.is_none(), "{:?}", result.error);
        assert_eq!(result.tested_combinations, 1);
        assert!(hint_proof(&result.hint.expect("singleton deduction")).verify(&grid));
    }
}

#[test]
fn high_beam_search_without_deductions_stops_at_the_evaluation_budget() {
    let grid = Grid::new_classic();
    let before = serde_json::to_value(&grid).unwrap();
    let options = ArithmeticSearchOptions {
        max_sources: 2,
        max_weight: 2,
        beam_width: 512,
        max_combinations: 1_000,
    };
    let result = Solver::new().search_arithmetic(&grid, &options);
    assert!(result.error.is_none(), "{:?}", result.error);
    assert!(result.hint.is_none());
    assert_eq!(result.tested_combinations, options.max_combinations);
    assert!(result.budget_exhausted);
    assert!(result.beam_pruned);
    assert_eq!(serde_json::to_value(&grid).unwrap(), before);
}

#[test]
fn custom_and_unrestored_grids_cannot_reuse_classic_arithmetic_proofs() {
    struct NamedFake(&'static str);
    impl crate::Constraint for NamedFake {
        fn validate(&self, _: &[[Option<u8>; 9]; 9], _: Position, _: u8) -> bool {
            true
        }

        fn affected_cells(&self, _: Position) -> Vec<Position> {
            Vec::new()
        }

        fn name(&self) -> &'static str {
            self.0
        }
    }

    let original = singleton_grid();
    let proof = ArithmeticProof::from_terms(&original, singleton_terms(), 80, 9, true).unwrap();
    let mut raw: Grid = serde_json::from_value(serde_json::to_value(&original).unwrap()).unwrap();
    assert_eq!(raw.variant(), original.variant());
    assert_arithmetic_rejects_grid(&raw, "installed standard Sudoku constraints");
    assert!(!proof.verify(&raw));
    assert!(!proof.verify(&raw.clone()));
    assert!(!proof.verify(&raw.deep_clone()));
    raw.restore_constraints();
    assert!(proof.verify(&raw));
    assert!(hint_proof(&Solver::new().get_arithmetic_hint(&raw).unwrap()).verify(&raw));

    for constraints in [
        Vec::<crate::ConstraintBox>::new(),
        vec![
            Box::new(NamedFake("Row")) as crate::ConstraintBox,
            Box::new(NamedFake("Column")),
            Box::new(NamedFake("Box")),
        ],
        vec![
            Box::new(crate::RowConstraint) as crate::ConstraintBox,
            Box::new(crate::ColumnConstraint),
            Box::new(crate::BoxConstraint),
        ],
    ] {
        let mut custom = Grid::new_with_constraints(constraints);
        custom
            .cell_mut(pos(80))
            .set_candidates(BitSet::from_slice(&[9]));
        assert_eq!(custom.variant(), original.variant());
        assert_arithmetic_rejects_grid(&custom, "installed standard Sudoku constraints");
        assert!(!proof.verify(&custom));
        assert!(!proof.verify(&custom.clone()));
        assert!(!proof.verify(&custom.deep_clone()));
        // Restoration installs the real built-in rules, even for a custom grid.
        custom.restore_constraints();
        assert!(proof.verify(&custom));
        assert!(hint_proof(&Solver::new().get_arithmetic_hint(&custom).unwrap()).verify(&custom));
    }
}

#[test]
fn built_in_variant_constraints_accept_standard_arithmetic_consequences() {
    for mut grid in [
        Grid::new_classic(),
        Grid::new_x_sudoku(),
        Grid::new_killer(vec![crate::KillerCageConstraint::new(vec![pos(80)], 9)]),
    ] {
        grid.cell_mut(pos(80))
            .set_candidates(BitSet::from_slice(&[9]));
        let proof = ArithmeticProof::from_terms(&grid, singleton_terms(), 80, 9, true).unwrap();
        assert!(proof.verify(&grid));
        let mut raw: Grid = serde_json::from_value(serde_json::to_value(&grid).unwrap()).unwrap();
        assert_arithmetic_rejects_grid(&raw, "installed standard Sudoku constraints");
        raw.restore_constraints();
        assert!(proof.verify(&raw));
        assert!(hint_proof(&Solver::new().get_arithmetic_hint(&raw).unwrap()).verify(&raw));
    }
}

#[test]
fn invalid_local_values_and_candidate_masks_are_rejected_without_panicking() {
    for raw in [0, 1, 1 << 10, 1 << 15, u16::MAX] {
        let mut grid = singleton_grid();
        grid.cell_mut(pos(0)).set_candidates(BitSet::from_raw(raw));
        assert_arithmetic_rejects_grid(&grid, "empty or invalid candidate mask");
    }
    for value in [0, 10, u8::MAX] {
        let mut grid = singleton_grid();
        grid.cell_mut(pos(0)).set_value(Some(value));
        assert_arithmetic_rejects_grid(&grid, "Placed digit is outside 1..9");
    }
}

#[test]
fn duplicate_placed_values_and_absent_sector_digits_are_rejected() {
    // Exercise a row, a column, and a box with no common row or column.
    for second in [1, 9, 10] {
        let mut grid = singleton_grid();
        grid.cell_mut(pos(0)).set_value(Some(7));
        grid.cell_mut(pos(second)).set_value(Some(7));
        assert_arithmetic_rejects_grid(&grid, "Duplicate placed digit");
    }
    for sector in [0, 9, 18] {
        let mut grid = singleton_grid();
        restrict_one(&mut grid, sector, &[]);
        assert_arithmetic_rejects_grid(&grid, "no remaining position for a digit");
    }
}
