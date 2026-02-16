//! Unified Fish Engine: ONE `search_fish()` parameterized by SectorConstraint.
//!
//! Replaces: find_x_wing, find_swordfish, find_jellyfish, find_finned_fish_generic,
//! find_franken_fish, find_siamese_fish, find_mutant_fish.
//!
//! Core idea: For each digit, enumerate base-sector combos of size n, find cover-sector
//! combos of size n. Fins = base_cells \ cover_cells. Classification by
//! (size, has_fins, sector_types) -> Technique variant.

use super::explain::{ExplanationData, Finding, InferenceResult, ProofCertificate};
use super::fabric::{
    idx_to_pos, sector_cells, CandidateFabric, SECTOR_BOX_BASE, SECTOR_COL_BASE, SECTOR_ROW_BASE,
};
use super::types::Technique;

/// Controls which sector types are allowed in base/cover sets.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum SectorConstraint {
    /// Rows as bases, cols as covers (or vice versa)
    Basic,
    /// Basic + boxes allowed in cover sets (Franken)
    Franken,
    /// All three sector types in both base and cover (Mutant)
    Mutant,
}

fn sector_name(sector: usize) -> String {
    if sector < 9 {
        format!("r{}", sector + 1)
    } else if sector < 18 {
        format!("c{}", sector - 9 + 1)
    } else {
        format!("b{}", sector - 18 + 1)
    }
}

fn sector_type(sector: usize) -> u8 {
    if sector < 9 {
        0
    }
    // row
    else if sector < 18 {
        1
    }
    // col
    else {
        2
    } // box
}

/// Get a bitmask of cells (linear indices) that are in a given sector and have a candidate.
fn sector_candidate_mask(fab: &CandidateFabric, sector: usize, digit: u8) -> u128 {
    let cells = sector_cells(sector);
    let mut mask = 0u128;
    for &c in &cells {
        if fab.values[c].is_none() && fab.cell_cands[c].contains(digit) {
            mask |= 1u128 << c;
        }
    }
    mask
}

/// Search for fish patterns of a given size, sector constraint, and digit.
fn search_fish_for_digit(
    fab: &CandidateFabric,
    digit: u8,
    size: usize,
    constraint: SectorConstraint,
) -> Option<Finding> {
    // Determine which sectors can be bases and covers
    let (base_sectors, cover_sectors) = match constraint {
        SectorConstraint::Basic => {
            // Row bases -> col covers, then col bases -> row covers
            let rows: Vec<usize> = (SECTOR_ROW_BASE..SECTOR_ROW_BASE + 9).collect();
            let cols: Vec<usize> = (SECTOR_COL_BASE..SECTOR_COL_BASE + 9).collect();
            // Try both orientations
            return search_fish_oriented(fab, digit, size, &rows, &cols, constraint)
                .or_else(|| search_fish_oriented(fab, digit, size, &cols, &rows, constraint));
        }
        SectorConstraint::Franken => {
            // Lines as bases, lines+boxes as covers
            let lines: Vec<usize> = (0..18).collect();
            let all: Vec<usize> = (0..27).collect();
            return search_fish_oriented(fab, digit, size, &lines, &all, constraint);
        }
        SectorConstraint::Mutant => {
            let all: Vec<usize> = (0..27).collect();
            (all.clone(), all)
        }
    };

    search_fish_oriented(fab, digit, size, &base_sectors, &cover_sectors, constraint)
}

fn search_fish_oriented(
    fab: &CandidateFabric,
    digit: u8,
    size: usize,
    base_pool: &[usize],
    cover_pool: &[usize],
    constraint: SectorConstraint,
) -> Option<Finding> {
    // Filter base sectors: must have at least 2 candidates for this digit
    let eligible_bases: Vec<usize> = base_pool
        .iter()
        .filter(|&&s| fab.sector_cand_count(s, digit) >= 2)
        .copied()
        .collect();

    if eligible_bases.len() < size {
        return None;
    }

    // Enumerate base combos
    for base_combo in combinations(&eligible_bases, size) {
        let mut base_cells = 0u128;
        for &bs in &base_combo {
            base_cells |= sector_candidate_mask(fab, bs, digit);
        }

        if base_cells.count_ones() < size as u32 {
            continue;
        }

        // Find cover sectors that overlap with base cells
        let eligible_covers: Vec<usize> = cover_pool
            .iter()
            .filter(|&&s| {
                !base_combo.contains(&s) && (sector_candidate_mask(fab, s, digit) & base_cells) != 0
            })
            .copied()
            .collect();

        if eligible_covers.len() < size {
            continue;
        }

        // For mutant/franken: ensure base and cover don't share same sector type pattern
        // For basic: base and cover must be different types
        for cover_combo in combinations(&eligible_covers, size) {
            // Validate sector type constraints
            if !validate_fish_types(&base_combo, &cover_combo, constraint) {
                continue;
            }

            let mut cover_cells = 0u128;
            for &cs in &cover_combo {
                cover_cells |= sector_candidate_mask(fab, cs, digit);
            }

            let fins = base_cells & !cover_cells;
            let eliminations = cover_cells & !base_cells;

            if fins.count_ones() == 0 {
                // Basic fish (no fins) - eliminate from cover_cells \ base_cells
                if let Some(f) = make_fish_finding(
                    fab,
                    digit,
                    size,
                    &base_combo,
                    &cover_combo,
                    &[],
                    eliminations,
                    constraint,
                ) {
                    return Some(f);
                }
            } else {
                // Finned fish: all fins must share one box
                let fin_cells: Vec<usize> = (0..81).filter(|&i| fins & (1u128 << i) != 0).collect();
                let fin_box = idx_to_pos(fin_cells[0]).box_index();
                if fin_cells
                    .iter()
                    .all(|&c| idx_to_pos(c).box_index() == fin_box)
                {
                    // Eliminations restricted to cells in cover that are also in fin box
                    let fin_box_mask = sector_candidate_mask(fab, SECTOR_BOX_BASE + fin_box, digit);
                    let restricted_elims = eliminations & fin_box_mask;
                    if let Some(f) = make_fish_finding(
                        fab,
                        digit,
                        size,
                        &base_combo,
                        &cover_combo,
                        &fin_cells,
                        restricted_elims,
                        constraint,
                    ) {
                        return Some(f);
                    }
                }
            }
        }
    }
    None
}

fn validate_fish_types(bases: &[usize], covers: &[usize], constraint: SectorConstraint) -> bool {
    match constraint {
        SectorConstraint::Basic => {
            // All bases same type, all covers different type
            let bt = sector_type(bases[0]);
            let ct = sector_type(covers[0]);
            bt != ct
                && bt != 2 && ct != 2  // No boxes in basic
                && bases.iter().all(|&s| sector_type(s) == bt)
                && covers.iter().all(|&s| sector_type(s) == ct)
        }
        SectorConstraint::Franken => {
            // Bases are lines (row or col), covers can include boxes
            // At least one base and one cover must be different types
            let has_box_cover = covers.iter().any(|&s| sector_type(s) == 2);
            let bases_are_lines = bases.iter().all(|&s| sector_type(s) != 2);
            bases_are_lines && has_box_cover
        }
        SectorConstraint::Mutant => {
            // At least 3 different sector types across all bases+covers
            let mut types = std::collections::HashSet::new();
            for &s in bases.iter().chain(covers.iter()) {
                types.insert(sector_type(s));
            }
            types.len() >= 3
                || (types.len() >= 2
                    && bases
                        .iter()
                        .chain(covers.iter())
                        .any(|&s| sector_type(s) == 2))
        }
    }
}

fn classify_fish(size: usize, has_fins: bool, constraint: SectorConstraint) -> Technique {
    match (constraint, size, has_fins) {
        (SectorConstraint::Basic, 2, false) => Technique::XWing,
        (SectorConstraint::Basic, 2, true) => Technique::FinnedXWing,
        (SectorConstraint::Basic, 3, false) => Technique::Swordfish,
        (SectorConstraint::Basic, 3, true) => Technique::FinnedSwordfish,
        (SectorConstraint::Basic, 4, false) => Technique::Jellyfish,
        (SectorConstraint::Basic, 4, true) => Technique::FinnedJellyfish,
        (SectorConstraint::Franken, _, _) => Technique::FrankenFish,
        (SectorConstraint::Mutant, _, _) => Technique::MutantFish,
        _ => Technique::Jellyfish, // fallback
    }
}

#[allow(clippy::too_many_arguments)]
fn make_fish_finding(
    fab: &CandidateFabric,
    digit: u8,
    size: usize,
    bases: &[usize],
    covers: &[usize],
    fins: &[usize],
    elim_mask: u128,
    constraint: SectorConstraint,
) -> Option<Finding> {
    if elim_mask == 0 {
        return None;
    }

    // Find the first elimination cell
    for cell in 0..81 {
        if elim_mask & (1u128 << cell) != 0
            && fab.values[cell].is_none()
            && fab.cell_cands[cell].contains(digit)
        {
            // Don't eliminate from cells that are in base sectors
            let in_base = bases.iter().any(|&s| {
                let sc = sector_cells(s);
                sc.contains(&cell)
            });
            if in_base {
                continue;
            }

            let technique = classify_fish(size, !fins.is_empty(), constraint);
            let variant: String = match constraint {
                SectorConstraint::Basic => {
                    if fins.is_empty() {
                        "Basic".into()
                    } else {
                        "Finned".into()
                    }
                }
                SectorConstraint::Franken => "Franken".into(),
                SectorConstraint::Mutant => "Mutant".into(),
            };

            let base_names: Vec<String> = bases.iter().map(|&s| sector_name(s)).collect();
            let cover_names: Vec<String> = covers.iter().map(|&s| sector_name(s)).collect();

            // Collect involved cells (all candidates for this digit in base sectors)
            let mut involved = Vec::new();
            for &bs in bases {
                for &c in &sector_cells(bs) {
                    if fab.values[c].is_none()
                        && fab.cell_cands[c].contains(digit)
                        && !involved.contains(&c)
                    {
                        involved.push(c);
                    }
                }
            }

            return Some(Finding {
                technique,
                inference: InferenceResult::Elimination {
                    cell,
                    values: vec![digit],
                },
                involved_cells: involved,
                explanation: ExplanationData::Fish {
                    size,
                    digit,
                    base_sectors: base_names,
                    cover_sectors: cover_names,
                    fins: fins.to_vec(),
                    variant,
                },
                proof: Some(ProofCertificate::Fish {
                    digit,
                    base_sectors: bases.to_vec(),
                    cover_sectors: covers.to_vec(),
                    fins: fins.to_vec(),
                }),
            });
        }
    }
    None
}

// ==================== Siamese Fish ====================

/// Siamese Fish: two overlapping fish patterns that share fins.
pub fn find_siamese_fish(fab: &CandidateFabric) -> Option<Finding> {
    // Siamese fish are detected as finned fish where two separate base combos
    // produce fin cells in the same box. We search for finned fish with an
    // additional constraint that there's a second valid fish pattern sharing the fin.
    for digit in 1..=9u8 {
        for size in 2..=4usize {
            // Get all row-based fish patterns with fins
            let rows: Vec<usize> = (SECTOR_ROW_BASE..SECTOR_ROW_BASE + 9)
                .filter(|&s| fab.sector_cand_count(s, digit) >= 2)
                .collect();
            let cols: Vec<usize> = (SECTOR_COL_BASE..SECTOR_COL_BASE + 9).collect();

            if rows.len() < size + 1 {
                continue;
            }

            // Find pairs of base combos that produce fins in the same box
            let all_combos = combinations(&rows, size);
            for i in 0..all_combos.len() {
                for j in (i + 1)..all_combos.len() {
                    let combo_a = &all_combos[i];
                    let combo_b = &all_combos[j];

                    // Check if these produce valid finned fish with shared fin box
                    if let Some(finding) =
                        check_siamese_pair(fab, digit, size, combo_a, combo_b, &cols)
                    {
                        return Some(finding);
                    }
                }
            }

            // Column-based
            let eligible_cols: Vec<usize> = (SECTOR_COL_BASE..SECTOR_COL_BASE + 9)
                .filter(|&s| fab.sector_cand_count(s, digit) >= 2)
                .collect();
            let row_covers: Vec<usize> = (SECTOR_ROW_BASE..SECTOR_ROW_BASE + 9).collect();

            if eligible_cols.len() < size + 1 {
                continue;
            }

            let all_col_combos = combinations(&eligible_cols, size);
            for i in 0..all_col_combos.len() {
                for j in (i + 1)..all_col_combos.len() {
                    if let Some(finding) = check_siamese_pair(
                        fab,
                        digit,
                        size,
                        &all_col_combos[i],
                        &all_col_combos[j],
                        &row_covers,
                    ) {
                        return Some(finding);
                    }
                }
            }
        }
    }
    None
}

fn check_siamese_pair(
    fab: &CandidateFabric,
    digit: u8,
    size: usize,
    combo_a: &[usize],
    combo_b: &[usize],
    cover_pool: &[usize],
) -> Option<Finding> {
    // Both combos must share at least one base sector (overlapping)
    let shared = combo_a.iter().filter(|s| combo_b.contains(s)).count();
    if shared == 0 || shared == size {
        return None;
    }

    for cover_combo in combinations(cover_pool, size) {
        let mut cover_cells_a = 0u128;
        let mut cover_cells_b = 0u128;
        let mut base_cells_a = 0u128;
        let mut base_cells_b = 0u128;

        for &cs in &cover_combo {
            let mask = sector_candidate_mask(fab, cs, digit);
            cover_cells_a |= mask;
            cover_cells_b |= mask;
        }
        for &bs in combo_a {
            base_cells_a |= sector_candidate_mask(fab, bs, digit);
        }
        for &bs in combo_b {
            base_cells_b |= sector_candidate_mask(fab, bs, digit);
        }

        let fins_a = base_cells_a & !cover_cells_a;
        let fins_b = base_cells_b & !cover_cells_b;

        if fins_a == 0 || fins_b == 0 {
            continue;
        }

        // Fins must be in the same box
        let fin_cells_a: Vec<usize> = (0..81).filter(|&i| fins_a & (1u128 << i) != 0).collect();
        let fin_cells_b: Vec<usize> = (0..81).filter(|&i| fins_b & (1u128 << i) != 0).collect();

        let box_a = idx_to_pos(fin_cells_a[0]).box_index();
        let box_b = idx_to_pos(fin_cells_b[0]).box_index();
        if box_a != box_b {
            continue;
        }
        if !fin_cells_a
            .iter()
            .all(|&c| idx_to_pos(c).box_index() == box_a)
        {
            continue;
        }
        if !fin_cells_b
            .iter()
            .all(|&c| idx_to_pos(c).box_index() == box_a)
        {
            continue;
        }

        // Combined eliminations: restricted to fin box
        let fin_box_mask = sector_candidate_mask(fab, SECTOR_BOX_BASE + box_a, digit);
        let elim_a = (cover_cells_a & !base_cells_a) & fin_box_mask;
        let elim_b = (cover_cells_b & !base_cells_b) & fin_box_mask;
        let combined = elim_a & elim_b;

        for cell in 0..81 {
            if combined & (1u128 << cell) != 0
                && fab.values[cell].is_none()
                && fab.cell_cands[cell].contains(digit)
            {
                let mut involved = Vec::new();
                for &bs in combo_a.iter().chain(combo_b.iter()) {
                    for &c in &sector_cells(bs) {
                        if fab.values[c].is_none()
                            && fab.cell_cands[c].contains(digit)
                            && !involved.contains(&c)
                        {
                            involved.push(c);
                        }
                    }
                }
                return Some(Finding {
                    technique: Technique::SiameseFish,
                    inference: InferenceResult::Elimination {
                        cell,
                        values: vec![digit],
                    },
                    involved_cells: involved,
                    explanation: ExplanationData::Fish {
                        size,
                        digit,
                        base_sectors: combo_a
                            .iter()
                            .chain(combo_b.iter())
                            .map(|&s| sector_name(s))
                            .collect(),
                        cover_sectors: cover_combo.iter().map(|&s| sector_name(s)).collect(),
                        fins: fin_cells_a
                            .iter()
                            .chain(fin_cells_b.iter())
                            .copied()
                            .collect(),
                        variant: "Siamese".into(),
                    },
                    proof: Some(ProofCertificate::Fish {
                        digit,
                        base_sectors: combo_a.iter().chain(combo_b.iter()).copied().collect(),
                        cover_sectors: cover_combo.clone(),
                        fins: fin_cells_a
                            .iter()
                            .chain(fin_cells_b.iter())
                            .copied()
                            .collect(),
                    }),
                });
            }
        }
    }
    None
}

// ==================== Size-1 fish (intersections) ====================

/// Long-form sector name for intersection explanations.
fn sector_name_long(sector: usize) -> String {
    if sector < 9 {
        format!("row {}", sector + 1)
    } else if sector < 18 {
        format!("column {}", sector - 9 + 1)
    } else {
        format!("box {}", sector - 18 + 1)
    }
}

/// Size-1 fish: pointing pair (box → line intersection).
///
/// When all candidates for a digit within a box are confined to a single
/// row or column, that digit can be eliminated from the rest of that row/col.
pub fn find_pointing_pair(fab: &CandidateFabric) -> Option<Finding> {
    for digit in 1..=9u8 {
        for box_idx in 0..9 {
            let base = SECTOR_BOX_BASE + box_idx;
            let base_mask = sector_candidate_mask(fab, base, digit);
            let base_count = base_mask.count_ones();
            if !(2..=3).contains(&base_count) {
                continue;
            }

            // Try each row and column as the cover
            for cover in 0..18 {
                let cover_mask = sector_candidate_mask(fab, cover, digit);
                // All base cells must lie within the cover sector
                if base_mask & cover_mask != base_mask {
                    continue;
                }
                // Eliminations: candidates in cover sector outside the base
                let elim_mask = cover_mask & !base_mask;
                if elim_mask == 0 {
                    continue;
                }
                for cell in 0..81usize {
                    if elim_mask & (1u128 << cell) != 0 {
                        let involved: Vec<usize> =
                            (0..81).filter(|&c| base_mask & (1u128 << c) != 0).collect();
                        return Some(Finding {
                            technique: Technique::PointingPair,
                            inference: InferenceResult::Elimination {
                                cell,
                                values: vec![digit],
                            },
                            involved_cells: involved,
                            explanation: ExplanationData::Intersection {
                                kind: "Pointing Pair",
                                digit,
                                from_sector: sector_name_long(base),
                                to_sector: sector_name_long(cover),
                            },
                            proof: Some(ProofCertificate::Fish {
                                digit,
                                base_sectors: vec![base],
                                cover_sectors: vec![cover],
                                fins: vec![],
                            }),
                        });
                    }
                }
            }
        }
    }
    None
}

/// Size-1 fish: box-line reduction (line → box intersection).
///
/// When all candidates for a digit within a row or column are confined to
/// a single box, that digit can be eliminated from the rest of that box.
pub fn find_box_line_reduction(fab: &CandidateFabric) -> Option<Finding> {
    for digit in 1..=9u8 {
        for line in 0..18usize {
            let base_mask = sector_candidate_mask(fab, line, digit);
            let base_count = base_mask.count_ones();
            if !(2..=3).contains(&base_count) {
                continue;
            }

            for box_idx in 0..9 {
                let cover = SECTOR_BOX_BASE + box_idx;
                let cover_mask = sector_candidate_mask(fab, cover, digit);
                // All base cells must lie within the cover sector
                if base_mask & cover_mask != base_mask {
                    continue;
                }
                // Eliminations: candidates in cover sector outside the base
                let elim_mask = cover_mask & !base_mask;
                if elim_mask == 0 {
                    continue;
                }
                for cell in 0..81usize {
                    if elim_mask & (1u128 << cell) != 0 {
                        let involved: Vec<usize> =
                            (0..81).filter(|&c| base_mask & (1u128 << c) != 0).collect();
                        return Some(Finding {
                            technique: Technique::BoxLineReduction,
                            inference: InferenceResult::Elimination {
                                cell,
                                values: vec![digit],
                            },
                            involved_cells: involved,
                            explanation: ExplanationData::Intersection {
                                kind: "Box/Line Reduction",
                                digit,
                                from_sector: sector_name_long(line),
                                to_sector: sector_name_long(cover),
                            },
                            proof: Some(ProofCertificate::Fish {
                                digit,
                                base_sectors: vec![line],
                                cover_sectors: vec![cover],
                                fins: vec![],
                            }),
                        });
                    }
                }
            }
        }
    }
    None
}

// ==================== Public API (size 2+) ====================

/// Search for basic fish (X-Wing, Swordfish, Jellyfish) and finned variants.
pub fn find_basic_fish(fab: &CandidateFabric, size: usize) -> Option<Finding> {
    for digit in 1..=9u8 {
        if let Some(f) = search_fish_for_digit(fab, digit, size, SectorConstraint::Basic) {
            return Some(f);
        }
    }
    None
}

/// Search for finned fish of given size.
pub fn find_finned_fish(fab: &CandidateFabric, size: usize) -> Option<Finding> {
    // The search_fish_for_digit already handles finned variants
    // We need to filter for finned-only results
    for digit in 1..=9u8 {
        if let Some(f) = search_fish_for_digit(fab, digit, size, SectorConstraint::Basic) {
            // Check if it's actually finned (the technique classification tells us)
            match f.technique {
                Technique::FinnedXWing
                | Technique::FinnedSwordfish
                | Technique::FinnedJellyfish => {
                    return Some(f);
                }
                _ => {}
            }
        }
    }
    None
}

/// Search for Franken fish (mixed line+box sectors).
pub fn find_franken_fish(fab: &CandidateFabric) -> Option<Finding> {
    for digit in 1..=9u8 {
        for size in 2..=4 {
            if let Some(f) = search_fish_for_digit(fab, digit, size, SectorConstraint::Franken) {
                return Some(f);
            }
        }
    }
    None
}

/// Search for Mutant fish (all three sector types).
pub fn find_mutant_fish(fab: &CandidateFabric) -> Option<Finding> {
    for digit in 1..=9u8 {
        for size in 2..=4 {
            if let Some(f) = search_fish_for_digit(fab, digit, size, SectorConstraint::Mutant) {
                return Some(f);
            }
        }
    }
    None
}

// ==================== Combination utility ====================

fn combinations(items: &[usize], k: usize) -> Vec<Vec<usize>> {
    let mut result = Vec::new();
    if k == 0 || k > items.len() {
        return result;
    }
    let mut indices: Vec<usize> = (0..k).collect();
    loop {
        result.push(indices.iter().map(|&i| items[i]).collect());
        let mut i = k;
        loop {
            if i == 0 {
                return result;
            }
            i -= 1;
            indices[i] += 1;
            if indices[i] <= items.len() - k + i {
                break;
            }
        }
        for j in (i + 1)..k {
            indices[j] = indices[j - 1] + 1;
        }
    }
}
