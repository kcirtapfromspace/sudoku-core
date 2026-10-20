//! Backtracking solver and propagation helpers.
//!
//! Contains: solve_recursive, count_solutions_recursive,
//! propagate_singles, propagate_full, has_contradiction.

use crate::{BitSet, Grid, Position};

/// Count classic Sudoku completions using only occupancy masks and empty cells.
/// Stored pencil marks and given flags do not restrict solution counting. The
/// caller must establish that the installed constraints are exactly classic.
pub(super) fn count_classic_solutions(grid: &Grid, limit: usize) -> usize {
    let mut state = ClassicCounter {
        rows: [0; 9],
        columns: [0; 9],
        boxes: [0; 9],
        empty: [Position::new(0, 0); 81],
    };
    let mut empty_count = 0;
    for pos in Position::all_9x9() {
        if let Some(value) = grid.get(pos) {
            if !(1..=9).contains(&value) {
                return 0;
            }
            let bit = 1 << value;
            let box_index = pos.box_index();
            if (state.rows[pos.row] | state.columns[pos.col] | state.boxes[box_index]) & bit != 0 {
                return 0;
            }
            state.rows[pos.row] |= bit;
            state.columns[pos.col] |= bit;
            state.boxes[box_index] |= bit;
        } else {
            state.empty[empty_count] = pos;
            empty_count += 1;
        }
    }
    state.count(empty_count, limit)
}

struct ClassicCounter {
    rows: [u16; 9],
    columns: [u16; 9],
    boxes: [u16; 9],
    empty: [Position; 81],
}

impl ClassicCounter {
    fn candidates(&self, pos: Position) -> u16 {
        BitSet::all_9().as_raw()
            & !(self.rows[pos.row] | self.columns[pos.col] | self.boxes[pos.box_index()])
    }

    fn count(&mut self, remaining: usize, limit: usize) -> usize {
        if remaining == 0 {
            return 1;
        }

        // MRV also propagates forced placements by trying singleton masks first.
        let mut best = 0;
        let mut candidates = self.candidates(self.empty[0]);
        for index in 1..remaining {
            if candidates.count_ones() <= 1 {
                break;
            }
            let mask = self.candidates(self.empty[index]);
            if mask.count_ones() < candidates.count_ones() {
                best = index;
                candidates = mask;
            }
        }
        if candidates == 0 {
            return 0;
        }

        // Remove the chosen cell from the active prefix, restoring the ordering
        // and all occupancy masks on return, including when the cap is reached.
        self.empty.swap(best, remaining - 1);
        let pos = self.empty[remaining - 1];
        let box_index = pos.box_index();
        let mut count = 0;
        while candidates != 0 && count < limit {
            let bit = 1 << candidates.trailing_zeros();
            candidates &= candidates - 1;
            self.rows[pos.row] |= bit;
            self.columns[pos.col] |= bit;
            self.boxes[box_index] |= bit;
            count += self.count(remaining - 1, limit - count);
            self.rows[pos.row] &= !bit;
            self.columns[pos.col] &= !bit;
            self.boxes[box_index] &= !bit;
        }
        self.empty.swap(best, remaining - 1);
        count
    }
}

/// Check if a grid has a contradiction.
pub fn has_contradiction(grid: &Grid) -> bool {
    for pos in grid.empty_positions() {
        if grid.get_candidates(pos).is_empty() {
            return true;
        }
    }
    for i in 0..9 {
        let mut row_seen = [false; 10];
        let mut col_seen = [false; 10];
        let mut box_seen = [false; 10];
        for j in 0..9 {
            if let Some(v) = grid.get(Position::new(i, j)) {
                if row_seen[v as usize] {
                    return true;
                }
                row_seen[v as usize] = true;
            }
            if let Some(v) = grid.get(Position::new(j, i)) {
                if col_seen[v as usize] {
                    return true;
                }
                col_seen[v as usize] = true;
            }
            let box_row = (i / 3) * 3 + j / 3;
            let box_col = (i % 3) * 3 + j % 3;
            if let Some(v) = grid.get(Position::new(box_row, box_col)) {
                if box_seen[v as usize] {
                    return true;
                }
                box_seen[v as usize] = true;
            }
        }
    }
    false
}

/// Propagate naked+hidden singles from an assumption until no more progress.
pub fn propagate_singles(grid: &Grid, pos: Position, val: u8) -> (Grid, bool) {
    let mut g = grid.deep_clone();
    g.set_cell_unchecked(pos, Some(val));
    g.recalculate_candidates();

    for _ in 0..200 {
        if has_contradiction(&g) {
            return (g, true);
        }
        if g.is_complete() {
            return (g, false);
        }
        let mut progress = false;
        // Naked singles
        for p in g.empty_positions() {
            if let Some(v) = g.get_candidates(p).single_value() {
                g.set_cell_unchecked(p, Some(v));
                g.recalculate_candidates();
                progress = true;
                break;
            }
        }
        if !progress {
            // Hidden singles
            'outer: for unit in 0..27 {
                let positions: Vec<Position> = if unit < 9 {
                    (0..9).map(|c| Position::new(unit, c)).collect()
                } else if unit < 18 {
                    (0..9).map(|r| Position::new(r, unit - 9)).collect()
                } else {
                    let bi = unit - 18;
                    let br = (bi / 3) * 3;
                    let bc = (bi % 3) * 3;
                    (0..9)
                        .map(|i| Position::new(br + i / 3, bc + i % 3))
                        .collect()
                };
                for value in 1..=9u8 {
                    let mut candidates: Vec<Position> = Vec::new();
                    for &p in &positions {
                        if g.cell(p).is_empty() && g.get_candidates(p).contains(value) {
                            candidates.push(p);
                        }
                    }
                    if candidates.len() == 1 {
                        g.set_cell_unchecked(candidates[0], Some(value));
                        g.recalculate_candidates();
                        progress = true;
                        break 'outer;
                    }
                }
            }
        }
        if !progress {
            break;
        }
    }
    let contradiction = has_contradiction(&g);
    (g, contradiction)
}

/// Solve a grid using backtracking with MRV heuristic.
pub fn solve_recursive(grid: &mut Grid) -> bool {
    // Apply singles first
    apply_naked_singles(grid);
    apply_hidden_singles(grid);

    if grid.is_complete() {
        return true;
    }

    let empty_positions = grid.empty_positions();
    if empty_positions.is_empty() {
        return false;
    }

    let best_pos = empty_positions
        .into_iter()
        .min_by_key(|&pos| grid.get_candidates(pos).count())
        .unwrap();

    let candidates = grid.get_candidates(best_pos);
    if candidates.is_empty() {
        return false;
    }

    for value in candidates.iter() {
        let mut test_grid = grid.deep_clone();
        test_grid.set_cell_unchecked(best_pos, Some(value));
        test_grid.recalculate_candidates();

        if test_grid.validate().is_valid && solve_recursive(&mut test_grid) {
            for row in 0..9 {
                for col in 0..9 {
                    let pos = Position::new(row, col);
                    grid.set_cell_unchecked(pos, test_grid.get(pos));
                }
            }
            return true;
        }
    }

    false
}

/// Count solutions up to a limit.
pub fn count_solutions_recursive(grid: &mut Grid, count: &mut usize, limit: usize) {
    if *count >= limit {
        return;
    }

    apply_naked_singles(grid);
    apply_hidden_singles(grid);

    if grid.is_complete() {
        *count += 1;
        return;
    }

    let empty_positions = grid.empty_positions();
    if empty_positions.is_empty() {
        return;
    }

    let best_pos = empty_positions
        .into_iter()
        .min_by_key(|&pos| grid.get_candidates(pos).count())
        .unwrap();

    let candidates = grid.get_candidates(best_pos);
    if candidates.is_empty() {
        return;
    }

    for value in candidates.iter() {
        if *count >= limit {
            return;
        }

        let mut test_grid = grid.deep_clone();
        test_grid.set_cell_unchecked(best_pos, Some(value));
        test_grid.recalculate_candidates();

        if test_grid.validate().is_valid {
            count_solutions_recursive(&mut test_grid, count, limit);
        }
    }
}

// ==================== Inline simple technique appliers ====================

fn apply_naked_singles(grid: &mut Grid) {
    // Preserve the first rebuild, which discards incoming pencil eliminations.
    // After that, only new values are added, so updating their peers is enough.
    let mut rebuilt = false;
    loop {
        let mut progress = false;
        for pos in grid.empty_positions() {
            if let Some(v) = grid.get_candidates(pos).single_value() {
                grid.set_cell_unchecked(pos, Some(v));
                if rebuilt {
                    grid.update_candidates_after_move(pos, v);
                } else {
                    grid.recalculate_candidates();
                    rebuilt = true;
                }
                progress = true;
                break;
            }
        }
        if !progress {
            break;
        }
    }
}

fn apply_hidden_singles(grid: &mut Grid) {
    // Match the first-rebuild behavior of the naked-single pass above.
    let mut rebuilt = false;
    loop {
        let mut progress = false;
        'outer: for unit in 0..27 {
            let positions: Vec<Position> = if unit < 9 {
                (0..9).map(|c| Position::new(unit, c)).collect()
            } else if unit < 18 {
                (0..9).map(|r| Position::new(r, unit - 9)).collect()
            } else {
                let bi = unit - 18;
                let br = (bi / 3) * 3;
                let bc = (bi % 3) * 3;
                (0..9)
                    .map(|i| Position::new(br + i / 3, bc + i % 3))
                    .collect()
            };
            for value in 1..=9u8 {
                let mut candidates: Vec<Position> = Vec::new();
                for &p in &positions {
                    if grid.cell(p).is_empty() && grid.get_candidates(p).contains(value) {
                        candidates.push(p);
                    }
                }
                if candidates.len() == 1 {
                    grid.set_cell_unchecked(candidates[0], Some(value));
                    if rebuilt {
                        grid.update_candidates_after_move(candidates[0], value);
                    } else {
                        grid.recalculate_candidates();
                        rebuilt = true;
                    }
                    progress = true;
                    break 'outer;
                }
            }
        }
        if !progress {
            break;
        }
    }
}

/// Generate a backtracking Finding for the hint system.
pub fn find_backtracking_hint(grid: &Grid) -> Option<super::explain::Finding> {
    let mut g = grid.deep_clone();
    if solve_recursive(&mut g) {
        for pos in grid.empty_positions() {
            if let Some(val) = g.get(pos) {
                let cell = pos.row * 9 + pos.col;
                return Some(super::explain::Finding {
                    technique: super::types::Technique::Backtracking,
                    inference: super::explain::InferenceResult::Placement { cell, value: val },
                    involved_cells: vec![cell],
                    explanation: super::explain::ExplanationData::Backtracking { cell, value: val },
                    proof: Some(super::explain::ProofCertificate::Backtracking),
                });
            }
        }
    }
    None
}

#[cfg(test)]
#[path = "backtrack_tests.rs"]
mod tests;
