//! Verus-verified sudoku validity and solver soundness (Tier 2).
//!
//! This file gives the engine's core rules a single spec-level definition of
//! "valid sudoku grid" (`valid`) and proves the production algorithms against it:
//!
//! - `validate_row` / `validate_col` / `validate_box` mirror the classic
//!   `Constraint::validate` impls in `src/constraint.rs` (same loops, same
//!   `(row / 3) * 3` index math) and are proven equivalent to their specs.
//! - `lemma_place_preserves_valid` is the placement soundness theorem: a value
//!   accepted by all three constraint checks keeps the grid valid.
//! - `has_duplicates` mirrors the duplicate scan of `has_contradiction` in
//!   `src/solver/backtrack.rs` (the 27-unit row/col/box decomposition with the
//!   exact `(u / 3) * 3 + k / 3` box arithmetic used there) and is proven to
//!   return true **iff** the board violates `valid`.
//! - `solve_from` / `solve` is a backtracking solver proven *sound*: if it
//!   returns true, the board is complete, valid, and preserves every given.
//!   (Completeness — "false means unsolvable" — is intentionally out of scope.)
//!
//! The board model is the flat `[Option<u8>; 81]` value grid (index = row*9+col,
//! as in the solver's `Finding.cell`); candidates/pencil marks are out of scope.
//! Note `has_duplicates` *requires* `wf` (values in 1..=9): the production scan
//! indexes `seen[v as usize]` and would panic on a corrupt value — Verus makes
//! that hidden invariant an explicit precondition.
//!
//! Verify with: `scripts/verify.sh`.

use vstd::prelude::*;

verus! {

// ==================== Board model ====================

pub type Board = [Option<u8>; 81];

/// Cell accessor on the abstract board.
pub open spec fn cell(s: Seq<Option<u8>>, r: int, c: int) -> Option<u8> {
    s[r * 9 + c]
}

pub open spec fn range_ok(o: Option<u8>) -> bool {
    o matches Some(v) ==> 1 <= v <= 9
}

/// Well-formed board: 81 cells, every filled value in 1..=9
/// (the implicit invariant of `Grid`/`Cell` in the production engine).
pub open spec fn wf(s: Seq<Option<u8>>) -> bool {
    &&& s.len() == 81
    &&& forall|i: int| 0 <= i < 81 ==> range_ok(#[trigger] s[i])
}

pub open spec fn same_box(r1: int, c1: int, r2: int, c2: int) -> bool {
    r1 / 3 == r2 / 3 && c1 / 3 == c2 / 3
}

/// Two distinct positions constrain each other in classic sudoku.
pub open spec fn sees(r1: int, c1: int, r2: int, c2: int) -> bool {
    r1 == r2 || c1 == c2 || same_box(r1, c1, r2, c2)
}

/// THE validity spec: no two distinct seeing cells hold the same value.
pub open spec fn valid(s: Seq<Option<u8>>) -> bool {
    forall|r1: int, c1: int, r2: int, c2: int|
        0 <= r1 < 9 && 0 <= c1 < 9 && 0 <= r2 < 9 && 0 <= c2 < 9 && !(r1 == r2 && c1 == c2)
            && sees(r1, c1, r2, c2) && (#[trigger] cell(s, r1, c1)) is Some ==> cell(s, r1, c1)
            != #[trigger] cell(s, r2, c2)
}

pub open spec fn complete(s: Seq<Option<u8>>) -> bool {
    forall|i: int| 0 <= i < 81 ==> (#[trigger] s[i]) is Some
}

/// Every filled cell of `orig` (the givens) is unchanged in `s`.
pub open spec fn extends(orig: Seq<Option<u8>>, s: Seq<Option<u8>>) -> bool {
    forall|i: int| 0 <= i < 81 && (#[trigger] orig[i]) is Some ==> s[i] == orig[i]
}

// ==================== Arithmetic lemmas ====================

/// Division shape: (3m + d) / 3 == m and (3m + d) % 3 == d for 0 <= d < 3.
proof fn lemma_div3_shape(m: int, d: int)
    requires
        0 <= m,
        0 <= d < 3,
    ensures
        (3 * m + d) / 3 == m,
        (3 * m + d) % 3 == d,
{
}

// ==================== Constraint checks (mirror src/constraint.rs) ====================

pub open spec fn spec_row_ok(s: Seq<Option<u8>>, r: int, c: int, v: u8) -> bool {
    forall|cc: int| 0 <= cc < 9 && cc != c ==> #[trigger] cell(s, r, cc) != Some(v)
}

pub open spec fn spec_col_ok(s: Seq<Option<u8>>, r: int, c: int, v: u8) -> bool {
    forall|rr: int| 0 <= rr < 9 && rr != r ==> #[trigger] cell(s, rr, c) != Some(v)
}

pub open spec fn spec_box_ok(s: Seq<Option<u8>>, r: int, c: int, v: u8) -> bool {
    forall|rr: int, cc: int|
        0 <= rr < 9 && 0 <= cc < 9 && same_box(rr, cc, r, c) && !(rr == r && cc == c)
            ==> #[trigger] cell(s, rr, cc) != Some(v)
}

/// Placing `v` at `(r, c)` violates no classic constraint (the conjunction the
/// production `Grid` checks by looping over its `ConstraintBox`es).
pub open spec fn spec_legal(s: Seq<Option<u8>>, r: int, c: int, v: u8) -> bool {
    spec_row_ok(s, r, c, v) && spec_col_ok(s, r, c, v) && spec_box_ok(s, r, c, v)
}

/// Mirrors `RowConstraint::validate`.
pub fn validate_row(cells: &Board, r: usize, c: usize, v: u8) -> (ok: bool)
    requires
        r < 9,
        c < 9,
    ensures
        ok == spec_row_ok(cells@, r as int, c as int, v),
{
    let mut col: usize = 0;
    while col < 9
        invariant
            r < 9,
            c < 9,
            col <= 9,
            forall|cc: int|
                0 <= cc < col && cc != c ==> #[trigger] cell(cells@, r as int, cc) != Some(v),
        decreases 9 - col,
    {
        if col != c {
            match cells[r * 9 + col] {
                Some(x) => if x == v {
                    proof {
                        assert(cell(cells@, r as int, col as int) == Some(v));
                    }
                    return false;
                },
                None => {},
            }
        }
        col += 1;
    }
    true
}

/// Mirrors `ColumnConstraint::validate`.
pub fn validate_col(cells: &Board, r: usize, c: usize, v: u8) -> (ok: bool)
    requires
        r < 9,
        c < 9,
    ensures
        ok == spec_col_ok(cells@, r as int, c as int, v),
{
    let mut row: usize = 0;
    while row < 9
        invariant
            r < 9,
            c < 9,
            row <= 9,
            forall|rr: int|
                0 <= rr < row && rr != r ==> #[trigger] cell(cells@, rr, c as int) != Some(v),
        decreases 9 - row,
    {
        if row != r {
            match cells[row * 9 + c] {
                Some(x) => if x == v {
                    proof {
                        assert(cell(cells@, row as int, c as int) == Some(v));
                    }
                    return false;
                },
                None => {},
            }
        }
        row += 1;
    }
    true
}

/// A position inside the 3x3 block starting at `(r/3)*3` shares `r`'s box row.
proof fn lemma_in_block_same_third(r: int, rr: int)
    requires
        0 <= r < 9,
        (r / 3) * 3 <= rr < (r / 3) * 3 + 3,
    ensures
        rr / 3 == r / 3,
        0 <= rr < 9,
{
    lemma_div3_shape(r / 3, rr - (r / 3) * 3);
}

/// Conversely, sharing `r`'s box row places `rr` inside that 3x3 block.
proof fn lemma_same_third_in_block(r: int, rr: int)
    requires
        0 <= r < 9,
        0 <= rr < 9,
        rr / 3 == r / 3,
    ensures
        (r / 3) * 3 <= rr < (r / 3) * 3 + 3,
{
    assert(rr == 3 * (rr / 3) + rr % 3);
}

/// Mirrors `BoxConstraint::validate` (same `(pos.row / 3) * 3` arithmetic).
pub fn validate_box(cells: &Board, r: usize, c: usize, v: u8) -> (ok: bool)
    requires
        r < 9,
        c < 9,
    ensures
        ok == spec_box_ok(cells@, r as int, c as int, v),
{
    let box_row = (r / 3) * 3;
    let box_col = (c / 3) * 3;
    let mut row = box_row;
    while row < box_row + 3
        invariant
            r < 9,
            c < 9,
            box_row == (r / 3) * 3,
            box_col == (c / 3) * 3,
            box_row <= row <= box_row + 3,
            forall|rr: int, cc: int|
                box_row <= rr < row && box_col <= cc < box_col + 3 && !(rr == r && cc == c)
                    ==> #[trigger] cell(cells@, rr, cc) != Some(v),
        decreases box_row + 3 - row,
    {
        let mut col = box_col;
        while col < box_col + 3
            invariant
                r < 9,
                c < 9,
                box_row == (r / 3) * 3,
                box_col == (c / 3) * 3,
                box_row <= row < box_row + 3,
                box_col <= col <= box_col + 3,
                forall|rr: int, cc: int|
                    box_row <= rr < row && box_col <= cc < box_col + 3 && !(rr == r && cc == c)
                        ==> #[trigger] cell(cells@, rr, cc) != Some(v),
                forall|cc: int|
                    box_col <= cc < col && !(row == r && cc == c) ==> #[trigger] cell(
                        cells@,
                        row as int,
                        cc,
                    ) != Some(v),
            decreases box_col + 3 - col,
        {
            proof {
                lemma_in_block_same_third(r as int, row as int);
                lemma_in_block_same_third(c as int, col as int);
            }
            if row != r || col != c {
                match cells[row * 9 + col] {
                    Some(x) => if x == v {
                        proof {
                            assert(same_box(row as int, col as int, r as int, c as int));
                            assert(cell(cells@, row as int, col as int) == Some(v));
                        }
                        return false;
                    },
                    None => {},
                }
            }
            col += 1;
        }
        row += 1;
    }
    proof {
        assert forall|rr: int, cc: int|
            0 <= rr < 9 && 0 <= cc < 9 && same_box(rr, cc, r as int, c as int) && !(rr == r as int
                && cc == c as int) implies #[trigger] cell(cells@, rr, cc) != Some(v) by {
            lemma_same_third_in_block(r as int, rr);
            lemma_same_third_in_block(c as int, cc);
        }
    }
    true
}

/// The conjunction the production engine evaluates over its constraint set.
pub fn is_legal(cells: &Board, r: usize, c: usize, v: u8) -> (ok: bool)
    requires
        r < 9,
        c < 9,
    ensures
        ok == spec_legal(cells@, r as int, c as int, v),
{
    validate_row(cells, r, c, v) && validate_col(cells, r, c, v) && validate_box(cells, r, c, v)
}

// ==================== Placement soundness ====================

/// THE placement theorem: on a valid board, a placement passing the three
/// constraint checks yields a valid board. This is what justifies every
/// `set_cell` the solver performs after `Constraint::validate` approves it.
proof fn lemma_place_preserves_valid(s: Seq<Option<u8>>, r: int, c: int, v: u8)
    requires
        s.len() == 81,
        valid(s),
        0 <= r < 9,
        0 <= c < 9,
        spec_legal(s, r, c, v),
    ensures
        valid(s.update(r * 9 + c, Some(v))),
{
    let t = s.update(r * 9 + c, Some(v));
    assert forall|r1: int, c1: int, r2: int, c2: int|
        0 <= r1 < 9 && 0 <= c1 < 9 && 0 <= r2 < 9 && 0 <= c2 < 9 && !(r1 == r2 && c1 == c2)
            && sees(r1, c1, r2, c2) && (#[trigger] cell(t, r1, c1)) is Some implies cell(t, r1, c1)
        != #[trigger] cell(t, r2, c2) by {
        if r1 == r && c1 == c {
            // p1 is the new placement; p2 is untouched.
            assert(r2 * 9 + c2 != r * 9 + c);
            assert(cell(t, r2, c2) == cell(s, r2, c2));
            assert(cell(t, r1, c1) == Some(v));
            if r2 == r1 {
                assert(cell(s, r, c2) != Some(v));
            } else if c2 == c1 {
                assert(cell(s, r2, c) != Some(v));
            } else {
                assert(same_box(r2, c2, r, c));
                assert(cell(s, r2, c2) != Some(v));
            }
        } else if r2 == r && c2 == c {
            // p2 is the new placement; p1 is untouched.
            assert(r1 * 9 + c1 != r * 9 + c);
            assert(cell(t, r1, c1) == cell(s, r1, c1));
            assert(cell(t, r2, c2) == Some(v));
            if cell(s, r1, c1) == Some(v) {
                if r1 == r2 {
                    assert(cell(s, r, c1) != Some(v));
                } else if c1 == c2 {
                    assert(cell(s, r1, c) != Some(v));
                } else {
                    assert(same_box(r1, c1, r, c));
                    assert(cell(s, r1, c1) != Some(v));
                }
            }
        } else {
            // Neither endpoint is the new placement: valid(s) applies verbatim.
            assert(r1 * 9 + c1 != r * 9 + c && r2 * 9 + c2 != r * 9 + c);
            assert(cell(t, r1, c1) == cell(s, r1, c1));
            assert(cell(t, r2, c2) == cell(s, r2, c2));
        }
    }
}

// ==================== The 27-unit decomposition ====================
// Exactly the unit indexing of `src/solver/backtrack.rs` (propagate_singles /
// apply_hidden_singles): unit 0-8 = rows, 9-17 = columns, 18-26 = boxes with
// box slot math `(bi / 3) * 3 + k / 3` / `(bi % 3) * 3 + k % 3`.

pub open spec fn unit_row(u: int, k: int) -> int {
    if u < 9 {
        u
    } else if u < 18 {
        k
    } else {
        ((u - 18) / 3) * 3 + k / 3
    }
}

pub open spec fn unit_col(u: int, k: int) -> int {
    if u < 9 {
        k
    } else if u < 18 {
        u - 9
    } else {
        ((u - 18) % 3) * 3 + k % 3
    }
}

pub open spec fn unit_cell(s: Seq<Option<u8>>, u: int, k: int) -> Option<u8> {
    cell(s, unit_row(u, k), unit_col(u, k))
}

/// No two distinct slots of unit `u` hold the same value.
pub open spec fn unit_dupfree(s: Seq<Option<u8>>, u: int) -> bool {
    forall|k1: int, k2: int|
        0 <= k1 < k2 < 9 && (#[trigger] unit_cell(s, u, k1)) is Some ==> unit_cell(s, u, k1)
            != #[trigger] unit_cell(s, u, k2)
}

proof fn lemma_unit_pos_bounds(u: int, k: int)
    requires
        0 <= u < 27,
        0 <= k < 9,
    ensures
        0 <= unit_row(u, k) < 9,
        0 <= unit_col(u, k) < 9,
{
}

/// Distinct slots of one unit are distinct, mutually-seeing positions.
proof fn lemma_unit_slots_see(u: int, k1: int, k2: int)
    requires
        0 <= u < 27,
        0 <= k1 < 9,
        0 <= k2 < 9,
        k1 != k2,
    ensures
        !(unit_row(u, k1) == unit_row(u, k2) && unit_col(u, k1) == unit_col(u, k2)),
        sees(unit_row(u, k1), unit_col(u, k1), unit_row(u, k2), unit_col(u, k2)),
{
    if u >= 18 {
        let b = u - 18;
        assert(k1 == 3 * (k1 / 3) + k1 % 3);
        assert(k2 == 3 * (k2 / 3) + k2 % 3);
        lemma_div3_shape(b / 3, k1 / 3);
        lemma_div3_shape(b / 3, k2 / 3);
        assert(b == 3 * (b / 3) + b % 3);
        assert(unit_row(u, k1) / 3 == b / 3 && unit_row(u, k2) / 3 == b / 3);
        lemma_div3_shape(b % 3, k1 % 3);
        lemma_div3_shape(b % 3, k2 % 3);
        assert(unit_col(u, k1) / 3 == b % 3 && unit_col(u, k2) / 3 == b % 3);
        assert(same_box(unit_row(u, k1), unit_col(u, k1), unit_row(u, k2), unit_col(u, k2)));
    }
}

/// Every distinct seeing pair lives in some common unit at distinct slots.
proof fn lemma_seeing_pair_in_unit(r1: int, c1: int, r2: int, c2: int)
    requires
        0 <= r1 < 9,
        0 <= c1 < 9,
        0 <= r2 < 9,
        0 <= c2 < 9,
        !(r1 == r2 && c1 == c2),
        sees(r1, c1, r2, c2),
    ensures
        exists|u: int, k1: int, k2: int|
            0 <= u < 27 && 0 <= k1 < 9 && 0 <= k2 < 9 && k1 != k2 && #[trigger] unit_row(u, k1)
                == r1 && unit_col(u, k1) == c1 && #[trigger] unit_row(u, k2) == r2 && unit_col(
                u,
                k2,
            ) == c2,
{
    if r1 == r2 {
        assert(unit_row(r1, c1) == r1 && unit_col(r1, c1) == c1 && unit_row(r1, c2) == r2
            && unit_col(r1, c2) == c2);
    } else if c1 == c2 {
        assert(unit_row(9 + c1, r1) == r1 && unit_col(9 + c1, r1) == c1 && unit_row(9 + c1, r2)
            == r2 && unit_col(9 + c1, r2) == c2);
    } else {
        let b = (r1 / 3) * 3 + c1 / 3;
        let u = 18 + b;
        let k1 = (r1 % 3) * 3 + c1 % 3;
        let k2 = (r2 % 3) * 3 + c2 % 3;
        assert(same_box(r1, c1, r2, c2));
        lemma_div3_shape(r1 / 3, c1 / 3);
        assert(b / 3 == r1 / 3 && b % 3 == c1 / 3);
        lemma_div3_shape(r1 % 3, c1 % 3);
        lemma_div3_shape(r2 % 3, c2 % 3);
        assert(k1 / 3 == r1 % 3 && k1 % 3 == c1 % 3);
        assert(k2 / 3 == r2 % 3 && k2 % 3 == c2 % 3);
        assert(r1 == 3 * (r1 / 3) + r1 % 3 && c1 == 3 * (c1 / 3) + c1 % 3);
        assert(r2 == 3 * (r2 / 3) + r2 % 3 && c2 == 3 * (c2 / 3) + c2 % 3);
        assert(unit_row(u, k1) == r1 && unit_col(u, k1) == c1);
        assert(unit_row(u, k2) == r2 && unit_col(u, k2) == c2);
        assert(k1 != k2);
    }
}

/// A duplicate inside any unit falsifies `valid`.
proof fn lemma_unit_dup_invalid(s: Seq<Option<u8>>, u: int, k1: int, k2: int)
    requires
        s.len() == 81,
        0 <= u < 27,
        0 <= k1 < 9,
        0 <= k2 < 9,
        k1 != k2,
        unit_cell(s, u, k1) is Some,
        unit_cell(s, u, k1) == unit_cell(s, u, k2),
    ensures
        !valid(s),
{
    lemma_unit_pos_bounds(u, k1);
    lemma_unit_pos_bounds(u, k2);
    lemma_unit_slots_see(u, k1, k2);
    if valid(s) {
        assert(cell(s, unit_row(u, k1), unit_col(u, k1)) != cell(
            s,
            unit_row(u, k2),
            unit_col(u, k2),
        ));
    }
}

/// If every unit is duplicate-free, the board is valid.
proof fn lemma_all_units_dupfree_valid(s: Seq<Option<u8>>)
    requires
        s.len() == 81,
        forall|u: int| 0 <= u < 27 ==> #[trigger] unit_dupfree(s, u),
    ensures
        valid(s),
{
    assert forall|r1: int, c1: int, r2: int, c2: int|
        0 <= r1 < 9 && 0 <= c1 < 9 && 0 <= r2 < 9 && 0 <= c2 < 9 && !(r1 == r2 && c1 == c2)
            && sees(r1, c1, r2, c2) && (#[trigger] cell(s, r1, c1)) is Some implies cell(s, r1, c1)
        != #[trigger] cell(s, r2, c2) by {
        lemma_seeing_pair_in_unit(r1, c1, r2, c2);
        let (u, k1, k2) = choose|u: int, k1: int, k2: int|
            0 <= u < 27 && 0 <= k1 < 9 && 0 <= k2 < 9 && k1 != k2 && #[trigger] unit_row(u, k1)
                == r1 && unit_col(u, k1) == c1 && #[trigger] unit_row(u, k2) == r2 && unit_col(
                u,
                k2,
            ) == c2;
        assert(unit_dupfree(s, u));
        assert(unit_cell(s, u, k1) == cell(s, r1, c1));
        assert(unit_cell(s, u, k2) == cell(s, r2, c2));
        if cell(s, r1, c1) == cell(s, r2, c2) {
            if k1 < k2 {
                assert(unit_cell(s, u, k1) != unit_cell(s, u, k2));
            } else {
                assert(unit_cell(s, u, k2) != unit_cell(s, u, k1));
            }
        }
    }
}

// ==================== Duplicate scan (mirrors has_contradiction) ====================

/// Exec computation of a unit slot's row, identical to the index arithmetic in
/// `src/solver/backtrack.rs`.
fn unit_row_exec(u: usize, k: usize) -> (r: usize)
    requires
        u < 27,
        k < 9,
    ensures
        r as int == unit_row(u as int, k as int),
        r < 9,
{
    proof {
        lemma_unit_pos_bounds(u as int, k as int);
    }
    if u < 9 {
        u
    } else if u < 18 {
        k
    } else {
        ((u - 18) / 3) * 3 + k / 3
    }
}

fn unit_col_exec(u: usize, k: usize) -> (c: usize)
    requires
        u < 27,
        k < 9,
    ensures
        c as int == unit_col(u as int, k as int),
        c < 9,
{
    proof {
        lemma_unit_pos_bounds(u as int, k as int);
    }
    if u < 9 {
        k
    } else if u < 18 {
        u - 9
    } else {
        ((u - 18) % 3) * 3 + k % 3
    }
}

/// Scan one unit with a seen-array, as `has_contradiction` does per row/col/box.
/// Note the `wf` precondition: the production code's `seen[v as usize]` would
/// panic on a value > 9 — that latent invariant is explicit here.
fn unit_has_dup(cells: &Board, u: usize) -> (dup: bool)
    requires
        u < 27,
        wf(cells@),
    ensures
        dup == !unit_dupfree(cells@, u as int),
{
    let mut seen = [false; 10];
    let mut k: usize = 0;
    while k < 9
        invariant
            u < 27,
            k <= 9,
            wf(cells@),
            seen@.len() == 10,
            forall|vv: u8|
                1 <= vv <= 9 ==> ((#[trigger] seen@[vv as int]) <==> exists|k2: int|
                    0 <= k2 < k && unit_cell(cells@, u as int, k2) == Some(vv)),
            forall|k1: int, k2: int|
                0 <= k1 < k2 < k && (#[trigger] unit_cell(cells@, u as int, k1)) is Some
                    ==> unit_cell(cells@, u as int, k1) != #[trigger] unit_cell(
                    cells@,
                    u as int,
                    k2,
                ),
        decreases 9 - k,
    {
        let r = unit_row_exec(u, k);
        let c = unit_col_exec(u, k);
        let ghost s = cells@;
        match cells[r * 9 + c] {
            Some(v) => {
                proof {
                    assert(range_ok(s[r as int * 9 + c as int]));
                    assert(unit_cell(s, u as int, k as int) == Some(v));
                }
                if seen[v as usize] {
                    proof {
                        let k2 = choose|k2: int|
                            0 <= k2 < k && unit_cell(s, u as int, k2) == Some(v);
                        if unit_dupfree(s, u as int) {
                            assert(unit_cell(s, u as int, k2) != unit_cell(s, u as int, k as int));
                        }
                    }
                    return true;
                }
                seen[v as usize] = true;
                proof {
                    let ghost k_new = k + 1;
                    assert forall|vv: u8|
                        1 <= vv <= 9 implies ((#[trigger] seen@[vv as int]) <==> exists|k2: int|
                        0 <= k2 < k_new && unit_cell(s, u as int, k2) == Some(vv)) by {
                        if vv == v {
                            assert(unit_cell(s, u as int, k as int) == Some(vv));
                        } else {
                            if seen@[vv as int] {
                                let k2 = choose|k2: int|
                                    0 <= k2 < k && unit_cell(s, u as int, k2) == Some(vv);
                                assert(0 <= k2 < k_new && unit_cell(s, u as int, k2) == Some(vv));
                            }
                            if exists|k2: int|
                                0 <= k2 < k_new && unit_cell(s, u as int, k2) == Some(vv) {
                                let k2 = choose|k2: int|
                                    0 <= k2 < k_new && unit_cell(s, u as int, k2) == Some(vv);
                                assert(k2 < k);
                            }
                        }
                    }
                }
            },
            None => {
                proof {
                    assert(unit_cell(s, u as int, k as int) is None);
                }
            },
        }
        k += 1;
    }
    false
}

/// Mirrors the duplicate scan of `has_contradiction`, proven exact:
/// returns true iff the board violates the validity spec.
pub fn has_duplicates(cells: &Board) -> (dup: bool)
    requires
        wf(cells@),
    ensures
        dup == !valid(cells@),
{
    let mut u: usize = 0;
    while u < 27
        invariant
            u <= 27,
            wf(cells@),
            forall|uu: int| 0 <= uu < u ==> #[trigger] unit_dupfree(cells@, uu),
        decreases 27 - u,
    {
        if unit_has_dup(cells, u) {
            proof {
                let s = cells@;
                let (k1, k2) = choose|k1: int, k2: int|
                    0 <= k1 < k2 < 9 && (#[trigger] unit_cell(s, u as int, k1)) is Some
                        && unit_cell(s, u as int, k1) == #[trigger] unit_cell(s, u as int, k2);
                lemma_unit_dup_invalid(s, u as int, k1, k2);
            }
            return true;
        }
        u += 1;
    }
    proof {
        lemma_all_units_dupfree_valid(cells@);
    }
    false
}

// ==================== Backtracking solver, proven sound ====================

pub open spec fn filled_from(s: Seq<Option<u8>>, i: int) -> bool {
    forall|j: int| i <= j < 81 ==> (#[trigger] s[j]) is Some
}

/// Backtracking over cells `i..81` (the kernel of `solve_recursive`).
/// SOUNDNESS: a `true` return means the board is valid and filled from `i`,
/// with all previously-filled cells (the givens among them) untouched.
/// A `false` return restores the board exactly.
pub fn solve_from(b: &mut Board, i: usize) -> (solved: bool)
    requires
        wf(old(b)@),
        valid(old(b)@),
        i <= 81,
    ensures
        wf(final(b)@),
        extends(old(b)@, final(b)@),
        forall|j: int| 0 <= j < i ==> final(b)@[j] == old(b)@[j],
        solved ==> valid(final(b)@) && filled_from(final(b)@, i as int),
        !solved ==> final(b)@ == old(b)@,
    decreases 81 - i,
{
    if i == 81 {
        return true;
    }
    match b[i] {
        Some(_) => solve_from(b, i + 1),
        None => {
            let ghost s0 = b@;
            let r = i / 9;
            let c = i % 9;
            proof {
                assert(i as int == (i as int / 9) * 9 + (i as int % 9));
            }
            let mut v: u8 = 1;
            while v <= 9
                invariant
                    wf(b@),
                    valid(b@),
                    b@ == s0,
                    s0 == old(b)@,
                    i < 81,
                    r == i / 9,
                    c == i % 9,
                    r < 9,
                    c < 9,
                    i as int == r as int * 9 + c as int,
                    s0[i as int] is None,
                    1 <= v <= 10,
                decreases 10 - v,
            {
                if is_legal(b, r, c, v) {
                    b[i] = Some(v);
                    proof {
                        lemma_place_preserves_valid(s0, r as int, c as int, v);
                        assert(b@ == s0.update(i as int, Some(v)));
                        assert forall|j: int| 0 <= j < 81 implies range_ok(#[trigger] b@[j]) by {
                            if j != i {
                                assert(b@[j] == s0[j]);
                            }
                        }
                    }
                    let solved = solve_from(b, i + 1);
                    if solved {
                        proof {
                            let sfin = b@;
                            let smid = s0.update(i as int, Some(v));
                            assert forall|j: int|
                                0 <= j < 81 && (#[trigger] s0[j]) is Some implies sfin[j] == s0[j]
                            by {
                                assert(j != i);
                                assert(smid[j] == s0[j]);
                            }
                            assert(smid[i as int] is Some);
                            assert(sfin[i as int] == smid[i as int]);
                            assert forall|j: int| 0 <= j < i implies sfin[j] == s0[j] by {
                                assert(smid[j] == s0[j]);
                            }
                        }
                        return true;
                    }
                    b[i] = None;
                    proof {
                        assert(b@ =~= s0);
                    }
                }
                v += 1;
            }
            false
        },
    }
}

/// Top-level solver. SOUNDNESS THEOREM: `true` means the final board is a
/// complete, valid sudoku solution that preserves every given of the input.
pub fn solve(b: &mut Board) -> (solved: bool)
    requires
        wf(old(b)@),
        valid(old(b)@),
    ensures
        wf(final(b)@),
        solved ==> valid(final(b)@) && complete(final(b)@) && extends(old(b)@, final(b)@),
        !solved ==> final(b)@ == old(b)@,
{
    solve_from(b, 0)
}

} // verus!
