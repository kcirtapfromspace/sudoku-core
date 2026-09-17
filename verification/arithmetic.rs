//! Verus-verified soundness of the Arithmetic Counting certificate checker
//! (Tier 3) — the machine-checked version of the theorem stated in
//! `docs/arithmetic-counting.md`.
//!
//! The production checker (`src/solver/arithmetic.rs`) reconstructs weighted
//! exactly-one equations into a coefficient vector `a` and a right-hand side
//! `rhs`, then reports a contradiction for the assumption "target variable
//! takes value `assumed`" using one of two terminals:
//!
//! - `contradiction` (IntervalGcd): residual outside `[L, U]`, or gcd of the
//!   remaining coefficients does not divide the residual;
//! - `residue_contradiction` (Residue): the residual's remainder mod `q` is
//!   not in the bitset of reachable subset-sum remainders.
//!
//! This file mirrors those functions (`gcd`, `contradiction`, `check_bounds`,
//! `add_residue_coefficient`, `reachable_residues`, `residue_contradiction`)
//! and proves the SOUNDNESS THEOREM for each terminal:
//!
//!     a reported contradiction really is one — no Boolean assignment `x`
//!     with `x[target] == assumed` satisfies  Σ a_j·x_j == rhs.
//!
//! Since every completion of the candidate state satisfies the reconstructed
//! equations (they are sums of exactly-one requirements), a verified
//! contradiction proves the advertised placement or elimination in every
//! completion — with no uniqueness assumption, exactly as the doc argues.
//! The equation *reconstruction* from the grid snapshot (sudoku-side) is the
//! remaining unverified step, documented as the trusted interface.
//!
//! Mirror deviations, each value-preserving and noted inline: Rust's signed
//! `%`/`rem_euclid` have no Verus spec, so remainders are computed through
//! provably-in-range unsigned casts; the residue fold's early exit is a loop
//! condition instead of a `break`; iterator pipelines are index loops.
//!
//! Verify with: `scripts/verify.sh`.

use vstd::arithmetic::div_mod::{
    lemma_add_mod_noop, lemma_fundamental_div_mod, lemma_mod_multiples_vanish, lemma_small_mod,
};
use vstd::prelude::*;

verus! {

// ==================== Divisibility ====================

pub open spec fn divides(d: int, n: int) -> bool {
    exists|k: int| n == #[trigger] (d * k)
}

proof fn lemma_divides_refl(a: int)
    ensures
        divides(a, a),
        divides(a, 0),
{
    assert(a == a * 1) by (nonlinear_arith);
    assert(0 == a * 0) by (nonlinear_arith);
}

proof fn lemma_divides_add(d: int, a: int, b: int)
    requires
        divides(d, a),
        divides(d, b),
    ensures
        divides(d, a + b),
{
    let k1 = choose|k: int| a == #[trigger] (d * k);
    let k2 = choose|k: int| b == #[trigger] (d * k);
    assert(a + b == d * (k1 + k2)) by (nonlinear_arith)
        requires
            a == d * k1,
            b == d * k2,
    ;
}

proof fn lemma_divides_mul(d: int, a: int, m: int)
    requires
        divides(d, a),
    ensures
        divides(d, a * m),
{
    let k = choose|k: int| a == #[trigger] (d * k);
    assert(a * m == d * (k * m)) by (nonlinear_arith)
        requires
            a == d * k,
    ;
}

proof fn lemma_divides_neg(d: int, a: int)
    requires
        divides(d, a),
    ensures
        divides(d, -a),
{
    let k = choose|k: int| a == #[trigger] (d * k);
    assert(-a == d * (-k)) by (nonlinear_arith)
        requires
            a == d * k,
    ;
}

/// "Zero divides only zero" — the g == 0 case of the theorem.
proof fn lemma_divides_zero_only(a: int)
    requires
        divides(0, a),
    ensures
        a == 0,
{
    let k = choose|k: int| a == #[trigger] (0 * k);
    assert(0 * k == 0) by (nonlinear_arith);
}

proof fn lemma_divides_iff_mod_zero(d: int, n: int)
    requires
        d > 0,
    ensures
        divides(d, n) <==> n % d == 0,
{
    if divides(d, n) {
        let k = choose|k: int| n == #[trigger] (d * k);
        lemma_mod_multiples_vanish(k, 0, d);
        assert(d * k + 0 == n);
        lemma_small_mod(0, d as nat);
    }
    if n % d == 0 {
        lemma_fundamental_div_mod(n, d);
        assert(n == d * (n / d));
    }
}

// ==================== Boolean-weighted sums ====================

/// Σ_{j < n} coeffs[j]·x[j] — the weighted combination a·x from the doc.
pub open spec fn bsum(coeffs: Seq<i16>, x: Seq<bool>, n: int) -> int
    decreases n,
{
    if n <= 0 {
        0
    } else {
        bsum(coeffs, x, n - 1) + if x[n - 1] {
            coeffs[n - 1] as int
        } else {
            0
        }
    }
}

/// The same sum with the target variable's term removed.
pub open spec fn bsum_excl(coeffs: Seq<i16>, x: Seq<bool>, t: int, n: int) -> int
    decreases n,
{
    if n <= 0 {
        0
    } else {
        bsum_excl(coeffs, x, t, n - 1) + if n - 1 != t && x[n - 1] {
            coeffs[n - 1] as int
        } else {
            0
        }
    }
}

/// Below the excluded index the two sums agree.
proof fn lemma_bsum_excl_above(coeffs: Seq<i16>, x: Seq<bool>, t: int, n: int)
    requires
        t >= n,
    ensures
        bsum_excl(coeffs, x, t, n) == bsum(coeffs, x, n),
    decreases n,
{
    if n > 0 {
        lemma_bsum_excl_above(coeffs, x, t, n - 1);
    }
}

proof fn lemma_bsum_split(coeffs: Seq<i16>, x: Seq<bool>, t: int, n: int)
    requires
        0 <= t < n,
    ensures
        bsum(coeffs, x, n) == bsum_excl(coeffs, x, t, n) + if x[t] {
            coeffs[t] as int
        } else {
            0
        },
    decreases n,
{
    if t < n - 1 {
        lemma_bsum_split(coeffs, x, t, n - 1);
    } else {
        lemma_bsum_excl_above(coeffs, x, t, n - 1);
    }
}

/// L = Σ_{j<n, j≠t} min(0, a_j) — the doc's lower endpoint.
pub open spec fn spec_lo(coeffs: Seq<i16>, t: int, n: int) -> int
    decreases n,
{
    if n <= 0 {
        0
    } else {
        spec_lo(coeffs, t, n - 1) + if n - 1 != t && coeffs[n - 1] < 0 {
            coeffs[n - 1] as int
        } else {
            0
        }
    }
}

/// U = Σ_{j<n, j≠t} max(0, a_j) — the doc's upper endpoint.
pub open spec fn spec_hi(coeffs: Seq<i16>, t: int, n: int) -> int
    decreases n,
{
    if n <= 0 {
        0
    } else {
        spec_hi(coeffs, t, n - 1) + if n - 1 != t && coeffs[n - 1] > 0 {
            coeffs[n - 1] as int
        } else {
            0
        }
    }
}

/// Interval argument: each Boolean contributes zero or its coefficient, so
/// every assignment's remaining sum lies in [L, U].
proof fn lemma_interval(coeffs: Seq<i16>, x: Seq<bool>, t: int, n: int)
    ensures
        spec_lo(coeffs, t, n) <= bsum_excl(coeffs, x, t, n) <= spec_hi(coeffs, t, n),
    decreases n,
{
    if n > 0 {
        lemma_interval(coeffs, x, t, n - 1);
    }
}

/// Divisibility argument: g | a_j for all remaining j forces g | Σ a_j·x_j.
proof fn lemma_divides_bsum_excl(d: int, coeffs: Seq<i16>, x: Seq<bool>, t: int, n: int)
    requires
        forall|j: int| 0 <= j < n && j != t ==> divides(d, #[trigger] coeffs[j] as int),
    ensures
        divides(d, bsum_excl(coeffs, x, t, n)),
    decreases n,
{
    if n <= 0 {
        lemma_divides_refl(d);
    } else {
        lemma_divides_bsum_excl(d, coeffs, x, t, n - 1);
        if n - 1 != t && x[n - 1] {
            lemma_divides_add(d, bsum_excl(coeffs, x, t, n - 1), coeffs[n - 1] as int);
        } else {
            lemma_divides_refl(d);
            lemma_divides_add(d, bsum_excl(coeffs, x, t, n - 1), 0);
        }
    }
}

// ==================== gcd (mirrors src/solver/arithmetic.rs) ====================

/// Euclid's algorithm as in the production checker. The `%` runs through u32
/// (identical value for the nonnegative operands the checker supplies).
fn gcd(a: i32, b: i32) -> (g: i32)
    requires
        a >= 0,
        b >= 0,
    ensures
        g >= 0,
        divides(g as int, a as int),
        divides(g as int, b as int),
{
    let ghost a0 = a as int;
    let ghost b0 = b as int;
    let mut a = a;
    let mut b = b;
    proof {
        lemma_divides_refl(a0);
        lemma_divides_refl(b0);
    }
    while b != 0
        invariant
            a >= 0,
            b >= 0,
            forall|d: int|
                #[trigger] divides(d, a as int) && divides(d, b as int) ==> divides(d, a0)
                    && divides(d, b0),
        decreases b,
    {
        let r = ((a as u32) % (b as u32)) as i32;
        proof {
            lemma_fundamental_div_mod(a as int, b as int);
            assert forall|d: int|
                #[trigger] divides(d, b as int) && divides(d, r as int) implies divides(d, a0)
                && divides(d, b0) by {
                lemma_divides_mul(d, b as int, a as int / b as int);
                lemma_divides_add(d, (b as int) * (a as int / b as int), r as int);
                assert(divides(d, a as int));
            }
        }
        a = b;
        b = r;
    }
    proof {
        lemma_divides_refl(a as int);
    }
    a
}

// ==================== IntervalGcd terminal ====================

/// Mirrors `ArithmeticCheck` (the reachable-residues listing is omitted; it is
/// presentation data, not part of the contradiction argument).
pub enum Check {
    Interval { residual: i32, lower: i32, upper: i32 },
    Divisibility { residual: i32, divisor: i32 },
    Residue { residual: i32, modulus: u8, required_residue: u8 },
}

/// Mirrors `check_bounds`: reports which contradiction the numbers exhibit.
/// The divisibility test runs `|residual| % divisor` through u32 — the same
/// boolean as the production `residual % divisor != 0`.
fn check_bounds(residual: i32, lower: i32, upper: i32, divisor: i32) -> (r: Option<Check>)
    requires
        divisor >= 0,
        residual > i32::MIN,
    ensures
        r is Some ==> residual < lower || residual > upper || (divisor == 0 && residual != 0) || (
        divisor > 0 && !divides(divisor as int, residual as int)),
{
    if residual < lower || residual > upper {
        Some(Check::Interval { residual, lower, upper })
    } else {
        let rem_nonzero = if divisor == 0 {
            false
        } else if residual >= 0 {
            (residual as u32) % (divisor as u32) != 0
        } else {
            ((-residual) as u32) % (divisor as u32) != 0
        };
        proof {
            if divisor > 0 {
                lemma_divides_iff_mod_zero(divisor as int, residual as int);
                if residual < 0 {
                    lemma_divides_iff_mod_zero(divisor as int, -(residual as int));
                    if divides(divisor as int, residual as int) {
                        lemma_divides_neg(divisor as int, residual as int);
                    }
                    if divides(divisor as int, -(residual as int)) {
                        lemma_divides_neg(divisor as int, -(residual as int));
                    }
                }
            }
        }
        if (divisor == 0 && residual != 0) || (divisor != 0 && rem_nonzero) {
            Some(Check::Divisibility { residual, divisor })
        } else {
            None
        }
    }
}

/// Mirrors `contradiction`, proven sound: a `Some` result means NO Boolean
/// assignment with `x[target] == assumed` satisfies Σ coeffs[j]·x[j] == rhs.
/// (In the engine, `assumed = !value` — refuting it proves the placement or
/// elimination `value` in every completion.)
pub fn contradiction(coeffs: &[i16], rhs: i32, target: usize, assumed: bool) -> (r: Option<Check>)
    requires
        target < coeffs@.len(),
        coeffs@.len() <= 729,
        -100_000 <= rhs <= 100_000,
    ensures
        r is Some ==> forall|x: Seq<bool>|
            x.len() == coeffs@.len() && x[target as int] == assumed ==> #[trigger] bsum(
                coeffs@,
                x,
                coeffs@.len() as int,
            ) != rhs as int,
{
    let n = coeffs.len();
    let residual = rhs - (coeffs[target] as i32) * (if assumed {
        1
    } else {
        0
    });
    let mut lower: i32 = 0;
    let mut upper: i32 = 0;
    let mut divisor: i32 = 0;
    let mut var: usize = 0;
    proof {
        lemma_divides_refl(0);
    }
    while var < n
        invariant
            n == coeffs@.len(),
            n <= 729,
            target < n,
            var <= n,
            divisor >= 0,
            lower == spec_lo(coeffs@, target as int, var as int),
            upper == spec_hi(coeffs@, target as int, var as int),
            -32768 * (var as int) <= lower <= 0,
            0 <= upper <= 32767 * (var as int),
            forall|j: int|
                0 <= j < var && j != target ==> divides(
                    divisor as int,
                    #[trigger] coeffs@[j] as int,
                ),
        decreases n - var,
    {
        if var != target {
            let coefficient = coeffs[var] as i32;
            // production: lower += coefficient.min(0); upper += coefficient.max(0);
            if coefficient < 0 {
                lower += coefficient;
            } else {
                upper += coefficient;
            }
            let magnitude = if coefficient < 0 {
                -coefficient
            } else {
                coefficient
            };
            let old_divisor = divisor;
            divisor = gcd(old_divisor, magnitude);
            proof {
                assert forall|j: int|
                    0 <= j < var + 1 && j != target implies divides(
                        divisor as int,
                        #[trigger] coeffs@[j] as int,
                    ) by {
                    if 0 <= j < var && j != target {
                        let k = choose|k: int| coeffs@[j] as int == #[trigger] ((old_divisor as int) * k);
                        let m = choose|k: int| old_divisor as int == #[trigger] ((divisor as int) * k);
                        assert(coeffs@[j] as int == (divisor as int) * (m * k))
                            by (nonlinear_arith)
                            requires
                                coeffs@[j] as int == (old_divisor as int) * k,
                                old_divisor as int == (divisor as int) * m,
                        ;
                    }
                    if j == var {
                        if coefficient < 0 {
                            lemma_divides_neg(divisor as int, magnitude as int);
                        }
                    }
                }
            }
        }
        var += 1;
    }
    let result = check_bounds(residual, lower, upper, divisor);
    proof {
        if result is Some {
            let s = coeffs@;
            let t = target as int;
            let len = s.len() as int;
            assert forall|x: Seq<bool>|
                x.len() == s.len() && x[t] == assumed implies #[trigger] bsum(s, x, len)
                != rhs as int by {
                lemma_bsum_split(s, x, t, len);
                lemma_interval(s, x, t, len);
                if bsum(s, x, len) == rhs as int {
                    // remaining sum equals the residual...
                    assert(bsum_excl(s, x, t, len) == residual as int);
                    // ...which the terminal proved impossible.
                    if divisor == 0 {
                        if residual != 0 {
                            assert forall|j: int| 0 <= j < len && j != t implies #[trigger] s[j]
                                as int == 0 by {
                                lemma_divides_zero_only(s[j] as int);
                            }
                            lemma_divides_bsum_excl(0, s, x, t, len);
                            lemma_divides_zero_only(bsum_excl(s, x, t, len));
                        }
                    } else {
                        lemma_divides_bsum_excl(divisor as int, s, x, t, len);
                        lemma_divides_iff_mod_zero(divisor as int, residual as int);
                    }
                }
            }
        }
    }
    result
}

// ==================== Residue terminal ====================

/// Bit `r` of a residue bitset.
pub open spec fn resbit(m: u16, r: u8) -> bool {
    r < 16 && ((m >> (r as u16)) & 1) == 1
}

/// All bits at and above the modulus are clear.
/// Euclidean remainder of a signed value, replacing `i32::rem_euclid` (which
/// has no Verus spec). Same value: nonnegative remainder in [0, q).
fn euclid_rem(a: i32, q: u32) -> (r: u32)
    requires
        0 < q <= 0xFFFF,
        -1_000_000 <= a <= 1_000_000,
    ensures
        r as int == (a as int) % (q as int),
        r < q,
{
    if a >= 0 {
        (a as u32) % q
    } else {
        let m = ((-a) as u32) % q;
        proof {
            let ai = a as int;
            let qi = q as int;
            lemma_fundamental_div_mod(-ai, qi);
            let k = (-ai) / qi;
            assert(-ai == qi * k + m);
            if m == 0 {
                // a == q·(−k), so a mod q == 0
                assert(ai == qi * (-k) + 0) by (nonlinear_arith)
                    requires
                        -ai == qi * k + 0,
                ;
                lemma_mod_multiples_vanish(-k, 0, qi);
                lemma_small_mod(0, qi as nat);
            } else {
                // a == q·(−k−1) + (q − m) with 0 < q − m < q
                assert(ai == qi * (-k - 1) + (qi - m as int)) by (nonlinear_arith)
                    requires
                        -ai == qi * k + m as int,
                ;
                lemma_mod_multiples_vanish(-k - 1, qi - m as int, qi);
                lemma_small_mod((qi - m as int) as nat, qi as nat);
            }
        }
        if m == 0 {
            0
        } else {
            q - m
        }
    }
}

/// Mirrors `add_residue_coefficient`: rotating the bitset by the coefficient
/// implements choosing that Boolean variable once. Proven: the result keeps
/// every reachable residue and adds each residue shifted by the coefficient.
fn add_residue_coefficient(reachable: u16, coefficient: i16, modulus: u8) -> (out: u16)
    requires
        2 <= modulus <= 16,
    ensures
        forall|r: u8| r < modulus && resbit(reachable, r) ==> resbit(out, r),
        forall|r: u8|
            r < modulus && resbit(reachable, r) ==> resbit(
                out,
                (((r as int + coefficient as int) % (modulus as int)) as u8),
            ),
{
    let shift = euclid_rem(coefficient as i32, modulus as u32);
    let bits = reachable as u32;
    proof {
        assert(2 <= modulus <= 16 ==> (1u32 << (modulus as u32)) >= 1) by (bit_vector);
    }
    let mask = (1u32 << (modulus as u32)) - 1;
    let rot = (bits << shift | bits >> (modulus as u32 - shift)) & mask;
    proof {
        // the OR of two 16-bit-bounded words stays 16-bit, so the cast is exact
        assert(2 <= modulus as u32 <= 16 ==> sub(1u32 << (modulus as u32), 1) <= 0xFFFF)
            by (bit_vector);
        assert(forall|x: u32, m: u32| (#[trigger] (x & m)) <= m) by (bit_vector);
        assert(forall|x: u32, y: u32|
            x <= 0xFFFF && y <= 0xFFFF ==> (#[trigger] (x | y)) <= 0xFFFF) by (bit_vector);
        assert(mask == sub(1u32 << (modulus as u32), 1));
    }
    let out = (bits | rot) as u16;
    proof {
        let q = modulus as u32;
        let s = shift;
        // bridge the exec values to the truncating operators the bit-vector proofs use
        assert(modulus as u32 - shift == sub(q, s));
        assert(out == ((reachable as u32) | ((((reachable as u32) << s) | ((reachable as u32)
            >> sub(q, s))) & sub(1u32 << q, 1))) as u16);
        // superset: old residues survive the OR
        assert forall|r: u8| r < modulus && resbit(reachable, r) implies resbit(out, r) by {
            assert((r < 16 && ((reachable >> (r as u16)) & 1) == 1) ==> (((((reachable as u32) | ((
            ((reachable as u32) << s) | ((reachable as u32) >> sub(q, s))) & sub(1u32 << q, 1)))
                as u16) >> (r as u16)) & 1) == 1) by (bit_vector);
        }
        // rotation: each old residue r also lands on (r + coefficient) mod q
        assert forall|r: u8|
            r < modulus && resbit(reachable, r) implies resbit(
            out,
            (((r as int + coefficient as int) % (modulus as int)) as u8),
        ) by {
            let ri = r as int;
            let ci = coefficient as int;
            let qi = modulus as int;
            let si = s as int;
            // (r + c) mod q == (r + s) mod q, since s == c mod q
            lemma_add_mod_noop(ri, ci, qi);
            lemma_add_mod_noop(ri, si, qi);
            lemma_small_mod(ri as nat, qi as nat);
            lemma_small_mod((si % qi) as nat, qi as nat);
            assert((ri + ci) % qi == (ri + si) % qi);
            if ri + si < qi {
                lemma_small_mod((ri + si) as nat, qi as nat);
                assert((add(r as u32, s) < q && r < 16 && s < q && 2 <= q <= 16 && ((reachable >> (
                r as u16)) & 1) == 1) ==> (((((reachable as u32) | ((((reachable as u32) << s) | ((
                reachable as u32) >> sub(q, s))) & sub(1u32 << q, 1))) as u16) >> (add(
                    r as u32,
                    s,
                ) as u16)) & 1) == 1) by (bit_vector);
                assert(add(r as u32, s) as int == ri + si);
            } else {
                lemma_mod_multiples_vanish(1, ri + si - qi, qi);
                lemma_small_mod((ri + si - qi) as nat, qi as nat);
                assert((ri + si) % qi == ri + si - qi);
                assert((add(r as u32, s) >= q && (r as u32) < q && s < q && 2 <= q <= 16 && ((
                reachable >> (r as u16)) & 1) == 1) ==> (((((reachable as u32) | ((((reachable
                    as u32) << s) | ((reachable as u32) >> sub(q, s))) & sub(1u32 << q, 1)))
                    as u16) >> (sub(add(r as u32, s), q) as u16)) & 1) == 1) by (bit_vector);
                assert(add(r as u32, s) as int == ri + si);
                assert(sub(add(r as u32, s), q) as int == ri + si - qi);
            }
        }
    }
    out
}

/// Mirrors `reachable_residues` over the remaining coefficients (index loop
/// with a loop-condition early exit instead of the iterator + `break`; the
/// fixed point makes the output identical). Proven: the bitset covers the
/// residue of EVERY Boolean assignment's remaining sum — the completeness
/// that makes an absent residue a genuine contradiction.
fn reachable_residues(coeffs: &[i16], target: usize, modulus: u8) -> (reachable: u16)
    requires
        target < coeffs@.len(),
        2 <= modulus <= 16,
    ensures
        forall|x: Seq<bool>|
            x.len() == coeffs@.len() ==> resbit(
                reachable,
                ((#[trigger] bsum_excl(coeffs@, x, target as int, coeffs@.len() as int)) % (
                modulus as int)) as u8,
            ),
{
    let n = coeffs.len();
    let full = {
        proof {
            assert(2 <= modulus <= 16 ==> (1u32 << (modulus as u32)) >= 1) by (bit_vector);
        }
        ((1u32 << (modulus as u32)) - 1) as u16
    };
    let mut reachable: u16 = 1;
    let mut j: usize = 0;
    proof {
        assert(((1u16 >> 0u16) & 1) == 1) by (bit_vector);
        lemma_small_mod(0, modulus as nat);
    }
    while j < n && reachable != full
        invariant
            n == coeffs@.len(),
            target < n,
            2 <= modulus <= 16,
            j <= n,
            full == sub((1u32 << (modulus as u32)), 1) as u16,
            forall|x: Seq<bool>|
                x.len() == n ==> resbit(
                    reachable,
                    ((#[trigger] bsum_excl(coeffs@, x, target as int, j as int)) % (
                    modulus as int)) as u8,
                ),
        decreases n - j,
    {
        let ghost old_reachable = reachable;
        let ghost old_j = j as int;
        if j != target && coeffs[j] != 0 {
            reachable = add_residue_coefficient(reachable, coeffs[j], modulus);
        }
        j += 1;
        proof {
            let s = coeffs@;
            let t = target as int;
            let qi = modulus as int;
            assert forall|x: Seq<bool>| x.len() == n implies resbit(
                reachable,
                ((#[trigger] bsum_excl(s, x, t, old_j + 1)) % qi) as u8,
            ) by {
                let prev = bsum_excl(s, x, t, old_j);
                let prev_bit = (prev % qi) as u8;
                assert(resbit(old_reachable, prev_bit));
                assert(0 <= prev % qi < qi);
                if old_j != t && x[old_j] && s[old_j] != 0 {
                    // this variable chosen: residue advances by the coefficient
                    let ci = s[old_j] as int;
                    lemma_add_mod_noop(prev, ci, qi);
                    lemma_add_mod_noop(prev % qi, ci, qi);
                    lemma_small_mod((prev % qi) as nat, qi as nat);
                    assert((prev + ci) % qi == ((prev % qi) + ci) % qi);
                    assert(resbit(reachable, (((prev_bit as int + ci) % qi) as u8)));
                    assert(bsum_excl(s, x, t, old_j + 1) == prev + ci);
                } else {
                    // skipped (target, zero coefficient, or unchosen): sum unchanged
                    assert(bsum_excl(s, x, t, old_j + 1) == prev);
                    if j <= n && old_j != t && x[old_j] && s[old_j] == 0 {
                        assert(bsum_excl(s, x, t, old_j + 1) == prev + 0);
                    }
                }
            }
        }
    }
    proof {
        // Early exit: a full bitset covers every residue.
        if reachable == full {
            let q = modulus;
            assert forall|x: Seq<bool>| x.len() == n implies resbit(
                reachable,
                ((#[trigger] bsum_excl(coeffs@, x, target as int, n as int)) % (q as int)) as u8,
            ) by {
                let rr = ((bsum_excl(coeffs@, x, target as int, n as int)) % (q as int)) as u8;
                assert(0 <= bsum_excl(coeffs@, x, target as int, n as int) % (q as int) < q as int);
                assert((2 <= q <= 16 && rr < q) ==> (((sub((1u32 << (q as u32)), 1) as u16) >> (rr
                    as u16)) & 1) == 1) by (bit_vector);
            }
        }
    }
    reachable
}

/// Mirrors `residue_contradiction`, proven sound: a `Some` result means NO
/// Boolean assignment with `x[target] == assumed` satisfies Σ coeffs[j]·x[j]
/// == rhs. (The rejected-both-endpoints refusal is preserved as in
/// production; it only makes the checker more conservative.)
pub fn residue_contradiction(
    coeffs: &[i16],
    rhs: i32,
    target: usize,
    assumed: bool,
    modulus: u8,
) -> (r: Option<Check>)
    requires
        target < coeffs@.len(),
        coeffs@.len() <= 729,
        -100_000 <= rhs <= 100_000,
        2 <= modulus <= 16,
    ensures
        r is Some ==> forall|x: Seq<bool>|
            x.len() == coeffs@.len() && x[target as int] == assumed ==> #[trigger] bsum(
                coeffs@,
                x,
                coeffs@.len() as int,
            ) != rhs as int,
{
    let reachable = reachable_residues(coeffs, target, modulus);
    let residual = rhs - (coeffs[target] as i32) * (if assumed {
        1
    } else {
        0
    });
    let other_residual = rhs - (coeffs[target] as i32) * (if !assumed {
        1
    } else {
        0
    });
    let required_residue = euclid_rem(residual, modulus as u32) as u8;
    let other = euclid_rem(other_residual, modulus as u32) as u8;
    // Version two refuses certificates that reject both candidate endpoints.
    if reachable & (1u16 << (required_residue as u16)) != 0 || reachable & (1u16 << (other
        as u16)) == 0 {
        return None;
    }
    proof {
        let s = coeffs@;
        let t = target as int;
        let len = s.len() as int;
        let qi = modulus as int;
        assert((required_residue < 16) ==> ((reachable & (1u16 << (required_residue as u16)))
            != 0) == (((reachable >> (required_residue as u16)) & 1) == 1)) by (bit_vector);
        assert forall|x: Seq<bool>|
            x.len() == s.len() && x[t] == assumed implies #[trigger] bsum(s, x, len)
            != rhs as int by {
            lemma_bsum_split(s, x, t, len);
            if bsum(s, x, len) == rhs as int {
                assert(bsum_excl(s, x, t, len) == residual as int);
                assert(resbit(reachable, ((residual as int) % qi) as u8));
            }
        }
    }
    Some(Check::Residue { residual, modulus, required_residue })
}

// ==================== The doc's worked example, end to end ====================

/// The "five-equation parity" Guardian from docs/arithmetic-counting.md:
/// a+b=1; b+c=1; c+d=1; d+e=1; e+a+t=1 sum (unit weights) to
/// 2a+2b+2c+2d+2e+t = 5. Assuming t = 0 leaves gcd 2 ∤ 5 — so the checker
/// reports a contradiction, and the verifier confirms what that MEANS:
/// every Boolean solution of the combined equation has t = 1.
fn demo_five_equation_parity() {
    let coeffs: [i16; 6] = [2, 2, 2, 2, 2, 1];
    let cs: &[i16] = coeffs.as_slice();
    let result = contradiction(cs, 5, 5, false);
    if result.is_some() {
        proof {
            // machine-checked consequence: no assignment with t = false sums to 5
            assert(cs@.len() == 6);
            assert(forall|x: Seq<bool>|
                x.len() == 6 && x[5] == false ==> #[trigger] bsum(cs@, x, 6) != 5);
        }
    }
}

} // verus!
