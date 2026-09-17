//! Verus-verified model of `src/bitset.rs`.
//!
//! Every operation of the production `BitSet` is mirrored here with a formal
//! contract, proven against a mathematical model of the set:
//!
//! - `bit(bits, v)` is the membership model (bit `v` of the `u16`).
//! - `popcount(bits, n)` is the mathematical cardinality of the low `n` bits.
//!
//! Verified properties include:
//! - every shift is in bounds (the production code will panic in debug builds
//!   on `insert(16)` etc.; here that call is a *compile-time* error),
//! - `insert`/`remove`/`toggle`/`contains` agree with the set model,
//! - `union`/`intersection`/`difference` obey set algebra,
//! - `count` equals the model cardinality,
//! - `single_value` returns `Some(v)` iff the set is exactly `{v}`,
//! - `from_slice` builds exactly the set of the slice's elements.
//!
//! Verify with: `scripts/verify.sh` (runs `verus --crate-type=lib` on this file).

use vstd::prelude::*;

verus! {

/// Membership model: is bit `v` set in `bits`? Values 16 and above are never members.
pub open spec fn bit(bits: u16, v: u8) -> bool {
    v < 16 && ((bits >> (v as u16)) & 1) == 1
}

/// Cardinality model: number of set bits among the low `n` bits.
pub open spec fn popcount(bits: u16, n: nat) -> nat
    decreases n,
{
    if n == 0 {
        0
    } else {
        popcount(bits, (n - 1) as nat) + if bit(bits, (n - 1) as u8) {
            1nat
        } else {
            0nat
        }
    }
}

/// If no bit below `n` is set, the popcount of the low `n` bits is zero.
proof fn lemma_popcount_zero(bits: u16, n: nat)
    requires
        n <= 16,
        forall|v: u8| v < n ==> !bit(bits, v),
    ensures
        popcount(bits, n) == 0,
    decreases n,
{
    if n > 0 {
        assert(!bit(bits, (n - 1) as u8));
        lemma_popcount_zero(bits, (n - 1) as nat);
    }
}

/// One known set bit below `n` forces popcount of the low `n` bits >= 1.
proof fn lemma_popcount_one(bits: u16, u: u8, n: nat)
    requires
        u < n,
        n <= 16,
        bit(bits, u),
    ensures
        popcount(bits, n) >= 1,
    decreases n,
{
    if u < n - 1 {
        lemma_popcount_one(bits, u, (n - 1) as nat);
    }
}

/// Two distinct known set bits below `n` force popcount of the low `n` bits >= 2.
proof fn lemma_popcount_two(bits: u16, u: u8, w: u8, n: nat)
    requires
        u < w,
        w < n,
        n <= 16,
        bit(bits, u),
        bit(bits, w),
    ensures
        popcount(bits, n) >= 2,
    decreases n,
{
    if w == n - 1 {
        lemma_popcount_one(bits, u, (n - 1) as nat);
    } else {
        lemma_popcount_two(bits, u, w, (n - 1) as nat);
    }
}

/// A nonzero word has at least one set bit.
proof fn lemma_nonzero_has_bit(b: u16)
    requires
        b != 0,
    ensures
        exists|v: u8| #[trigger] bit(b, v),
{
    assert(b == 0 || ((b >> 0u16) & 1) == 1 || ((b >> 1u16) & 1) == 1 || ((b >> 2u16) & 1) == 1
        || ((b >> 3u16) & 1) == 1 || ((b >> 4u16) & 1) == 1 || ((b >> 5u16) & 1) == 1
        || ((b >> 6u16) & 1) == 1 || ((b >> 7u16) & 1) == 1 || ((b >> 8u16) & 1) == 1
        || ((b >> 9u16) & 1) == 1 || ((b >> 10u16) & 1) == 1 || ((b >> 11u16) & 1) == 1
        || ((b >> 12u16) & 1) == 1 || ((b >> 13u16) & 1) == 1 || ((b >> 14u16) & 1) == 1
        || ((b >> 15u16) & 1) == 1) by (bit_vector);
    assert(bit(b, 0) || bit(b, 1) || bit(b, 2) || bit(b, 3) || bit(b, 4) || bit(b, 5)
        || bit(b, 6) || bit(b, 7) || bit(b, 8) || bit(b, 9) || bit(b, 10) || bit(b, 11)
        || bit(b, 12) || bit(b, 13) || bit(b, 14) || bit(b, 15));
}

/// A compact bitset for candidate values, mirroring `src/bitset.rs`.
/// Bit `i` of the `u16` represents value `i` (1-indexed for sudoku candidates,
/// so bit 1 = value 1; bit 0 is unused by convention).
#[derive(Clone, Copy)]
pub struct BitSet(pub u16);

impl BitSet {
    /// Membership in the abstract set.
    pub open spec fn has(self, v: u8) -> bool {
        bit(self.0, v)
    }

    /// Abstract cardinality.
    pub open spec fn card(self) -> nat {
        popcount(self.0, 16)
    }

    /// Create an empty bitset.
    pub fn empty() -> (r: Self)
        ensures
            forall|v: u8| !r.has(v),
            r.card() == 0,
    {
        let r = BitSet(0);
        proof {
            assert forall|v: u8| !#[trigger] bit(0u16, v) by {
                assert(v < 16 ==> ((0u16 >> (v as u16)) & 1) != 1) by (bit_vector);
            }
            lemma_popcount_zero(0, 16);
        }
        r
    }

    /// Create a bitset with all values 1-9 set.
    pub fn all_9() -> (r: Self)
        ensures
            forall|v: u8| #[trigger] r.has(v) <==> (1 <= v && v <= 9),
    {
        let r = BitSet(0b1111111110);
        proof {
            assert forall|v: u8| #[trigger] r.has(v) <==> (1 <= v && v <= 9) by {
                assert(v < 16 ==> ((((0b1111111110u16) >> (v as u16)) & 1) == 1) == (1 <= v && v
                    <= 9)) by (bit_vector);
            }
        }
        r
    }

    /// Create a bitset with all values 1-n set.
    /// The production code `assert!`s at runtime; here the bound is a static precondition.
    pub fn all(n: u8) -> (r: Self)
        requires
            n < 16,
        ensures
            forall|v: u8| #[trigger] r.has(v) <==> (1 <= v && v <= n),
    {
        let ones = 1u16 << (n as u16);
        proof {
            assert(n < 16 ==> (1u16 << (n as u16)) >= 1) by (bit_vector);
        }
        let raw = (ones - 1) << 1u16;
        let r = BitSet(raw);
        proof {
            assert forall|v: u8| #[trigger] r.has(v) <==> (1 <= v && v <= n) by {
                assert((n < 16 && v < 16) ==> ((((sub(1u16 << (n as u16), 1) << 1u16) >> (
                v as u16)) & 1) == 1) == (1 <= v && v <= n)) by (bit_vector);
                assert(raw == sub(1u16 << (n as u16), 1) << 1u16);
            }
        }
        r
    }

    /// Create a bitset with a single value.
    /// Production code does `1 << value` unchecked; `value >= 16` is a shift
    /// overflow there. Here it is a rejected call.
    pub fn single(value: u8) -> (r: Self)
        requires
            value < 16,
        ensures
            forall|v: u8| #[trigger] r.has(v) <==> v == value,
    {
        let r = BitSet(1u16 << (value as u16));
        proof {
            assert forall|v: u8| #[trigger] r.has(v) <==> v == value by {
                assert((value < 16 && v < 16) ==> ((((1u16 << (value as u16)) >> (v as u16)) & 1)
                    == 1) == (v == value)) by (bit_vector);
            }
        }
        r
    }

    /// Check if the bitset contains a value.
    pub fn contains(&self, value: u8) -> (r: bool)
        requires
            value < 16,
        ensures
            r == self.has(value),
    {
        (self.0 >> (value as u16)) & 1 == 1
    }

    /// Insert a value into the bitset.
    pub fn insert(&mut self, value: u8)
        requires
            value < 16,
        ensures
            forall|v: u8| #[trigger] final(self).has(v) == (v == value || old(self).has(v)),
    {
        let ghost ob = self.0;
        self.0 = self.0 | (1u16 << (value as u16));
        proof {
            assert forall|v: u8| #[trigger] self.has(v) == (v == value || old(self).has(v)) by {
                assert((value < 16 && v < 16) ==> ((((ob | (1u16 << (value as u16))) >> (
                v as u16)) & 1) == 1) == (v == value || (((ob >> (v as u16)) & 1) == 1)))
                    by (bit_vector);
            }
        }
    }

    /// Remove a value from the bitset.
    pub fn remove(&mut self, value: u8)
        requires
            value < 16,
        ensures
            forall|v: u8| #[trigger] final(self).has(v) == (v != value && old(self).has(v)),
    {
        let ghost ob = self.0;
        self.0 = self.0 & !(1u16 << (value as u16));
        proof {
            assert forall|v: u8| #[trigger] self.has(v) == (v != value && old(self).has(v)) by {
                assert((value < 16 && v < 16) ==> ((((ob & !(1u16 << (value as u16))) >> (
                v as u16)) & 1) == 1) == (v != value && (((ob >> (v as u16)) & 1) == 1)))
                    by (bit_vector);
            }
        }
    }

    /// Toggle a value in the bitset.
    pub fn toggle(&mut self, value: u8)
        requires
            value < 16,
        ensures
            forall|v: u8|
                #[trigger] final(self).has(v) == if v == value {
                    !old(self).has(v)
                } else {
                    old(self).has(v)
                },
    {
        let ghost ob = self.0;
        self.0 = self.0 ^ (1u16 << (value as u16));
        proof {
            assert forall|v: u8|
                #[trigger] self.has(v) == if v == value {
                    !old(self).has(v)
                } else {
                    old(self).has(v)
                } by {
                assert((value < 16 && v < 16) ==> ((((ob ^ (1u16 << (value as u16))) >> (
                v as u16)) & 1) == 1) == (if v == value {
                    !(((ob >> (v as u16)) & 1) == 1)
                } else {
                    ((ob >> (v as u16)) & 1) == 1
                })) by (bit_vector);
            }
        }
    }

    /// Get the number of values in the bitset.
    /// (The production `count_ones` intrinsic is replaced by a verified scan;
    /// both compute the popcount.)
    pub fn count(&self) -> (r: u32)
        ensures
            r as nat == self.card(),
    {
        let mut acc: u32 = 0;
        let mut v: u8 = 0;
        while v < 16
            invariant
                v <= 16,
                acc <= v,
                acc as nat == popcount(self.0, v as nat),
            decreases 16 - v,
        {
            if (self.0 >> (v as u16)) & 1 == 1 {
                acc = acc + 1;
            }
            v = v + 1;
        }
        acc
    }

    /// Check if the bitset is empty.
    pub fn is_empty(&self) -> (r: bool)
        ensures
            r <==> self.card() == 0,
    {
        let r = self.0 == 0;
        proof {
            let b = self.0;
            if b == 0 {
                assert forall|v: u8| !#[trigger] bit(self.0, v) by {
                    assert(b == 0 ==> (v < 16 ==> ((b >> (v as u16)) & 1) != 1)) by (bit_vector);
                }
                lemma_popcount_zero(self.0, 16);
            } else {
                lemma_nonzero_has_bit(self.0);
                let v = choose|v: u8| bit(self.0, v);
                lemma_popcount_one(self.0, v, 16);
            }
        }
        r
    }

    /// Get the union of two bitsets.
    pub fn union(&self, other: &BitSet) -> (r: BitSet)
        ensures
            forall|v: u8| #[trigger] r.has(v) == (self.has(v) || other.has(v)),
    {
        let r = BitSet(self.0 | other.0);
        proof {
            let a = self.0;
            let b = other.0;
            assert forall|v: u8| #[trigger] r.has(v) == (self.has(v) || other.has(v)) by {
                assert(v < 16 ==> ((((a | b) >> (v as u16)) & 1) == 1) == ((((a >> (v as u16)) & 1)
                    == 1) || (((b >> (v as u16)) & 1) == 1))) by (bit_vector);
            }
        }
        r
    }

    /// Get the intersection of two bitsets.
    pub fn intersection(&self, other: &BitSet) -> (r: BitSet)
        ensures
            forall|v: u8| #[trigger] r.has(v) == (self.has(v) && other.has(v)),
    {
        let r = BitSet(self.0 & other.0);
        proof {
            let a = self.0;
            let b = other.0;
            assert forall|v: u8| #[trigger] r.has(v) == (self.has(v) && other.has(v)) by {
                assert(v < 16 ==> ((((a & b) >> (v as u16)) & 1) == 1) == ((((a >> (v as u16)) & 1)
                    == 1) && (((b >> (v as u16)) & 1) == 1))) by (bit_vector);
            }
        }
        r
    }

    /// Get the difference (self - other).
    pub fn difference(&self, other: &BitSet) -> (r: BitSet)
        ensures
            forall|v: u8| #[trigger] r.has(v) == (self.has(v) && !other.has(v)),
    {
        let r = BitSet(self.0 & !other.0);
        proof {
            let a = self.0;
            let b = other.0;
            assert forall|v: u8| #[trigger] r.has(v) == (self.has(v) && !other.has(v)) by {
                assert(v < 16 ==> ((((a & !b) >> (v as u16)) & 1) == 1) == ((((a >> (v as u16)) & 1)
                    == 1) && !(((b >> (v as u16)) & 1) == 1))) by (bit_vector);
            }
        }
        r
    }

    /// Get the only value if the bitset has exactly one element.
    /// Full functional correctness: `Some(v)` iff the set is exactly `{v}`.
    pub fn single_value(&self) -> (r: Option<u8>)
        ensures
            match r {
                Some(v) => self.has(v) && self.card() == 1 && (forall|u: u8|
                    #[trigger] self.has(u) ==> u == v),
                None => self.card() != 1,
            },
    {
        let c = self.count();
        if c != 1 {
            return None;
        }
        let mut v: u8 = 0;
        while v < 16
            invariant
                v <= 16,
                self.card() == 1,
                forall|u: u8| u < v ==> !self.has(u),
            decreases 16 - v,
        {
            if (self.0 >> (v as u16)) & 1 == 1 {
                proof {
                    assert(self.has(v));
                    assert forall|u: u8| #[trigger] self.has(u) implies u == v by {
                        if u != v {
                            if u < v {
                                assert(!self.has(u));
                            } else {
                                lemma_popcount_two(self.0, v, u, 16);
                            }
                        }
                    }
                }
                return Some(v);
            }
            v = v + 1;
        }
        proof {
            assert forall|u: u8| !#[trigger] bit(self.0, u) by {
                if u < 16 {
                    assert(!self.has(u));
                }
            }
            lemma_popcount_zero(self.0, 16);
            assert(false);
        }
        None
    }

    /// Create from a slice of values.
    /// Production code panics (debug) on any element >= 16; here the slice
    /// contents are a static precondition and the result is proven to be
    /// exactly the set of the slice's elements.
    pub fn from_slice(values: &[u8]) -> (r: Self)
        requires
            forall|i: int| 0 <= i < values@.len() ==> values@[i] < 16,
        ensures
            forall|v: u8| #[trigger] r.has(v) <==> values@.contains(v),
    {
        let mut set = BitSet::empty();
        let mut k: usize = 0;
        while k < values.len()
            invariant
                k <= values@.len(),
                forall|i: int| 0 <= i < values@.len() ==> values@[i] < 16,
                forall|v: u8|
                    #[trigger] set.has(v) <==> exists|j: int| 0 <= j < k && values@[j] == v,
            decreases values@.len() - k,
        {
            let ghost oldset = set;
            let val = values[k];
            set.insert(val);
            k = k + 1;
            proof {
                assert forall|v: u8|
                    #[trigger] set.has(v) <==> exists|j: int|
                        0 <= j < k && values@[j] == v by {
                    if set.has(v) {
                        if v == val {
                            assert(0 <= k - 1 < k && values@[k - 1] == v);
                        } else {
                            assert(oldset.has(v));
                            let j = choose|j: int| 0 <= j < k - 1 && values@[j] == v;
                            assert(0 <= j < k && values@[j] == v);
                        }
                    } else {
                        assert forall|j: int| 0 <= j < k implies values@[j] != v by {
                            if j == k - 1 {
                                assert(set.has(val));
                            } else {
                                if values@[j] == v {
                                    assert(oldset.has(v));
                                    assert(set.has(v));
                                }
                            }
                        }
                    }
                }
            }
        }
        proof {
            assert forall|v: u8| #[trigger] set.has(v) <==> values@.contains(v) by {
                if set.has(v) {
                    let j = choose|j: int| 0 <= j < k && values@[j] == v;
                    assert(values@[j] == v);
                } else {
                    if values@.contains(v) {
                        let j = choose|i: int| 0 <= i < values@.len() && values@[i] == v;
                        assert(0 <= j < k && values@[j] == v);
                        assert(set.has(v));
                    }
                }
            }
        }
        set
    }

    /// Get the raw u16 representation.
    pub fn as_raw(&self) -> (r: u16)
        ensures
            r == self.0,
    {
        self.0
    }

    /// Create from a raw u16 value.
    pub fn from_raw(raw: u16) -> (r: Self)
        ensures
            r.0 == raw,
    {
        Self(raw)
    }
}

/// End-to-end client demo: the verifier checks these facts statically —
/// there is no runtime assertion here at all.
fn demo() {
    let mut s = BitSet::empty();
    s.insert(3);
    assert(s.has(3));
    assert(!s.has(4));
    s.insert(7);
    s.remove(3);
    assert(!s.has(3));
    assert(s.has(7));

    let nine = BitSet::all_9();
    assert(nine.has(9));
    assert(!nine.has(0));
    assert(!nine.has(10));

    let only5 = BitSet::single(5);
    let both = nine.intersection(&only5);
    assert(both.has(5));
    assert(!both.has(6));

    // s.insert(20);          // ← REJECTED at verification time: `value < 16` fails.
    // BitSet::single(16);    // ← REJECTED: shift overflow is impossible by contract.
}

} // verus!
