use serde::{Deserialize, Serialize};

/// A compact bitset for storing candidate values (1-9 for standard Sudoku)
/// Uses a u16 internally where bit i represents value i (1-indexed, so bit 1 = value 1)
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Serialize, Deserialize)]
pub struct BitSet(u16);

impl BitSet {
    /// Create an empty bitset
    pub const fn empty() -> Self {
        Self(0)
    }

    /// Create a bitset with all values 1-9 set
    pub const fn all_9() -> Self {
        Self(0b1111111110) // bits 1-9 set
    }

    /// Create a bitset with all values 1-n set
    pub fn all(n: u8) -> Self {
        Self(((1u16 << (n + 1)) - 1) & !1) // bits 1-n set
    }

    /// Create a bitset with a single value
    pub const fn single(value: u8) -> Self {
        Self(1 << value)
    }

    /// Check if the bitset contains a value
    pub const fn contains(&self, value: u8) -> bool {
        (self.0 & (1 << value)) != 0
    }

    /// Insert a value into the bitset
    pub fn insert(&mut self, value: u8) {
        self.0 |= 1 << value;
    }

    /// Remove a value from the bitset
    pub fn remove(&mut self, value: u8) {
        self.0 &= !(1 << value);
    }

    /// Toggle a value in the bitset
    pub fn toggle(&mut self, value: u8) {
        self.0 ^= 1 << value;
    }

    /// Get the number of values in the bitset
    pub const fn count(&self) -> u32 {
        self.0.count_ones()
    }

    /// Check if the bitset is empty
    pub const fn is_empty(&self) -> bool {
        self.0 == 0
    }

    /// Get the union of two bitsets
    pub const fn union(&self, other: &BitSet) -> BitSet {
        BitSet(self.0 | other.0)
    }

    /// Get the intersection of two bitsets
    pub const fn intersection(&self, other: &BitSet) -> BitSet {
        BitSet(self.0 & other.0)
    }

    /// Get the difference (self - other)
    pub const fn difference(&self, other: &BitSet) -> BitSet {
        BitSet(self.0 & !other.0)
    }

    /// Iterate over all values in the bitset
    pub fn iter(&self) -> BitSetIter {
        BitSetIter {
            bits: self.0,
            current: 1,
        }
    }

    /// Get the only value if the bitset has exactly one element
    pub fn single_value(&self) -> Option<u8> {
        if self.count() == 1 {
            Some(self.0.trailing_zeros() as u8)
        } else {
            None
        }
    }

    /// Get all values as a Vec
    pub fn to_vec(&self) -> Vec<u8> {
        self.iter().collect()
    }

    /// Create from a slice of values
    pub fn from_slice(values: &[u8]) -> Self {
        let mut set = BitSet::empty();
        for &v in values {
            set.insert(v);
        }
        set
    }

    /// Get the raw u16 representation
    pub const fn as_raw(&self) -> u16 {
        self.0
    }

    /// Create from a raw u16 value
    pub const fn from_raw(raw: u16) -> Self {
        Self(raw)
    }
}

/// Iterator over values in a BitSet
pub struct BitSetIter {
    bits: u16,
    current: u8,
}

impl Iterator for BitSetIter {
    type Item = u8;

    fn next(&mut self) -> Option<Self::Item> {
        while self.current < 16 {
            let value = self.current;
            self.current += 1;
            if (self.bits & (1u16 << value)) != 0 {
                return Some(value);
            }
        }
        None
    }
}

impl std::fmt::Display for BitSet {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let values: Vec<String> = self.iter().map(|v| v.to_string()).collect();
        write!(f, "{{{}}}", values.join(", "))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_empty() {
        let set = BitSet::empty();
        assert!(set.is_empty());
        assert_eq!(set.count(), 0);
    }

    #[test]
    fn test_all_9() {
        let set = BitSet::all_9();
        assert_eq!(set.count(), 9);
        for i in 1..=9 {
            assert!(set.contains(i));
        }
        assert!(!set.contains(0));
        assert!(!set.contains(10));
    }

    #[test]
    fn test_insert_remove() {
        let mut set = BitSet::empty();
        set.insert(5);
        assert!(set.contains(5));
        assert_eq!(set.count(), 1);

        set.insert(3);
        assert!(set.contains(3));
        assert_eq!(set.count(), 2);

        set.remove(5);
        assert!(!set.contains(5));
        assert_eq!(set.count(), 1);
    }

    #[test]
    fn test_single_value() {
        let mut set = BitSet::empty();
        set.insert(7);
        assert_eq!(set.single_value(), Some(7));

        set.insert(3);
        assert_eq!(set.single_value(), None);
    }

    #[test]
    fn test_iter() {
        let set = BitSet::from_slice(&[1, 3, 5, 9]);
        let values: Vec<u8> = set.iter().collect();
        assert_eq!(values, vec![1, 3, 5, 9]);
    }

    #[test]
    fn test_operations() {
        let a = BitSet::from_slice(&[1, 2, 3]);
        let b = BitSet::from_slice(&[2, 3, 4]);

        let union = a.union(&b);
        assert_eq!(union.to_vec(), vec![1, 2, 3, 4]);

        let intersection = a.intersection(&b);
        assert_eq!(intersection.to_vec(), vec![2, 3]);

        let difference = a.difference(&b);
        assert_eq!(difference.to_vec(), vec![1]);
    }
}
