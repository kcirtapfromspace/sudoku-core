//! Exactly-one partition relaxation for weighted Sudoku equations.
//!
//! The checker rebuilds all groups, including zero-coefficient alternatives.
//! Search uses prefix/suffix sums; verification recomputes both assumptions
//! directly. See docs/partition-residues.md for the soundness argument.

use super::*;

/// Each variant partitions all 729 candidate variables into 81 disjoint native
/// exactly-one requirements. Arbitrary overlapping groups are not accepted.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum ArithmeticPartition {
    Cells,
    RowDigits,
    ColumnDigits,
    BoxDigits,
}

impl ArithmeticPartition {
    pub(super) fn group(self, variable: usize) -> usize {
        let (cell, digit) = (variable / 9, variable % 9);
        match self {
            Self::Cells => cell,
            Self::RowDigits => cell / 9 * 9 + digit,
            Self::ColumnDigits => cell % 9 * 9 + digit,
            Self::BoxDigits => (cell / 27 * 3 + cell % 9 / 3) * 9 + digit,
        }
    }

    pub(super) fn label(self) -> &'static str {
        match self {
            Self::Cells => "cell",
            Self::RowDigits => "row/digit requirement",
            Self::ColumnDigits => "column/digit requirement",
            Self::BoxDigits => "box/digit requirement",
        }
    }
}

fn groups(snapshot: &Snapshot, partition: ArithmeticPartition) -> Vec<Vec<usize>> {
    let mut groups = vec![Vec::new(); 81];
    for cell in 0..81 {
        for digit in 0..9 {
            if snapshot.domains[cell] & (1 << digit) != 0 {
                let variable = cell * 9 + digit;
                groups[partition.group(variable)].push(variable);
            }
        }
    }
    groups
}

fn residue(coefficient: i16, modulus: u8) -> usize {
    i32::from(coefficient).rem_euclid(i32::from(modulus)) as usize
}

// Minkowski addition of residue sets. In particular, an empty choice set
// annihilates the result; it must never be treated as an optional zero choice.
fn sum_residues(left: u16, right: u16, modulus: u8) -> u16 {
    let mut result = 0;
    for a in 0..modulus {
        if left & (1 << a) != 0 {
            for b in 0..modulus {
                if right & (1 << b) != 0 {
                    result |= 1 << ((a + b) % modulus);
                }
            }
        }
    }
    result
}

fn assumption_residues(
    groups: &[Vec<usize>],
    coefficients: &[i16; VARIABLE_COUNT],
    target: usize,
    assumed: bool,
    partition: ArithmeticPartition,
    modulus: u8,
) -> u16 {
    let mut reachable = 1;
    for (index, group) in groups.iter().enumerate() {
        let mut choices = 0;
        for &variable in group {
            let allowed = if index == partition.group(target) {
                (variable == target) == assumed
            } else {
                true
            };
            if allowed {
                choices |= 1 << residue(coefficients[variable], modulus);
            }
        }
        // Even after saturation, a later group can become empty under x_t=0.
        reachable = sum_residues(reachable, choices, modulus);
    }
    reachable
}

pub(super) fn check(
    snapshot: &Snapshot,
    coefficients: &[i16; VARIABLE_COUNT],
    rhs: i32,
    target: usize,
    assumed: bool,
    partition: ArithmeticPartition,
    modulus: u8,
) -> Option<ArithmeticCheck> {
    let groups = groups(snapshot, partition);
    let required_residue = rhs.rem_euclid(i32::from(modulus)) as u8;
    let required_bit = 1 << required_residue;
    let rejected = assumption_residues(&groups, coefficients, target, assumed, partition, modulus);
    let accepted = assumption_residues(&groups, coefficients, target, !assumed, partition, modulus);
    // Refuse arbitrary conclusions from a relaxation rejecting both endpoints.
    if rejected & required_bit != 0 || accepted & required_bit == 0 {
        return None;
    }
    Some(ArithmeticCheck::PartitionResidue {
        rhs,
        modulus,
        partition,
        required_residue,
        reachable_residues: (0..modulus).filter(|r| rejected & (1 << r) != 0).collect(),
    })
}

pub(super) fn finding(
    node: &Node,
    snapshot: &Snapshot,
    partition: ArithmeticPartition,
    modulus: u8,
) -> Option<Finding> {
    let groups = groups(snapshot, partition);
    let choices: Vec<u16> = groups
        .iter()
        .map(|group| {
            group.iter().fold(0, |set, &var| {
                set | (1 << residue(node.coefficients[var], modulus))
            })
        })
        .collect();
    let mut prefix = [1; 82];
    let mut suffix = [1; 82];
    for i in 0..81 {
        prefix[i + 1] = sum_residues(prefix[i], choices[i], modulus);
    }
    for i in (0..81).rev() {
        suffix[i] = sum_residues(choices[i], suffix[i + 1], modulus);
    }
    let required_bit = 1 << node.rhs.rem_euclid(i32::from(modulus));
    for (index, group) in groups.iter().enumerate() {
        // With one residue every choice is indistinguishable, except removing
        // the sole candidate: the empty set can certify a forced singleton.
        if choices[index].count_ones() == 1 && group.len() != 1 {
            continue;
        }
        let others = sum_residues(prefix[index], suffix[index + 1], modulus);
        let mut counts = [0u8; MODULUS_LIMIT as usize];
        for &var in group {
            counts[residue(node.coefficients[var], modulus)] += 1;
        }
        for &var in group {
            if snapshot.placed[var / 9] != 0 {
                continue;
            }
            let r = residue(node.coefficients[var], modulus);
            let without = if counts[r] == 1 {
                choices[index] & !(1 << r)
            } else {
                choices[index]
            };
            let zero_rejected = sum_residues(others, without, modulus) & required_bit == 0;
            let one_rejected = sum_residues(others, 1 << r, modulus) & required_bit == 0;
            if zero_rejected == one_rejected {
                continue;
            }
            let proof = ArithmeticProof {
                version: PARTITION_VERSION,
                terminal: ArithmeticTerminal::PartitionResidue { modulus, partition },
                state_hash: snapshot.hash.clone(),
                terms: node
                    .terms
                    .iter()
                    .map(|&(id, weight)| ArithmeticTerm {
                        requirement: ArithmeticRequirement::from_id(id),
                        weight,
                    })
                    .collect(),
                cell: var / 9,
                digit: (var % 9 + 1) as u8,
                value: zero_rejected,
            };
            // Prefix/suffix search results never authorize a hint on their own.
            if let Some(check) = proof.check_snapshot(snapshot) {
                return Some(make_finding(snapshot, proof, check));
            }
        }
    }
    None
}

#[cfg(test)]
#[path = "arithmetic_partition_tests.rs"]
mod tests;
