//! Bounded search for replayable weighted exactly-one deductions.
//!
//! Candidate masks are premises, as for the other logical engines. The checker
//! proves that every completion respecting those masks respects the deduction;
//! it neither establishes satisfiability nor authenticates earlier eliminations.

use super::explain::{ExplanationData, Finding, InferenceResult, ProofCertificate};
use super::fabric::{idx_to_pos, sector_cells, CandidateFabric};
use super::{Hint, Technique};
use crate::Grid;
use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};
use std::cmp::Ordering;
use std::collections::{BinaryHeap, HashSet};

const LEGACY_VERSION: u8 = 1;
const RESIDUE_VERSION: u8 = 2;
const SOURCE_LIMIT: usize = 6;
const WEIGHT_LIMIT: i8 = 2;
const RESIDUE_WEIGHT_LIMIT: i8 = 8;
const MODULUS_LIMIT: u8 = 16;
const VARIABLE_COUNT: usize = 729;
// Bound materialized children independently of the evaluation budget.
const FRONTIER_LIMIT: usize = 8192;

/// One of the 324 native exactly-one requirements of a classic Sudoku.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub enum ArithmeticRequirement {
    Cell {
        cell: usize,
    },
    /// Sectors 0..8 are rows, 9..17 columns, and 18..26 boxes; digits are 1..9.
    SectorDigit {
        sector: usize,
        digit: u8,
    },
}

impl ArithmeticRequirement {
    fn id(&self) -> Option<usize> {
        match *self {
            Self::Cell { cell } if cell < 81 => Some(cell),
            Self::SectorDigit { sector, digit } if sector < 27 && (1..=9).contains(&digit) => {
                Some(81 + sector * 9 + usize::from(digit - 1))
            }
            _ => None,
        }
    }

    fn from_id(id: usize) -> Self {
        if id < 81 {
            Self::Cell { cell: id }
        } else {
            Self::SectorDigit {
                sector: (id - 81) / 9,
                digit: ((id - 81) % 9 + 1) as u8,
            }
        }
    }

    fn label(&self) -> String {
        match *self {
            Self::Cell { cell } => format!("cell r{}c{}", cell / 9 + 1, cell % 9 + 1),
            Self::SectorDigit { sector, digit } => {
                let (kind, number) = if sector < 9 {
                    ("row", sector + 1)
                } else if sector < 18 {
                    ("column", sector - 8)
                } else {
                    ("box", sector - 17)
                };
                format!("digit {digit} in {kind} {number}")
            }
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct ArithmeticTerm {
    pub requirement: ArithmeticRequirement,
    pub weight: i8,
}

/// The independently reconstructed contradiction required by a certificate.
#[derive(Debug, Clone, Default, PartialEq, Eq, Serialize, Deserialize)]
pub enum ArithmeticTerminal {
    /// Version-one interval and gcd semantics, also the default for old JSON.
    #[default]
    IntervalGcd,
    /// Version-two Boolean attainable remainders, with a modulus from 2 to 16.
    Residue { modulus: u8 },
}

impl ArithmeticTerminal {
    fn is_legacy(&self) -> bool {
        matches!(self, Self::IntervalGcd)
    }
}

/// A certificate bound to the exact values and candidate masks of one grid.
///
/// Both versions support at most six distinct native requirements. Version one
/// accepts nonzero weights of magnitude at most two; version two accepts at most
/// eight to support parity-tail compression. These are proof format limits,
/// not a claim that these certificates express every Sudoku deduction.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct ArithmeticProof {
    pub version: u8,
    #[serde(default, skip_serializing_if = "ArithmeticTerminal::is_legacy")]
    pub terminal: ArithmeticTerminal,
    pub state_hash: String,
    pub terms: Vec<ArithmeticTerm>,
    pub cell: usize,
    pub digit: u8,
    /// True means placement; false means elimination.
    pub value: bool,
}

/// The contradiction obtained by assuming the opposite of a certified result.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub enum ArithmeticCheck {
    Interval {
        residual: i32,
        lower: i32,
        upper: i32,
    },
    Divisibility {
        residual: i32,
        divisor: i32,
    },
    Residue {
        residual: i32,
        modulus: u8,
        required_residue: u8,
        /// Recomputed from the source equations, never supplied by the proof.
        reachable_residues: Vec<u8>,
    },
}

impl ArithmeticProof {
    /// Construct a certificate only when its deduction can be independently checked.
    pub fn from_terms(
        grid: &Grid,
        terms: Vec<ArithmeticTerm>,
        cell: usize,
        digit: u8,
        value: bool,
    ) -> Option<Self> {
        Self::from_terms_with_terminal(
            grid,
            terms,
            cell,
            digit,
            value,
            ArithmeticTerminal::IntervalGcd,
        )
    }

    /// Construct a bounded residue certificate, independently checked on this grid.
    pub fn from_residue_terms(
        grid: &Grid,
        terms: Vec<ArithmeticTerm>,
        cell: usize,
        digit: u8,
        value: bool,
        modulus: u8,
    ) -> Option<Self> {
        Self::from_terms_with_terminal(
            grid,
            terms,
            cell,
            digit,
            value,
            ArithmeticTerminal::Residue { modulus },
        )
    }

    pub fn from_terms_with_terminal(
        grid: &Grid,
        terms: Vec<ArithmeticTerm>,
        cell: usize,
        digit: u8,
        value: bool,
        terminal: ArithmeticTerminal,
    ) -> Option<Self> {
        let snapshot = Snapshot::from_grid(grid).ok()?;
        let proof = Self {
            version: if terminal.is_legacy() {
                LEGACY_VERSION
            } else {
                RESIDUE_VERSION
            },
            terminal,
            state_hash: snapshot.hash.clone(),
            terms,
            cell,
            digit,
            value,
        };
        proof.check_snapshot(&snapshot)?;
        Some(proof)
    }

    pub fn verify(&self, grid: &Grid) -> bool {
        self.check(grid).is_some()
    }

    /// Reconstruct the source equations from this grid, then check the contradiction.
    pub fn check(&self, grid: &Grid) -> Option<ArithmeticCheck> {
        self.check_snapshot(&Snapshot::from_grid(grid).ok()?)
    }

    fn check_snapshot(&self, snapshot: &Snapshot) -> Option<ArithmeticCheck> {
        let weight_limit = match (&self.terminal, self.version) {
            (ArithmeticTerminal::IntervalGcd, LEGACY_VERSION) => WEIGHT_LIMIT,
            (ArithmeticTerminal::Residue { modulus }, RESIDUE_VERSION)
                if (2..=MODULUS_LIMIT).contains(modulus) =>
            {
                RESIDUE_WEIGHT_LIMIT
            }
            _ => return None,
        };
        if self.state_hash != snapshot.hash
            || self.cell >= 81
            || !(1..=9).contains(&self.digit)
            || snapshot.placed[self.cell] != 0
            || snapshot.domains[self.cell] & (1 << (self.digit - 1)) == 0
            || self.terms.is_empty()
            || self.terms.len() > SOURCE_LIMIT
        {
            return None;
        }
        let mut seen = [false; 324];
        let mut coefficients = [0i16; VARIABLE_COUNT];
        let mut rhs = 0;
        for term in &self.terms {
            if term.weight == 0 || !(-weight_limit..=weight_limit).contains(&term.weight) {
                return None;
            }
            let id = term.requirement.id()?;
            if seen[id] {
                return None;
            }
            seen[id] = true;
            // Rebuild from original domains, independent of search-node arithmetic.
            for var in snapshot.requirement_variables(id) {
                coefficients[var] += i16::from(term.weight);
            }
            rhs += i32::from(term.weight);
        }
        let target = self.cell * 9 + usize::from(self.digit - 1);
        match self.terminal {
            ArithmeticTerminal::IntervalGcd => {
                contradiction(&coefficients, rhs, target, !self.value)
            }
            ArithmeticTerminal::Residue { modulus } => {
                residue_contradiction(&coefficients, rhs, target, !self.value, modulus)
            }
        }
    }

    /// Compress a checked unit-weight parity placement and one distinct incident
    /// exactly-one rule into direct residue eliminations. The union must contain
    /// at most six sources; tail cardinalities 2..9 give moduli 2..16 and weights
    /// at most eight. Invalid premises or unsupported shapes return no proofs.
    pub fn compile_parity_tail(&self, grid: &Grid, tail: ArithmeticRequirement) -> Vec<Self> {
        let Ok(snapshot) = Snapshot::from_grid(grid) else {
            return Vec::new();
        };
        self.compile_tail_snapshot(&snapshot, tail)
    }

    fn compile_tail_snapshot(&self, snapshot: &Snapshot, tail: ArithmeticRequirement) -> Vec<Self> {
        if !self.value
            || self.terms.len() >= SOURCE_LIMIT
            || self.terms.iter().any(|term| !matches!(term.weight, -1 | 1))
            || !matches!(self.check_snapshot(snapshot),
                Some(ArithmeticCheck::Divisibility { residual, divisor })
                    if residual % 2 != 0 && divisor != 0 && divisor % 2 == 0)
        {
            return Vec::new();
        }
        let Some(id) = tail.id() else {
            return Vec::new();
        };
        if self.terms.iter().any(|term| term.requirement == tail) {
            return Vec::new();
        }
        let guardian = self.cell * 9 + usize::from(self.digit - 1);
        let variables = snapshot.requirement_variables(id);
        if !(2..=9).contains(&variables.len()) || !variables.contains(&guardian) {
            return Vec::new();
        }
        let scale = (variables.len() - 1) as i8;
        let mut terms: Vec<_> = self
            .terms
            .iter()
            .map(|term| ArithmeticTerm {
                requirement: term.requirement.clone(),
                weight: term.weight * scale,
            })
            .collect();
        terms.push(ArithmeticTerm {
            requirement: tail,
            weight: 1,
        });
        variables
            .into_iter()
            .filter(|&var| var != guardian)
            .filter_map(|var| {
                let proof = Self {
                    version: RESIDUE_VERSION,
                    terminal: ArithmeticTerminal::Residue {
                        modulus: 2 * scale as u8,
                    },
                    state_hash: snapshot.hash.clone(),
                    terms: terms.clone(),
                    cell: var / 9,
                    digit: (var % 9 + 1) as u8,
                    value: false,
                };
                proof.check_snapshot(snapshot).map(|_| proof)
            })
            .collect()
    }
}

/// Limits for an intentionally incomplete, deterministic connected beam search.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct ArithmeticSearchOptions {
    pub max_sources: usize,
    pub max_weight: i8,
    pub beam_width: usize,
    pub max_combinations: usize,
}

impl Default for ArithmeticSearchOptions {
    fn default() -> Self {
        Self {
            max_sources: 6,
            max_weight: 2,
            beam_width: 64,
            max_combinations: 30_000,
        }
    }
}

#[derive(Debug, Clone)]
pub struct ArithmeticSearchResult {
    pub hint: Option<Hint>,
    pub tested_combinations: usize,
    /// The combination limit prevented further evaluations (including zero budget).
    pub budget_exhausted: bool,
    /// A beam or per-depth allocation discarded unexplored combinations.
    pub beam_pruned: bool,
    /// Invalid options, unsupported constraint topology, or invalid local state.
    pub error: Option<String>,
}

struct Snapshot {
    domains: [u16; 81],
    placed: [u8; 81],
    hash: String,
}

impl Snapshot {
    fn from_grid(grid: &Grid) -> Result<Self, String> {
        if !grid.has_standard_sudoku_constraints() {
            return Err(
                "Arithmetic Counting requires installed standard Sudoku constraints".into(),
            );
        }
        let mut snapshot = Self {
            domains: [0; 81],
            placed: [0; 81],
            hash: String::new(),
        };
        for cell in 0..81 {
            let pos = idx_to_pos(cell);
            if let Some(value) = grid.get(pos) {
                if !(1..=9).contains(&value) {
                    return Err("Placed digit is outside 1..9".into());
                }
                snapshot.placed[cell] = value;
                snapshot.domains[cell] = 1 << (value - 1);
            } else {
                let raw = grid.get_candidates(pos).as_raw();
                if raw == 0 || raw & !0b11_1111_1110 != 0 {
                    return Err("Empty cell has an empty or invalid candidate mask".into());
                }
                snapshot.domains[cell] = raw >> 1;
            }
        }
        for sector in 0..27 {
            let mut placed_digits = 0u16;
            for cell in sector_cells(sector) {
                if snapshot.placed[cell] != 0 {
                    let bit = snapshot.domains[cell];
                    if placed_digits & bit != 0 {
                        return Err("Duplicate placed digit in a standard sector".into());
                    }
                    placed_digits |= bit;
                }
            }
        }
        for id in 81..324 {
            if snapshot.requirement_variables(id).is_empty() {
                return Err("A standard sector has no remaining position for a digit".into());
            }
        }
        let mut hasher = Sha256::new();
        hasher.update(b"sudoku-arithmetic-v1-classic-324");
        for cell in 0..81 {
            hasher.update(snapshot.domains[cell].to_le_bytes());
            hasher.update([snapshot.placed[cell]]);
        }
        snapshot.hash = format!("{:x}", hasher.finalize());
        Ok(snapshot)
    }

    fn requirement_variables(&self, id: usize) -> Vec<usize> {
        if id < 81 {
            (0..9)
                .filter(|digit| self.domains[id] & (1 << digit) != 0)
                .map(|digit| id * 9 + digit)
                .collect()
        } else {
            let sector = (id - 81) / 9;
            let digit = (id - 81) % 9;
            sector_cells(sector)
                .into_iter()
                .filter(|&cell| self.domains[cell] & (1 << digit) != 0)
                .map(|cell| cell * 9 + digit)
                .collect()
        }
    }
}

fn gcd(mut a: i32, mut b: i32) -> i32 {
    while b != 0 {
        (a, b) = (b, a % b);
    }
    a
}

fn contradiction(
    coefficients: &[i16; VARIABLE_COUNT],
    rhs: i32,
    target: usize,
    assumed: bool,
) -> Option<ArithmeticCheck> {
    let residual = rhs - i32::from(coefficients[target]) * i32::from(assumed);
    let (mut lower, mut upper, mut divisor) = (0, 0, 0);
    for (var, &coefficient) in coefficients.iter().enumerate() {
        if var == target {
            continue;
        }
        let coefficient = i32::from(coefficient);
        lower += coefficient.min(0);
        upper += coefficient.max(0);
        divisor = gcd(divisor, coefficient.abs());
    }
    check_bounds(residual, lower, upper, divisor)
}

fn check_bounds(residual: i32, lower: i32, upper: i32, divisor: i32) -> Option<ArithmeticCheck> {
    if residual < lower || residual > upper {
        Some(ArithmeticCheck::Interval {
            residual,
            lower,
            upper,
        })
    } else if (divisor == 0 && residual != 0) || (divisor != 0 && residual % divisor != 0) {
        Some(ArithmeticCheck::Divisibility { residual, divisor })
    } else {
        None
    }
}

// Rotating a bitset by a coefficient implements choosing that Boolean variable
// once. Use the previous set, preserve duplicate coefficients, and reduce signed
// coefficients with Euclidean remainders. u32 intermediates handle modulus 16.
fn add_residue_coefficient(reachable: u16, coefficient: i16, modulus: u8) -> u16 {
    let shift = i32::from(coefficient).rem_euclid(i32::from(modulus)) as u32;
    let bits = u32::from(reachable);
    let mask = (1u32 << modulus) - 1;
    (bits | ((bits << shift | bits >> (u32::from(modulus) - shift)) & mask)) as u16
}

fn reachable_residues(coefficients: impl Iterator<Item = i16>, modulus: u8) -> u16 {
    debug_assert!((2..=MODULUS_LIMIT).contains(&modulus));
    let full = ((1u32 << modulus) - 1) as u16;
    let mut reachable = 1;
    for coefficient in coefficients {
        reachable = add_residue_coefficient(reachable, coefficient, modulus);
        if reachable == full {
            break;
        }
    }
    reachable
}

fn residue_contradiction(
    coefficients: &[i16; VARIABLE_COUNT],
    rhs: i32,
    target: usize,
    assumed: bool,
    modulus: u8,
) -> Option<ArithmeticCheck> {
    let reachable = reachable_residues(
        coefficients
            .iter()
            .enumerate()
            .filter(|&(var, &coefficient)| var != target && coefficient != 0)
            .map(|(_, &coefficient)| coefficient),
        modulus,
    );
    let residual = rhs - i32::from(coefficients[target]) * i32::from(assumed);
    let required_residue = residual.rem_euclid(i32::from(modulus)) as u8;
    let other = (rhs - i32::from(coefficients[target]) * i32::from(!assumed))
        .rem_euclid(i32::from(modulus)) as u8;
    // Version two refuses certificates that reject both candidate endpoints.
    if reachable & (1 << required_residue) != 0 || reachable & (1 << other) == 0 {
        return None;
    }
    Some(ArithmeticCheck::Residue {
        residual,
        modulus,
        required_residue,
        reachable_residues: (0..modulus).filter(|r| reachable & (1 << r) != 0).collect(),
    })
}

#[derive(Clone)]
struct Node {
    terms: Vec<(usize, i8)>,
    coefficients: [i16; VARIABLE_COUNT],
    rhs: i32,
}

// The heap keeps the worst retained rank at its root, allowing children to be
// discarded during construction rather than after a potentially large allocation.
struct RankedNode {
    rank: (usize, usize, usize, Vec<(usize, i8)>),
    node: Node,
}

impl PartialEq for RankedNode {
    fn eq(&self, other: &Self) -> bool {
        self.rank == other.rank
    }
}
impl Eq for RankedNode {}
impl PartialOrd for RankedNode {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for RankedNode {
    fn cmp(&self, other: &Self) -> Ordering {
        self.rank.cmp(&other.rank)
    }
}

impl Node {
    fn rank(&self) -> (usize, usize, usize, Vec<(usize, i8)>) {
        let mut odd = 0;
        let mut support = 0;
        let mut mass = 0;
        for &coefficient in &self.coefficients {
            odd += usize::from(coefficient % 2 != 0);
            support += usize::from(coefficient != 0);
            mass += usize::from(coefficient.unsigned_abs());
        }
        (odd, support, mass, self.terms.clone())
    }

    fn finding(&self, snapshot: &Snapshot) -> Option<Finding> {
        let support: Vec<_> = self
            .coefficients
            .iter()
            .enumerate()
            .filter(|(_, coefficient)| **coefficient != 0)
            .collect();
        let mut prefix = vec![0; support.len() + 1];
        let mut suffix = vec![0; support.len() + 1];
        let (mut lower, mut upper) = (0, 0);
        for (i, &(_, &coefficient)) in support.iter().enumerate() {
            let c = i32::from(coefficient);
            lower += c.min(0);
            upper += c.max(0);
            prefix[i + 1] = gcd(prefix[i], c.abs());
        }
        for (i, &(_, &coefficient)) in support.iter().enumerate().rev() {
            suffix[i] = gcd(suffix[i + 1], i32::from(coefficient).abs());
        }
        for (i, &(var, &coefficient)) in support.iter().enumerate() {
            let cell = var / 9;
            if snapshot.placed[cell] != 0 {
                continue;
            }
            let c = i32::from(coefficient);
            let lo = lower - c.min(0);
            let hi = upper - c.max(0);
            let divisor = gcd(prefix[i], suffix[i + 1]);
            let zero_rejected = check_bounds(self.rhs, lo, hi, divisor).is_some();
            let one_rejected = check_bounds(self.rhs - c, lo, hi, divisor).is_some();
            // Never report a target for which this combination rejects both values.
            if zero_rejected == one_rejected {
                continue;
            }
            let proof = ArithmeticProof {
                version: LEGACY_VERSION,
                terminal: ArithmeticTerminal::IntervalGcd,
                state_hash: snapshot.hash.clone(),
                terms: self
                    .terms
                    .iter()
                    .map(|&(id, weight)| ArithmeticTerm {
                        requirement: ArithmeticRequirement::from_id(id),
                        weight,
                    })
                    .collect(),
                cell,
                digit: (var % 9 + 1) as u8,
                value: zero_rejected,
            };
            // Reconstruct independently; search calculations alone never authorize a hint.
            if let Some(check) = proof.check_snapshot(snapshot) {
                return Some(make_finding(snapshot, proof, check));
            }
        }
        // Bounds and gcd remain the cheap first pass. Removing any one of two
        // coefficients equal modulo q leaves the same remainder set, so cache
        // the DP by coefficient residue rather than recomputing per candidate.
        for modulus in 2..=MODULUS_LIMIT {
            let mut without_residue = [None; MODULUS_LIMIT as usize];
            for &(var, &coefficient) in &support {
                let cell = var / 9;
                if snapshot.placed[cell] != 0 {
                    continue;
                }
                let residue = i32::from(coefficient).rem_euclid(i32::from(modulus)) as usize;
                if residue == 0 {
                    continue;
                }
                let reachable = *without_residue[residue].get_or_insert_with(|| {
                    reachable_residues(
                        support
                            .iter()
                            .filter(|&&(other, _)| other != var)
                            .map(|&(_, &c)| c),
                        modulus,
                    )
                });
                let zero = self.rhs.rem_euclid(i32::from(modulus)) as u8;
                let one = (self.rhs - i32::from(coefficient)).rem_euclid(i32::from(modulus)) as u8;
                let zero_rejected = reachable & (1 << zero) == 0;
                let one_rejected = reachable & (1 << one) == 0;
                if zero_rejected == one_rejected {
                    continue;
                }
                let proof = ArithmeticProof {
                    version: RESIDUE_VERSION,
                    terminal: ArithmeticTerminal::Residue { modulus },
                    state_hash: snapshot.hash.clone(),
                    terms: self
                        .terms
                        .iter()
                        .map(|&(id, weight)| ArithmeticTerm {
                            requirement: ArithmeticRequirement::from_id(id),
                            weight,
                        })
                        .collect(),
                    cell,
                    digit: (var % 9 + 1) as u8,
                    value: zero_rejected,
                };
                if let Some(check) = proof.check_snapshot(snapshot) {
                    return Some(make_finding(snapshot, proof, check));
                }
            }
        }
        None
    }
}

fn make_finding(snapshot: &Snapshot, proof: ArithmeticProof, check: ArithmeticCheck) -> Finding {
    let sources = proof
        .terms
        .iter()
        .map(|term| format!("{:+} × ({})", term.weight, term.requirement.label()))
        .collect::<Vec<_>>()
        .join("; ");
    let target = format!(
        "r{}c{}={}",
        proof.cell / 9 + 1,
        proof.cell % 9 + 1,
        proof.digit
    );
    let assumption = if proof.value {
        format!("{target} is false")
    } else {
        format!("{target} is true")
    };
    let reason = match check {
        ArithmeticCheck::Interval {
            residual,
            lower,
            upper,
        } => {
            format!("the remaining sum would be {residual}, outside its possible interval [{lower}, {upper}]")
        }
        ArithmeticCheck::Divisibility { residual, divisor } => {
            format!("the remaining coefficients are divisible by {divisor}, but their required sum {residual} is not")
        }
        ArithmeticCheck::Residue {
            residual,
            modulus,
            required_residue,
            reachable_residues,
        } => {
            let reachable = reachable_residues
                .iter()
                .map(u8::to_string)
                .collect::<Vec<_>>()
                .join(", ");
            format!("the remaining sum would be {residual}, which has remainder {required_residue} modulo {modulus}, but the remaining candidates can produce only remainders {{{reachable}}}")
        }
    };
    let conclusion = if proof.value {
        format!(
            "Place {} in r{}c{}",
            proof.digit,
            proof.cell / 9 + 1,
            proof.cell % 9 + 1
        )
    } else {
        format!(
            "Exclude {} from r{}c{}",
            proof.digit,
            proof.cell / 9 + 1,
            proof.cell % 9 + 1
        )
    };
    let mut involved = [false; 81];
    for term in &proof.terms {
        for var in snapshot.requirement_variables(term.requirement.id().expect("checked source")) {
            involved[var / 9] = true;
        }
    }
    Finding {
        technique: Technique::ArithmeticCounting,
        inference: if proof.value {
            InferenceResult::Placement { cell: proof.cell, value: proof.digit }
        } else {
            InferenceResult::Elimination { cell: proof.cell, values: vec![proof.digit] }
        },
        involved_cells: (0..81).filter(|&cell| involved[cell]).collect(),
        explanation: ExplanationData::Raw(format!(
            "Combine these exactly-one requirements: {sources}. If {assumption}, {reason}. {conclusion}."
        )),
        proof: Some(ProofCertificate::Arithmetic(proof)),
    }
}

/// Search the current stored masks without recalculating candidates or changing the grid.
pub fn search(grid: &Grid, options: &ArithmeticSearchOptions) -> ArithmeticSearchResult {
    let (finding, mut result) = search_inner(grid, options);
    result.hint = finding.map(|finding| finding.to_hint());
    result
}

pub(crate) fn find(grid: &Grid, _fab: &CandidateFabric) -> Option<Finding> {
    search_inner(grid, &ArithmeticSearchOptions::default()).0
}

fn search_inner(
    grid: &Grid,
    options: &ArithmeticSearchOptions,
) -> (Option<Finding>, ArithmeticSearchResult) {
    let mut result = ArithmeticSearchResult {
        hint: None,
        tested_combinations: 0,
        budget_exhausted: false,
        beam_pruned: false,
        error: None,
    };
    if !(1..=SOURCE_LIMIT).contains(&options.max_sources)
        || !(1..=WEIGHT_LIMIT).contains(&options.max_weight)
        || !(1..=512).contains(&options.beam_width)
        || options.max_combinations > 1_000_000
    {
        result.error = Some(
            "Expected sources 1..6, weight 1..2, beam width 1..512, and combinations 0..1000000"
                .into(),
        );
        return (None, result);
    }
    let snapshot = match Snapshot::from_grid(grid) {
        Ok(snapshot) => snapshot,
        Err(error) => {
            result.error = Some(error);
            return (None, result);
        }
    };
    if options.max_combinations == 0 {
        result.budget_exhausted = true;
        return (None, result);
    }
    let sources: Vec<_> = (0..324)
        .map(|id| snapshot.requirement_variables(id))
        .collect();
    let mut incidence: Vec<Vec<usize>> = vec![Vec::new(); VARIABLE_COUNT];
    let mut frontier = Vec::new();
    for (id, vars) in sources.iter().enumerate() {
        if vars.iter().all(|&var| snapshot.placed[var / 9] != 0) {
            continue;
        }
        for &var in vars {
            incidence[var].push(id);
        }
        let mut coefficients = [0; VARIABLE_COUNT];
        for &var in vars {
            coefficients[var] = 1;
        }
        frontier.push(Node {
            terms: vec![(id, 1)],
            coefficients,
            rhs: 1,
        });
    }
    let weights: &[i8] = if options.max_weight == 1 {
        &[1, -1]
    } else {
        &[1, -1, 2, -2]
    };
    for depth in 1..=options.max_sources {
        if frontier.is_empty() {
            break;
        }
        frontier.sort_by_cached_key(Node::rank);
        let remaining = options.max_combinations - result.tested_combinations;
        // Reserve work for larger source groups instead of consuming the whole
        // budget on pairs. Allocation limits are reported as pruning.
        let allocation = (remaining / (options.max_sources - depth + 1))
            .max(1)
            .min(remaining);
        let mut evaluated = Vec::new();
        for (evaluated_count, node) in frontier.into_iter().enumerate() {
            if evaluated_count == allocation {
                result.beam_pruned = true;
                break;
            }
            result.tested_combinations += 1;
            if let Some(finding) = node.finding(&snapshot) {
                // A checked parity placement can supply a small, motivated
                // extension without widening the beam or increasing its weight
                // range. Each distinct tail combination uses the same budget.
                if let Some(ProofCertificate::Arithmetic(proof)) = &finding.proof {
                    if proof.value && proof.terms.len() < options.max_sources {
                        let guardian = proof.cell * 9 + usize::from(proof.digit - 1);
                        for &id in &incidence[guardian] {
                            if result.tested_combinations == options.max_combinations {
                                break;
                            }
                            let tail = ArithmeticRequirement::from_id(id);
                            if sources[id].len() < 2
                                || sources[id].len() - 1 > options.max_weight as usize
                                || proof.terms.iter().any(|term| term.requirement == tail)
                                || proof
                                    .terms
                                    .iter()
                                    .any(|term| !matches!(term.weight, -1 | 1))
                            {
                                continue;
                            }
                            result.tested_combinations += 1;
                            if let Some(compiled) = proof
                                .compile_tail_snapshot(&snapshot, tail)
                                .into_iter()
                                .next()
                            {
                                if let Some(check) = compiled.check_snapshot(&snapshot) {
                                    return (
                                        Some(make_finding(&snapshot, compiled, check)),
                                        result,
                                    );
                                }
                            }
                        }
                    }
                }
                return (Some(finding), result);
            }
            if evaluated.len() < options.beam_width {
                evaluated.push(node);
            } else if depth < options.max_sources {
                result.beam_pruned = true;
            }
        }
        if result.tested_combinations == options.max_combinations {
            result.budget_exhausted = true;
            break;
        }
        if depth == options.max_sources {
            break;
        }
        // At most beam_width * 6 * 9 * 4 * 4 extensions are considered.
        // Retain at most FRONTIER_LIMIT children and their deduplication keys;
        // the configured evaluation budget never implies a large allocation.
        let mut seen = HashSet::new();
        let mut children: BinaryHeap<RankedNode> = BinaryHeap::new();
        for node in evaluated {
            let mut adjacent = [false; 324];
            for &(id, _) in &node.terms {
                for &var in &sources[id] {
                    for &neighbor in &incidence[var] {
                        adjacent[neighbor] = true;
                    }
                }
            }
            for &(id, _) in &node.terms {
                adjacent[id] = false;
            }
            for (id, is_adjacent) in adjacent.into_iter().enumerate() {
                if !is_adjacent {
                    continue;
                }
                for &weight in weights {
                    let mut terms = node.terms.clone();
                    terms.push((id, weight));
                    terms.sort_unstable();
                    let sign = if terms[0].1 < 0 { -1 } else { 1 };
                    for (_, w) in &mut terms {
                        *w *= sign;
                    }
                    if seen.contains(&terms) {
                        continue;
                    }
                    let mut coefficients = node.coefficients;
                    for &var in &sources[id] {
                        coefficients[var] += i16::from(weight);
                    }
                    if sign == -1 {
                        for coefficient in &mut coefficients {
                            *coefficient = -*coefficient;
                        }
                    }
                    let child = Node {
                        terms,
                        coefficients,
                        rhs: (node.rhs + i32::from(weight)) * i32::from(sign),
                    };
                    let rank = child.rank();
                    if children.len() == FRONTIER_LIMIT {
                        result.beam_pruned = true;
                        if rank >= children.peek().expect("full frontier").rank {
                            continue;
                        }
                        let evicted = children.pop().expect("full frontier");
                        seen.remove(&evicted.node.terms);
                    }
                    seen.insert(child.terms.clone());
                    children.push(RankedNode { rank, node: child });
                }
            }
        }
        frontier = children.into_iter().map(|child| child.node).collect();
    }
    (None, result)
}

#[cfg(test)]
#[path = "arithmetic_tests.rs"]
mod tests;
