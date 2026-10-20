//! Exact matching support for one crossing row, column and box.
//!
//! Search uses subset matching DP. Certificates instead cover every boundary
//! signature with a directly checked Hall-deficient subset; the checker never
//! invokes the search matcher. Candidate domains are conditional premises.
use super::candidate_state::CandidateState;
use super::explain::{ExplanationData, Finding, InferenceResult, ProofCertificate};
use super::fabric::{idx_to_pos, sector_cells};
use super::{Hint, Technique};
use crate::Grid;
use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};

const ALL: u16 = 0x1ff;
const VERSION: u8 = 1;
const MAX_BLOCKERS: usize = 366;

/// A Hall-deficient subset inside one membership atom, under the target assumption.
/// Atoms 1/2 are the two row/box and column/box arm cells; atoms 3/4 are
/// the six row-private and column-private cells; atom 5 is the four box-private
/// cells. Cells in each atom are ordered by their row-major index.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct HallObstruction {
    pub atom: u8,
    /// The exact digit set assigned to this atom by a boundary signature.
    pub digits: u16,
    /// Local cell indices, bit zero naming the first cell in the atom.
    pub cells: u8,
}

/// A non-modular elimination certificate for three native houses.
///
/// `witness` assigns all 21 union cells in ascending row-major order, showing
/// that the unforced relaxation remains feasible and avoids the target. The
/// blockers refute every signature under the target=true assumption. This
/// proves neither global satisfiability nor the validity of earlier deletions.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ThreeHouseProof {
    pub version: u8,
    pub state_hash: String,
    /// Intersection of the selected row, column and box (0..80).
    pub center: usize,
    pub cell: usize,
    pub digit: u8,
    pub witness: Vec<u8>,
    pub blockers: Vec<HallObstruction>,
}

/// Counts reconstructed by the checker, rather than trusted from a proof.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ThreeHouseCheck {
    pub signatures: usize,
    pub hall_obstructions: usize,
}

/// Deterministic prefix of the 81 crossing triples to inspect (0..=81).
#[derive(Debug, Clone)]
pub struct ThreeHouseSearchOptions {
    pub max_centers: usize,
}
impl Default for ThreeHouseSearchOptions {
    fn default() -> Self {
        Self { max_centers: 81 }
    }
}

#[derive(Debug, Clone)]
pub struct ThreeHouseSearchResult {
    pub hint: Option<Hint>,
    pub centers_checked: usize,
    /// True only if all 81 centers were checked and no deduction was found.
    pub exhausted: bool,
    pub error: Option<String>,
}

struct Layout {
    houses: [usize; 3],
    atoms: [Vec<usize>; 6],
    cells: Vec<usize>,
}
impl Layout {
    fn new(center: usize) -> Option<Self> {
        if center >= 81 {
            return None;
        }
        let (r, c) = (center / 9, center % 9);
        let houses = [r, 9 + c, 18 + r / 3 * 3 + c / 3];
        let mut atoms: [Vec<usize>; 6] = Default::default();
        let mut cells = Vec::with_capacity(21);
        for cell in 0..81 {
            let (rr, cc) = (cell / 9, cell % 9);
            let (row, col, bx) = (rr == r, cc == c, rr / 3 == r / 3 && cc / 3 == c / 3);
            let atom = match (row, col, bx) {
                (true, true, true) => 0,
                (true, false, true) => 1,
                (false, true, true) => 2,
                (true, false, false) => 3,
                (false, true, false) => 4,
                (false, false, true) => 5,
                _ => continue,
            };
            atoms[atom].push(cell);
            cells.push(cell);
        }
        Some(Self {
            houses,
            atoms,
            cells,
        })
    }
}

fn state_hash(state: &CandidateState) -> String {
    let mut hash = Sha256::new();
    hash.update(b"sudoku-three-house-v1-native");
    for (&value, &domain) in state.values.iter().zip(&state.domains) {
        hash.update([value]);
        hash.update(domain.to_le_bytes());
    }
    format!("{:x}", hash.finalize())
}

/// Visit all signatures allowed by the center domain. Returning false stops.
/// Digit masks indexed by atom: d, A, B, U\(d+A), U\(d+B), U\(d+A+B).
fn signatures(center: u16, mut visit: impl FnMut([u16; 6]) -> bool) -> bool {
    for digit in 0..9 {
        let d = 1 << digit;
        if center & d == 0 {
            continue;
        }
        for a in 1u16..=ALL {
            if a.count_ones() != 2 || a & d != 0 {
                continue;
            }
            for b in 1u16..=ALL {
                if b.count_ones() != 2 || b & (d | a) != 0 {
                    continue;
                }
                if !visit([d, a, b, ALL ^ (d | a), ALL ^ (d | b), ALL ^ (d | a | b)]) {
                    return false;
                }
            }
        }
    }
    true
}

impl ThreeHouseProof {
    /// Produce an independently checked certificate for a requested elimination.
    pub fn for_elimination(grid: &Grid, center: usize, cell: usize, digit: u8) -> Option<Self> {
        let state = CandidateState::capture(grid).ok()?;
        let layout = Layout::new(center)?;
        build_proof(&state, &layout, center, cell, digit)
    }

    pub fn verify(&self, grid: &Grid) -> bool {
        self.check(grid).is_some()
    }

    /// Check Hall cardinalities and complete signature coverage without matching DP.
    pub fn check(&self, grid: &Grid) -> Option<ThreeHouseCheck> {
        let state = CandidateState::capture(grid).ok()?;
        self.check_state(&state)
    }

    fn check_state(&self, state: &CandidateState) -> Option<ThreeHouseCheck> {
        let layout = Layout::new(self.center)?;
        if self.version != VERSION
            || self.state_hash != state_hash(state)
            || self.cell >= 81
            || !(1..=9).contains(&self.digit)
            || state.values[self.cell] != 0
            || !layout.cells.contains(&self.cell)
            || state.domains[self.cell] & (1 << (self.digit - 1)) == 0
            || self.witness.len() != 21
            || self.blockers.len() > MAX_BLOCKERS
        {
            return None;
        }
        // Independently validate the explicit base witness using native houses.
        let mut assignment = [0u8; 81];
        for (&cell, &digit) in layout.cells.iter().zip(&self.witness) {
            if !(1..=9).contains(&digit) || state.domains[cell] & (1 << (digit - 1)) == 0 {
                return None;
            }
            assignment[cell] = digit;
        }
        if assignment[self.cell] == self.digit {
            return None;
        }
        for sector in layout.houses {
            let mut seen = 0u16;
            for cell in sector_cells(sector) {
                let digit = assignment[cell];
                if digit == 0 || seen & (1 << (digit - 1)) != 0 {
                    return None;
                }
                seen |= 1 << (digit - 1);
            }
            if seen != ALL {
                return None;
            }
        }
        let mut domains = state.domains.clone();
        domains[self.cell] = 1 << (self.digit - 1);
        let mut blocked = [[false; 512]; 6];
        for block in &self.blockers {
            let atom = usize::from(block.atom);
            if !(1..=5).contains(&atom) {
                return None;
            }
            let cells = &layout.atoms[atom];
            if block.digits & !ALL != 0
                || block.digits.count_ones() as usize != cells.len()
                || block.cells == 0
                || usize::from(block.cells) >= (1 << cells.len())
                || blocked[atom][usize::from(block.digits)]
            {
                return None;
            }
            let mut neighbors = 0u16;
            for (i, &cell) in cells.iter().enumerate() {
                if block.cells & (1 << i) != 0 {
                    neighbors |= domains[cell] & block.digits;
                }
            }
            if neighbors.count_ones() >= block.cells.count_ones() {
                return None;
            }
            blocked[atom][usize::from(block.digits)] = true;
        }
        let mut count = 0;
        let covered = signatures(domains[self.center], |sets| {
            count += 1;
            (1..6).any(|atom| blocked[atom][usize::from(sets[atom])])
        });
        if !covered || count == 0 {
            return None;
        }
        Some(ThreeHouseCheck {
            signatures: count,
            hall_obstructions: self.blockers.len(),
        })
    }
}

#[derive(Clone)]
struct AtomMatch {
    support: Vec<u16>,
    witness: Vec<u8>,
}

// Search-only subset DP: memoized existence followed by traversal of successful
// transitions. The certificate checker above uses neither this DP nor its cache.
fn matching(domains: &[u16], digits: u16) -> Option<AtomMatch> {
    fn feasible(domains: &[u16], remaining: u16, memo: &mut [i8; 512]) -> bool {
        let slot = usize::from(remaining);
        if memo[slot] != 0 {
            return memo[slot] > 0;
        }
        let i = domains.len() - remaining.count_ones() as usize;
        let mut choices = domains[i] & remaining;
        while choices != 0 {
            let bit = 1 << choices.trailing_zeros();
            choices ^= bit;
            if feasible(domains, remaining ^ bit, memo) {
                memo[slot] = 1;
                return true;
            }
        }
        memo[slot] = -1;
        false
    }
    fn collect(
        domains: &[u16],
        remaining: u16,
        memo: &mut [i8; 512],
        visited: &mut [bool; 512],
        support: &mut [u16],
    ) {
        if remaining == 0 || visited[usize::from(remaining)] {
            return;
        }
        visited[usize::from(remaining)] = true;
        let i = domains.len() - remaining.count_ones() as usize;
        let mut choices = domains[i] & remaining;
        while choices != 0 {
            let bit = 1 << choices.trailing_zeros();
            choices ^= bit;
            if feasible(domains, remaining ^ bit, memo) {
                support[i] |= bit;
                collect(domains, remaining ^ bit, memo, visited, support);
            }
        }
    }
    let mut memo = [0i8; 512];
    memo[0] = 1;
    if !feasible(domains, digits, &mut memo) {
        return None;
    }
    let mut support = vec![0; domains.len()];
    collect(domains, digits, &mut memo, &mut [false; 512], &mut support);
    let mut remaining = digits;
    let mut witness = Vec::with_capacity(domains.len());
    for domain in domains {
        let mut choices = domain & remaining;
        while choices != 0 {
            let bit = 1 << choices.trailing_zeros();
            choices ^= bit;
            if feasible(domains, remaining ^ bit, &mut memo) {
                witness.push(bit.trailing_zeros() as u8 + 1);
                remaining ^= bit;
                break;
            }
        }
    }
    Some(AtomMatch { support, witness })
}

struct Matches<'a> {
    layout: &'a Layout,
    domains: &'a [u16],
    cache: Vec<Vec<Option<Option<AtomMatch>>>>,
}
impl<'a> Matches<'a> {
    fn new(layout: &'a Layout, domains: &'a [u16]) -> Self {
        Self {
            layout,
            domains,
            cache: vec![vec![None; 512]; 6],
        }
    }
    fn get(&mut self, atom: usize, digits: u16) -> Option<&AtomMatch> {
        let entry = &mut self.cache[atom][usize::from(digits)];
        entry
            .get_or_insert_with(|| {
                let domains: Vec<_> = self.layout.atoms[atom]
                    .iter()
                    .map(|&cell| self.domains[cell])
                    .collect();
                matching(&domains, digits)
            })
            .as_ref()
    }
    fn witness(&mut self, sets: [u16; 6]) -> Option<Vec<u8>> {
        let mut values = [0; 81];
        for (atom, digits) in sets.into_iter().enumerate() {
            let witness = self.get(atom, digits)?.witness.clone();
            for (&cell, digit) in self.layout.atoms[atom].iter().zip(witness) {
                values[cell] = digit;
            }
        }
        Some(self.layout.cells.iter().map(|&cell| values[cell]).collect())
    }
}

fn hall_subset(cells: &[usize], domains: &[u16], digits: u16) -> Option<u8> {
    for subset in 1u8..(1 << cells.len()) {
        let neighbors = cells.iter().enumerate().fold(0u16, |mask, (i, &cell)| {
            mask | if subset & (1 << i) != 0 {
                domains[cell] & digits
            } else {
                0
            }
        });
        if neighbors.count_ones() < subset.count_ones() {
            return Some(subset);
        }
    }
    None
}

fn build_proof(
    state: &CandidateState,
    layout: &Layout,
    center: usize,
    cell: usize,
    digit: u8,
) -> Option<ThreeHouseProof> {
    if cell >= 81
        || !(1..=9).contains(&digit)
        || state.values[cell] != 0
        || !layout.cells.contains(&cell)
        || state.domains[cell] & (1 << (digit - 1)) == 0
    {
        return None;
    }
    let mut unforced = Matches::new(layout, &state.domains);
    let mut witness = None;
    signatures(state.domains[center], |sets| {
        witness = unforced.witness(sets);
        witness.is_none()
    });
    let witness = witness?;
    let mut assumed = state.domains.clone();
    assumed[cell] = 1 << (digit - 1);
    let mut matcher = Matches::new(layout, &assumed);
    let mut recorded = [[false; 512]; 6];
    let mut blockers = Vec::new();
    let refuted = signatures(assumed[center], |sets| {
        for atom in 1..6 {
            let digits = sets[atom];
            if matcher.get(atom, digits).is_none() {
                if !recorded[atom][usize::from(digits)] {
                    let Some(cells) = hall_subset(&layout.atoms[atom], &assumed, digits) else {
                        return false;
                    };
                    blockers.push(HallObstruction {
                        atom: atom as u8,
                        digits,
                        cells,
                    });
                    recorded[atom][usize::from(digits)] = true;
                }
                return true;
            }
        }
        false
    });
    if !refuted {
        return None;
    }
    let proof = ThreeHouseProof {
        version: VERSION,
        state_hash: state_hash(state),
        center,
        cell,
        digit,
        witness,
        blockers,
    };
    proof.check_state(state)?;
    Some(proof)
}

fn finding(proof: ThreeHouseProof) -> Finding {
    let center = proof.center;
    let cell = proof.cell;
    let digit = proof.digit;
    let layout = Layout::new(center).expect("checked center");
    Finding {
        technique: Technique::ThreeHouseMatching,
        inference: InferenceResult::Elimination { cell, values: vec![digit] },
        involved_cells: layout.cells,
        explanation: ExplanationData::Raw(format!(
            "If r{}c{} were {}, row {}, column {} and box {} could not all be filled: every boundary signature has a Hall-deficient cell set. Therefore remove {} from r{}c{}.",
            cell / 9 + 1, cell % 9 + 1, digit, center / 9 + 1, center % 9 + 1,
            center / 27 * 3 + center % 9 / 3 + 1, digit, cell / 9 + 1, cell % 9 + 1)),
        proof: Some(ProofCertificate::ThreeHouse(proof)),
    }
}

pub(super) fn find(grid: &Grid) -> Option<Finding> {
    search(grid, &ThreeHouseSearchOptions::default())
        .hint
        .and_then(|hint| match hint.proof? {
            ProofCertificate::ThreeHouse(proof) => Some(finding(proof)),
            _ => None,
        })
}

pub(super) fn search(grid: &Grid, options: &ThreeHouseSearchOptions) -> ThreeHouseSearchResult {
    let mut result = ThreeHouseSearchResult {
        hint: None,
        centers_checked: 0,
        exhausted: false,
        error: None,
    };
    if options.max_centers > 81 {
        result.error = Some("Three-house center budget must be in 0..=81".into());
        return result;
    }
    let state = match CandidateState::capture(grid) {
        Ok(state) => state,
        Err(error) => {
            result.error = Some(error);
            return result;
        }
    };
    for center in 0..options.max_centers {
        let layout = Layout::new(center).expect("bounded center");
        let mut matches = Matches::new(&layout, &state.domains);
        let mut support = [0u16; 81];
        let mut feasible = false;
        signatures(state.domains[center], |sets| {
            if (0..6).any(|atom| matches.get(atom, sets[atom]).is_none()) {
                return true;
            }
            feasible = true;
            for (atom, digits) in sets.into_iter().enumerate() {
                let supported = &matches.get(atom, digits).expect("feasible atom").support;
                for (&cell, &mask) in layout.atoms[atom].iter().zip(supported) {
                    support[cell] |= mask;
                }
            }
            true
        });
        result.centers_checked += 1;
        if !feasible {
            continue;
        }
        for &cell in &layout.cells {
            if state.values[cell] != 0 {
                continue;
            }
            let missing = state.domains[cell] & !support[cell];
            if missing == 0 {
                continue;
            }
            let digit = missing.trailing_zeros() as u8 + 1;
            if let Some(proof) = build_proof(&state, &layout, center, cell, digit) {
                result.hint = Some(finding(proof).to_hint());
                return result;
            }
        }
    }
    result.exhausted = options.max_centers == 81;
    result
}

/// Portable one-step elimination, with exact candidate premises.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ThreeHouseReplay {
    pub version: u8,
    pub state: CandidateState,
    pub proof: ThreeHouseProof,
}
impl ThreeHouseReplay {
    pub fn capture(grid: &Grid, proof: &ThreeHouseProof) -> Result<Self, String> {
        let replay = Self {
            version: VERSION,
            state: CandidateState::capture(grid)?,
            proof: proof.clone(),
        };
        replay.check()?;
        Ok(replay)
    }
    pub fn check(&self) -> Result<ThreeHouseCheck, String> {
        if self.version != VERSION {
            return Err("Unsupported three-house replay version".into());
        }
        // Validate serialized state before indexing its vectors in the checker.
        self.state.to_grid()?;
        self.proof
            .check_state(&self.state)
            .ok_or_else(|| "Three-house certificate does not verify".into())
    }
    pub fn verify(&self) -> bool {
        self.check().is_ok()
    }
    /// Atomic application to exactly the recorded native state; masks are preserved.
    pub fn apply(&self, grid: &mut Grid) -> Result<(), String> {
        self.check()?;
        if CandidateState::capture(grid)? != self.state {
            return Err("Three-house replay does not match the current state".into());
        }
        let mut next = grid.deep_clone();
        next.cell_mut(idx_to_pos(self.proof.cell))
            .remove_candidate(self.proof.digit);
        if !next.validate().is_valid {
            return Err("Three-house replay would violate an installed constraint".into());
        }
        CandidateState::capture(&next)?;
        *grid = next;
        Ok(())
    }
}

#[cfg(test)]
#[path = "three_house_tests.rs"]
mod tests;
