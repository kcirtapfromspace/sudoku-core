//! Declarative technique witnesses. This checker never calls a technique detector.
#![allow(deprecated)]
use crate::{
    CandidateState, ExactCoverProof, ExactPurpose, Grid, Hint, ProofCertificate, Technique,
};
use serde::{Deserialize, Serialize};
use std::collections::HashSet;

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PatternSet {
    pub cells: Vec<usize>,
    pub digits: Vec<u8>,
    pub house: usize,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PatternFish {
    pub digit: u8,
    pub bases: Vec<usize>,
    pub covers: Vec<usize>,
    pub fins: Vec<usize>,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PatternStep {
    pub candidate: u16,
    pub on: bool,
    /// Positive: exactly-one requirement. Negative: incompatible true candidate.
    pub reason: u16,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PatternBranch {
    pub assumption: u16,
    pub steps: Vec<PatternStep>,
    pub conflict: Option<u16>,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub enum PatternWitness {
    Subset {
        house: usize,
        cells: Vec<usize>,
        digits: Vec<u8>,
    },
    Fish {
        parts: Vec<PatternFish>,
    },
    Wing {
        cells: Vec<usize>,
        digits: Vec<u8>,
    },
    Als {
        sets: Vec<PatternSet>,
        links: Vec<u8>,
        z: Option<u8>,
    },
    Chain {
        nodes: Vec<u16>,
        strong: Vec<bool>,
        colors: Vec<bool>,
    },
    Trade {
        cells: Vec<usize>,
        digits: Vec<u8>,
    },
    Bug,
    Aligned {
        cells: Vec<usize>,
        exclusions: Vec<Vec<usize>>,
    },
    Forcing {
        source: super::solver::explain::ForcingSource,
        fish: Option<PatternFish>,
        branches: Vec<PatternBranch>,
    },
    Native {
        proof: Box<ProofCertificate>,
    },
    Search {
        proof: ExactCoverProof,
    },
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PatternUniqueness {
    pub givens: String,
    pub proof: ExactCoverProof,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct NamedPatternProof {
    pub version: u8,
    pub state_hash: String,
    pub technique: Technique,
    pub conclusion: ExactPurpose,
    pub witness: PatternWitness,
    pub uniqueness: Option<PatternUniqueness>,
}

fn peers(a: usize, b: usize) -> bool {
    a != b && (a / 9 == b / 9 || a % 9 == b % 9 || (a / 27 == b / 27 && a % 9 / 3 == b % 9 / 3))
}
fn house(h: usize) -> Vec<usize> {
    (0..81)
        .filter(|&c| {
            if h < 9 {
                c / 9 == h
            } else if h < 18 {
                c % 9 == h - 9
            } else {
                h < 27 && c / 27 * 3 + c % 9 / 3 == h - 18
            }
        })
        .collect()
}
fn mask(ds: &[u8]) -> Option<u16> {
    let mut m = 0;
    for &d in ds {
        if !(1..=9).contains(&d) || m & (1 << (d - 1)) != 0 {
            return None;
        }
        m |= 1 << (d - 1);
    }
    Some(m)
}
fn unique(xs: &[usize], limit: usize) -> bool {
    !xs.is_empty()
        && xs.iter().all(|&c| c < limit)
        && xs.iter().copied().collect::<HashSet<_>>().len() == xs.len()
}
struct Context<'a> {
    s: CandidateState,
    hint: &'a Hint,
    cell: usize,
    ds: Vec<u8>,
    placement: bool,
}
impl Context<'_> {
    fn live(&self, c: usize, d: u8) -> bool {
        c < 81
            && (1..=9).contains(&d)
            && self.s.values[c] == 0
            && self.s.domains[c] & (1 << (d - 1)) != 0
    }
    fn locs(&self, h: usize, d: u8) -> Vec<usize> {
        house(h).into_iter().filter(|&c| self.live(c, d)).collect()
    }
    fn union(&self, cs: &[usize]) -> u16 {
        cs.iter().fold(0, |m, &c| m | self.s.domains[c])
    }
    fn elim(&self, d: u8) -> bool {
        !self.placement && self.ds == [d]
    }
    fn sees_all(&self, cs: &[usize], d: u8) -> bool {
        let locs: Vec<_> = cs.iter().copied().filter(|&c| self.live(c, d)).collect();
        !locs.is_empty() && locs.iter().all(|&c| peers(c, self.cell))
    }
    fn conflict(a: u16, b: u16) -> bool {
        a != b
            && (a / 9 == b / 9 || (a % 9 == b % 9 && peers(usize::from(a / 9), usize::from(b / 9))))
    }
    fn link(&self, a: u16, b: u16, strong: bool) -> bool {
        if a >= 729
            || b >= 729
            || !self.live(usize::from(a / 9), (a % 9 + 1) as u8)
            || !self.live(usize::from(b / 9), (b % 9 + 1) as u8)
            || !Self::conflict(a, b)
        {
            return false;
        }
        !strong
            || if a / 9 == b / 9 {
                self.s.domains[usize::from(a / 9)].count_ones() == 2
            } else {
                (0..27).any(|h| {
                    self.locs(h, (a % 9 + 1) as u8)
                        .iter()
                        .map(|&c| c as u16 * 9 + a % 9)
                        .collect::<HashSet<_>>()
                        == HashSet::from([a, b])
                })
            }
    }
    fn fish(&self, f: &PatternFish, endpoint: bool) -> bool {
        if !(1..=9).contains(&f.digit)
            || !unique(&f.bases, 27)
            || !unique(&f.covers, 27)
            || f.bases.len() != f.covers.len()
            || f.bases.len() > 4
            || f.bases.iter().any(|h| f.covers.contains(h))
        {
            return false;
        }
        let independent = |hs: &[usize]| {
            let mut seen = HashSet::new();
            hs.iter().all(|&h| {
                let ls = self.locs(h, f.digit);
                !ls.is_empty() && ls.into_iter().all(|c| seen.insert(c))
            })
        };
        if !independent(&f.bases) || !independent(&f.covers) {
            return false;
        }
        let base: HashSet<_> = f
            .bases
            .iter()
            .flat_map(|&h| self.locs(h, f.digit))
            .collect();
        let cover: HashSet<_> = f
            .covers
            .iter()
            .flat_map(|&h| self.locs(h, f.digit))
            .collect();
        let fins: HashSet<_> = base.difference(&cover).copied().collect();
        fins == f.fins.iter().copied().collect()
            && f.fins.len() == fins.len()
            && f.fins.iter().all(|&c| c < 81)
            && (!endpoint
                || (self.elim(f.digit)
                    && cover.contains(&self.cell)
                    && !base.contains(&self.cell)
                    && fins.iter().all(|&c| peers(c, self.cell))))
    }
    fn als(&self, a: &PatternSet) -> bool {
        a.house < 27
            && unique(&a.cells, 81)
            && a.cells.len() <= 8
            && a.cells
                .iter()
                .all(|&c| self.s.values[c] == 0 && house(a.house).contains(&c))
            && mask(&a.digits) == Some(self.union(&a.cells))
            && a.digits.len() == a.cells.len() + 1
    }
    fn rcc(&self, a: &PatternSet, b: &PatternSet, d: u8) -> bool {
        a.digits.contains(&d)
            && b.digits.contains(&d)
            && a.cells.iter().filter(|&&c| self.live(c, d)).all(|&c| {
                b.cells
                    .iter()
                    .filter(|&&o| self.live(o, d))
                    .all(|&o| peers(c, o))
            })
    }
    fn requirement(&self, r: u16) -> Vec<u16> {
        if r < 81 {
            (1..=9)
                .filter(|&d| self.s.domains[usize::from(r)] & (1 << (d - 1)) != 0)
                .map(|d| r * 9 + d - 1)
                .collect()
        } else if r < 324 {
            let h = usize::from((r - 81) / 9);
            let d = (r - 81) % 9 + 1;
            house(h)
                .into_iter()
                .filter(|&c| self.s.domains[c] & (1 << (d - 1)) != 0)
                .map(|c| c as u16 * 9 + d - 1)
                .collect()
        } else {
            vec![]
        }
    }
    fn branch(&self, b: &PatternBranch, digit_only: Option<u8>) -> bool {
        if b.assumption >= 729
            || !self.live(usize::from(b.assumption / 9), (b.assumption % 9 + 1) as u8)
            || b.steps.len() > 729
        {
            return false;
        }
        let mut assigned = [0i8; 729];
        for c in 0..81 {
            for d in 1..=9 {
                let n = c * 9 + d - 1;
                if self.s.domains[c] & (1 << (d - 1)) == 0 {
                    assigned[n] = -1;
                } else if self.s.values[c] != 0 {
                    assigned[n] = 1;
                }
            }
        }
        assigned[usize::from(b.assumption)] = 1;
        for step in &b.steps {
            let n = usize::from(step.candidate);
            if n >= 729
                || assigned[n] != 0
                || digit_only.is_some_and(|d| n % 9 + 1 != usize::from(d))
            {
                return false;
            }
            if step.on {
                if step.reason >= 324 || (digit_only.is_some() && step.reason < 81) {
                    return false;
                }
                let ls = self.requirement(step.reason);
                if !ls.contains(&step.candidate)
                    || ls
                        .iter()
                        .any(|&o| o != step.candidate && assigned[usize::from(o)] != -1)
                {
                    return false;
                }
            } else if step.reason >= 729
                || assigned[usize::from(step.reason)] != 1
                || !Self::conflict(step.candidate, step.reason)
            {
                return false;
            }
            assigned[n] = if step.on { 1 } else { -1 };
        }
        if let Some(r) = b.conflict {
            return r < 324
                && (digit_only.is_none() || r >= 81)
                && self
                    .requirement(r)
                    .iter()
                    .all(|&n| assigned[usize::from(n)] == -1);
        }
        if self.placement {
            assigned[self.cell * 9 + usize::from(self.ds[0] - 1)] == 1
        } else {
            self.ds
                .iter()
                .all(|&d| assigned[self.cell * 9 + usize::from(d - 1)] == -1)
        }
    }
}

impl NamedPatternProof {
    pub fn verify(&self, grid: &Grid, hint: &Hint) -> bool {
        self.check(grid, hint).is_ok()
    }
    pub fn check(&self, grid: &Grid, hint: &Hint) -> Result<(), String> {
        let s = CandidateState::capture(grid)?;
        if self.version != 1
            || self.technique != hint.technique
            || self.state_hash != crate::exact_cover::hash(&s)
            || Some(self.conclusion.clone()) != crate::exact_cover::purpose(hint)
        {
            return Err("Named pattern has wrong state, method or conclusion".into());
        }
        let ExactPurpose::Deduction {
            cell,
            digits,
            placement,
        } = &self.conclusion
        else {
            return Err("Pattern lacks conclusion".into());
        };
        if *cell >= 81
            || s.values[*cell] != 0
            || mask(digits).is_none_or(|m| m == 0 || m & !s.domains[*cell] != 0)
            || (*placement && digits.len() != 1)
        {
            return Err("Invalid pattern endpoint".into());
        }
        for c in 0..81 {
            if s.values[c] == 0
                && (0..81).any(|o| {
                    s.values[o] != 0 && peers(c, o) && s.domains[c] & (1 << (s.values[o] - 1)) != 0
                })
            {
                return Err("Pattern contains a candidate excluded by a placed peer".into());
            }
        }
        let unique = matches!(
            self.technique,
            Technique::UniqueRectangle
                | Technique::HiddenRectangle
                | Technique::AvoidableRectangle
                | Technique::ExtendedUniqueRectangle
                | Technique::BivalueUniversalGrave
        );
        if unique != self.uniqueness.is_some() {
            return Err("Wrong uniqueness premise for named method".into());
        }
        let x = Context {
            s,
            hint,
            cell: *cell,
            ds: digits.clone(),
            placement: *placement,
        };
        if self.validate(&x, grid) {
            Ok(())
        } else {
            Err(format!("Named {:?} pattern rejected", hint.technique))
        }
    }
}
impl NamedPatternProof {
    fn validate(&self, x: &Context<'_>, grid: &Grid) -> bool {
        use PatternWitness::*;
        use Technique::*;
        let t = self.technique;
        match &self.witness {
            Subset {
                house: h,
                cells,
                digits,
            } => {
                if *h >= 27
                    || !unique(cells, 81)
                    || cells
                        .iter()
                        .any(|&c| x.s.values[c] != 0 || !house(*h).contains(&c))
                {
                    return false;
                }
                let Some(m) = mask(digits) else {
                    return false;
                };
                let (n, hidden) = match t {
                    NakedSingle => (1, false),
                    HiddenSingle => (1, true),
                    NakedPair => (2, false),
                    HiddenPair => (2, true),
                    NakedTriple => (3, false),
                    HiddenTriple => (3, true),
                    NakedQuad => (4, false),
                    HiddenQuad => (4, true),
                    _ => return false,
                };
                if cells.len() != n || digits.len() != n {
                    return false;
                }
                if n == 1 {
                    return x.placement
                        && cells == &[x.cell]
                        && digits == &x.ds
                        && if hidden {
                            x.locs(*h, digits[0]) == *cells
                        } else {
                            x.s.domains[x.cell] == m
                        };
                }
                !x.placement
                    && if hidden {
                        cells.contains(&x.cell)
                            && x.ds.iter().all(|&d| m & (1 << (d - 1)) == 0)
                            && digits.iter().all(|&d| {
                                let ls = x.locs(*h, d);
                                !ls.is_empty() && ls.iter().all(|c| cells.contains(c))
                            })
                    } else {
                        x.union(cells) == m
                            && !cells.contains(&x.cell)
                            && house(*h).contains(&x.cell)
                            && x.ds.iter().all(|&d| digits.contains(&d))
                    }
            }
            Fish { parts } => {
                if parts.is_empty() || parts.len() > 2 {
                    return false;
                }
                if t == SiameseFish {
                    if parts.len() != 2 {
                        return false;
                    }
                    let (a, b) = (&parts[0], &parts[1]);
                    return x.fish(a, true)
                        && x.fish(b, false)
                        && a.digit == b.digit
                        && a.bases == b.bases
                        && a.covers.len() == b.covers.len()
                        && a.covers.iter().filter(|h| !b.covers.contains(h)).count() == 1
                        && !a.fins.is_empty()
                        && !b.fins.is_empty()
                        && (0..81).any(|c| {
                            c != x.cell
                                && x.live(c, b.digit)
                                && b.covers.iter().any(|&h| house(h).contains(&c))
                                && !b.bases.iter().any(|&h| house(h).contains(&c))
                                && b.fins.iter().all(|&f| peers(c, f))
                        });
                }
                if parts.len() != 1 || !x.fish(&parts[0], true) {
                    return false;
                }
                let f = &parts[0];
                let n = f.bases.len();
                let rows = |hs: &[usize]| hs.iter().all(|&h| h < 9);
                let cols = |hs: &[usize]| hs.iter().all(|&h| (9..18).contains(&h));
                let basic =
                    (rows(&f.bases) && cols(&f.covers)) || (cols(&f.bases) && rows(&f.covers));
                let mixed = |hs: &[usize]| {
                    hs.iter().any(|&h| h < 9) && hs.iter().any(|&h| (9..18).contains(&h))
                };
                match t {
                    PointingPair => {
                        n == 1
                            && f.bases[0] >= 18
                            && f.covers[0] < 18
                            && f.fins.is_empty()
                            && (2..=3).contains(&x.locs(f.bases[0], f.digit).len())
                    }
                    BoxLineReduction => {
                        n == 1 && f.bases[0] < 18 && f.covers[0] >= 18 && f.fins.is_empty()
                    }
                    XWing | Swordfish | Jellyfish => {
                        basic
                            && f.fins.is_empty()
                            && n == match t {
                                XWing => 2,
                                Swordfish => 3,
                                _ => 4,
                            }
                    }
                    FinnedXWing | FinnedSwordfish | FinnedJellyfish => {
                        basic
                            && !f.fins.is_empty()
                            && f.fins.iter().all(|&c| {
                                c / 27 * 3 + c % 9 / 3 == f.fins[0] / 27 * 3 + f.fins[0] % 9 / 3
                            })
                            && n == match t {
                                FinnedXWing => 2,
                                FinnedSwordfish => 3,
                                _ => 4,
                            }
                    }
                    FrankenFish => {
                        n >= 2
                            && f.bases.iter().chain(&f.covers).any(|&h| h >= 18)
                            && !mixed(&f.bases)
                            && !mixed(&f.covers)
                    }
                    MutantFish => n >= 2 && (mixed(&f.bases) || mixed(&f.covers)),
                    EmptyRectangle => {
                        n == 2
                            && f.bases.iter().filter(|&&h| h >= 18).count() == 1
                            && f.bases.iter().filter(|&&h| h < 18).count() == 1
                            && !f.fins.is_empty()
                            && f.covers.iter().all(|&h| h < 18)
                            && {
                                let bx = *f.bases.iter().find(|&&h| h >= 18).unwrap();
                                let ls = x.locs(bx, f.digit);
                                (0..9).any(|r| {
                                    (0..9).any(|c| {
                                        ls.iter().all(|&o| o / 9 == r || o % 9 == c)
                                            && ls.iter().any(|&o| o / 9 == r && o % 9 != c)
                                            && ls.iter().any(|&o| o % 9 == c && o / 9 != r)
                                    })
                                })
                            }
                    }
                    _ => false,
                }
            }
            Wing { cells, digits } => {
                if !unique(cells, 81)
                    || mask(digits).is_none()
                    || !x.elim(*digits.last().unwrap_or(&0))
                    || cells.contains(&x.cell)
                    || cells.iter().any(|&c| x.s.values[c] != 0)
                {
                    return false;
                }
                let z = *digits.last().unwrap();
                if t == WXYZWing {
                    return cells.len() == 4
                        && digits.len() == 4
                        && x.union(cells) == mask(digits).unwrap()
                        && x.sees_all(cells, z)
                        && digits.iter().filter(|&&d| d != z).all(|&d| {
                            let ls: Vec<_> =
                                cells.iter().copied().filter(|&c| x.live(c, d)).collect();
                            ls.iter()
                                .all(|&a| ls.iter().all(|&b| a == b || peers(a, b)))
                        })
                        && {
                            let ls: Vec<_> =
                                cells.iter().copied().filter(|&c| x.live(c, z)).collect();
                            ls.iter()
                                .any(|&a| ls.iter().any(|&b| a != b && !peers(a, b)))
                        };
                }
                if !matches!(t, XYWing | XYZWing) || cells.len() != 3 || digits.len() != 3 {
                    return false;
                }
                let (a, b, p) = (cells[1], cells[2], cells[0]);
                let (d, e) = (digits[0], digits[1]);
                x.s.domains[a] == (1 << (d - 1)) | (1 << (z - 1))
                    && x.s.domains[b] == (1 << (e - 1)) | (1 << (z - 1))
                    && x.s.domains[p]
                        == (1 << (d - 1))
                            | (1 << (e - 1))
                            | if t == XYZWing { 1 << (z - 1) } else { 0 }
                    && peers(p, a)
                    && peers(p, b)
                    && peers(x.cell, a)
                    && peers(x.cell, b)
                    && (t == XYWing || peers(x.cell, p))
            }
            Als { sets, links, z } => {
                if sets.is_empty() || sets.len() > 8 {
                    return false;
                }
                let all: Vec<_> = sets.iter().flat_map(|a| a.cells.iter().copied()).collect();
                if !unique(&all, 81) {
                    return false;
                }
                if t == SueDeCoq {
                    if sets.len() != 3 || z.is_some() || links.len() != 4 {
                        return false;
                    }
                    let (i, b, l) = (&sets[0], &sets[1], &sets[2]);
                    return i.cells.len() == 2
                        && i.digits.len() == 4
                        && i.house < 18
                        && b.house >= 18
                        && l.house == i.house
                        && x.als(b)
                        && x.als(l)
                        && mask(&i.digits) == Some(x.union(&i.cells))
                        && i.cells.iter().all(|&c| {
                            x.s.values[c] == 0
                                && house(i.house).contains(&c)
                                && house(b.house).contains(&c)
                        })
                        && b.cells.iter().all(|c| !house(i.house).contains(c))
                        && l.cells.iter().all(|c| !house(b.house).contains(c))
                        && mask(&b.digits).unwrap() & mask(&l.digits).unwrap() == 0
                        && mask(&links[..2]) == mask(&b.digits)
                        && mask(&links[2..]) == mask(&l.digits)
                        && x.union(&all) == mask(&i.digits).unwrap()
                        && !x.placement
                        && !all.contains(&x.cell)
                        && x.ds.len() == 1
                        && local_refutes(x, &all);
                }
                if t == DeathBlossom {
                    if sets.len() < 3
                        || sets[0].cells.len() != 1
                        || links.len() + 1 != sets.len()
                        || mask(links) != Some(x.s.domains[sets[0].cells[0]])
                        || sets[0].house >= 27
                        || !house(sets[0].house).contains(&sets[0].cells[0])
                        || mask(&sets[0].digits) != Some(x.s.domains[sets[0].cells[0]])
                    {
                        return false;
                    }
                    let Some(z) = z else {
                        return false;
                    };
                    return x.elim(*z)
                        && !all.contains(&x.cell)
                        && sets.iter().skip(1).enumerate().all(|(i, a)| {
                            x.als(a)
                                && links[i] != *z
                                && a.digits.contains(&links[i])
                                && a.digits.contains(z)
                                && a.cells
                                    .iter()
                                    .filter(|&&c| x.live(c, links[i]))
                                    .all(|&c| peers(c, sets[0].cells[0]))
                                && x.sees_all(&a.cells, *z)
                        });
                }
                let Some(z) = z else {
                    return false;
                };
                let count_ok = match t {
                    AlsXz => sets.len() == 2,
                    AlsXyWing => sets.len() == 3,
                    AlsChain => sets.len() >= 4,
                    _ => false,
                };
                count_ok
                    && links.len() + 1 == sets.len()
                    && x.elim(*z)
                    && !all.contains(&x.cell)
                    && sets.iter().all(|a| x.als(a))
                    && links.iter().enumerate().all(|(i, &d)| {
                        d != *z && (i == 0 || links[i - 1] != d) && x.rcc(&sets[i], &sets[i + 1], d)
                    })
                    && x.sees_all(&sets[0].cells, *z)
                    && x.sees_all(&sets.last().unwrap().cells, *z)
            }
            Chain {
                nodes,
                strong,
                colors,
            } => {
                if nodes.len() < 2
                    || nodes.len() > 729
                    || nodes
                        .iter()
                        .any(|&n| n >= 729 || !x.live(usize::from(n / 9), (n % 9 + 1) as u8))
                    || nodes.iter().copied().collect::<HashSet<_>>().len() != nodes.len()
                {
                    return false;
                }
                if t == ThreeDMedusa {
                    if colors.len() != nodes.len()
                        || !strong.is_empty()
                        || nodes.iter().map(|n| n % 9).collect::<HashSet<_>>().len() < 2
                    {
                        return false;
                    }
                    let mut seen = HashSet::from([0]);
                    loop {
                        let old = seen.len();
                        for i in 0..nodes.len() {
                            for j in 0..i {
                                if x.link(nodes[i], nodes[j], true) {
                                    if colors[i] == colors[j] {
                                        return false;
                                    }
                                    if seen.contains(&i) || seen.contains(&j) {
                                        seen.insert(i);
                                        seen.insert(j);
                                    }
                                }
                            }
                        }
                        if seen.len() == old {
                            break;
                        }
                    }
                    if seen.len() != nodes.len() {
                        return false;
                    }
                    let viable = |color: bool| {
                        !(0..nodes.len()).any(|i| {
                            (0..i).any(|j| {
                                colors[i] == color
                                    && colors[j] == color
                                    && Context::conflict(nodes[i], nodes[j])
                            })
                        })
                    };
                    let mut count = 0;
                    for color in [false, true] {
                        if !viable(color) {
                            continue;
                        }
                        count += 1;
                        for &d in &x.ds {
                            let n = (x.cell * 9 + usize::from(d - 1)) as u16;
                            let proved = if x.placement {
                                (0..nodes.len()).any(|i| colors[i] == color && nodes[i] == n)
                            } else {
                                (0..nodes.len())
                                    .any(|i| colors[i] == color && Context::conflict(nodes[i], n))
                            };
                            if !proved {
                                return false;
                            }
                        }
                    }
                    return count > 0;
                }
                if !matches!(t, XChain | AIC | WWing)
                    || !colors.is_empty()
                    || strong.len() + 1 != nodes.len()
                    || strong
                        .iter()
                        .enumerate()
                        .any(|(i, &s)| s != (i % 2 == 0) || !x.link(nodes[i], nodes[i + 1], s))
                    || strong.len() % 2 != 1
                    || x.placement
                {
                    return false;
                }
                if t == XChain && nodes.iter().any(|n| n % 9 != nodes[0] % 9) {
                    return false;
                }
                if t == WWing
                    && (nodes.len() != 6
                        || nodes[0] / 9 != nodes[1] / 9
                        || nodes[4] / 9 != nodes[5] / 9
                        || nodes[0] % 9 != nodes[5] % 9
                        || nodes[1] % 9 != nodes[2] % 9
                        || nodes[2] % 9 != nodes[3] % 9
                        || nodes[3] % 9 != nodes[4] % 9
                        || x.s.domains[usize::from(nodes[0] / 9)]
                            != x.s.domains[usize::from(nodes[5] / 9)])
                {
                    return false;
                }
                x.ds.iter().all(|&d| {
                    let n = (x.cell * 9 + usize::from(d - 1)) as u16;
                    Context::conflict(n, nodes[0]) && Context::conflict(n, *nodes.last().unwrap())
                })
            }
            Trade { cells, digits } => self.trade(x, grid, cells, digits),
            Bug => {
                t == BivalueUniversalGrave && x.placement && self.unique_original(x, grid) && {
                    let mut ds = x.s.domains.clone();
                    ds[x.cell] &= !(1 << (x.ds[0] - 1));
                    (0..81).all(|c| x.s.values[c] != 0 || ds[c].count_ones() == 2)
                        && (0..27).all(|h| {
                            (1..=9).all(|d| {
                                let cs = house(h);
                                let placed = cs.iter().any(|&c| x.s.values[c] == d);
                                cs.iter()
                                    .filter(|&&c| x.s.values[c] == 0 && ds[c] & (1 << (d - 1)) != 0)
                                    .count()
                                    == if placed { 0 } else { 2 }
                            })
                        })
                }
            }
            Aligned { cells, exclusions } => {
                let n = match t {
                    AlignedPairExclusion => 2,
                    AlignedTripletExclusion => 3,
                    _ => return false,
                };
                cells.len() == n
                    && unique(cells, 81)
                    && cells.contains(&x.cell)
                    && !x.placement
                    && cells.iter().all(|&c| {
                        x.s.values[c] == 0 && cells.iter().all(|&o| c == o || peers(c, o))
                    })
                    && !exclusions.is_empty()
                    && exclusions.len() <= 32
                    && exclusions.iter().all(|cs| {
                        unique(cs, 81)
                            && cs.len() <= 4
                            && cs.iter().all(|&c| {
                                x.s.values[c] == 0
                                    && !cells.contains(&c)
                                    && cs.iter().all(|&o| c == o || peers(c, o))
                            })
                    })
                    && aligned_refutes(x, cells, exclusions)
            }
            Forcing {
                source,
                fish,
                branches,
            } => {
                use super::solver::explain::ForcingSource;
                if branches.is_empty() || branches.len() > 9 {
                    return false;
                }
                let assumptions: Vec<_> = branches.iter().map(|b| b.assumption).collect();
                if assumptions.iter().copied().collect::<HashSet<_>>().len() != assumptions.len() {
                    return false;
                }
                let mut expected = match (t, source) {
                    (
                        NishioForcingChain | DynamicForcingChain,
                        ForcingSource::Nishio { cell, digit },
                    ) if *cell == x.cell && x.elim(*digit) => {
                        vec![(*cell * 9 + usize::from(*digit - 1)) as u16]
                    }
                    (CellForcingChain | DynamicForcingChain, ForcingSource::Cell(c))
                        if *c < 81 && x.s.values[*c] == 0 && x.s.domains[*c].count_ones() >= 2 =>
                    {
                        x.requirement(*c as u16)
                    }
                    (RegionForcingChain, ForcingSource::Region { sector, digit })
                        if *sector < 27 && (1..=9).contains(digit) =>
                    {
                        x.requirement(81 + (*sector as u16) * 9 + u16::from(*digit - 1))
                    }
                    (KrakenFish, ForcingSource::Region { sector, digit }) => {
                        let Some(f) = fish else {
                            return false;
                        };
                        if !x.fish(f, false)
                            || f.bases.len() < 2
                            || *sector != f.bases[0]
                            || *digit != f.digit
                            || f.fins.is_empty()
                            || !x.elim(f.digit)
                            || f.bases.iter().any(|&h| house(h).contains(&x.cell))
                            || !f.covers.iter().any(|&h| house(h).contains(&x.cell))
                        {
                            return false;
                        }
                        f.fins
                            .iter()
                            .map(|&c| (c * 9 + usize::from(f.digit - 1)) as u16)
                            .collect()
                    }
                    _ => return false,
                };
                expected.sort_unstable();
                let mut actual = assumptions;
                actual.sort_unstable();
                expected == actual
                    && (t == KrakenFish || fish.is_none())
                    && branches.iter().all(|b| {
                        x.branch(
                            b,
                            if t == NishioForcingChain {
                                Some(x.ds[0])
                            } else {
                                None
                            },
                        )
                    })
            }
            Native { proof } => match (t, proof.as_ref()) {
                (ArithmeticCounting, ProofCertificate::Arithmetic(p)) => {
                    x.ds == [p.digit]
                        && x.cell == p.cell
                        && x.placement == p.value
                        && p.verify(grid)
                }
                (ThreeHouseMatching, ProofCertificate::ThreeHouse(p)) => {
                    x.elim(p.digit) && x.cell == p.cell && p.verify(grid)
                }
                (GlobalColoring, ProofCertificate::GlobalColoring(p)) => {
                    x.ds == [p.digit]
                        && x.cell == p.cell
                        && x.placement == p.value
                        && p.verify(grid)
                }
                _ => false,
            },
            Search { proof } => {
                t == Backtracking && proof.matches_hint(x.hint) && proof.verify(grid)
            }
        }
    }
    fn unique_original(&self, x: &Context<'_>, grid: &Grid) -> bool {
        let Some(u) = &self.uniqueness else {
            return false;
        };
        let Some(g) = Grid::from_string(&u.givens) else {
            return false;
        };
        if u.proof.purpose != ExactPurpose::Unique || !u.proof.verify(&g) {
            return false;
        }
        (0..81).all(|c| {
            let p = crate::Position::new(c / 9, c % 9);
            g.get(p)
                .is_none_or(|d| grid.cell(p).is_given() && x.s.values[c] == d)
                && grid.cell(p).is_given() == g.get(p).is_some()
                && (x.s.values[c] == 0 || x.s.values[c] == u.proof.witness[c])
                && x.s.domains[c] & (1 << (u.proof.witness[c] - 1)) != 0
        })
    }
    fn trade(&self, x: &Context<'_>, grid: &Grid, cells: &[usize], digits: &[u8]) -> bool {
        use Technique::*;
        if !matches!(
            self.technique,
            UniqueRectangle | HiddenRectangle | AvoidableRectangle | ExtendedUniqueRectangle
        ) || !unique(cells, 81)
            || !cells.contains(&x.cell)
            || x.placement
            || mask(digits).is_none()
            || !self.unique_original(x, grid)
            || cells
                .iter()
                .any(|&c| grid.cell(crate::Position::new(c / 9, c % 9)).is_given())
        {
            return false;
        }
        let rows: HashSet<_> = cells.iter().map(|c| c / 9).collect();
        let cols: HashSet<_> = cells.iter().map(|c| c % 9).collect();
        let boxes: HashSet<_> = cells.iter().map(|c| c / 27 * 3 + c % 9 / 3).collect();
        if self.technique == ExtendedUniqueRectangle {
            if cells.len() != 6
                || digits.len() != 3
                || boxes.len() != 3
                || !((rows.len() == 2 && cols.len() == 3) || (rows.len() == 3 && cols.len() == 2))
                || cells
                    .iter()
                    .filter(|&&c| x.s.domains[c] & !mask(digits).unwrap() != 0)
                    .count()
                    != 1
            {
                return false;
            }
        } else if cells.len() != 4
            || digits.len() != 2
            || rows.len() != 2
            || cols.len() != 2
            || boxes.len() != 2
        {
            return false;
        }
        if self.technique == AvoidableRectangle {
            if cells.iter().filter(|&&c| x.s.values[c] != 0).count() != 3 {
                return false;
            }
        } else if cells.iter().any(|&c| x.s.values[c] != 0) {
            return false;
        }
        if self.technique == HiddenRectangle
            && !cells
                .iter()
                .any(|&c| x.s.domains[c] == mask(digits).unwrap())
        {
            return false;
        }
        let core = mask(digits).unwrap();
        if self.technique == Technique::AvoidableRectangle {
            if cells
                .iter()
                .any(|&c| x.s.values[c] != 0 && !digits.contains(&x.s.values[c]))
            {
                return false;
            }
        } else if self.technique != Technique::ExtendedUniqueRectangle
            && (cells.iter().any(|&c| x.s.domains[c] & core != core)
                || !cells.iter().any(|&c| x.s.domains[c] == core))
        {
            return false;
        }
        let mut any = false;
        let mut ok = true;
        assignments(x, cells, &mut |vs| {
            if !opposite_holds(x, cells, vs) || !local_houses(x, cells, vs) {
                return;
            }
            any = true;
            let swapped = (rows.len() == 2 && swap_valid(cells, vs, true))
                || (cols.len() == 2 && swap_valid(cells, vs, false));
            if !swapped {
                ok = false;
            }
        });
        any && ok
    }
}

fn opposite_holds(x: &Context<'_>, cs: &[usize], vs: &[u8]) -> bool {
    cs.iter()
        .position(|&c| c == x.cell)
        .is_none_or(|i| x.ds.contains(&vs[i]) != x.placement)
}
fn assignments(x: &Context<'_>, cs: &[usize], visit: &mut impl FnMut(&[u8])) {
    fn go(x: &Context<'_>, cs: &[usize], vs: &mut Vec<u8>, visit: &mut impl FnMut(&[u8])) {
        if vs.len() == cs.len() {
            visit(vs);
            return;
        }
        let i = vs.len();
        for d in 1..=9 {
            if x.s.domains[cs[i]] & (1 << (d - 1)) != 0
                && (0..i).all(|j| vs[j] != d || !peers(cs[i], cs[j]))
            {
                vs.push(d);
                go(x, cs, vs, visit);
                vs.pop();
            }
        }
    }
    go(x, cs, &mut vec![], visit);
}
fn local_houses(x: &Context<'_>, cs: &[usize], vs: &[u8]) -> bool {
    (0..27).all(|h| {
        (1..=9).all(|d| {
            if house(h).iter().any(|&c| x.s.values[c] == d) {
                return true;
            }
            let ls = x.locs(h, d);
            !ls.iter().all(|c| cs.contains(c))
                || cs
                    .iter()
                    .enumerate()
                    .any(|(i, c)| ls.contains(c) && vs[i] == d)
        })
    })
}
fn local_refutes(x: &Context<'_>, cs: &[usize]) -> bool {
    let mut any = false;
    assignments(x, cs, &mut |vs| {
        if !opposite_holds(x, cs, vs) {
            return;
        }
        if vs
            .iter()
            .enumerate()
            .all(|(i, &d)| cs[i] == x.cell || !peers(cs[i], x.cell) || !x.ds.contains(&d))
        {
            any = true;
        }
    });
    !any
}
fn aligned_refutes(x: &Context<'_>, cs: &[usize], ex: &[Vec<usize>]) -> bool {
    let mut bad = false;
    assignments(x, cs, &mut |vs| {
        if !opposite_holds(x, cs, vs) {
            return;
        }
        let excluded = ex.iter().any(|es| {
            let mut m = 0u16;
            for &c in es {
                let mut ds = x.s.domains[c];
                for (i, &o) in cs.iter().enumerate() {
                    if peers(c, o) {
                        ds &= !(1 << (vs[i] - 1));
                    }
                }
                m |= ds;
            }
            m.count_ones() < (es.len() as u32)
        });
        if !excluded {
            bad = true;
        }
    });
    !bad
}
fn swap_valid(cs: &[usize], vs: &[u8], rows: bool) -> bool {
    let mut swapped = vec![];
    for &c in cs {
        let Some(i) = cs.iter().position(|&o| {
            if rows {
                o % 9 == c % 9 && o / 9 != c / 9
            } else {
                o / 9 == c / 9 && o % 9 != c % 9
            }
        }) else {
            return false;
        };
        swapped.push(vs[i]);
    }
    swapped != vs
        && (0..27).all(|h| {
            let hs = house(h);
            let mut a: Vec<_> = cs
                .iter()
                .enumerate()
                .filter(|(_, c)| hs.contains(c))
                .map(|(i, _)| vs[i])
                .collect();
            let mut b: Vec<_> = cs
                .iter()
                .enumerate()
                .filter(|(_, c)| hs.contains(c))
                .map(|(i, _)| swapped[i])
                .collect();
            a.sort_unstable();
            b.sort_unstable();
            a == b
        })
}
