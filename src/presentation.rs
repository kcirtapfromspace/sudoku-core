//! Renderer-neutral grid scenes. Construction checks evidence; rendering does not.
//!
//! Scenes preserve candidate masks and never modify the source grid. Exact native
//! certificates are checked against both their premises and the hint endpoint.
//! Other technique descriptions remain explicitly unverified metadata.

use crate::solver::fabric::sector_cells;
use crate::{
    ArithmeticRequirement, ColoringNode, ForcingLeaf, ForcingProof, ForcingStep, Grid, Hint,
    HintType, LinkType, Polarity, Position, ProofCertificate, Technique,
};
use serde::Serialize;

// This limits display output, not proof checking. Reserve a frame for the result.
const MAX_TRACE_FRAMES: usize = 255;

/// A portable view of one board cell. Candidate digits use 1..=9.
#[derive(Debug, Clone, PartialEq, Eq, Serialize)]
pub struct SceneCell {
    pub value: Option<u8>,
    pub given: bool,
    pub candidates: Vec<u8>,
}

/// Display roles are independent of any application's colors.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum HighlightRole {
    Pattern,
    Base,
    Cover,
    Fin,
    Floor,
    Roof,
    On,
    Off,
    Assumption,
    Elimination,
    Placement,
    Conflict,
}

/// A cell or candidate annotation. `digit: None` means the whole cell.
#[derive(Debug, Clone, PartialEq, Eq, Serialize)]
pub struct SceneHighlight {
    pub cell: usize,
    pub digit: Option<u8>,
    pub role: HighlightRole,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
pub struct SceneCandidate {
    pub cell: usize,
    pub digit: u8,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum SceneLinkKind {
    Strong,
    Weak,
}

/// Structural links from detector metadata; a link does not certify a deduction.
#[derive(Debug, Clone, PartialEq, Eq, Serialize)]
pub struct SceneLink {
    pub from: SceneCandidate,
    pub to: SceneCandidate,
    pub kind: SceneLinkKind,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum FrameKind {
    Pattern,
    Assumption,
    Propagation,
    Contradiction,
    Conclusion,
    Validation,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize)]
pub struct SceneFrame {
    pub kind: FrameKind,
    pub title: String,
    pub description: String,
    /// Exactly 81 cells in row-major order, including branch-local domains.
    pub cells: Vec<SceneCell>,
    pub highlights: Vec<SceneHighlight>,
    pub links: Vec<SceneLink>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum EvidenceStatus {
    CheckedCertificate,
    CheckedDeduction,
    CheckedPatternAndDeduction,
    MetadataOnly,
    Rejected,
    RulesSatisfied,
    RuleViolation,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize)]
pub struct SceneEvidence {
    pub status: EvidenceStatus,
    pub message: String,
}

/// Version-one scene format for native, web and mobile renderers.
///
/// Serialized scenes are presentation output, not proof certificates. Verify
/// the original certificate in Rust before trusting a received scene's status.
#[derive(Debug, Clone, Serialize)]
pub struct GridScene {
    pub version: u8,
    pub technique: Option<Technique>,
    pub title: String,
    pub evidence: SceneEvidence,
    pub frames: Vec<SceneFrame>,
    /// True if a checked tree exceeds the display limit. Checking is never cut short.
    pub trace_truncated: bool,
}

fn pos(cell: usize) -> Position {
    Position::new(cell / 9, cell % 9)
}

fn snapshot(grid: &Grid) -> Result<Vec<SceneCell>, String> {
    (0..81)
        .map(|cell| {
            let c = grid.cell(pos(cell));
            if c.value().is_some_and(|v| !(1..=9).contains(&v))
                || c.candidates().as_raw() & !1022 != 0
            {
                return Err("Scene contains a digit outside 1..9".into());
            }
            Ok(SceneCell {
                value: c.value(),
                given: c.is_given(),
                candidates: c
                    .value()
                    .map_or_else(|| c.candidates().iter().collect(), |v| vec![v]),
            })
        })
        .collect()
}

fn mark(frame: &mut SceneFrame, cell: usize, digit: Option<u8>, role: HighlightRole) {
    // Detector metadata can come from callers. Never export invalid coordinates.
    if cell < 81 && digit.is_none_or(|d| (1..=9).contains(&d)) {
        let h = SceneHighlight { cell, digit, role };
        if !frame.highlights.contains(&h) {
            frame.highlights.push(h);
        }
    }
}

fn house(frame: &mut SceneFrame, sector: usize, digit: Option<u8>, role: HighlightRole) {
    if sector < 27 {
        for cell in sector_cells(sector) {
            if digit.is_none_or(|d| frame.cells[cell].candidates.contains(&d)) {
                mark(frame, cell, digit, role);
            }
        }
    }
}

fn triple(frame: &mut SceneFrame, center: usize, role: HighlightRole) {
    if center < 81 {
        for sector in [center / 9, 9 + center % 9, 18 + pos(center).box_index()] {
            house(frame, sector, None, role);
        }
    }
}

fn endpoint(hint: &Hint) -> Result<(usize, Vec<u8>, bool), String> {
    let (p, digits, value) = match &hint.hint_type {
        HintType::SetValue { pos, value } => (*pos, vec![*value], true),
        HintType::EliminateCandidates { pos, values } => (*pos, values.clone(), false),
    };
    if p.row >= 9
        || p.col >= 9
        || digits.is_empty()
        || digits.iter().any(|d| !(1..=9).contains(d))
        || digits
            .iter()
            .enumerate()
            .any(|(i, d)| digits[..i].contains(d))
        || hint.involved_cells.iter().any(|p| p.row >= 9 || p.col >= 9)
    {
        return Err("Invalid hint position or candidate digits".into());
    }
    Ok((p.row * 9 + p.col, digits, value))
}

fn evidence(grid: &Grid, hint: &Hint, cell: usize, digits: &[u8], value: bool) -> SceneEvidence {
    if let Some(ProofCertificate::Certified {
        pattern,
        deduction,
        named,
    }) = &hint.proof
    {
        let checked = pattern
            .as_deref()
            .is_none_or(|p| !matches!(p, ProofCertificate::Certified { .. }))
            && deduction.matches_hint(hint)
            && deduction.verify(grid)
            && named.as_deref().is_none_or(|p| p.verify(grid, hint));
        return SceneEvidence {
            status: if !checked {EvidenceStatus::Rejected} else if named.is_some() {EvidenceStatus::CheckedPatternAndDeduction} else {EvidenceStatus::CheckedDeduction},
            message: if !checked {"Certificate rejected: the state, named pattern or deduction does not verify."}
                else if named.is_some() {"Independent checkers verified the named pattern and its exact deduction from these candidates. Complete example bundles also authenticate the original puzzle and setup premises."}
                else {"The independent exact-cover checker verified this deduction. The technique name and pattern remain detector descriptions."}.into(),
        };
    }
    let matches = |c, d, v| c == cell && digits == [d] && v == value;
    let checked = match &hint.proof {
        Some(ProofCertificate::Arithmetic(p)) => {
            Some(matches(p.cell, p.digit, p.value) && p.verify(grid))
        }
        Some(ProofCertificate::ThreeHouse(p)) => {
            Some(matches(p.cell, p.digit, false) && p.verify(grid))
        }
        Some(ProofCertificate::GlobalColoring(p)) => {
            Some(matches(p.cell, p.digit, p.value) && p.verify(grid))
        }
        Some(ProofCertificate::ForcingProof(p)) => {
            Some(matches(p.cell, p.digit, p.value) && p.verify(grid))
        }
        Some(ProofCertificate::Guardian(p)) => {
            Some(matches(p.cell, p.digit, p.value) && p.verify(grid))
        }
        _ => None,
    };
    let (status, message) = match checked {
        Some(true) => (EvidenceStatus::CheckedCertificate,
            "The native-rule checker accepted the exact state and conclusion. Earlier imported eliminations remain premises."),
        Some(false) => (EvidenceStatus::Rejected,
            "The certificate failed verification or does not prove this hint's exact conclusion."),
        None => (EvidenceStatus::MetadataOnly,
            "Pattern metadata only. This scene does not independently verify the deduction or any uniqueness assumption."),
    };
    SceneEvidence {
        status,
        message: message.into(),
    }
}

impl GridScene {
    /// Render any technique through common cell/endpoint data, adding supplied
    /// structural metadata and checked native proof traces when available.
    pub fn from_hint(grid: &Grid, hint: &Hint) -> Result<Self, String> {
        let cells = snapshot(grid)?;
        let (cell, digits, value) = endpoint(hint)?;
        if cells[cell].value.is_some() || digits.iter().any(|d| !cells[cell].candidates.contains(d))
        {
            return Err("Hint endpoint is not a live candidate in this grid".into());
        }
        let mut scene = Self {
            version: 1,
            technique: Some(hint.technique),
            title: hint.technique.to_string(),
            evidence: evidence(grid, hint, cell, &digits, value),
            frames: vec![],
            trace_truncated: false,
        };
        if !grid.validate().is_valid
            || cells
                .iter()
                .any(|c| c.value.is_none() && c.candidates.is_empty())
        {
            scene.evidence = SceneEvidence {
                status: EvidenceStatus::Rejected,
                message: "The source grid has a rule violation or an empty candidate domain."
                    .into(),
            };
        }
        let mut after = grid.clone();
        if value {
            after.place_solver_value(pos(cell), digits[0]);
        } else {
            for &d in &digits {
                after.cell_mut(pos(cell)).remove_candidate(d);
            }
        }
        if !after.validate().is_valid
            || (0..81)
                .any(|c| after.get(pos(c)).is_none() && after.get_candidates(pos(c)).is_empty())
        {
            scene.evidence = SceneEvidence { status: EvidenceStatus::Rejected,
                message: "The proposed conclusion violates an installed rule or empties a candidate domain.".into() };
        }
        let mut pattern = SceneFrame {
            kind: FrameKind::Pattern,
            title: "Inspect the pattern".into(),
            description: hint.explanation.clone(),
            cells,
            highlights: vec![],
            links: vec![],
        };
        for p in &hint.involved_cells {
            mark(
                &mut pattern,
                p.row * 9 + p.col,
                None,
                HighlightRole::Pattern,
            );
        }
        if scene.evidence.status != EvidenceStatus::Rejected {
            annotate(&mut pattern, hint.proof.as_ref(), hint.technique);
        }
        for &d in &digits {
            mark(
                &mut pattern,
                cell,
                Some(d),
                if value {
                    HighlightRole::Placement
                } else {
                    HighlightRole::Elimination
                },
            );
        }
        scene.frames.push(pattern);
        if scene.evidence.status == EvidenceStatus::Rejected {
            let mut rejected = scene.frames[0].clone();
            rejected.kind = FrameKind::Validation;
            rejected.title = "Reject this evidence".into();
            rejected.description = scene.evidence.message.clone();
            for &d in &digits {
                mark(&mut rejected, cell, Some(d), HighlightRole::Conflict);
            }
            scene.frames.push(rejected);
            return Ok(scene);
        }
        if matches!(
            scene.evidence.status,
            EvidenceStatus::CheckedCertificate
                | EvidenceStatus::CheckedDeduction
                | EvidenceStatus::CheckedPatternAndDeduction
        ) {
            scene.proof_frames(grid, hint.proof.as_ref().unwrap(), cell, digits[0], value);
        }
        let mut conclusion = scene.frames[0].clone();
        conclusion.kind = FrameKind::Conclusion;
        conclusion.title = if matches!(
            scene.evidence.status,
            EvidenceStatus::CheckedCertificate
                | EvidenceStatus::CheckedDeduction
                | EvidenceStatus::CheckedPatternAndDeduction
        ) {
            "Apply the checked conclusion"
        } else {
            "Preview the proposed conclusion"
        }
        .into();
        conclusion.description = scene.evidence.message.clone();
        conclusion.cells = snapshot(&after)?;
        scene.frames.push(conclusion);
        Ok(scene)
    }

    /// Show installed-rule validation on the original grid, including conflicts.
    /// A valid partial grid is not a proof of solvability or uniqueness.
    pub fn from_validation(grid: &Grid) -> Result<Self, String> {
        let cells = snapshot(grid)?;
        let result = grid.validate();
        let mut frame = SceneFrame {
            kind: FrameKind::Validation,
            title: if result.is_valid {
                "No placed-value conflict"
            } else {
                "Rule violation"
            }
            .into(),
            description: if result.is_valid {
                "The installed rules accept the placed values. Solvability, candidate notes and uniqueness are not checked.".into()
            } else {
                result.violations.join(" ")
            },
            cells,
            highlights: vec![],
            links: vec![],
        };
        for p in result.invalid_cells {
            mark(&mut frame, p.row * 9 + p.col, None, HighlightRole::Conflict);
        }
        Ok(Self {
            version: 1,
            technique: None,
            title: "Board validation".into(),
            evidence: SceneEvidence {
                status: if result.is_valid {
                    EvidenceStatus::RulesSatisfied
                } else {
                    EvidenceStatus::RuleViolation
                },
                message: frame.description.clone(),
            },
            frames: vec![frame],
            trace_truncated: false,
        })
    }

    fn push_trace(&mut self, frame: SceneFrame) -> bool {
        if self.frames.len() >= MAX_TRACE_FRAMES {
            self.trace_truncated = true;
            return false;
        }
        self.frames.push(frame);
        true
    }

    fn proof_frames(
        &mut self,
        grid: &Grid,
        proof: &ProofCertificate,
        cell: usize,
        digit: u8,
        value: bool,
    ) {
        if let ProofCertificate::Certified { deduction, .. } = proof {
            self.exact_frames(deduction);
            return;
        }
        let mut assumption = self.frames[0].clone();
        assumption.kind = FrameKind::Assumption;
        assumption.title = "Assume the opposite conclusion".into();
        assumption.description = format!(
            "Assume r{}c{} {} {digit}. The certificate must refute this case.",
            cell / 9 + 1,
            cell % 9 + 1,
            if value { "is not" } else { "is" }
        );
        assumption.highlights.clear();
        assumption.links.clear();
        restrict(&mut assumption.cells, cell, digit, !value);
        mark(
            &mut assumption,
            cell,
            Some(digit),
            HighlightRole::Assumption,
        );
        self.push_trace(assumption.clone());
        match proof {
            ProofCertificate::ForcingProof(p) => self.forcing_frames(
                p,
                p.nodes.len() - 1,
                assumption,
                assumption_label(cell, digit, value),
            ),
            ProofCertificate::GlobalColoring(p) => self.coloring_frames(
                &p.cells,
                &p.nodes,
                p.nodes.len() - 1,
                assumption,
                assumption_label(cell, digit, value),
            ),
            ProofCertificate::Arithmetic(p) => {
                assumption.kind = FrameKind::Contradiction;
                assumption.title = "The source equations reject the assumption".into();
                assumption.description = match p.check(grid).unwrap() {
                    crate::ArithmeticCheck::Interval { residual, lower, upper } => format!("The equations require {residual}, outside the attainable interval [{lower}, {upper}]."),
                    crate::ArithmeticCheck::Divisibility { residual, divisor } => format!("The equations require {residual}, which is not divisible by {divisor}."),
                    crate::ArithmeticCheck::Residue { required_residue, modulus, reachable_residues, .. } => format!("The required remainder {required_residue} modulo {modulus} is absent from {reachable_residues:?}."),
                };
                mark(&mut assumption, cell, Some(digit), HighlightRole::Conflict);
                self.push_trace(assumption);
            }
            ProofCertificate::ThreeHouse(p) => {
                assumption.kind = FrameKind::Contradiction;
                assumption.title = "No boundary signature remains".into();
                assumption.description = format!("The checked cover has {} Hall obstructions. Every assignment to these three houses fails under the assumption.", p.blockers.len());
                triple(&mut assumption, p.center, HighlightRole::Pattern);
                mark(&mut assumption, cell, Some(digit), HighlightRole::Conflict);
                self.push_trace(assumption.clone());
                for blocker in &p.blockers {
                    let mut frame = assumption.clone();
                    frame.highlights.clear();
                    frame.title = "Inspect a Hall obstruction".into();
                    // Match the certificate's row-major atom convention. Each
                    // blocker covers the boundary signatures with this digit set.
                    let atom: Vec<_> = (0..81)
                        .filter(|&c| {
                            let row = c / 9 == p.center / 9;
                            let col = c % 9 == p.center % 9;
                            let same_box = pos(c).box_index() == pos(p.center).box_index();
                            match blocker.atom {
                                1 => row && same_box && c != p.center,
                                2 => col && same_box && c != p.center,
                                3 => row && !same_box,
                                4 => col && !same_box,
                                5 => same_box && !row && !col,
                                _ => false,
                            }
                        })
                        .collect();
                    let mut union = 0u16;
                    for (i, &c) in atom.iter().enumerate() {
                        frame.cells[c]
                            .candidates
                            .retain(|d| blocker.digits & (1 << (d - 1)) != 0);
                        mark(&mut frame, c, None, HighlightRole::Pattern);
                        if blocker.cells & (1 << i) != 0 {
                            for &d in &frame.cells[c].candidates {
                                union |= 1 << (d - 1);
                            }
                            mark(&mut frame, c, None, HighlightRole::Conflict);
                        }
                    }
                    frame.description = format!("In boundary cases assigning {:?} to this atom, {} highlighted cells have only {} available digits. This blocker rules out those cases.", digits(blocker.digits), blocker.cells.count_ones(), union.count_ones());
                    if !self.push_trace(frame) {
                        break;
                    }
                }
            }
            ProofCertificate::Guardian(p) => {
                assumption.kind = FrameKind::Contradiction;
                assumption.title = "The retained OR relation rejects the assumption".into();
                assumption.description = "The checked guardian relation requires a surviving literal that conflicts with this assumption.".into();
                for g in &p.relation.guardians {
                    mark(
                        &mut assumption,
                        g.cell,
                        Some(g.digit),
                        HighlightRole::Pattern,
                    );
                }
                mark(&mut assumption, cell, Some(digit), HighlightRole::Conflict);
                self.push_trace(assumption);
            }
            _ => {}
        }
    }

    fn exact_frames(&mut self, proof: &crate::ExactCoverProof) {
        let crate::ExactPurpose::Deduction {
            cell,
            digits,
            placement,
        } = &proof.purpose
        else {
            return;
        };
        let mut frame = self.frames[0].clone();
        frame.kind = FrameKind::Assumption;
        frame.title = "Independently assume the opposite".into();
        frame.description = format!(
            "Check r{}c{} {} {:?}. The independent proof has {} nodes.",
            cell / 9 + 1,
            cell % 9 + 1,
            if *placement { "outside" } else { "inside" },
            digits,
            proof.nodes.len()
        );
        frame.highlights.clear();
        frame.links.clear();
        frame.cells[*cell]
            .candidates
            .retain(|d| digits.contains(d) != *placement);
        for &d in digits {
            mark(&mut frame, *cell, Some(d), HighlightRole::Assumption);
        }
        self.push_trace(frame.clone());
        self.exact_node(&proof.nodes, proof.nodes.len() - 1, frame);
    }

    fn exact_node(&mut self, nodes: &[crate::ExactNode], id: usize, mut state: SceneFrame) {
        if self.frames.len() >= MAX_TRACE_FRAMES {
            self.trace_truncated = true;
            return;
        }
        match &nodes[id] {
            crate::ExactNode::Conflict { requirement } => {
                state.kind = FrameKind::Contradiction;
                state.title = "A native requirement has no support".into();
                state
                    .highlights
                    .retain(|h| h.role == HighlightRole::Assumption);
                let r = usize::from(*requirement);
                if r < 81 {
                    state.description = format!(
                        "Cell r{}c{} has no possible digit under these assumptions.",
                        r / 9 + 1,
                        r % 9 + 1
                    );
                    mark(&mut state, r, None, HighlightRole::Conflict);
                } else {
                    let sector = (r - 81) / 9;
                    let digit = ((r - 81) % 9 + 1) as u8;
                    let (kind, number) = if sector < 9 {
                        ("row", sector + 1)
                    } else if sector < 18 {
                        ("column", sector - 8)
                    } else {
                        ("box", sector - 17)
                    };
                    state.description = format!("Digit {digit} has no possible cell in {kind} {number} under these assumptions.");
                    house(&mut state, sector, None, HighlightRole::Conflict);
                }
                self.push_trace(state);
            }
            crate::ExactNode::Split { branches, .. } => {
                for branch in branches {
                    let mut next = state.clone();
                    let cell = usize::from(branch.candidate) / 9;
                    let digit = (branch.candidate % 9 + 1) as u8;
                    next.kind = FrameKind::Assumption;
                    next.title = "Check an independent proof branch".into();
                    next.description = format!("Assume r{}c{}={digit}. The checker requires a child for every available candidate in this native requirement.", cell / 9 + 1, cell % 9 + 1);
                    next.highlights
                        .retain(|h| h.role == HighlightRole::Assumption);
                    next.cells[cell].value = Some(digit);
                    next.cells[cell].given = false;
                    next.cells[cell].candidates = vec![digit];
                    mark(&mut next, cell, Some(digit), HighlightRole::Assumption);
                    for other in 0..81 {
                        if next.cells[other].value.is_none()
                            && peers(cell, other)
                            && next.cells[other].candidates.contains(&digit)
                        {
                            restrict(&mut next.cells, other, digit, false);
                            mark(&mut next, other, Some(digit), HighlightRole::Elimination);
                        }
                    }
                    if !self.push_trace(next.clone()) {
                        return;
                    }
                    self.exact_node(nodes, branch.child, next);
                    if self.trace_truncated {
                        return;
                    }
                }
            }
            crate::ExactNode::Solution => {}
        }
    }

    fn forcing_frames(
        &mut self,
        proof: &ForcingProof,
        id: usize,
        mut state: SceneFrame,
        path: String,
    ) {
        if self.frames.len() >= MAX_TRACE_FRAMES {
            self.trace_truncated = true;
            return;
        }
        let node = &proof.nodes[id];
        for step in &node.steps {
            let before = state.cells.clone();
            state.highlights.clear();
            let changes = match step {
                ForcingStep::House { sector, edits } => {
                    house(&mut state, *sector, None, HighlightRole::Pattern);
                    edits
                }
                ForcingStep::ThreeHouse { center, edits } => {
                    triple(&mut state, *center, HighlightRole::Pattern);
                    edits
                }
            };
            for edit in changes {
                state.cells[edit.cell].candidates = digits(edit.after);
                for d in &before[edit.cell].candidates {
                    if !state.cells[edit.cell].candidates.contains(d) {
                        mark(&mut state, edit.cell, Some(*d), HighlightRole::Elimination);
                    }
                }
            }
            state.kind = FrameKind::Propagation;
            state.title = "Propagate the branch".into();
            state.description = format!("{path} Apply the recorded native rule. Highlighted candidates are removed in this branch.").trim().into();
            if !self.push_trace(state.clone()) {
                return;
            }
        }
        if let Some(leaf) = &node.leaf {
            state.highlights.clear();
            state.kind = FrameKind::Contradiction;
            state.title = "This branch has no assignment".into();
            state.description = match leaf {
                ForcingLeaf::Hall { sector, cells } => {
                    house(&mut state, *sector, None, HighlightRole::Pattern);
                    let mut union = vec![];
                    for &c in cells {
                        for &d in &state.cells[c].candidates {
                            if !union.contains(&d) {
                                union.push(d);
                            }
                        }
                        mark(&mut state, c, None, HighlightRole::Conflict);
                    }
                    format!(
                        "{path} Hall obstruction: {} cells have only {} available digits.",
                        cells.len(),
                        union.len()
                    )
                }
                ForcingLeaf::ThreeHouse { center, blockers } => {
                    triple(&mut state, *center, HighlightRole::Conflict);
                    format!("{path} All boundary signatures fail. The checker accepted {} Hall obstructions.", blockers.len())
                }
            };
            self.push_trace(state);
        } else {
            let cell = node.cell.unwrap();
            for branch in &node.branches {
                let mut next = state.clone();
                next.kind = FrameKind::Assumption;
                next.title = "Check the next branch".into();
                let label = format!("r{}c{}={}", cell / 9 + 1, cell % 9 + 1, branch.digit);
                let path = if path.is_empty() {
                    label
                } else {
                    format!("{path} → {label}")
                };
                next.description =
                    format!("{path}. Every live candidate at this split has a checked child.");
                next.highlights.clear();
                restrict(&mut next.cells, cell, branch.digit, true);
                mark(
                    &mut next,
                    cell,
                    Some(branch.digit),
                    HighlightRole::Assumption,
                );
                if !self.push_trace(next.clone()) {
                    return;
                }
                self.forcing_frames(proof, branch.child, next, path);
                if self.trace_truncated {
                    return;
                }
            }
        }
    }

    fn coloring_frames(
        &mut self,
        selected: &[usize],
        nodes: &[ColoringNode],
        id: usize,
        mut state: SceneFrame,
        path: String,
    ) {
        if self.frames.len() >= MAX_TRACE_FRAMES {
            self.trace_truncated = true;
            return;
        }
        match &nodes[id] {
            ColoringNode::Conflict { cell } => {
                state.kind = FrameKind::Contradiction;
                state.title = "No candidate remains".into();
                state.description = format!(
                    "{path}. The selected peer graph leaves r{}c{} with an empty domain.",
                    cell / 9 + 1,
                    cell % 9 + 1
                );
                mark(&mut state, *cell, None, HighlightRole::Conflict);
                self.push_trace(state);
            }
            ColoringNode::Split { cell, branches } => {
                for branch in branches {
                    let mut next = state.clone();
                    next.kind = FrameKind::Assumption;
                    next.title = "Check a graph-coloring branch".into();
                    let label = format!("r{}c{}={}", cell / 9 + 1, cell % 9 + 1, branch.digit);
                    let path = if path.is_empty() {
                        label
                    } else {
                        format!("{path} → {label}")
                    };
                    next.description =
                        format!("{path}. Remove this digit from selected peer cells.");
                    next.highlights.clear();
                    restrict(&mut next.cells, *cell, branch.digit, true);
                    mark(
                        &mut next,
                        *cell,
                        Some(branch.digit),
                        HighlightRole::Assumption,
                    );
                    for &other in selected {
                        if peers(*cell, other)
                            && next.cells[other].candidates.contains(&branch.digit)
                        {
                            restrict(&mut next.cells, other, branch.digit, false);
                            mark(
                                &mut next,
                                other,
                                Some(branch.digit),
                                HighlightRole::Elimination,
                            );
                        }
                    }
                    if !self.push_trace(next.clone()) {
                        return;
                    }
                    self.coloring_frames(selected, nodes, branch.child, next, path);
                    if self.trace_truncated {
                        return;
                    }
                }
            }
        }
    }
}

fn assumption_label(cell: usize, digit: u8, conclusion: bool) -> String {
    format!(
        "r{}c{}{}{}",
        cell / 9 + 1,
        cell % 9 + 1,
        if conclusion { "≠" } else { "=" },
        digit
    )
}

fn digits(mask: u16) -> Vec<u8> {
    (1..=9).filter(|d| mask & (1 << (d - 1)) != 0).collect()
}

fn restrict(cells: &mut [SceneCell], cell: usize, digit: u8, value: bool) {
    cells[cell].candidates.retain(|&d| (d == digit) == value);
}

fn peers(a: usize, b: usize) -> bool {
    a != b && (a / 9 == b / 9 || a % 9 == b % 9 || pos(a).box_index() == pos(b).box_index())
}

// Independent local relation check; no inference-graph producer is called.
fn structural_link(
    cells: &[SceneCell],
    a: SceneCandidate,
    b: SceneCandidate,
    kind: SceneLinkKind,
) -> bool {
    if a == b
        || a.cell >= 81
        || b.cell >= 81
        || !(1..=9).contains(&a.digit)
        || !(1..=9).contains(&b.digit)
        || cells[a.cell].value.is_some()
        || cells[b.cell].value.is_some()
        || !cells[a.cell].candidates.contains(&a.digit)
        || !cells[b.cell].candidates.contains(&b.digit)
    {
        return false;
    }
    if a.cell == b.cell {
        return kind == SceneLinkKind::Weak || cells[a.cell].candidates.len() == 2;
    }
    if a.digit != b.digit || !peers(a.cell, b.cell) {
        return false;
    }
    if kind == SceneLinkKind::Weak {
        return true;
    }
    (0..27).any(|sector| {
        let members = sector_cells(sector);
        members.contains(&a.cell)
            && members.contains(&b.cell)
            && !members.iter().any(|&c| cells[c].value == Some(a.digit))
            && members
                .iter()
                .filter(|&&c| cells[c].value.is_none() && cells[c].candidates.contains(&a.digit))
                .count()
                == 2
    })
}

fn annotate(frame: &mut SceneFrame, proof: Option<&ProofCertificate>, technique: Technique) {
    match proof {
        Some(ProofCertificate::Certified { pattern, named, .. }) => {
            if let Some(p) = named {
                annotate_named(frame, &p.witness, technique);
            } else {
                annotate(frame, pattern.as_deref(), technique);
            }
        }
        Some(ProofCertificate::Fish {
            digit,
            base_sectors,
            cover_sectors,
            fins,
        }) => {
            for &s in base_sectors {
                house(frame, s, Some(*digit), HighlightRole::Base);
            }
            for &s in cover_sectors {
                house(frame, s, Some(*digit), HighlightRole::Cover);
            }
            for &c in fins {
                mark(frame, c, Some(*digit), HighlightRole::Fin);
            }
        }
        Some(ProofCertificate::Als { als_chain, .. }) => {
            for als in als_chain {
                for &c in &als.cells {
                    mark(frame, c, None, HighlightRole::Pattern);
                    for &d in &als.candidates {
                        if c < 81 && frame.cells[c].candidates.contains(&d) {
                            mark(frame, c, Some(d), HighlightRole::On);
                        }
                    }
                }
            }
        }
        Some(ProofCertificate::Aic { chain, link_types }) => {
            for &(c, d, polarity) in chain {
                mark(
                    frame,
                    c,
                    Some(d),
                    if polarity == Polarity::On {
                        HighlightRole::On
                    } else {
                        HighlightRole::Off
                    },
                );
            }
            #[allow(deprecated)]
            if technique == Technique::ThreeDMedusa {
                // Legacy Medusa metadata lists a component in hash-map order,
                // not a path. Reconstruct its actual strong edges from premises.
                let mut nodes: Vec<_> = chain.iter().map(|&(c, d, _)| (c, d)).collect();
                nodes.sort_unstable();
                nodes.dedup();
                for (i, &(c, d)) in nodes.iter().enumerate() {
                    for &(other, digit) in &nodes[i + 1..] {
                        let from = SceneCandidate { cell: c, digit: d };
                        let to = SceneCandidate { cell: other, digit };
                        if structural_link(&frame.cells, from, to, SceneLinkKind::Strong) {
                            frame.links.push(SceneLink {
                                from,
                                to,
                                kind: SceneLinkKind::Strong,
                            });
                        }
                    }
                }
                return;
            }
            for (pair, kind) in chain.windows(2).zip(link_types) {
                if pair.iter().all(|(c, d, _)| *c < 81 && (1..=9).contains(d)) {
                    let link = SceneLink {
                        from: SceneCandidate {
                            cell: pair[0].0,
                            digit: pair[0].1,
                        },
                        to: SceneCandidate {
                            cell: pair[1].0,
                            digit: pair[1].1,
                        },
                        kind: if *kind == LinkType::Strong {
                            SceneLinkKind::Strong
                        } else {
                            SceneLinkKind::Weak
                        },
                    };
                    if structural_link(&frame.cells, link.from, link.to, link.kind) {
                        frame.links.push(link);
                    }
                }
            }
        }
        Some(ProofCertificate::Uniqueness {
            floor_cells,
            roof_cells,
            ..
        }) => {
            for &c in floor_cells {
                mark(frame, c, None, HighlightRole::Floor);
            }
            for &c in roof_cells {
                mark(frame, c, None, HighlightRole::Roof);
            }
        }
        Some(ProofCertificate::Arithmetic(p)) => {
            for term in &p.terms {
                match term.requirement {
                    ArithmeticRequirement::Cell { cell } => {
                        mark(frame, cell, None, HighlightRole::Pattern)
                    }
                    ArithmeticRequirement::SectorDigit { sector, digit } => {
                        house(frame, sector, Some(digit), HighlightRole::Pattern)
                    }
                }
            }
        }
        Some(ProofCertificate::ThreeHouse(p)) => triple(frame, p.center, HighlightRole::Pattern),
        Some(ProofCertificate::GlobalColoring(p)) => {
            for &c in &p.cells {
                mark(frame, c, None, HighlightRole::Pattern);
            }
        }
        Some(ProofCertificate::Guardian(p)) => {
            for &c in &p.relation.cells {
                mark(frame, c, None, HighlightRole::Pattern);
            }
            for g in &p.relation.guardians {
                mark(frame, g.cell, Some(g.digit), HighlightRole::On);
            }
        }
        _ => {}
    }
}

fn annotate_named(frame: &mut SceneFrame, w: &crate::PatternWitness, technique: Technique) {
    use crate::PatternWitness::*;
    let cells = |frame: &mut SceneFrame, cs: &[usize], role| {
        for &c in cs {
            mark(frame, c, None, role);
        }
    };
    let fish = |frame: &mut SceneFrame, f: &crate::PatternFish| {
        for &h in &f.bases {
            house(frame, h, Some(f.digit), HighlightRole::Base);
        }
        for &h in &f.covers {
            house(frame, h, Some(f.digit), HighlightRole::Cover);
        }
        for &c in &f.fins {
            mark(frame, c, Some(f.digit), HighlightRole::Fin);
        }
    };
    match w {
        Subset {
            cells: cs, digits, ..
        } => {
            cells(frame, cs, HighlightRole::Pattern);
            for &c in cs {
                for &d in digits {
                    if c < 81 && frame.cells[c].candidates.contains(&d) {
                        mark(frame, c, Some(d), HighlightRole::On);
                    }
                }
            }
        }
        Fish { parts } => {
            for f in parts {
                fish(frame, f);
            }
        }
        Wing { cells: cs, digits } => {
            cells(frame, cs, HighlightRole::Cover);
            if let Some(&p) = cs.first() {
                mark(frame, p, None, HighlightRole::Base);
            }
            for &c in cs {
                if c < 81 {
                    for &d in digits {
                        if frame.cells[c].candidates.contains(&d) {
                            mark(frame, c, Some(d), HighlightRole::On);
                        }
                    }
                }
            }
        }
        Als { sets, .. } => {
            for a in sets {
                cells(frame, &a.cells, HighlightRole::Pattern);
            }
        }
        Chain {
            nodes,
            strong,
            colors,
        } => {
            for (i, &n) in nodes.iter().enumerate() {
                mark(
                    frame,
                    usize::from(n / 9),
                    Some((n % 9 + 1) as u8),
                    if colors.get(i).copied().unwrap_or(i % 2 == 1) {
                        HighlightRole::On
                    } else {
                        HighlightRole::Off
                    },
                );
            }
            for i in 0..nodes.len() {
                for j in i + 1..nodes.len() {
                    let a = SceneCandidate {
                        cell: usize::from(nodes[i] / 9),
                        digit: (nodes[i] % 9 + 1) as u8,
                    };
                    let b = SceneCandidate {
                        cell: usize::from(nodes[j] / 9),
                        digit: (nodes[j] % 9 + 1) as u8,
                    };
                    let kind = if colors.is_empty() {
                        if j != i + 1 {
                            continue;
                        }
                        if strong.get(i) == Some(&true) {
                            SceneLinkKind::Strong
                        } else {
                            SceneLinkKind::Weak
                        }
                    } else {
                        SceneLinkKind::Strong
                    };
                    if structural_link(&frame.cells, a, b, kind) {
                        frame.links.push(SceneLink {
                            from: a,
                            to: b,
                            kind,
                        });
                    }
                }
            }
        }
        Trade { cells: cs, digits } => {
            for &c in cs {
                if c < 81 {
                    let roof = frame.cells[c]
                        .candidates
                        .iter()
                        .any(|d| !digits.contains(d));
                    mark(
                        frame,
                        c,
                        None,
                        if roof {
                            HighlightRole::Roof
                        } else {
                            HighlightRole::Floor
                        },
                    );
                }
            }
        }
        Aligned {
            cells: cs,
            exclusions,
        } => {
            cells(frame, cs, HighlightRole::Base);
            for cs in exclusions {
                cells(frame, cs, HighlightRole::Cover);
            }
        }
        Forcing {
            fish: f, branches, ..
        } => {
            if let Some(f) = f {
                fish(frame, f);
            }
            for b in branches {
                mark(
                    frame,
                    usize::from(b.assumption / 9),
                    Some((b.assumption % 9 + 1) as u8),
                    HighlightRole::Assumption,
                );
            }
        }
        Native { proof } => annotate(frame, Some(proof), technique),
        Bug | Search { .. } => {}
    }
}
