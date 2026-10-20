//! Solver-owned masks and learned relations, isolated from player pencil marks.
use super::explain::ProofCertificate;
use super::global_coloring::{self, GuardianProof, GuardianPropagation};
use super::{
    CandidateState, ForcingSearchOptions, GlobalColoringSearchOptions, Hint, HintType,
    PropagationPolicy, Solver, Technique,
};
use crate::{Grid, Position};

#[derive(Debug, Clone)]
pub struct SolverState {
    grid: Grid,
    guardians: Vec<GuardianProof>,
    history: Vec<Grid>,
}

pub(super) fn validate(grid: &Grid) -> Result<(), String> {
    let state = CandidateState::capture(grid)?;
    if !grid.validate().is_valid {
        return Err("Solver state violates an installed constraint".into());
    }
    for p in grid.empty_positions() {
        if grid.get_candidates(p).as_raw() & !grid.compute_candidates(p).as_raw() != 0 {
            return Err("Solver candidates conflict with a placed value".into());
        }
    }
    for sector in 0..27 {
        if super::three_house::house_domain_support(&state.domains, sector).is_none() {
            return Err("Solver candidates have a Hall contradiction in a native house".into());
        }
    }
    Ok(())
}
impl SolverState {
    /// Start from placed values, rebuilding solver candidates independently of
    /// any player pencil marks stored on the supplied grid.
    pub fn new(grid: &Grid) -> Result<Self, String> {
        let mut working = grid.deep_clone();
        working.recalculate_candidates();
        Self::from_candidates(&working)
    }
    /// Explicit opt-in for already derived, retained candidate premises.
    /// Validation checks native consistency, not the provenance of eliminations.
    pub fn from_candidates(grid: &Grid) -> Result<Self, String> {
        validate(grid)?;
        Ok(Self {
            grid: grid.deep_clone(),
            guardians: vec![],
            history: vec![],
        })
    }
    pub fn grid(&self) -> &Grid {
        &self.grid
    }
    pub fn guardians(&self) -> &[GuardianProof] {
        &self.guardians
    }
    /// Accept only an independently verified relation bound to this exact state.
    pub fn retain_guardian(&mut self, proof: GuardianProof) -> Result<(), String> {
        if CandidateState::capture(&self.grid)? != proof.premise || !proof.verify() {
            return Err("Guardian relation does not match or verify".into());
        }
        if !self
            .guardians
            .iter()
            .any(|p| p.cells == proof.cells && p.guardians == proof.guardians)
        {
            self.guardians.push(proof);
        }
        Ok(())
    }
    pub fn learn_guardian(
        &mut self,
        options: &GlobalColoringSearchOptions,
    ) -> Result<Option<GuardianProof>, String> {
        if options.max_patterns > 16384 || options.max_nodes > 1_000_000 {
            return Err("Guardian search limits exceed supported bounds".into());
        }
        let proof = global_coloring::find_guardian(&self.grid, options);
        if let Some(proof) = &proof {
            self.retain_guardian(proof.clone())?;
        }
        Ok(proof)
    }
    fn relation_hint(&self) -> Result<Option<Hint>, String> {
        let state = CandidateState::capture(&self.grid)?;
        for relation in &self.guardians {
            if relation
                .guardians
                .iter()
                .all(|g| state.domains[g.cell] & (1 << (g.digit - 1)) == 0)
            {
                return Err("All literals of a retained guardian relation are false".into());
            }
            if let Some(f) = GuardianPropagation::finding(relation, &state) {
                return Ok(Some(f.to_hint()));
            }
        }
        Ok(None)
    }
    pub fn get_hint(&self, solver: &Solver) -> Result<Option<Hint>, String> {
        self.get_hint_with_options(
            solver,
            PropagationPolicy::Full,
            &ForcingSearchOptions::default(),
        )
    }
    pub fn get_hint_with_options(
        &self,
        solver: &Solver,
        policy: PropagationPolicy,
        forcing: &ForcingSearchOptions,
    ) -> Result<Option<Hint>, String> {
        self.select_hint(solver, policy, forcing, policy == PropagationPolicy::Full)
            .map(|(hint, _)| hint)
    }

    // Queries and steps use identical selection. A query returns self-contained
    // guardian evidence; a step additionally stores any newly verified relation.
    fn select_hint(
        &self,
        solver: &Solver,
        policy: PropagationPolicy,
        forcing: &ForcingSearchOptions,
        backtracking: bool,
    ) -> Result<(Option<Hint>, Option<GuardianProof>), String> {
        forcing.validate()?;
        let mut hint = self.relation_hint()?.or_else(|| {
            solver
                .dispatch(&self.grid, policy, forcing)
                .map(|f| f.to_hint())
        });
        let mut learned = None;
        if policy == PropagationPolicy::Full
            && hint
                .as_ref()
                .is_none_or(|h| h.technique >= Technique::NishioForcingChain)
        {
            learned =
                global_coloring::find_guardian(&self.grid, &GlobalColoringSearchOptions::default());
            if let Some(relation) = &learned {
                if let Some(f) =
                    GuardianPropagation::finding(relation, &CandidateState::capture(&self.grid)?)
                {
                    hint = Some(f.to_hint());
                }
            }
        }
        if hint.is_none() && backtracking {
            hint = super::backtrack::find_backtracking_hint(&self.grid).map(|f| f.to_hint());
        }
        Ok((hint, learned))
    }

    /// Apply one complete-ladder deduction while retaining masks and relations.
    pub fn step(&mut self, solver: &Solver) -> Result<Option<Hint>, String> {
        self.step_with_options(
            solver,
            PropagationPolicy::Full,
            &ForcingSearchOptions::default(),
        )
    }
    /// Apply a deduction under the selected dispatch policy and forcing budgets.
    pub fn step_with_options(
        &mut self,
        solver: &Solver,
        policy: PropagationPolicy,
        forcing: &ForcingSearchOptions,
    ) -> Result<Option<Hint>, String> {
        self.step_inner(solver, policy, forcing, policy == PropagationPolicy::Full)
    }
    pub(super) fn step_logical(&mut self, solver: &Solver) -> Result<Option<Hint>, String> {
        self.step_inner(
            solver,
            PropagationPolicy::Full,
            &ForcingSearchOptions::default(),
            false,
        )
    }
    fn step_inner(
        &mut self,
        solver: &Solver,
        policy: PropagationPolicy,
        forcing: &ForcingSearchOptions,
        backtracking: bool,
    ) -> Result<Option<Hint>, String> {
        let (hint, learned) = self.select_hint(solver, policy, forcing, backtracking)?;
        if let Some(proof) = learned {
            self.retain_guardian(proof)?;
        }
        let Some(hint) = hint else {
            return Ok(None);
        };
        let mut next = self.grid.deep_clone();
        match &hint.hint_type {
            HintType::SetValue { pos, value } => next.place_solver_value(*pos, *value),
            HintType::EliminateCandidates { pos, values } => {
                for &v in values {
                    next.cell_mut(*pos).remove_candidate(v);
                }
            }
        }
        validate(&next)?;
        self.validate_relations(&next)?;
        if let Some(ProofCertificate::GlobalColoring(p)) = &hint.proof {
            if let Some(relation) = GuardianProof::from_coloring(&self.grid, p) {
                self.retain_guardian(relation)?;
            }
        }
        self.history.push(self.grid.deep_clone());
        self.grid = next;
        Ok(Some(hint))
    }
    /// Apply an explicit checked forcing deduction to this exact state.
    pub fn apply_forcing(&mut self, proof: &super::ForcingProof) -> Result<(), String> {
        let mut next = self.grid.deep_clone();
        proof.apply(&mut next)?;
        self.validate_relations(&next)?;
        self.history.push(self.grid.deep_clone());
        self.grid = next;
        Ok(())
    }

    fn validate_relations(&self, grid: &Grid) -> Result<(), String> {
        let state = CandidateState::capture(grid)?;
        if self.guardians.iter().any(|r| {
            r.guardians
                .iter()
                .all(|g| state.domains[g.cell] & (1 << (g.digit - 1)) == 0)
        }) {
            return Err("Deduction would falsify a retained guardian relation".into());
        }
        Ok(())
    }

    /// Edits rebuild all solver masks and invalidate every learned relation.
    pub fn edit(&mut self, pos: Position, value: Option<u8>) -> Result<(), String> {
        let mut next = self.grid.deep_clone();
        match value {
            Some(v) => next.set_cell(pos, v),
            None => next.clear_cell(pos),
        }
        .map_err(|e| e.to_string())?;
        next.recalculate_candidates();
        validate(&next)?;
        self.history.push(self.grid.deep_clone());
        self.grid = next;
        self.guardians.clear();
        Ok(())
    }
    /// Restore the preceding candidate snapshot and invalidate learned relations.
    pub fn undo(&mut self) -> bool {
        self.guardians.clear();
        if let Some(grid) = self.history.pop() {
            self.grid = grid;
            true
        } else {
            false
        }
    }
}
