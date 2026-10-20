//! Sudoku Core Engine
//!
//! A cross-platform Sudoku engine supporting puzzle generation, solving, hints,
//! and variant rules (X-Sudoku, Killer Sudoku, etc.)

mod bitset;
mod cell;
mod constraint;
mod diversity;
mod exact_cover;
mod generator;
mod grid;
mod hash;
mod position;
mod presentation;
mod puzzle_id;
mod solver;
mod technique_examples;

pub use bitset::BitSet;
pub use cell::Cell;
pub use constraint::{
    BoxConstraint, ColumnConstraint, Constraint, ConstraintBox, DiagonalConstraint,
    KillerCageConstraint, RowConstraint, ThermoConstraint,
};
pub use diversity::{DifficultyStats, DiversityAnalyzer, DiversityReport, TheoreticalEstimates};
pub use exact_cover::{
    ExactBranch, ExactCoverCheck, ExactCoverProof, ExactCoverSearchOptions, ExactCoverSearchResult,
    ExactNode, ExactPurpose,
};
pub use generator::{Generator, GeneratorConfig, SymmetryType};
pub use grid::{Grid, MoveError, ValidationResult};
pub use hash::{canonical_puzzle_hash, canonical_puzzle_hash_str};
pub use position::Position;
pub use presentation::{
    EvidenceStatus, FrameKind, GridScene, HighlightRole, SceneCandidate, SceneCell, SceneEvidence,
    SceneFrame, SceneHighlight, SceneLink, SceneLinkKind,
};
pub use puzzle_id::PuzzleId;
pub use solver::{
    AlsProofDescriptor, ArithmeticCheck, ArithmeticProof, ArithmeticReplay, ArithmeticRequirement,
    ArithmeticSearchOptions, ArithmeticSearchResult, ArithmeticState, ArithmeticTerm,
    ArithmeticTerminal, CandidateState, ColoringBranch, ColoringNode, Difficulty, DomainEdit,
    ForcingBranch, ForcingCheck, ForcingLeaf, ForcingNode, ForcingProof, ForcingPropagation,
    ForcingSearchOptions, ForcingSearchResult, ForcingSource, ForcingStep, GlobalColoringCheck,
    GlobalColoringProof, GlobalColoringReplay, GlobalColoringSearchOptions,
    GlobalColoringSearchResult, GuardianLiteral, GuardianProof, GuardianPropagation,
    HallObstruction, Hint, HintType, LinkType, Polarity, ProofCertificate, PropagationPolicy,
    Solver, SolverState, Technique, ThreeHouseCheck, ThreeHouseProof, ThreeHouseReplay,
    ThreeHouseSearchOptions, ThreeHouseSearchResult,
};
pub use technique_examples::{TechniqueExample, TechniqueExampleCorpus, TechniqueExamplePuzzle};

mod named_patterns;
pub use named_patterns::{
    NamedPatternProof, PatternBranch, PatternFish, PatternSet, PatternStep, PatternUniqueness,
    PatternWitness,
};
