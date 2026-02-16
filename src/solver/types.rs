use crate::Position;
use serde::{Deserialize, Serialize};

/// Difficulty level of a puzzle
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub enum Difficulty {
    Beginner,
    Easy,
    Medium,
    Intermediate,
    Hard,
    Expert,
    Master,
    Extreme,
}

impl Difficulty {
    /// Get the maximum technique level allowed for this difficulty
    pub fn max_technique(&self) -> Technique {
        match self {
            Difficulty::Beginner => Technique::NakedSingle,
            Difficulty::Easy => Technique::NakedSingle,
            Difficulty::Medium => Technique::HiddenSingle,
            Difficulty::Intermediate => Technique::HiddenTriple,
            Difficulty::Hard => Technique::BoxLineReduction,
            Difficulty::Expert => Technique::HiddenRectangle,
            Difficulty::Master => Technique::BivalueUniversalGrave,
            Difficulty::Extreme => Technique::Backtracking,
        }
    }

    /// Check if this is a secret/locked difficulty
    pub fn is_secret(&self) -> bool {
        matches!(self, Difficulty::Master | Difficulty::Extreme)
    }

    /// Get all standard (non-secret) difficulties
    pub fn standard_levels() -> &'static [Difficulty] {
        &[
            Difficulty::Beginner,
            Difficulty::Easy,
            Difficulty::Medium,
            Difficulty::Intermediate,
            Difficulty::Hard,
            Difficulty::Expert,
        ]
    }

    /// Get all difficulties including secret ones
    pub fn all_levels() -> &'static [Difficulty] {
        &[
            Difficulty::Beginner,
            Difficulty::Easy,
            Difficulty::Medium,
            Difficulty::Intermediate,
            Difficulty::Hard,
            Difficulty::Expert,
            Difficulty::Master,
            Difficulty::Extreme,
        ]
    }

    /// SE rating range for this difficulty tier (min, max)
    pub fn se_range(&self) -> (f32, f32) {
        match self {
            Difficulty::Beginner => (1.5, 2.0),
            Difficulty::Easy => (2.0, 2.5),
            Difficulty::Medium => (2.5, 3.4),
            Difficulty::Intermediate => (3.4, 3.8),
            Difficulty::Hard => (3.8, 4.5),
            Difficulty::Expert => (4.5, 5.5),
            Difficulty::Master => (5.5, 7.0),
            Difficulty::Extreme => (7.0, 11.0),
        }
    }

    /// Default SE target (midpoint of range)
    pub fn se_default(&self) -> f32 {
        let (lo, hi) = self.se_range();
        (lo + hi) / 2.0
    }

    /// Brief description of techniques used at this tier
    pub fn technique_hint(&self) -> &'static str {
        match self {
            Difficulty::Beginner => "Hidden singles",
            Difficulty::Easy => "Naked singles",
            Difficulty::Medium => "Pairs & triples",
            Difficulty::Intermediate => "Hidden triples",
            Difficulty::Hard => "Box/line reduction",
            Difficulty::Expert => "Fish & rectangles",
            Difficulty::Master => "Wings & chains",
            Difficulty::Extreme => "Advanced techniques",
        }
    }
}

impl std::fmt::Display for Difficulty {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Difficulty::Beginner => write!(f, "Beginner"),
            Difficulty::Easy => write!(f, "Easy"),
            Difficulty::Medium => write!(f, "Medium"),
            Difficulty::Intermediate => write!(f, "Intermediate"),
            Difficulty::Hard => write!(f, "Hard"),
            Difficulty::Expert => write!(f, "Expert"),
            Difficulty::Master => write!(f, "Master"),
            Difficulty::Extreme => write!(f, "Extreme"),
        }
    }
}

/// Solving technique used (ordered by difficulty)
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub enum Technique {
    // Basic (Beginner/Easy)
    NakedSingle,
    HiddenSingle,

    // Intermediate
    NakedPair,
    HiddenPair,
    NakedTriple,
    HiddenTriple,

    // Hard
    PointingPair,
    BoxLineReduction,

    // Expert (fish family + quads + rectangles)
    XWing,
    FinnedXWing,
    Swordfish,
    FinnedSwordfish,
    Jellyfish,
    FinnedJellyfish,
    NakedQuad,
    HiddenQuad,
    EmptyRectangle,
    AvoidableRectangle,
    UniqueRectangle,
    HiddenRectangle,

    // Master (wings + chains + ALS-XZ + advanced patterns)
    XYWing,
    XYZWing,
    WXYZWing,
    WWing,
    XChain,
    /// Legacy: 3D Medusa is strong-link-only AIC coloring, subsumed by general AIC.
    /// Retained for SE rating compatibility (SE 5.0).
    ThreeDMedusa,
    SueDeCoq,
    AIC,
    FrankenFish,
    SiameseFish,
    AlsXz,
    ExtendedUniqueRectangle,
    BivalueUniversalGrave,

    // Extreme (ALS chains + advanced fish + forcing chains)
    AlsXyWing,
    AlsChain,
    MutantFish,
    /// Legacy: APE is subsumed by ALS chains. Retained for SE rating compatibility.
    AlignedPairExclusion,
    /// Legacy: ATE is subsumed by ALS chains. Retained for SE rating compatibility.
    AlignedTripletExclusion,
    DeathBlossom,
    NishioForcingChain,
    KrakenFish,
    RegionForcingChain,
    CellForcingChain,
    DynamicForcingChain,
    Backtracking,
}

impl Technique {
    /// Get the Sudoku Explainer (SE) numerical rating for this technique.
    /// This is the community-standard difficulty scale.
    pub fn se_rating(&self) -> f32 {
        match self {
            Technique::HiddenSingle => 1.5,
            Technique::NakedSingle => 2.3,
            Technique::PointingPair => 2.6,
            Technique::BoxLineReduction => 2.8,
            Technique::NakedPair => 3.0,
            Technique::XWing => 3.2,
            Technique::FinnedXWing => 3.4,
            Technique::HiddenPair => 3.4,
            Technique::NakedTriple => 3.6,
            Technique::Swordfish => 3.8,
            Technique::HiddenTriple => 3.8,
            Technique::FinnedSwordfish => 4.0,
            Technique::XYWing => 4.2,
            Technique::XYZWing => 4.4,
            Technique::WWing => 4.4,
            Technique::XChain => 4.5,
            Technique::EmptyRectangle => 4.6,
            Technique::AvoidableRectangle => 4.6,
            Technique::UniqueRectangle => 4.6,
            Technique::WXYZWing => 4.6,
            Technique::HiddenRectangle => 4.7,
            Technique::NakedQuad => 5.0,
            Technique::ThreeDMedusa => 5.0,
            Technique::SueDeCoq => 5.0,
            Technique::Jellyfish => 5.2,
            Technique::FinnedJellyfish => 5.4,
            Technique::HiddenQuad => 5.4,
            Technique::AlsXz => 5.5,
            Technique::FrankenFish => 5.5,
            Technique::SiameseFish => 5.5,
            Technique::ExtendedUniqueRectangle => 5.5,
            Technique::BivalueUniversalGrave => 5.6,
            Technique::AIC => 6.0,
            Technique::AlignedPairExclusion => 6.2,
            Technique::MutantFish => 6.5,
            Technique::AlsXyWing => 7.0,
            Technique::AlsChain => 7.5,
            Technique::AlignedTripletExclusion => 7.5,
            Technique::NishioForcingChain => 7.5,
            Technique::KrakenFish => 8.0,
            Technique::CellForcingChain => 8.3,
            Technique::DeathBlossom => 8.5,
            Technique::RegionForcingChain => 8.5,
            Technique::DynamicForcingChain => 9.3,
            Technique::Backtracking => 11.0,
        }
    }
}

impl std::fmt::Display for Technique {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Technique::NakedSingle => write!(f, "Naked Single"),
            Technique::HiddenSingle => write!(f, "Hidden Single"),
            Technique::NakedPair => write!(f, "Naked Pair"),
            Technique::HiddenPair => write!(f, "Hidden Pair"),
            Technique::NakedTriple => write!(f, "Naked Triple"),
            Technique::HiddenTriple => write!(f, "Hidden Triple"),
            Technique::PointingPair => write!(f, "Pointing Pair"),
            Technique::BoxLineReduction => write!(f, "Box/Line Reduction"),
            Technique::XWing => write!(f, "X-Wing"),
            Technique::FinnedXWing => write!(f, "Finned X-Wing"),
            Technique::Swordfish => write!(f, "Swordfish"),
            Technique::FinnedSwordfish => write!(f, "Finned Swordfish"),
            Technique::Jellyfish => write!(f, "Jellyfish"),
            Technique::FinnedJellyfish => write!(f, "Finned Jellyfish"),
            Technique::NakedQuad => write!(f, "Naked Quad"),
            Technique::HiddenQuad => write!(f, "Hidden Quad"),
            Technique::EmptyRectangle => write!(f, "Empty Rectangle"),
            Technique::AvoidableRectangle => write!(f, "Avoidable Rectangle"),
            Technique::XYWing => write!(f, "XY-Wing"),
            Technique::XYZWing => write!(f, "XYZ-Wing"),
            Technique::WXYZWing => write!(f, "WXYZ-Wing"),
            Technique::WWing => write!(f, "W-Wing"),
            Technique::XChain => write!(f, "X-Chain"),
            Technique::ThreeDMedusa => write!(f, "3D Medusa"),
            Technique::SueDeCoq => write!(f, "Sue de Coq"),
            Technique::AIC => write!(f, "AIC"),
            Technique::FrankenFish => write!(f, "Franken Fish"),
            Technique::SiameseFish => write!(f, "Siamese Fish"),
            Technique::AlsXz => write!(f, "ALS-XZ"),
            Technique::AlsXyWing => write!(f, "ALS-XY-Wing"),
            Technique::AlsChain => write!(f, "ALS Chain"),
            Technique::UniqueRectangle => write!(f, "Unique Rectangle"),
            Technique::HiddenRectangle => write!(f, "Hidden Rectangle"),
            Technique::ExtendedUniqueRectangle => write!(f, "Extended Unique Rectangle"),
            Technique::MutantFish => write!(f, "Mutant Fish"),
            Technique::AlignedPairExclusion => write!(f, "Aligned Pair Exclusion"),
            Technique::AlignedTripletExclusion => write!(f, "Aligned Triplet Exclusion"),
            Technique::BivalueUniversalGrave => write!(f, "BUG+1"),
            Technique::DeathBlossom => write!(f, "Death Blossom"),
            Technique::NishioForcingChain => write!(f, "Nishio Forcing Chain"),
            Technique::KrakenFish => write!(f, "Kraken Fish"),
            Technique::RegionForcingChain => write!(f, "Region Forcing Chain"),
            Technique::CellForcingChain => write!(f, "Cell Forcing Chain"),
            Technique::DynamicForcingChain => write!(f, "Dynamic Forcing Chain"),
            Technique::Backtracking => write!(f, "Backtracking"),
        }
    }
}

/// Type of hint provided
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum HintType {
    /// Place this value in this cell
    SetValue { pos: Position, value: u8 },
    /// Remove these candidates from this cell
    EliminateCandidates { pos: Position, values: Vec<u8> },
}

/// A hint for the player
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Hint {
    /// The technique used to find this hint
    pub technique: Technique,
    /// The type of hint
    pub hint_type: HintType,
    /// Explanation of the hint
    pub explanation: String,
    /// Cells involved in the reasoning
    pub involved_cells: Vec<Position>,
    /// Structural proof certificate (not serialized — WASM rendering only)
    #[serde(skip)]
    pub proof: Option<super::explain::ProofCertificate>,
}
