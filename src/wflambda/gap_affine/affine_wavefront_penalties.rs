#[derive(Clone, Copy)]
pub struct AffinePenalties {
    pub a_match: i32, // (Penalty representation; usually M <= 0) // match is a rust_keyword
    pub mismatch: i32, // (Penalty representation; usually X > 0)
    pub gap_opening: i32, // (Penalty representation; usually O > 0)
    pub gap_extension: i32, // (Penalty representation; usually E > 0)
}

/*
 * Wavefront Strategy
 */
pub enum WavefrontsPenaltiesStrategy {
    MatchZero,        // WavefrontsPenaltiesMatchZero
    ForceZeroMatch,   // WavefrontsPenaltiesForceZeroMatch
    ShiftedPenalties, // WavefrontsPenaltiesShiftedPenalties
    OddPairPenalties, // WavefrontsPenaltiesOddPairPenalties
}

/*
 * Wavefront Penalties
 */
pub struct AffineWavefrontsPenalties {
    base_penalties: AffinePenalties, // Input base Gap-Affine penalties
    wavefront_penalties: AffinePenalties, // Wavefront Gap-Affine penalties
    penalties_strategy: WavefrontsPenaltiesStrategy, // Penalties adaptation strategy
}

impl AffineWavefrontsPenalties {
    pub fn init(
        penalties: &AffinePenalties,
        mut penalties_strategy: WavefrontsPenaltiesStrategy,
    ) -> Self {
        // TODO: is this right?
        let base_penalties = *penalties;
        let mut wavefront_penalties = *penalties;

        if penalties.a_match == 0 {
            penalties_strategy = WavefrontsPenaltiesStrategy::MatchZero;
        } else {
            penalties_strategy = penalties_strategy;
        };

        match penalties_strategy {
            WavefrontsPenaltiesStrategy::MatchZero
            | WavefrontsPenaltiesStrategy::ForceZeroMatch => {
                affine_penalties_mzero(penalties, &mut wavefront_penalties);
            }
            WavefrontsPenaltiesStrategy::ShiftedPenalties => {
                affine_penalties_shift(penalties, &mut wavefront_penalties, false);
            }
            // WavefrontsPenaltiesStrategy::OddPairPenalties
            _ => {
                affine_penalties_shift(penalties, &mut wavefront_penalties, true);
            }
        };

        Self {
            base_penalties,
            wavefront_penalties,
            penalties_strategy,
        }
    }
}

fn affine_penalties_mzero(
    base_penalties: &AffinePenalties,
    shifted_penalties: &mut AffinePenalties,
) {
    // Check base penalties
    if base_penalties.a_match > 0 {
        eprintln!(
            "Match score must be negative or zero (M={})",
            base_penalties.a_match
        );
        panic!();
    }

    if base_penalties.mismatch <= 0
        || base_penalties.gap_opening <= 0
        || base_penalties.gap_extension <= 0
    {
        eprintln!(
            "Mismatch/Gap scores must be strictly positive (X={},O={},E={})\n",
            base_penalties.mismatch, base_penalties.gap_opening, base_penalties.gap_extension
        );
        panic!();
    }

    // Copy base penalties
    *shifted_penalties = *base_penalties;
    // Zero match score
    shifted_penalties.a_match = 0;
}

fn affine_penalties_shift(
    base_penalties: &AffinePenalties,
    shifted_penalties: &mut AffinePenalties,
    pair_odd_heuristic: bool,
) {
    // Check base penalties
    if base_penalties.a_match > 0 {
        eprintln!(
            "Match score must be negative or zero (M={})",
            base_penalties.a_match
        );
        panic!();
    }

    if base_penalties.mismatch <= 0
        || base_penalties.gap_opening <= 0
        || base_penalties.gap_extension <= 0
    {
        eprintln!(
            "Mismatch/Gap scores must be strictly positive (X={},O={},E={})\n",
            base_penalties.mismatch, base_penalties.gap_opening, base_penalties.gap_extension
        );
        panic!();
    }

    // Copy base penalties
    *shifted_penalties = *base_penalties;

    // Shift to zero match score
    shifted_penalties.a_match = 0;
    shifted_penalties.mismatch -= base_penalties.a_match;
    shifted_penalties.gap_opening -= base_penalties.a_match;
    shifted_penalties.gap_extension -= base_penalties.a_match;

    // Odd/Pair shift heuristic
    if pair_odd_heuristic {
        let is_mismatch_pair: bool = (shifted_penalties.mismatch % 2) == 0;
        let is_gap_opening_pair: bool = (shifted_penalties.gap_opening % 2) == 0;
        let is_gap_extension_pair: bool = (shifted_penalties.gap_extension % 2) == 0;
        let total_odd: i32 =
            !is_mismatch_pair as i32 + !is_gap_opening_pair as i32 + !is_gap_extension_pair as i32;
        let total_pair: i32 =
            is_mismatch_pair as i32 + is_gap_opening_pair as i32 + is_gap_extension_pair as i32;
        if total_odd > total_pair {
            // Shift all to odd
            if is_mismatch_pair {
                shifted_penalties.mismatch += 1;
            }
            if is_gap_opening_pair {
                shifted_penalties.gap_opening += 1;
            }
            if is_gap_extension_pair {
                shifted_penalties.gap_extension += 1;
            }
        } else {
            // Shift all to pair
            if !is_mismatch_pair {
                shifted_penalties.mismatch += 1;
            }
            if !is_gap_opening_pair {
                shifted_penalties.gap_opening += 1;
            }
            if !is_gap_extension_pair {
                shifted_penalties.gap_extension += 1;
            }
        }
    }
}
