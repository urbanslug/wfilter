pub type AwfOffset = usize;

struct AffinePenalties {
    a_match: i32,       // (Penalty representation; usually M <= 0) // match is a rust_keyword
    mismatch: i32,      // (Penalty representation; usually X > 0)
    gap_opening: i32,   // (Penalty representation; usually O > 0)
    gap_extension: i32, // (Penalty representation; usually E > 0)
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
struct AffineWavefrontsPenalties {
    base_penalties: AffinePenalties,                 // Input base Gap-Affine penalties
    wavefront_penalties: AffinePenalties,            // Wavefront Gap-Affine penalties
    penalties_strategy: WavefrontsPenaltiesStrategy, // Penalties adaptation strategy
}

impl AffineWavefrontsPenalties {
    fn init(&mut self,
            penalties: &AffinePenalties,
            penalties_strategy: WavefrontsPenaltiesStrategy) {

        self.base_penalties = *penalties;

        if penalties.a_match == 0 {
            self.penalties_strategy = WavefrontsPenaltiesStrategy::MatchZero;
        } else {
            self.penalties_strategy = penalties_strategy;
        };

        match (self.penalties_strategy) {
            WavefrontsPenaltiesStrategy::MatchZero | WavefrontsPenaltiesStrategy::ForceZeroMatch => {
                affine_penalties_mzero(penalties, &self.wavefront_penalties);
            },
            WavefrontsPenaltiesStrategy::ShiftedPenalties => {
                affine_penalties_shift(penalties, &(self.wavefront_penalties), false);
            },
            // WavefrontsPenaltiesStrategy::OddPairPenalties
            _ => affine_penalties_shift(penalties, &(self.wavefront_penalties), true);
        };
    }
}

fn affine_penalties_mzero(base_penalties: AffinePenalties, shifted_penalties: &AffinePenalties) {
    // Check base penalties
    if base_penalties.match > 0 {
        eprintln!("Match score must be negative or zero (M={})",base_penalties.match);
        panic!();
    }

    if base_penalties.mismatch <= 0 || base_penalties.gap_opening <= 0 || base_penalties.gap_extension <= 0 {
        eprintln!("Mismatch/Gap scores must be strictly positive (X={},O={},E={})\n",
                 base_penalties.mismatch,
                  base_penalties.gap_opening,
                  base_penalties.gap_extension);
        panic!();
    }

    // Copy base penalties
    *shifted_penalties = *base_penalties;
    // Zero match score
    shifted_penalties.match = 0;
}

fn affine_penalties_shift(base_penalties: AffinePenalties,
                          shifted_penalties: &AffinePenalties
                          pair_odd_heuristic: bool) {
    // Check base penalties
    if base_penalties.match > 0 {
        eprintln!("Match score must be negative or zero (M={})",base_penalties.match);
        panic!();
    }

    if (base_penalties.mismatch <= 0 ||
        base_penalties.gap_opening <= 0 ||
        base_penalties.gap_extension <= 0) {
        eprintln!("Mismatch/Gap scores must be strictly positive (X={},O={},E={})\n",
                  base_penalties.mismatch,
                  base_penalties.gap_opening,
                  base_penalties.gap_extension);
        panic!();
    }

    // Copy base penalties
    *shifted_penalties = *base_penalties;

    // Shift to zero match score
  shifted_penalties.match = 0;
  shifted_penalties.mismatch -= base_penalties.match;
  shifted_penalties.gap_opening -= base_penalties.match;
  shifted_penalties.gap_extension -= base_penalties.match;
  // Odd/Pair shift heuristic
  if (pair_odd_heuristic) {
    let is_mismatch_pair: bool = ((shifted_penalties.mismatch%2)==0);
    let is_gap_opening_pair: bool = ((shifted_penalties.gap_opening%2)==0);
    let is_gap_extension_pair: bool = ((shifted_penalties.gap_extension%2)==0);
    let total_odd: i32 = !is_mismatch_pair + !is_gap_opening_pair + !is_gap_extension_pair;
    let total_pair: i32 = is_mismatch_pair + is_gap_opening_pair + is_gap_extension_pair;
    if (total_odd > total_pair) {
      // Shift all to odd
      if (is_mismatch_pair) ++(shifted_penalties.mismatch);
      if (is_gap_opening_pair) ++(shifted_penalties.gap_opening);
      if (is_gap_extension_pair) ++(shifted_penalties.gap_extension);
    } else {
      // Shift all to pair
      if (!is_mismatch_pair) ++(shifted_penalties.mismatch);
      if (!is_gap_opening_pair) ++(shifted_penalties.gap_opening);
      if (!is_gap_extension_pair) ++(shifted_penalties.gap_extension);
    }
  }
}
