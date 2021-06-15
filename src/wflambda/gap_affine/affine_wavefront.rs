use super::super::edit::edit_cigar;
use super::affine_wavefront_penalties as penalties;
use std::cmp;

pub type AwfOffset = usize;

pub fn affine_lambda_wavefront_v(k: i32, offset: i32) -> i32 {
    offset-k
}

pub fn affine_lambda_wavefront_h(k: i32, offset: i32) -> i32 {
    offset
}

pub fn affine_lambda_wavefront_diagonal(h: i32, v: i32) -> i32 {
    h - v
}
pub fn affine_lambda_wavefront_offset(h: i32, v: i32) -> i32  {
    h
}



pub struct AffineWavefront {
    // Range
    pub null: bool,     // Is null interval?
    pub lo: i32,      // Effective lowest diagonal (inclusive)
    pub hi: i32,      // Effective highest diagonal (inclusive)
    pub lo_base: i32, // Lowest diagonal before reduction (inclusive)
    pub hi_base: i32, // Highest diagonal before reduction (inclusive)

    // Offsets
    pub offsets: Vec<AwfOffset>

                                 /* TODO: make this an option type
                                 // #ifdef AFFINE_LAMBDA_WAVEFRONT_DEBUG
                                 // offsets_base: AwfOffset, // Offsets increment
                                 // #endif
                                  */
}

impl AffineWavefront {
    pub fn new(null: bool, lo: i32, hi: i32, lo_base: i32, hi_base: i32) -> Self {
        let offsets: Vec<AwfOffset> = Vec::new();
        Self {null, lo, hi, lo_base, hi_base, offsets}
    }
}

type WavefrontsStats = ();

#[derive(PartialEq)]
pub enum WavefrontReductionType {
    WavefrontsReductionNone,
    WavefrontsReductionDynamic,
}

#[derive(PartialEq)]
pub struct AffineWavefrontsReduction {
    pub reduction_strategy: WavefrontReductionType, // Reduction strategy
    pub min_wavefront_length: i32,                  // Dynamic: Minimum wavefronts length to reduce
    pub max_distance_threshold: i32,                // Dynamic: Maximum distance between offsets allowed
}

impl AffineWavefrontsReduction {
    // TODO: remove this
    fn null() -> AffineWavefrontsReduction {
        Self {
            reduction_strategy: WavefrontReductionType::WavefrontsReductionNone,
            min_wavefront_length: 0,
            max_distance_threshold: 0,
        }
    }
}

pub struct AffineWavefronts {
    // Dimensions
    pub pattern_length: i32, // Pattern length
    pub text_length: i32,    // Text length
    pub num_wavefronts: i32, // Total number of allocatable wavefronts

    // Limits
    pub max_penalty: i32, // MAX(mismatch_penalty,single_gap_penalty)
    pub max_k: i32, // Maximum diagonal k (used for null-wf, display, and banding) TODO: make these option?
    pub min_k: i32, // Maximum diagonal k (used for null-wf, display, and banding) TODO: make these option?

    // Wavefronts
    // TODO: should they be 2D
    pub mwavefronts: Vec<AffineWavefront>, // M-wavefronts
    pub iwavefronts: Vec<AffineWavefront>, // I-wavefronts
    pub dwavefronts: Vec<AffineWavefront>, // D-wavefronts
    // TODO: make this an option type?
    pub wavefront_null: Option<AffineWavefront>, // Null wavefront (used to gain orthogonality)

    // Reduction
    pub reduction: AffineWavefrontsReduction, // Reduction parameters

    // Penalties
    pub penalties: penalties::AffineWavefrontsPenalties, // Penalties parameters

    // CIGAR
    pub edit_cigar: edit_cigar::EditCigar, // Alignment CIGAR

    /*
    // MM
    mm_allocator_t* mm_allocator;              // MM-Allocator
    wavefronts_mem: AffineWavefront*,          // MM-Slab for AffineWavefront (base)
    wavefronts_current: AffineWavefront*,      // MM-Slab for AffineWavefront (next)
     */
    // STATS
    pub wavefronts_stats: WavefrontsStats, // Stats

                                           /*
                                           // DEBUG
                                           TODO: make this an option type
                                           #ifdef AFFINE_LAMBDA_WAVEFRONT_DEBUG
                                             affine_table_t gap_affine_table;             // DP-Table encoded by the wavefronts
                                           #endif
                                           */
}

// affine_wavefronts_new
impl AffineWavefronts {
    pub fn new(
        pattern_length: i32,
        text_length: i32,
        penalties: &mut penalties::AffinePenalties,
        penalties_strategy: penalties::WavefrontsPenaltiesStrategy,
        // mm_allocator_t* const mm_allocator
    ) -> Self {
        // Allocate
        // affine_wavefronts_t* const affine_wavefronts = mm_allocator_alloc(mm_allocator,affine_wavefronts_t);

        // Dimensions
        let max_score_misms: i32 = cmp::min(pattern_length, text_length) * penalties.mismatch; // TODO: remove std
        let max_score_indel: i32 = penalties.gap_opening
            + i32::abs(pattern_length - text_length) * penalties.gap_extension;
        let num_wavefronts: i32 = max_score_misms + max_score_indel;

        // MM
        // affine_wavefronts.mm_allocator = mm_allocator;
        // Limits
        let single_gap_penalty: i32 = penalties.gap_opening + penalties.gap_extension;
        let max_penalty: i32 = cmp::max(penalties.mismatch, single_gap_penalty);

        // Penalties
        let affine_wavefronts_penalties =
            penalties::AffineWavefrontsPenalties::init(penalties, penalties_strategy);

        // Allocate wavefronts
        // affine_wavefronts_allocate_wavefront_components(affine_wavefronts);
        // affine_wavefronts_allocate_wavefront_null(affine_wavefronts);

        // CIGAR
        let edit_cigar = edit_cigar::edit_cigar_allocate(pattern_length, text_length);

        // STATS
        let wavefronts_stats = ();

        AffineWavefronts {
            pattern_length,
            text_length,
            num_wavefronts,

            max_penalty,
            max_k: 0,
            min_k: 0,

            penalties: affine_wavefronts_penalties,
            reduction: AffineWavefrontsReduction::null(),

            mwavefronts: Vec::new(),
            iwavefronts: Vec::new(),
            dwavefronts: Vec::new(),
            wavefront_null: Option::None,

            edit_cigar,

            wavefronts_stats,
        }
    }

    pub fn affine_wavefronts_new_complete(
        pattern_length: i32,
        text_length: i32,
        penalties: &mut penalties::AffinePenalties,
        wavefronts_stats: &WavefrontsStats,
        // mm_allocator_t* const mm_allocator
    ) -> AffineWavefronts {
        // Create new
        let mut affine_wavefronts = AffineWavefronts::new(
            pattern_length,
            text_length,
            penalties,
            penalties::WavefrontsPenaltiesStrategy::ForceZeroMatch,
            // mm_allocator
        );

        // Limits
        affine_wavefronts.max_k = text_length;
        affine_wavefronts.min_k = -pattern_length;
        // Reduction
        affine_wavefronts.reduction.reduction_strategy =
            WavefrontReductionType::WavefrontsReductionNone;
        // affine_wavefronts_reduction_set_none(&affine_wavefronts.reduction);
        // Stats
        affine_wavefronts.wavefronts_stats = *wavefronts_stats;

        // Return
        affine_wavefronts
    }
}


/*
 * Allocate individual wavefront
 */
pub fn affine_wavefronts_allocate_wavefront(
    affine_wavefronts: &mut AffineWavefronts,
    lo_base: i32,
    hi_base: i32) -> AffineWavefront {
    // TODO: implement?
    // Compute limits
    let wavefront_length: i32 = hi_base - lo_base + 2; // (+1) for k=0
    let wavefront = AffineWavefront::new(false, lo_base, hi_base, lo_base, hi_base);
    return wavefront;
}
