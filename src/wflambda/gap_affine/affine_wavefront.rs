use std::cmp;
use super::affine_wavefront_penalties as penalties;

pub struct AffineWavefront {
    // Range
    pub null: bool,   // Is null interval?
    pub lo: usize,  // Effective lowest diagonal (inclusive)
    pub hi: usize,  // Effective highest diagonal (inclusive)
    pub lo_base: usize, // Lowest diagonal before reduction (inclusive)
    pub hi_base: usize, // Highest diagonal before reduction (inclusive)

    // Offsets
    pub offsets: Vec<AwfOffset>, // Offsets

        /* TODO: make this an option type
// #ifdef AFFINE_LAMBDA_WAVEFRONT_DEBUG
// offsets_base: AwfOffset, // Offsets increment
// #endif
 */
}

struct EditCigar {
    operations: String,
    max_operations: i32,
    begin_offset: i32,
    end_offset: i32,
    score: i32,
}

struct WavefrontsStats();

enum WavefrontReductionType {
    WavefrontsReductionNone,
    WavefrontsReductionDynamic,
}

struct AffineWavefrontsReduction {
    reduction_strategy: WavefrontReductionType, // Reduction strategy
    min_wavefront_length: i32,                  // Dynamic: Minimum wavefronts length to reduce
    max_distance_threshold: i32, // Dynamic: Maximum distance between offsets allowed
}

pub struct AffineWavefronts {
    // Dimensions
    pub pattern_length: i32, // Pattern length
    pub text_length: i32,    // Text length
    pub num_wavefronts: i32, // Total number of allocatable wavefronts

    // Limits
    pub max_penalty: i32, // MAX(mismatch_penalty,single_gap_penalty)
    pub max_k: i32,       // Maximum diagonal k (used for null-wf, display, and banding)
    pub min_k: i32,       // Maximum diagonal k (used for null-wf, display, and banding)

    // Wavefronts
    pub mwavefronts: Vec<Vec<AffineWavefront>>, // M-wavefronts
    pub iwavefronts: Vec<Vec<AffineWavefront>>,     // I-wavefronts
    pub dwavefronts: Vec<Vec<AffineWavefront>>,     // D-wavefronts
    // TODO: make this an option type?
    pub wavefront_null: AffineWavefront, // Null wavefront (used to gain orthogonality)

    // Reduction
    pub reduction: AffineWavefrontsReduction, // Reduction parameters

    // Penalties
    pub penalties: AffineWavefrontsPenalties, // Penalties parameters

    // CIGAR
    pub edit_cigar: EditCigar, // Alignment CIGAR

    /*
    // MM
    mm_allocator_t* mm_allocator;              // MM-Allocator
    wavefronts_mem: AffineWavefront*,          // MM-Slab for AffineWavefront (base)
    wavefronts_current: AffineWavefront*,      // MM-Slab for AffineWavefront (next)
     */
    // STATS
    pub wavefronts_stats: WavefrontsStats,     // Stats

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
    pub fn new (
        pattern_length: i32,
        text_length: i32,
        penalties: &mut AffinePenalties,
        penalties_strategy: WavefrontsPenaltiesStrategy,
        // mm_allocator_t* const mm_allocator
    ) -> Self {
        // Allocate
        // affine_wavefronts_t* const affine_wavefronts = mm_allocator_alloc(mm_allocator,affine_wavefronts_t);

        // Dimensions
        let max_score_misms: i32 = cmp::min(pattern_length,text_length) * penalties.mismatch; // TODO: remove std
        let max_score_indel: i32 = penalties.gap_opening + i32::abs(pattern_length-text_length) * penalties.gap_extension;
        let num_wavefronts: i32 = max_score_misms + max_score_indel;

        // MM
        // affine_wavefronts.mm_allocator = mm_allocator;
        // Limits
        let single_gap_penalty: i32 = penalties.gap_opening + penalties.gap_extension;
        let max_penalty: i32 = cmp::max(penalties.mismatch,single_gap_penalty);

        // Penalties
        penalties.init(penalties, penalties_strategy);

        // Allocate wavefronts
        affine_wavefronts_allocate_wavefront_components(affine_wavefronts);
        affine_wavefronts_allocate_wavefront_null(affine_wavefronts);

        // CIGAR
        edit_cigar_allocate(&affine_wavefronts.edit_cigar,pattern_length,text_length,mm_allocator);

        // STATS
        let wavefronts_stats = ();

        AffineWavefronts {
            pattern_length,
            text_length,
            num_wavefronts,

            max_penalty,

            wavefronts_stats,
        }
    }

    pub fn affine_wavefronts_new_complete(
        pattern_length: i32,
        text_length: i32,
        penalties: &wavefront_types::AffinePenalties,
        wavefronts_stats: &wavefront_types::WavefrontStats,
        // mm_allocator_t* const mm_allocator
    ) -> wavefront_types::AffineWavefronts {
        // Create new
        let affine_wavefronts = new(
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
        affine_wavefronts_reduction_set_none(&affine_wavefronts.reduction);
        // Stats
        affine_wavefronts.wavefronts_stats = wavefronts_stats;

        // Return
        affine_wavefronts
    }

}

