#![allow(dead_code)]

/*
 * Offset size

//#define AFFINE_LAMBDA_WAVEFRONT_W8
//#define AFFINE_LAMBDA_WAVEFRONT_W16
#define AFFINE_LAMBDA_WAVEFRONT_W32

#ifdef AFFINE_LAMBDA_WAVEFRONT_W8
typedef int8_t awf_offset_t;
#else
#ifdef AFFINE_LAMBDA_WAVEFRONT_W16
typedef int16_t awf_offset_t;
#else // AFFINE_LAMBDA_WAVEFRONT_W32
typedef int32_t awf_offset_t;
#endif
#endif
 */

pub mod wavefront {

    pub type AwfOffset = usize;

    struct AffinePenalties {
        a_match: i32,     // (Penalty representation; usually M <= 0) // match is a rust_keyword
        mismatch: i32,    // (Penalty representation; usually X > 0)
        gap_opening: i32, // (Penalty representation; usually O > 0)
        gap_extension: i32, // (Penalty representation; usually E > 0)
    }

    /*
     * Wavefront Strategy
     */
    enum WavefrontsPenaltiesStrategy {
        WavefrontsPenaltiesMatchZero,
        WavefrontsPenaltiesForceZeroMatch,
        WavefrontsPenaltiesShiftedPenalties,
        WavefrontsPenaltiesOddPairPenalties,
    }

    /*
     * Wavefront Penalties
     */
    struct AffineWavefrontsPenalties {
        base_penalties: AffinePenalties, // Input base Gap-Affine penalties
        wavefront_penalties: AffinePenalties, // Wavefront Gap-Affine penalties
        penalties_strategy: WavefrontsPenaltiesStrategy, // Penalties adaptation strategy
    }

    
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

    enum WavefrontReductionType {
        WavefrontsReductionNone,
        WavefrontsReductionDynamic,
    }

    /*
     * Wavefront Penalties
     */
    struct AffineWavefrontsReduction {
        reduction_strategy: WavefrontReductionType, // Reduction strategy
        min_wavefront_length: i32,                  // Dynamic: Minimum wavefronts length to reduce
        max_distance_threshold: i32, // Dynamic: Maximum distance between offsets allowed
    }

    pub struct AffineWavefronts {
        // Dimensions
        pattern_length: i32, // Pattern length
        text_length: i32,    // Text length
        num_wavefronts: i32, // Total number of allocatable wavefronts

        // Limits
        max_penalty: i32, // MAX(mismatch_penalty,single_gap_penalty)
        max_k: i32,       // Maximum diagonal k (used for null-wf, display, and banding)
        min_k: i32,       // Maximum diagonal k (used for null-wf, display, and banding)

        // Wavefronts
        // reference to AffineWavefront
        // pub mwavefronts: Vec<AffineWavefront>, // M-wavefronts
        // iwavefronts: Vec<AffineWavefront>,     // I-wavefronts
        // dwavefronts: Vec<AffineWavefront>,     // D-wavefronts
        // make this an option type?
        wavefront_null: AffineWavefront, // Null wavefront (used to gain orthogonality)

        // Reduction
        reduction: AffineWavefrontsReduction, // Reduction parameters

        // Penalties
        penalties: AffineWavefrontsPenalties, // Penalties parameters

        // CIGAR
        edit_cigar: EditCigar, // Alignment CIGAR

        // MM-like
        // mm_allocator_t* mm_allocator;              // MM-Allocator
        // wavefronts_mem: Vec< &'a Vec<AffineWavefront>>,         // MM-Slab for AffineWavefront (base)
        // pub wavefronts_current: Vec< &'a Vec<AffineWavefront>>, // MM-Slab for AffineWavefront (next)

        /*
        // MM
        mm_allocator_t* mm_allocator;              // MM-Allocator
        wavefronts_mem: AffineWavefront*,          // MM-Slab for AffineWavefront (base)
        wavefronts_current: AffineWavefront*,      // MM-Slab for AffineWavefront (next)

        // STATS
        wavefronts_stats: wavefronts_stats_t*,     // Stats

        // DEBUG
        TODO: make this an option type
        #ifdef AFFINE_LAMBDA_WAVEFRONT_DEBUG
          affine_table_t gap_affine_table;             // DP-Table encoded by the wavefronts
        #endif
        */
    }
}
