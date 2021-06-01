#![allow(non_camel_case_types, non_snake_case, dead_code)]

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
type awf_offset_t = i32;

/*
 * Counters
 *   (from http://www.johndcook.com/standard_deviation.html)
 */
struct profiler_counter_t {
    total: u64,
    samples: u64,
    min: u64,
    max: u64,
    m_oldM: f64,
    m_newM: f64,
    m_oldS: f64,
    m_newS: f64,
}

struct profiler_timer_t {
    /* Timer */
    // struct timespec begin_timer;     // Timer begin // TODO: how to translate

    /* Total time & samples taken */
    time_ns: profiler_counter_t,
    accumulated: u64,
}

/*
 * Wavefront Stats
 */
struct wavefronts_stats_t {
    wf_score: profiler_counter_t,         // Score reached by WF-alignment
    wf_steps: profiler_counter_t,         // Step performed by WF-alignment
    wf_steps_null: profiler_counter_t,    // Avoided WF-alignment steps due to null wavefronts
    wf_steps_extra: profiler_counter_t, // Extra steps performed by WF-alignment searching for a better solution
    wf_operations: profiler_counter_t,  // Single cell WF-operations performed
    wf_extensions: profiler_counter_t,  // Single cell WF-extensions performed
    wf_reduction: profiler_counter_t,   // Calls to reduce wavefront
    wf_reduced_cells: profiler_counter_t, // Total cells reduced
    wf_null_used: profiler_counter_t,   // Total times a null-WF was used as padding
    wf_extend_inner_loop: profiler_counter_t, // Total times SIMD-extension had to re-iterate
    //wf_compute_kernel[4]: profiler_counter_t,  // Specialized WF computation kernel used // array of size 4
    wf_time_backtrace: profiler_timer_t, // Time spent doing backtrace
    wf_backtrace_paths: profiler_counter_t, // Total paths explored by the backtrace
    wf_backtrace_alg: profiler_counter_t, // Total alignments explored by the backtrace
}

struct affine_wavefront {
    // Range
    null: bool,   // Is null interval?
    lo: i32,      // Effective lowest diagonal (inclusive)
    hi: i32,      // Effective highest diagonal (inclusive)
    lo_base: i32, // Lowest diagonal before reduction (inclusive)
    hi_base: i32, // Highest diagonal before reduction (inclusive)

    // Offsets
    offsets: awf_offset_t, // Offsets

                           // TODO: make this an option type
                           // #ifdef AFFINE_LAMBDA_WAVEFRONT_DEBUG
                           // offsets_base: awf_offset_t, // Offsets increment
                           // #endif
}

enum wavefront_reduction_type {
    wavefronts_reduction_none,
    wavefronts_reduction_dynamic,
}

struct edit_cigar_t {
    operations: String,
    max_operations: i32,
    begin_offset: i32,
    end_offset: i32,
    score: i32,
}

struct affine_penalties_t {
    a_match: i32,       // (Penalty representation; usually M <= 0)
    mismatch: i32,      // (Penalty representation; usually X > 0)
    gap_opening: i32,   // (Penalty representation; usually O > 0)
    gap_extension: i32, // (Penalty representation; usually E > 0)
}

/*
 * Wavefront Strategy
 */
enum wavefronts_penalties_strategy {
    wavefronts_penalties_match_zero,
    wavefronts_penalties_force_zero_match,
    wavefronts_penalties_shifted_penalties,
    wavefronts_penalties_odd_pair_penalties,
}

/*
 * Wavefront Penalties
 */
struct affine_wavefronts_penalties_t {
    base_penalties: affine_penalties_t, // Input base Gap-Affine penalties
    wavefront_penalties: affine_penalties_t, // Wavefront Gap-Affine penalties
    penalties_strategy: wavefronts_penalties_strategy, // Penalties adaptation strategy
}

/*
 * Wavefront Penalties
 */
struct affine_wavefronts_reduction_t {
    reduction_strategy: wavefront_reduction_type, // Reduction strategy
    min_wavefront_length: i32,                    // Dynamic: Minimum wavefronts length to reduce
    max_distance_threshold: i32, // Dynamic: Maximum distance between offsets allowed
}

struct affine_wavefronts {
    // Dimensions
    pattern_length: i32, // Pattern length
    text_length: i32,    // Text length
    num_wavefronts: i32, // Total number of allocatable wavefronts

    // Limits
    max_penalty: i32, // MAX(mismatch_penalty,single_gap_penalty)
    max_k: i32,       // Maximum diagonal k (used for null-wf, display, and banding)
    min_k: i32,       // Maximum diagonal k (used for null-wf, display, and banding)

    // Wavefronts
    // mwavefronts: affine_wavefront**,            // M-wavefronts
    // iwavefronts: affine_wavefront**,            // I-wavefronts
    // dwavefronts: affine_wavefront**,            // D-wavefronts
    wavefront_null: affine_wavefront, // Null wavefront (used to gain orthogonality)

    // Reduction
    reduction: affine_wavefronts_reduction_t, // Reduction parameters

    // Penalties
    penalties: affine_wavefronts_penalties_t, // Penalties parameters

    // CIGAR
    edit_cigar: edit_cigar_t, // Alignment CIGAR

                              // MM
                              //mm_allocator_t* mm_allocator;                // MM-Allocator
                              //wavefronts_mem: affine_wavefront_t*,          // MM-Slab for affine_wavefront_t (base)
                              //wavefronts_current: affine_wavefront_t*,      // MM-Slab for affine_wavefront_t (next)

                              // STATS
                              //wavefronts_stats: wavefronts_stats_t*, // Stats

                              // DEBUG
                              // TODO: make this an option type
                              // #ifdef AFFINE_LAMBDA_WAVEFRONT_DEBUG
                              // affine_table_t gap_affine_table;             // DP-Table encoded by the wavefronts
                              // #endif
}
