use super::affine_wavefront;
use super::affine_wavefront_utils as utils;
use super::affine_wavefront_extend as extend;

/*
 * Computation using Wavefronts
 */
pub fn affine_wavefronts_align<T>(
    affine_wavefronts: &mut affine_wavefront::AffineWavefronts,
    match_lambda: T,
    traceback_lambda: T,
    pattern_length: i32,
    text_length: i32)
where T: Fn(i32, i32) -> bool {

    // Initialize wavefront
    utils::affine_wavefront_initialize(&mut affine_wavefronts);

  // Compute wavefronts for increasing score
    let score: i32 = 0;
    loop {
        // Exact extend s-wavefront
        extend::affine_wavefronts_extend_wavefront(
            affine_wavefronts,
            match_lambda,
            pattern_length,
            text_length,score);
        // Exit condition
        if affine_wavefront_end_reached(affine_wavefronts,pattern_length,text_length,score) {
            // Backtrace & check alignment reached
            affine_wavefronts_backtrace(
                affine_wavefronts,
                traceback_lambda,
                pattern_length,
                text_length,score);
            break;
        }

        // Update all wavefronts
        score += 1; // Increase score
        affine_wavefronts_compute_wavefront(affine_wavefronts,
                                            match_lambda,
                                            pattern_length,
                                            text_length,
                                            score);

        // DEBUG
        // affine_wavefronts_debug_step(affine_wavefronts,pattern,text,score);
        // WAVEFRONT_STATS_COUNTER_ADD(affine_wavefronts,wf_steps,1);
    }

    // DEBUG
    //affine_wavefronts_debug_step(affine_wavefronts,pattern,text,score);
    //WAVEFRONT_STATS_COUNTER_ADD(affine_wavefronts,wf_score,score); // STATS
}

