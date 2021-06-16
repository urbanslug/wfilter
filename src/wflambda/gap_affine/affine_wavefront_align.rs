use super::affine_wavefront;
use super::affine_wavefront_penalties as penalties;
use super::affine_wavefront_utils as utils;
use super::affine_wavefront_extend as extend;
use super::affine_wavefront_backtrace as backtrace;


/*
 * Fetch & allocate wavefronts
 */
fn affine_wavefronts_fetch_wavefronts(
    affine_wavefronts: &mut affine_wavefront::AffineWavefronts,
    wavefront_set: &mut affine_wavefront::AffineWavefrontSet,
    score: i32) {

    // Compute scores
    let wavefront_penalties: &penalties::AffinePenalties = &affine_wavefronts.penalties.wavefront_penalties;
    let mismatch_score: i32 = score - wavefront_penalties.mismatch;
    let gap_open_score: i32 = score - wavefront_penalties.gap_opening - wavefront_penalties.gap_extension;
    let gap_extend_score: i32 = score - wavefront_penalties.gap_extension;

    // Fetch wavefronts
    if let Some(m) =  utils::affine_wavefronts_get_source_mwavefront(affine_wavefronts,mismatch_score) {
        wavefront_set.in_mwavefront_sub = m;
    };

    if let Some(m) = utils::affine_wavefronts_get_source_mwavefront(affine_wavefronts,gap_open_score) {
        wavefront_set.in_mwavefront_gap = m;
    };
    if let Some(i) = utils::affine_wavefronts_get_source_iwavefront(affine_wavefronts,gap_extend_score) {
        wavefront_set.in_iwavefront_ext = i;
    };
    if let Some(d) = utils::affine_wavefronts_get_source_dwavefront(affine_wavefronts,gap_extend_score) {
        wavefront_set.in_dwavefront_ext = d;
    }
}


fn affine_wavefronts_compute_limits(
    affine_wavefronts: &mut affine_wavefront::AffineWavefronts,
    wavefront_set: &affine_wavefront::AffineWavefrontSet,
    score: i32,
    lo_effective: &i32,
    hi_effective: &i32) {
  // Set limits (min_lo)
  let lo: i32 = wavefront_set.in_mwavefront_sub.lo;
  if lo > wavefront_set.in_mwavefront_gap.lo {lo = wavefront_set.in_mwavefront_gap.lo;}
  if lo > wavefront_set.in_iwavefront_ext.lo {lo = wavefront_set.in_iwavefront_ext.lo;}
  if lo > wavefront_set.in_dwavefront_ext.lo {lo = wavefront_set.in_dwavefront_ext.lo;}
  lo -= 1;
  // Set limits (max_hi)
  let hi: i32 = wavefront_set.in_mwavefront_sub.hi;
  if hi < wavefront_set.in_mwavefront_gap.hi {hi = wavefront_set.in_mwavefront_gap.hi;}
  if hi < wavefront_set.in_iwavefront_ext.hi {hi = wavefront_set.in_iwavefront_ext.hi;}
  if hi < wavefront_set.in_dwavefront_ext.hi {hi = wavefront_set.in_dwavefront_ext.hi;}
  hi += 1;
  // Set effective limits values
  *hi_effective = hi;
  *lo_effective = lo;
}


/*
 * Compute wavefront
 */
fn affine_wavefronts_compute_wavefront<T>(
    affine_wavefronts: &mut affine_wavefront::AffineWavefronts,
    lambda: T,
    pattern_length: i32,
    text_length: i32,
    score: i32)
where  T: Fn(i32, i32) -> bool {

    // Select wavefronts
    let wavefront_set: affine_wavefront::AffineWavefrontSet;
    affine_wavefronts_fetch_wavefronts(affine_wavefronts,&mut wavefront_set,score);

    /*
    We are not doing stats

    // Check null wavefronts
    if (wavefront_set.in_mwavefront_sub->null &&
        wavefront_set.in_mwavefront_gap->null &&
        wavefront_set.in_iwavefront_ext->null &&
        wavefront_set.in_dwavefront_ext->null) {
        WAVEFRONT_STATS_COUNTER_ADD(affine_wavefronts,wf_steps_null,1);
        return;
    }
    WAVEFRONT_STATS_COUNTER_ADD(affine_wavefronts,wf_null_used,(wavefront_set.in_mwavefront_sub->null?1:0));
    WAVEFRONT_STATS_COUNTER_ADD(affine_wavefronts,wf_null_used,(wavefront_set.in_mwavefront_gap->null?1:0));
    WAVEFRONT_STATS_COUNTER_ADD(affine_wavefronts,wf_null_used,(wavefront_set.in_iwavefront_ext->null?1:0));
    WAVEFRONT_STATS_COUNTER_ADD(affine_wavefronts,wf_null_used,(wavefront_set.in_dwavefront_ext->null?1:0));
     */

    // Set limits
    let hi: i32;
    let lo: i32;
    affine_wavefronts_compute_limits(affine_wavefronts,&wavefront_set,score,&lo,&hi);
    // Allocate score-wavefronts
    affine_wavefronts_allocate_wavefronts(affine_wavefronts,&wavefront_set,score,lo,hi);
    // Compute WF
    // TODO: fix this bit shift and null
    let kernel: i32 = ((wavefront_set.out_iwavefront!=NULL) << 1) | (wavefront_set.out_dwavefront!=NULL);
    //WAVEFRONT_STATS_COUNTER_ADD(affine_wavefronts,wf_compute_kernel[kernel],1);
    

    match kernel {
        3 => { affine_wavefronts_compute_offsets_idm(affine_wavefronts,&wavefront_set,lo,hi); },
        2 => { affine_wavefronts_compute_offsets_im(affine_wavefronts,&wavefront_set,lo,hi); },
        1 => { affine_wavefronts_compute_offsets_dm(affine_wavefronts,&wavefront_set,lo,hi); },
        0 => { affine_wavefronts_compute_offsets_m(affine_wavefronts,&wavefront_set,lo,hi); },
    };

    /*
    // Account for WF operations performed
    WAVEFRONT_STATS_COUNTER_ADD(affine_wavefronts,wf_operations,hi-lo+1);
    // DEBUG
    #ifdef AFFINE_LAMBDA_WAVEFRONT_DEBUG
    // Copy offsets base before extension (for display purposes)
    affine_wavefront_t* const mwavefront = affine_wavefronts->mwavefronts[score];
    if (mwavefront!=NULL) {
        int k;
        for (k=mwavefront->lo;k<=mwavefront->hi;++k) {
            mwavefront->offsets_base[k] = mwavefront->offsets[k];
        }
    }
    #endif
    */
}

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
        if utils::affine_wavefront_end_reached(affine_wavefronts,pattern_length,text_length,score) {
            // Backtrace & check alignment reached
            backtrace::affine_wavefronts_backtrace(
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

