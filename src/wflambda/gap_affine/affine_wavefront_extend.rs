use super::affine_wavefront;
use std::cmp;
use super::affine_wavefront_utils as utils;

/*
 * Reduce wavefront
 */
fn affine_wavefronts_reduce_wavefront_offsets(
    affine_wavefronts: &mut affine_wavefront::AffineWavefronts,
    wavefront: &mut affine_wavefront::AffineWavefront,
    pattern_length: i32,
    text_length: i32,
    min_distance: i32,
    max_distance_threshold: i32,
    alignment_k: i32) {
    // Parameters
    let offsets: Vec<affine_wavefront::AwfOffset> = wavefront.offsets;

    // Reduce from bottom

    let top_limit: i32 = cmp::min(alignment_k-1,wavefront.hi);
    // let k: i32;
    for k in wavefront.lo..=top_limit {
        let distance: i32 = utils::affine_wavefronts_compute_distance(pattern_length,text_length,offsets[k as usize],k);
        if distance - min_distance  <= max_distance_threshold {break};
        wavefront.lo += 1;

    };

    // Reduce from top
    let bottom_limit: i32 = cmp::max(alignment_k+1,wavefront.lo);
    for k in (bottom_limit+1..=wavefront.hi).rev()  {
        let distance: i32 = utils::affine_wavefronts_compute_distance(pattern_length,text_length,offsets[k as usize],k);
        if distance - min_distance <= max_distance_threshold { break };
        wavefront.hi -= 1;
  }
  // Check hi/lo range
  if wavefront.lo > wavefront.hi {
    wavefront.null = true
  }
  // STATS
  // WAVEFRONT_STATS_COUNTER_ADD(affine_wavefronts,wf_reduced_cells, (wavefront->hi_base-wavefront->hi)+(wavefront->lo-wavefront->lo_base));
}

fn affine_wavefronts_extend_mwavefront_compute<T>(
    affine_wavefronts: &mut affine_wavefront::AffineWavefronts,
    lambda: T,
    pattern_length: i32,
    text_length: i32,
    score: i32)
where  T: Fn(i32, i32) -> bool {

    // Fetch m-wavefront
    //if let Some(mwavefront) = affine_wavefronts.mwavefronts.get(score as usize) {}
    let mwavefront: &affine_wavefront::AffineWavefront;
    match affine_wavefronts.mwavefronts.get(score as usize) {
        None => return,
        Some(m) => {mwavefront=m},
    };

    // Extend diagonally each wavefront point
    let offsets: &Vec<affine_wavefront::AwfOffset> = &mwavefront.offsets;

    for k in mwavefront.lo as usize ..= mwavefront.hi as usize {
        let offset: affine_wavefront::AwfOffset = offsets[k];
    }

    for k in mwavefront.lo as usize ..= mwavefront.hi as usize {
        // Exact extend
        let offset: affine_wavefront::AwfOffset = offsets[k];
        let v: i32 = affine_wavefront::affine_lambda_wavefront_v(k as i32, offset as i32);
        let h: i32 =  affine_wavefront::affine_lambda_wavefront_v(k as i32,offset as i32);
        while lambda(v , h) {
            offsets[k] += 1;
            v+=1;
            h+=1;
        }
  }
}


fn affine_wavefronts_reduce_wavefronts(
    affine_wavefronts: &mut affine_wavefront::AffineWavefronts,
    pattern_length: i32,
    text_length: i32,
    score: i32
) {

    // Parameters
    let min_wavefront_length: i32 = affine_wavefronts.reduction.min_wavefront_length;
    let max_distance_threshold: i32 = affine_wavefronts.reduction.max_distance_threshold;
    let alignment_k: i32 = affine_wavefront::affine_lambda_wavefront_diagonal(text_length,pattern_length);

    // Fetch m-wavefront
    let mwavefront: &affine_wavefront::AffineWavefront;
    match affine_wavefronts.mwavefronts.get(score as usize) {
        None => return,
        Some(m) => {mwavefront=m},
    };


    if (mwavefront.hi - mwavefront.lo + 1) < min_wavefront_length { return }
    // WAVEFRONT_STATS_COUNTER_ADD(affine_wavefronts,wf_reduction,1); // STATS

    // Compute min-distance
    let offsets: Vec<affine_wavefront::AwfOffset> = mwavefront.offsets;
    let min_distance: i32 = cmp::max(pattern_length,text_length);
    let k: i32;

    for k in mwavefront.lo ..= mwavefront.hi {
        let distance: i32 = utils::affine_wavefronts_compute_distance(pattern_length,text_length,offsets[k as usize],k);
        min_distance = cmp::min(min_distance,distance);
    }

    // Reduce m-wavefront
    affine_wavefronts_reduce_wavefront_offsets(
        affine_wavefronts,&mut mwavefront,pattern_length,text_length,
        min_distance,max_distance_threshold,alignment_k);

    // Reduce i-wavefront
    if let Some(iwavefront) = affine_wavefronts.iwavefronts.get(score as usize) {
        if mwavefront.lo > iwavefront.lo {iwavefront.lo = mwavefront.lo;}
        if mwavefront.hi < iwavefront.hi {iwavefront.hi = mwavefront.hi;}
        if iwavefront.lo > iwavefront.hi {iwavefront.null = true;}
    };

    // Reduce d-wavefront
    if let Some(dwavefront) = affine_wavefronts.dwavefronts.get(score as usize) {
        if mwavefront.lo > dwavefront.lo {dwavefront.lo = mwavefront.lo};
        if mwavefront.hi < dwavefront.hi {dwavefront.hi = mwavefront.hi};
        if dwavefront.lo > dwavefront.hi {dwavefront.null = true};
  }
}

/*
 * Gap-Affine Wavefront exact extension
 */
pub fn affine_wavefronts_extend_wavefront<T>(
    affine_wavefronts: &mut affine_wavefront::AffineWavefronts,
    lambda: T,
    pattern_length: i32,
    text_length: i32,
    score: i32)
where  T: Fn(i32, i32) -> bool {
    // Extend wavefront
    affine_wavefronts_extend_mwavefront_compute(
        affine_wavefronts,
        lambda,
        pattern_length,
        text_length,score);

    // Reduce wavefront dynamically
    if affine_wavefronts.reduction.reduction_strategy == affine_wavefront::WavefrontReductionType::WavefrontsReductionDynamic {
        affine_wavefronts_reduce_wavefronts(affine_wavefronts, pattern_length, text_length, score);
    }
}
