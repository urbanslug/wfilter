use super::types::wavefront::AffineWavefronts;

pub fn  affine_wavefronts_extend_mwavefront_compute <T>(
    affine_wavefronts: &AffineWavefronts,
    lambda: T,
    pattern_length: i32,
    text_length: i32,
    score: i32)  where T: Fn(i32, i32) -> bool  {

  // Fetch m-wavefront
    let mwavefront: &affine_wavefront_t = affine_wavefronts.mwavefronts[score];
  if (mwavefront==NULL) return;
  // Extend diagonally each wavefront point
  awf_offset_t* const offsets = mwavefront->offsets;
  int k;
  for (k=mwavefront->lo;k<=mwavefront->hi;++k) {
    // Exact extend
    const awf_offset_t offset = offsets[k];
    int v = AFFINE_LAMBDA_WAVEFRONT_V(k,offset);
    int h = AFFINE_LAMBDA_WAVEFRONT_H(k,offset);
    while (lambda(v++,h++)) {
      ++(offsets[k]);
/*
//#define AFFINE_LAMBDA_WAVEFRONT_SHOW
#ifdef AFFINE_LAMBDA_WAVEFRONT_SHOW
      std::cerr << v-1 << "\t" << h-1 << "\t" << score << "\t" << "aligned" << std::endl;
#endif
    }
#ifdef AFFINE_LAMBDA_WAVEFRONT_SHOW
    std::cerr << v << "\t" << h << "\t" << score << "\t" << "unaligned" << std::endl;
#endif
*/
    }

  // DEBUG
  affine_wavefronts_extend_mwavefront_epiloge(affine_wavefronts,score,pattern_length,text_length);
  }
}
/*
 * Gap-Affine Wavefront exact extension
 */
pub fn affine_wavefronts_extend_wavefront <T>(
    affine_wavefronts: &AffineWavefronts,
    lambda: T,
    pattern_length: i32,
    text_length: i32,
    score: i32) where T: Fn(i32, i32) -> bool {
  // Extend wavefront
  affine_wavefronts_extend_mwavefront_compute(
      affine_wavefronts,
      lambda,
      pattern_length,
      text_length,
      score);
    // Reduce wavefront dynamically
    /* TODO: enable
  if (affine_wavefronts->reduction.reduction_strategy == wavefronts_reduction_dynamic) {
    affine_wavefronts_reduce_wavefronts(
        affine_wavefronts,pattern_length,
        text_length,score);
  }
    */
}

