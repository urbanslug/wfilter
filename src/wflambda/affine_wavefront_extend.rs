use super::types::wavefront::AffineWavefronts;

#[allow(dead_code, unused_variables)]
pub fn affine_wavefronts_extend_mwavefront_compute<T>(
    affine_wavefronts: &AffineWavefronts,
    lambda: T,
    pattern_length: i32,
    text_length: i32,
    score: i32,
) where
    T: Fn(i32, i32) -> bool,
{
}
/*
 * Gap-Affine Wavefront exact extension
 */
pub fn affine_wavefronts_extend_wavefront<T>(
    affine_wavefronts: &AffineWavefronts,
    lambda: T,
    pattern_length: i32,
    text_length: i32,
    score: i32,
) where
    T: Fn(i32, i32) -> bool,
{
    // Extend wavefront
    affine_wavefronts_extend_mwavefront_compute(
        affine_wavefronts,
        lambda,
        pattern_length,
        text_length,
        score,
    );
    // Reduce wavefront dynamically
    /* TODO: enable
    if (affine_wavefronts->reduction.reduction_strategy == wavefronts_reduction_dynamic) {
      affine_wavefronts_reduce_wavefronts(
          affine_wavefronts,pattern_length,
          text_length,score);
    }
      */
}
