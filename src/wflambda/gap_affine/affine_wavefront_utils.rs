use super::affine_wavefront;
use std::cmp;

/*
 * Initial Conditions and finalization
 */
pub fn affine_wavefront_initialize(affine_wavefronts: &mut affine_wavefront::AffineWavefronts) {
    let lo_base = 0;
    let hi_base = 0;
    let mut wavefront = affine_wavefront::AffineWavefront::new(false, lo_base, hi_base, lo_base, hi_base);
    wavefront.offsets.push(0);
    //let wavefronts: Vec<affine_wavefront::AffineWavefront> = vec![wavefront];
    // TODO: make more idiomatic rust
    affine_wavefronts.mwavefronts.push(wavefront);
}


pub fn affine_wavefronts_compute_distance(
    pattern_length: i32,
    text_length: i32,
    offset: affine_wavefront::AwfOffset,
    k: i32) -> i32 {
    let v: i32 = affine_wavefront::affine_lambda_wavefront_v(k, offset as i32);
    let h: i32 = affine_wavefront::affine_lambda_wavefront_h(k, offset as i32);
    let left_v: i32 = pattern_length - v;
    let left_h: i32 = text_length - h;

    cmp::max(left_v,left_h)
}
