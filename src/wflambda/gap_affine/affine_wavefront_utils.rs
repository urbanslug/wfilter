use super::affine_wavefront;
use std::cmp;


/*
 * Initial Conditions and finalization
 */
pub fn affine_wavefront_initialize(
    affine_wavefronts: &mut affine_wavefront::AffineWavefronts) {
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

pub fn affine_wavefront_end_reached(
    affine_wavefronts: &mut affine_wavefront::AffineWavefronts,
    pattern_length: i32,
    text_length: i32,
    score: i32) -> bool {
    // Parameters
    let alignment_k: usize = affine_wavefront::affine_lambda_wavefront_diagonal(text_length,pattern_length) as usize;
    let alignment_offset: usize = affine_wavefront::affine_lambda_wavefront_offset(text_length,pattern_length) as usize;
    // Fetch wavefront and check termination

    // let mwavefront: &affine_wavefront::AffineWavefront;
    match affine_wavefronts.mwavefronts.get(score as usize) {
        None => false,
        Some(mwavefront) => {
            let offsets: &Vec<affine_wavefront::AwfOffset> = &mwavefront.offsets;
            if (mwavefront.lo as usize) <= alignment_k &&
                alignment_k <= (mwavefront.hi as usize) &&
                offsets[alignment_k] >= alignment_offset {
                    true
                } else {
                    false
                }
        }
    }

}


/*
 * Accessors
 */
pub fn affine_wavefronts_get_source_mwavefront(
    affine_wavefronts: &affine_wavefront::AffineWavefronts,
    score: i32) -> Option<&affine_wavefront::AffineWavefront> {

    if score < 0 {
        return None
    }

    match affine_wavefronts.mwavefronts.get(score as usize) {
        Some(mwavefront)  => { Some(&mwavefront) }
        _ => { None }
    }
}


pub fn affine_wavefronts_get_source_iwavefront(
    affine_wavefronts: &affine_wavefront::AffineWavefronts,
    score: i32) -> Option<&affine_wavefront::AffineWavefront> {

    if score < 0 {
        return None
    }

    match affine_wavefronts.iwavefronts.get(score as usize) {
        Some(iwavefront)  => { Some(&iwavefront) }
        _ => { None }
    }
}

pub fn affine_wavefronts_get_source_dwavefront(
    affine_wavefronts: &affine_wavefront::AffineWavefronts,
    score: i32) -> Option<&affine_wavefront::AffineWavefront> {

    if score < 0 {
        return None
    }

    match affine_wavefronts.dwavefronts.get(score as usize) {
        Some(dwavefront)  => { Some(&dwavefront) }
        _ => { None }
    }
}
