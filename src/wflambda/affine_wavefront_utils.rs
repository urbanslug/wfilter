use super::types::wavefront as wavefront_types;
use super::types::wavefront::{AffineWavefront, AffineWavefronts};

// Easier to write these macros as fns for now
// TODO: make these macros?
/*
 * Translate k and offset to coordinates h,v
 */
// #define AFFINE_LAMBDA_WAVEFRONT_V(k,offset) ((offset)-(k))
// #define AFFINE_LAMBDA_WAVEFRONT_H(k,offset) (offset)

fn affine_lambda_wavefront_diagonal(h: usize, v: usize) -> usize {
    (h) - (v)
}
fn affine_lambda_wavefront_offset(h: usize, v: usize) -> usize {
    h
}

const zero: usize = 0;

/*
 * Allocate individual wavefront
 */
fn affine_wavefronts_allocate_wavefront<'a>(
    affine_wavefronts: &'a AffineWavefronts,
    lo_base: usize,
    hi_base: usize) -> &'a mut AffineWavefront {
    // Compute limits
    let wavefront_length: usize = hi_base - lo_base + 2; // (+1) for k=0

    // Allocate wavefront
    let mut wavefront: &mut AffineWavefront = affine_wavefronts.wavefronts_current[0];
    // ++(affine_wavefronts.wavefronts_current); // Next

    // Configure offsets
    wavefront.null = false;
    wavefront.lo = lo_base;
    wavefront.hi = hi_base;
    wavefront.lo_base = lo_base;
    wavefront.hi_base = hi_base;

    /*
    // Allocate offsets

    let offsets_mem: wavefront_types::AwfOffset = mm_allocator_calloc(
        affine_wavefronts.mm_allocator,wavefront_length,awf_offset_t,false);
    let offsets: Vec<wavefront_types::AwfOffset> = offsets_mem - lo_base; // Center at k=0
    wavefront.offsets = offsets;

    // DEBUG
#ifdef AFFINE_LAMBDA_WAVEFRONT_DEBUG
  awf_offset_t* const offsets_base_mem = mm_allocator_calloc(
      affine_wavefronts->mm_allocator,wavefront_length,awf_offset_t,false);
  wavefront->offsets_base = offsets_base_mem - lo_base; // Center at k=0
#endif
  
     */

    // Return
    wavefront


}


/*
 * Initial Conditions and finalization
 */
pub fn affine_wavefront_initialize(affine_wavefronts: &mut AffineWavefronts) {
    affine_wavefronts.mwavefronts[zero] =
        affine_wavefronts_allocate_wavefront(affine_wavefronts, 0, 0);
    affine_wavefronts.mwavefronts[0].offsets[0] = 0;
}

pub fn affine_wavefront_end_reached(
    affine_wavefronts: &AffineWavefronts,
    pattern_length: usize,
    text_length: usize,
    score: usize,
) -> bool {
    // Parameters
    let alignment_k: usize = affine_lambda_wavefront_diagonal(text_length, pattern_length);
    let alignment_offset: usize = affine_lambda_wavefront_offset(text_length, pattern_length);

    // Fetch wavefront and check termination
    let mwavefront: &AffineWavefront = &affine_wavefronts.mwavefronts[score];
    // TODO how to handle this null?
    if true
    // TODO: mwavefront!=NULL
    {
        let offsets: &Vec<wavefront_types::AwfOffset> = &mwavefront.offsets;
        if mwavefront.lo <= alignment_k
            && alignment_k <= mwavefront.hi
            && offsets[alignment_k] >= alignment_offset
        {
            return true;
        }
    }
    return false;
}
