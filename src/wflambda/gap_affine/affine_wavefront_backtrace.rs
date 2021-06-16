use super::affine_wavefront;
use super::affine_wavefront_penalties as penalties;
use super::super::edit::edit_cigar;
use std::cmp;

/*
 * WF type
 */
#[derive(PartialEq)]
pub enum BacktraceWavefront {
    BacktraceWavefrontM = 0,
    BacktraceWavefrontI = 1,
    BacktraceWavefrontD = 2
}

fn affine_wavefronts_valid_location(
    k: i32,
    offset: affine_wavefront::AwfOffset,
    pattern_length: i32,
    text_length: i32) -> bool {
    // Locate offset (remember that backtrace is always +1 offset ahead)
    let v: i32 = affine_wavefront::affine_lambda_wavefront_v(k,offset as i32);
    let h: i32 = affine_wavefront::affine_lambda_wavefront_h(k,offset as i32);

    v > 0 && v <= pattern_length && h > 0 && h <= text_length
}

fn backtrace_wavefront_trace_deletion_extend_offset(
    affine_wavefronts: &mut affine_wavefront::AffineWavefronts,
    score: i32,
    k: i32,
    offset: affine_wavefront::AwfOffset) -> affine_wavefront::AwfOffset {

    if score < 0 { return affine_wavefront::AFFINE_LAMBDA_WAVEFRONT_OFFSET_NULL; }

    match affine_wavefronts.dwavefronts.get(score as usize) {
        Some(dwavefront) => {
            if dwavefront.lo_base <= k+1 && k+1 <= dwavefront.hi_base {
                dwavefront.offsets[k as usize +1]
            } else {
                affine_wavefront::AFFINE_LAMBDA_WAVEFRONT_OFFSET_NULL
            }
        },
        _ => { affine_wavefront::AFFINE_LAMBDA_WAVEFRONT_OFFSET_NULL }
    }

}

fn backtrace_wavefront_trace_insertion_open_offset(
    affine_wavefronts: &mut affine_wavefront::AffineWavefronts,
    score: i32,
    k: i32,
    offset: affine_wavefront::AwfOffset) -> affine_wavefront::AwfOffset  {

    if score < 0 { return affine_wavefront::AFFINE_LAMBDA_WAVEFRONT_OFFSET_NULL; }

    match affine_wavefronts.mwavefronts.get(score as usize) {
        Some(mwavefront) => {
            if mwavefront.lo_base <= k-1 && k-1 <= mwavefront.hi_base {
                mwavefront.offsets[k as usize - 1] + 1
            } else {
                affine_wavefront::AFFINE_LAMBDA_WAVEFRONT_OFFSET_NULL
            }
        },
        _ => {
            affine_wavefront::AFFINE_LAMBDA_WAVEFRONT_OFFSET_NULL
        }
    }
}

pub fn affine_wavefronts_backtrace_matches__check<T>(
    affine_wavefronts: &mut affine_wavefront::AffineWavefronts,
    lambda: T,
    k: i32,
    offset: affine_wavefront::AwfOffset,
    valid_location: bool,
    num_matches: i32,
    edit_cigar: &mut edit_cigar::EditCigar)
where T: Fn(i32, i32) -> bool {
    let i: i32;
    for i in 0..num_matches {
        // run our traceback function on the alignments
        let v: i32 = affine_wavefront::affine_lambda_wavefront_v(k,offset as i32);
        let h: i32 = affine_wavefront::affine_lambda_wavefront_h(k,offset as i32);

        if !valid_location { // Check inside table
            panic!("Backtrace error: Match outside DP-Table\n");
        } else if !lambda(v-1, h-1) { // Check match
            panic!("Backtrace error: Not a match traceback\n");
        }
        // Set Match
        edit_cigar.operations[(edit_cigar.begin_offset) -= 1] = 'M';
        // Update state
        offset -= 1;
    }
}

fn affine_wavefronts_offset_add_trailing_gap(
    edit_cigar: &mut edit_cigar::EditCigar,
    k: i32,
    alignment_k: i32) {
    // Parameters
    let operations: &str = &edit_cigar.operations[..];
    let op_sentinel: i32 = edit_cigar.begin_offset;

    // Add trailing gap
    let i: i32;
    if k < alignment_k {
        for i in k..alignment_k { operations[op_sentinel -= 1] = 'I'; }
    } else if k > alignment_k {
        for i in alignment_k..k { operations[op_sentinel -= 1] = 'D'; }
    }
    edit_cigar.begin_offset = op_sentinel;
}

/*
 * Backtrace (single solution)
 */
fn affine_wavefronts_backtrace<T>(
    affine_wavefronts: &mut affine_wavefront::AffineWavefronts,
    lambda: T,
    pattern_length: i32,
    text_length: i32,
    alignment_score: i32)
where T: Fn(i32, i32) -> bool {
    // STATS
    // WAVEFRONT_STATS_TIMER_START(affine_wavefronts,wf_time_backtrace);
    // Parameters
    let wavefront_penalties: &penalties::AffinePenalties = &affine_wavefronts.penalties.wavefront_penalties;
    let cigar: &mut edit_cigar::EditCigar = &mut affine_wavefronts.edit_cigar;
    let alignment_k: i32 = affine_wavefront::affine_lambda_wavefront_diagonal(text_length,pattern_length);
    // Compute starting location
    let score: i32 = alignment_score;
    let k: i32 = alignment_k;
    let offset: affine_wavefront::AwfOffset = affine_wavefronts.mwavefronts[alignment_score as usize].offsets[k as usize];
    let valid_location: bool = affine_wavefronts_valid_location(k,offset,pattern_length,text_length);
    // Trace the alignment back
    let backtrace_type: BacktraceWavefront = BacktraceWavefront::BacktraceWavefrontM;
    let v: i32 = affine_wavefront::affine_lambda_wavefront_v(k,offset as i32);
    let h: i32 =  affine_wavefront::affine_lambda_wavefront_h(k,offset as i32);
    while v > 0 && h > 0 && score > 0 {
        // Check location
        if !valid_location {
            valid_location = affine_wavefronts_valid_location(k,offset,pattern_length,text_length);
            if valid_location {
                affine_wavefronts_offset_add_trailing_gap(cigar,k,alignment_k);
            }
        }
        // Compute scores
        let gap_open_score: i32 = score - wavefront_penalties.gap_opening - wavefront_penalties.gap_extension;
        let gap_extend_score: i32 = score - wavefront_penalties.gap_extension;
        let mismatch_score: i32 = score - wavefront_penalties.mismatch;

        // Compute source offsets
        let del_ext: affine_wavefront::AwfOffset = if backtrace_type == BacktraceWavefront::BacktraceWavefrontI {
            affine_wavefront::AFFINE_LAMBDA_WAVEFRONT_OFFSET_NULL
        } else {
            backtrace_wavefront_trace_deletion_extend_offset(affine_wavefronts,gap_extend_score,k,offset)
        };

        let del_open: affine_wavefront::AwfOffset = if backtrace_type == BacktraceWavefront::BacktraceWavefrontI {
            affine_wavefront::AFFINE_LAMBDA_WAVEFRONT_OFFSET_NULL
        } else {
            backtrace_wavefront_trace_deletion_open_offset(affine_wavefronts,gap_open_score,k,offset)
        };

        let ins_ext: affine_wavefront::AwfOffset = if backtrace_type == BacktraceWavefront::BacktraceWavefrontD {
            AFFINE_LAMBDA_WAVEFRONT_OFFSET_NULL
        } else {
            backtrace_wavefront_trace_insertion_extend_offset(affine_wavefronts,gap_extend_score,k,offset)
        };

        let ins_open: affine_wavefront::AwfOffset = if backtrace_type == BacktraceWavefront::BacktraceWavefrontD {
            affine_wavefront::AFFINE_LAMBDA_WAVEFRONT_OFFSET_NULL
        } else {
            backtrace_wavefront_trace_insertion_open_offset(affine_wavefronts,gap_open_score,k,offset)
        };

        let misms: affine_wavefront::AwfOffset = if backtrace_type != BacktraceWavefront::BacktraceWavefrontM {
            affine_wavefront::AFFINE_LAMBDA_WAVEFRONT_OFFSET_NULL
        } else {
            backtrace_wavefront_trace_mismatch_offset(affine_wavefronts,mismatch_score,k,offset)
        };

        // Compute maximum offset
        let max_del: affine_wavefront::AwfOffset = cmp::max(del_ext,del_open);
        let max_ins: affine_wavefront::AwfOffset = cmp::max(ins_ext,ins_open);
        let max_all: affine_wavefront::AwfOffset = cmp::max(misms,cmp::max(max_ins,max_del));

        // Traceback Matches
        if backtrace_type == BacktraceWavefront::BacktraceWavefrontM {
            let num_matches: i32 = offset as i32 - max_all as i32;
            affine_wavefronts_backtrace_matches__check(
                affine_wavefronts,
                lambda,
                k,
                offset,
                valid_location,
                num_matches,
                cigar);
            offset = max_all;
        }

        // Traceback Operation
        if max_all == del_ext {
            // Add Deletion
            if valid_location { cigar.operations[cigar.begin_offset -= 1] = 'D' };
            // Update state
            score = gap_extend_score;
            k += 1;
            backtrace_type = BacktraceWavefront::BacktraceWavefrontD;
        } else if max_all == del_open {
            // Add Deletion
            if valid_location { cigar.operations[cigar.begin_offset -= 1] = 'D' };
            // Update state
            score = gap_open_score;
            k += 1;
            backtrace_type = BacktraceWavefront::BacktraceWavefrontM;
        } else if max_all == ins_ext {
            // Add Insertion
            if valid_location { cigar.operations[cigar.begin_offset -= 1] = 'I' };
            // Update state
            score = gap_extend_score;
            k -= 1;
            offset -= 1;
            backtrace_type = BacktraceWavefront::BacktraceWavefrontI;
        } else if max_all == ins_open {
            // Add Insertion
            if valid_location cigar.operations[cigar.begin_offset -= 1] = 'I';
            // Update state
            score = gap_open_score;
            k -= 1;
            offset -= 1;
            backtrace_type = BacktraceWavefront::BacktraceWavefrontM;
        } else if max_all == misms {
            // Add Mismatch
            if valid_location cigar.operations[cigar.begin_offset -= 1] = 'X';
            // Update state
            score = mismatch_score;
            offset -= 1;
        } else {
            panic!("Backtrace error: No link found during backtrace\n");
        }

        // Update coordinates
        v = affine_wavefront::affine_lambda_wavefront_v(k,offset as i32);
        h = affine_wavefront::affine_lambda_wavefront_h(k,offset as i32);
    }
    // Account for last operations
    if score == 0 {
        // Account for last stroke of matches
        affine_wavefronts_backtrace_matches__check(
            affine_wavefronts,
            lambda,
            k,
            offset,
            valid_location,
            offset as i32,
            cigar);
    } else {
        // Account for last stroke of insertion/deletion
        while v > 0 {
            cigar.operations[cigar.begin_offset -= 1] = 'D';
            v -= 1;
        };
        while h > 0 {
            cigar.operations[cigar.begin_offset -= 1] = 'I';
            h -= 1;
        };
    }
    cigar.begin_offset += 1; // Set CIGAR length
    // STATS
    //WAVEFRONT_STATS_TIMER_STOP(affine_wavefronts,wf_time_backtrace);
}
