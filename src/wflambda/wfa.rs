#![allow(dead_code, unused_variables)]
use std::cmp::{max};
use self::utils::*;
use self::types::*;
use self::backtrace_utils::*;

const VERBOSITY_LEVEL: usize = 0;
const NULL_OFFSET: isize = -10;

/*
k = (v-h) + (qlen-1)
 */

#[macro_use]
mod macros {
    macro_rules! abs {
        ($x: expr) => {{
            if $x >= 0 {
                $x
            } else {
                -$x
            }
        }};
    }

    #[allow(unused_macros)]
    macro_rules! v {
        ($offset: expr, $k: expr, $a_k: expr) => {{
            if $k <= $a_k {
                $offset
            } else {
                $offset + abs!($k - a_k )
            }
        }};
    }

    #[allow(unused_macros)]
    macro_rules! h {
        ($offset: expr, $k: expr, $a_k: expr) => {{
            if $k <= $a_k {
                $offset + abs!($k - a_k )
            } else {
                $offset
            }
        }};
    }
}

mod utils {
    use super::types::*;

    pub fn compute_k(k: isize, central_diagonal: usize) -> usize {
        (k + central_diagonal as isize) as usize
    }

    pub fn compute_kk(k: usize, central_diagonal: usize) -> isize {
        k as isize - central_diagonal as isize
    }

    pub fn v(offset: Offset, diagonal: usize, central_diagonal: usize) -> Offset {
        if diagonal <= central_diagonal {
            offset
        } else {
            offset + abs!(diagonal as isize - central_diagonal as isize ) as usize
        }
    }

    pub fn h(offset: Offset, diagonal: usize, central_diagonal: usize) -> Offset {
        if diagonal <= central_diagonal {
            offset + abs!(diagonal as isize - central_diagonal as isize ) as usize
        } else {
            offset
        }
    }

    pub fn compute_v(offset: Offset, diagonal: usize, central_diagonal: usize) -> Offset {
        v(offset, diagonal, central_diagonal)
    }

    pub fn compute_h(offset: Offset, diagonal: usize, central_diagonal: usize) -> Offset {
        h(offset, diagonal, central_diagonal)
    }

    // just make the cigar proper
    pub fn run_length_encode(cigar: &str, reverse: bool) -> String {

        let mut cigar = String::from(cigar);
        if reverse {
            cigar = cigar.chars().rev().collect::<String>();
        }

        let mut xcigar = String::new();

        // edge cases
        if cigar.len() == 0 {
            panic!("[run_length_encode] empty cigar");
        } else if cigar.len() == 1 {
            xcigar.push_str(&format!("{}{}", 1, cigar));
            xcigar
        } else {
            let mut chars = cigar.chars();
            let mut current: char = chars.next().unwrap();
            let mut count = 1;

            for c in chars {
                if c == current {
                    count += 1;
                } else {
                    xcigar.push_str(&format!("{}{}", count, current));
                    current = c;
                    count = 1;
                }
            }

            // last run
            xcigar.push_str(&format!("{}{}", count, current));

            xcigar
        }
    }
}

mod types {
    use std::cmp::max;

    pub type Offset = usize;
    pub type DpMatrix = Vec<Vec<Option<usize>>>;

    pub struct Alignment {
        pub score: usize,
        pub cigar: String,
    }

    #[derive(PartialEq, Debug)]
    pub enum Operation {
        Insertion,
        Deletion,
        MatchMismatch
    }

    #[derive(Copy, Clone)]
    pub struct Penalties {
        pub mismatch: usize,
        pub matches: usize,
        pub gap_open: usize,
        pub gap_extend: usize,
    }

    #[derive(Debug)]
    pub struct Wavefront {
        // TODO: make usize?
        pub hi: isize,
        pub lo: isize,
        pub offsets: Vec<Offset>,
    }

    impl Wavefront {
        pub fn new(offset_len: usize) -> Self {
            Self {
                hi: 0,
                lo: 0,
                offsets: vec![0; offset_len]
            }
        }
    }

    pub struct WavefrontSet {
        pub mwavefront: Wavefront,
        pub dwavefront: Wavefront,
        pub iwavefront: Wavefront,
    }

    impl WavefrontSet {
        pub fn new(offset_len: usize) -> Self {
            Self {
                mwavefront: Wavefront::new(offset_len),
                dwavefront: Wavefront::new(offset_len),
                iwavefront: Wavefront::new(offset_len),
            }
        }
    }

    pub struct Wavefronts<'a> {
        pub query: &'a [u8],
        pub text: &'a [u8],

        pub wavefronts: Vec<WavefrontSet>,

        pub diagonals: usize,

        pub central_diagonal: usize, // k_zero
        // max_diagonal: usize,
        pub min_diagonal: isize,
        pub a_offset: usize,

        pub penalties: Penalties,

        pub dp_matrix: Vec<Vec<Option<usize>>>,
    }

    impl<'a> Wavefronts<'a> {
        pub fn new(query: &'a[u8], text: &'a[u8], penalties: Penalties) -> Self {

            let qlen = query.len();
            let tlen = text.len();
            let diagonals = qlen+tlen - 1;

            let mut wavefronts = Vec::with_capacity(diagonals);
            (0..=diagonals).for_each(|_| {
                wavefronts.push(WavefrontSet::new(diagonals));
            });

            let dp_matrix: Vec<Vec<Option<usize>>> = vec![vec![None; qlen]; tlen];

            let a_offset = max(qlen, tlen);

            Self {
                query,
                text,
                wavefronts,
                diagonals,
                penalties,
                central_diagonal: qlen - 1,
                // max_diagonal: diagonals - 2,
                min_diagonal: -(qlen as isize -1),
                a_offset,
                dp_matrix,
            }
        }

        pub fn get_wavefront(&self, score: usize) -> Option<&WavefrontSet> {
            self.wavefronts.get(score)
        }

        pub fn get_wavefront_mut(&mut self, score: usize) -> Option<&mut WavefrontSet> {
            self.wavefronts.get_mut(score)
        }

        pub fn print(&self) {
            let dp_matrix = &self.dp_matrix;
            let query = self.query;
            let text = self.text;

            // print chars
            eprint!("\t");
            dp_matrix[0].iter().enumerate().for_each(|(j, _)| {
                eprint!("{}\t", j);
            });
            eprint!("\n");

            // print col nums
            eprint!("\t");
            dp_matrix[0].iter().enumerate().for_each(|(j, _)| {
                eprint!("{}\t", query[j] as char);
            });

            eprint!("\n");
            dp_matrix.iter().enumerate().for_each(|(i, row)| {
                eprint!("{} {}\t", i, text[i] as char);
                row.iter().for_each(|col| {
                    match col {
                        Some(o) => { eprint!("{}", o) },
                        None =>  { eprint!("*") },
                    }
                    eprint!("\t");
                });
                eprint!("\n");
            });

            eprint!("\n\n");
        }

        pub fn print_tsv(&self) {
            let dp_matrix = &self.dp_matrix;
            let query = self.query;
            let text = self.text;

            println!("text\tquery\tscore\tqbase\ttbase");
            dp_matrix.iter().enumerate().for_each(|(i, row)| {
                row.iter().enumerate().for_each(|(j, col)| {
                    match col {
                        Some(score) => {
                            println!("{}\t{}\t{}\t{}\t{}",
                                     i, j, score, query[j] as char, text[i] as char);
                        },
                        None =>  {},
                    }
                });
            });
        }
    }
}

fn wf_extend<T>(mwavefront: &mut Wavefront,
                match_lambda: T,
                central_diagonal: usize,
                score: usize,
                dp_matrix: &mut Vec<Vec<Option<usize>>>)
where
    T: Fn(usize, usize, usize, usize, &mut Vec<Vec<Option<usize>>>) -> bool,
{
    let lo = mwavefront.lo;
    let hi = mwavefront.hi;

    if VERBOSITY_LEVEL > 2 {
        eprintln!("[wf_extend] Extending wavefront with score {}", score);
        eprintln!("\tlo={}, hi={}", lo, hi);
    }

    for k in lo..=hi {
        let _lo = k;
        let k: usize = compute_k(k, central_diagonal);

        if  k >= mwavefront.offsets.len() {
            if VERBOSITY_LEVEL > 3 {
                eprintln!("[wf_extend] k={} is therefore out of scope. Skipping", k);
            }
            continue;
        }
        let offset = mwavefront.offsets[k];
        let mut v: usize = v(offset, k, central_diagonal);
        let mut h: usize = h(offset, k, central_diagonal);

        if VERBOSITY_LEVEL > 4 {
            eprintln!("\tk={} offset={}", k, offset);
            eprintln!("\tpre extend k={} offset={} ({},{})", k, offset, v, h);
        }

        while match_lambda(v, h, mwavefront.offsets[k], score, dp_matrix) {
            mwavefront.offsets[k] += 1;
            v += 1;
            h += 1;
        }

        if VERBOSITY_LEVEL > 4 {
            eprintln!("\tpost extend k={} offset={} ({},{})", k, mwavefront.offsets[k], v, h);
            eprintln!("");
        }

    }
}

fn wf_expand(wavefronts: &mut Wavefronts, score: usize) -> (isize, isize) {

    let s: isize = score as isize;

    let x: isize = wavefronts.penalties.mismatch as isize;
    let o: isize = wavefronts.penalties.gap_open as isize;
    let e: isize = wavefronts.penalties.gap_extend as isize;

    let m_hi =  wavefronts.get_wavefront(score).unwrap().mwavefront.hi;
    let i_hi =  wavefronts.get_wavefront(score).unwrap().iwavefront.hi;
    let d_hi =  wavefronts.get_wavefront(score).unwrap().dwavefront.hi;

    let m_lo =  wavefronts.get_wavefront(score).unwrap().mwavefront.lo;
    let d_lo =  wavefronts.get_wavefront(score).unwrap().dwavefront.lo;
    let i_lo =  wavefronts.get_wavefront(score).unwrap().iwavefront.lo;

    let hi = vec![
        if s-x   < 0 { m_hi } else { wavefronts.get_wavefront((s-x) as usize).unwrap().mwavefront.hi },
        if s-o-e < 0 { m_hi } else { wavefronts.get_wavefront((s-o-e) as usize).unwrap().mwavefront.hi },
        if s-e   < 0 { i_hi } else { wavefronts.get_wavefront((s-e) as usize).unwrap().iwavefront.hi },
        if s-e   < 0 { d_hi } else { wavefronts.get_wavefront((s-e) as usize).unwrap().dwavefront.hi},
    ].iter().max().unwrap() + 1;

    let lo = vec![
        if s-x   < 0 { m_lo } else { wavefronts.get_wavefront((s-x) as usize).unwrap().mwavefront.lo },
        if s-o-e < 0 { m_lo } else { wavefronts.get_wavefront((s-o-e) as usize).unwrap().mwavefront.lo },
        if s-e   < 0 { i_lo } else { wavefronts.get_wavefront((s-e) as usize).unwrap().iwavefront.lo },
        if s-e   < 0 { d_lo } else { wavefronts.get_wavefront((s-e) as usize).unwrap().dwavefront.lo},
    ].iter().min().unwrap() - 1;

    wavefronts.get_wavefront_mut(score).unwrap().mwavefront.hi = hi;
    wavefronts.get_wavefront_mut(score).unwrap().dwavefront.hi = hi;
    wavefronts.get_wavefront_mut(score).unwrap().iwavefront.hi = hi;

    wavefronts.get_wavefront_mut(score).unwrap().iwavefront.lo = lo;
    wavefronts.get_wavefront_mut(score).unwrap().mwavefront.lo = lo;
    wavefronts.get_wavefront_mut(score).unwrap().dwavefront.lo = lo;

    (lo, hi)
}

fn wf_next(wavefronts: &mut Wavefronts, score: usize) {
    if VERBOSITY_LEVEL > 2 {
        eprintln!("[wf_next] Computing wavefront for score {}", score);
    }

    let s: isize = score as isize;

    let x: isize = wavefronts.penalties.mismatch as isize;
    let o: isize = wavefronts.penalties.gap_open as isize;
    let e: isize = wavefronts.penalties.gap_extend as isize;

    let num_wavefronts = wavefronts.wavefronts.len();
    if num_wavefronts <= score {
        (num_wavefronts..=score).for_each(|_| {
            wavefronts.wavefronts.push(WavefrontSet::new(wavefronts.diagonals));
        });
    }

    let (lo, hi) = wf_expand(wavefronts, score);

    if VERBOSITY_LEVEL > 3 {
        eprintln!("\tk'\tk\tmmax\timax\tdmax");
    }

    for k in lo..=hi {
        if k < wavefronts.min_diagonal {
            if VERBOSITY_LEVEL > 3 {
                eprintln!("[wf_next] k={} is therefore out of scope. Skipping", k);
            }
            continue;
        }

        let k: usize = compute_k(k, wavefronts.central_diagonal);

        let imax = *vec![
            if s-o-e < 0 || k <= 0 || k+1 >= wavefronts.diagonals { 0 } else { wavefronts.get_wavefront((s-o-e) as usize).unwrap().mwavefront.offsets[(k-1) as usize] as isize },
            if s-e < 0   || k <= 0 || k+1 >= wavefronts.diagonals { 0 } else { wavefronts.get_wavefront((s-e) as usize).unwrap().iwavefront.offsets[(k-1) as usize] as isize },
        ].iter().max().unwrap();

        if k < wavefronts.diagonals {
            wavefronts.get_wavefront_mut(score).unwrap().iwavefront.offsets[k] = imax as usize;
        }

        let dmax = *vec![
            if s-o-e < 0 || k+1 >= wavefronts.diagonals { 0 } else { wavefronts.get_wavefront((s-o-e) as usize).unwrap().mwavefront.offsets[k+1] as isize },
            if s-e < 0   || k+1 >= wavefronts.diagonals { 0 } else { wavefronts.get_wavefront((s-e) as usize).unwrap().dwavefront.offsets[k+1] as isize},
        ].iter().max().unwrap();

        if k < wavefronts.diagonals  {
            wavefronts.get_wavefront_mut(score).unwrap().dwavefront.offsets[k] = dmax as usize;
        }

        let mmax = *vec![
            if s-x < 0 || k >= wavefronts.diagonals { 0 } else { (wavefronts.get_wavefront((s-x) as usize).unwrap().mwavefront.offsets[k] + 1 as usize) as isize},
            if            k >= wavefronts.diagonals { 0 } else { wavefronts.get_wavefront(score).unwrap().iwavefront.offsets[k] as isize },
            if            k >= wavefronts.diagonals { 0 } else { wavefronts.get_wavefront(score).unwrap().dwavefront.offsets[k] as isize },
        ].iter().max().unwrap();

        if k < wavefronts.diagonals {
            wavefronts.get_wavefront_mut(score).unwrap().mwavefront.offsets[k] = mmax as usize;
        }

        if VERBOSITY_LEVEL > 3 {
            let k_prime = k as isize - wavefronts.central_diagonal as isize;
            eprintln!("\t{}\t{}\t{}\t{}\t{}", k_prime, k, mmax, imax, dmax);
        }
    }
}

mod backtrace_utils {
    use super::utils::*;
    use super::types::*;

    const NULL_OFFSET: isize = -10;

    pub fn backtrace_matches_check(offset: &mut isize, cigar: &mut String, num_matches: isize, k: isize) {
        let k = compute_k(k, 5);

        (0..num_matches).for_each(|_| {
            let v = compute_v(*offset as usize, k, 5);
            let h = compute_h(*offset as usize, k, 5);

            cigar.push('M');
            *offset -= 1;
        });
    }

    pub fn backtrace_deletion_extend_offset(wavefronts: &Wavefronts, score: isize, k: isize) -> isize {
        if score < 0 { return NULL_OFFSET }

        let d = &wavefronts.get_wavefront(score as usize).unwrap().dwavefront;
        if d.lo <= k+1 && k+1 <= d.hi {
            let k = compute_k(k, wavefronts.central_diagonal);
            return d.offsets[k + 1] as isize
        } else {
            return NULL_OFFSET
        }
    }

    pub fn backtrace_deletion_open_offset(wavefronts: &Wavefronts, score: isize, k: isize) -> isize {
        if score < 0 { return NULL_OFFSET }

        let m = &wavefronts.get_wavefront(score as usize).unwrap().mwavefront;

        if m.lo <= k+1 && k+1 <= m.hi {
            let k = compute_k(k, wavefronts.central_diagonal);
            return m.offsets[k + 1] as isize;
        } else {
            return NULL_OFFSET
        }
    }

    pub fn backtrace_insertion_extend_offset(wavefronts: &Wavefronts, score: isize, k: isize) -> isize {
        if score < 0 { return NULL_OFFSET }


        let i = &wavefronts.get_wavefront(score as usize).unwrap().iwavefront;

        if i.lo <= k-1 && k-1 <= i.hi {
            let k = compute_k(k, wavefronts.central_diagonal);
            return i.offsets[k -1] as isize;
        } else {
            return NULL_OFFSET
        }
    }

    pub fn backtrace_insertion_open_offset(wavefronts: &Wavefronts, score: isize, k: isize) -> isize {
        if score < 0 { return NULL_OFFSET }

        let m = &wavefronts.get_wavefront(score as usize).unwrap().mwavefront;

        if m.lo <= k-1 && k-1 <= m.hi {
            let k = compute_k(k, wavefronts.central_diagonal);
            return m.offsets[k - 1]  as isize + 1;
        } else {
            return NULL_OFFSET
        }
    }

    pub fn backtrace_insertion_mismatch_offset(wavefronts: &Wavefronts, score: isize, k: isize) -> isize {
        if score < 0 { return NULL_OFFSET }

        let m = &wavefronts.get_wavefront(score as usize).unwrap().mwavefront;

        if m.lo <= k && k <= m.hi {
            let k = compute_k(k, wavefronts.central_diagonal);
            return m.offsets[k]  as isize + 1;
        } else {
            return NULL_OFFSET
        }
    }

    pub fn backtrace_mismatch_offset(wavefronts: &Wavefronts, score: isize, k: isize) -> isize {
        if score < 0 { return NULL_OFFSET }

        let m = &wavefronts.get_wavefront(score as usize).unwrap().mwavefront;

        if m.lo <= k && k <= m.hi {
            // eprintln!("\t[backtrace_mismatch_offset] score={} k={} offset={}", score, k, m.offsets[k as usize]);
            let k = compute_k(k, wavefronts.central_diagonal);
            return m.offsets[k]  as isize + 1;
        } else {
            return NULL_OFFSET
        }
    }
}

#[allow(unused_variables, unused_mut)]
fn backtrace(wavefronts: &mut Wavefronts, score: usize) -> String {
    if VERBOSITY_LEVEL > 2 {
        eprintln!("[backtrace]");
    }

    let a_k = wavefronts.central_diagonal;
    let a_offset = wavefronts.a_offset;
    let text = wavefronts.text;
    let query = wavefronts.query;

    let mut cigar = String::new();
    let mut q_aln = String::new();
    let mut t_aln = String::new();

    let wf = wavefronts.get_wavefront(score).unwrap();

    let mut score: isize = score as isize;
    let mut offset: isize = wf.mwavefront.offsets[a_k] as isize;

    let mut k = compute_kk(a_k, a_k);

    let mut v: usize = compute_v(offset as usize, k as usize, a_k);
    let mut h: usize = compute_h(offset as usize, k as usize, a_k);

    let mut backtrace_type = Operation::MatchMismatch;

    let x: isize = wavefronts.penalties.mismatch as isize;
    let o: isize = wavefronts.penalties.gap_open as isize;
    let e: isize = wavefronts.penalties.gap_extend as isize;

    while v > 0 && h > 0 && score > 0 {
        if VERBOSITY_LEVEL > 4 {
            eprintln!("\tbacktrace_type = {:?}", backtrace_type);
        }
        // compute scores
        let gap_open_score = score - o - e;
        let gap_extend_score = score - e;
        let mismatch_score = score - x;

        let del_ext: isize = if backtrace_type == Operation::Insertion {NULL_OFFSET} else {
            backtrace_deletion_extend_offset(&wavefronts, gap_extend_score, k) };
        let del_open: isize = if backtrace_type == Operation::Insertion { NULL_OFFSET } else {
            backtrace_deletion_open_offset(&wavefronts, gap_open_score, k) };
        let ins_ext: isize = if backtrace_type == Operation::Deletion { NULL_OFFSET } else {
            backtrace_insertion_extend_offset(&wavefronts, gap_extend_score, k)
        };
        let ins_open: isize = if backtrace_type == Operation::Deletion { NULL_OFFSET } else {
            backtrace_insertion_open_offset(&wavefronts, gap_open_score, k)
        };
        let misms: isize = if backtrace_type != Operation::MatchMismatch { NULL_OFFSET } else {
            backtrace_mismatch_offset(&wavefronts, mismatch_score, k)
        };

        // Compute maximum offset
        let max_all = *vec![del_ext, del_open, ins_ext, ins_open, misms].iter().max().unwrap();

        if VERBOSITY_LEVEL > 4 {
            eprintln!("\tscore={} offset={} k={} \
                       gap_open_score={} gap_extend_score={} mismatch_score={} \
                       max_all={} del_ext={} del_open={} ins_ext={} ins_open={} mims={} \
                       backtrace_type={:?}",
                      score, offset, k,
                      gap_open_score, gap_extend_score, mismatch_score,
                      max_all, del_ext, del_open, ins_ext, ins_open, misms,
                      backtrace_type
            );
        }

        // Traceback Matches
        if backtrace_type == Operation::MatchMismatch {
            let num_matches = offset - max_all;
            // eprintln!("num_matches={}", num_matches);
            backtrace_matches_check(&mut offset, &mut cigar, num_matches, k);
            offset = max_all;
        }

        if max_all == del_ext {
            // Add Deletion
            cigar.push('D');
            // Update state
            score = gap_extend_score;
            k += 1;
            backtrace_type = Operation::Deletion;
        } else if max_all == del_open {
            // Add Deletion
            cigar.push('D');
            // Update state
            score = gap_open_score;
            k += 1;
            backtrace_type = Operation::MatchMismatch;
        } else if max_all == ins_ext {
            // Add Insertion
            cigar.push('I');
            // Update state
            score = gap_extend_score;
            k -= 1;
            offset -= 1;
            backtrace_type = Operation::Insertion;
        } else if max_all == ins_open {
            // Add Insertion
            cigar.push('I');
            // Update state
            score = gap_open_score;
            k -= 1;
            offset -= 1;
            backtrace_type = Operation::MatchMismatch;
        } else if max_all == misms {
            // Add Mismatch
            cigar.push('X');

            // Update state
            score = mismatch_score;
            offset -= 1;
        } else {
            panic!("Backtrace error: No link found during backtrace\n")
        }

        // Update coordinates
        v = compute_v(offset as usize, k as usize, a_k);
        h = compute_h(offset as usize, k as usize, a_k);
    }

    if score == 0 {
        // backtrace matches check
        let num_matches = offset;
        backtrace_matches_check(&mut offset, &mut cigar, num_matches, k);
    } else {

        while v > 0 {
            cigar.push('D');
            v -= 1;
        }

        while h > 0 {
            cigar.push('I');
            h -= 1;
        }
    }

    cigar = run_length_encode(&cigar[..], true);
    cigar
}

fn wf_align(text: &[u8], query: &[u8], penalties: Penalties) -> Alignment {
    let mut wavefronts = Wavefronts::new(query, text, penalties);

    let qlen = query.len();
    let tlen = text.len();

    // central diagonal
    let a_k: usize = qlen - 1;

    // furthest offset along...
    let a_offset: usize = max(qlen, tlen);

    let mut score = 0;

    let exit_condition = |wavefronts: &Wavefronts, score: usize| {
        *wavefronts.get_wavefront(score).unwrap().mwavefront.offsets.get(a_k).unwrap() >= a_offset
    };

    let match_lambda =
        |v: usize, h: usize, offset: usize, score: usize, dp_matrix:  &mut DpMatrix| -> bool {
        if h < tlen && v < qlen && dp_matrix[h][v] == None {
            dp_matrix[h][v] = Some(offset);
        }

        h < tlen && v < qlen && text[h] == query[v]
    };

    loop {
        let dp_matrix = &mut wavefronts.dp_matrix;
        let m_s = &mut wavefronts.wavefronts.get_mut(score).unwrap().mwavefront;
        wf_extend(m_s, match_lambda, a_k, score, dp_matrix);

        if exit_condition(&wavefronts, score) {
            if VERBOSITY_LEVEL > 2 {
                eprintln!("Final state of the DP table");
                eprintln!("---------------------------");
                wavefronts.print();
            }

            if VERBOSITY_LEVEL > 0 {
                eprintln!("\tscore: {}\n\
                           \tcentral diagonal (a_k): {}\n\
                           \tmaximum offset (a_offset): {}\n\
                           \tquery: {}\n\
                           \ttext:  {}\n",
                          score,
                          a_k,
                          a_offset,
                          std::str::from_utf8(query).unwrap(),
                          std::str::from_utf8(text).unwrap());
                eprintln!("");
            }

            if VERBOSITY_LEVEL > 3 {
                wavefronts.print_tsv();
            }

            //wavefronts.print_tsv();
            let cigar = backtrace(&mut wavefronts, score);
            return Alignment {score, cigar};
        }

        if VERBOSITY_LEVEL > 2 {
            wavefronts.print();
        }

        score += 1;

        wf_next(&mut wavefronts, score);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    static PENALTIES: Penalties =  Penalties {
        mismatch: 4,
        matches: 0,
        gap_open: 6,
        gap_extend: 2,
    };

    mod backtrace {
        use super::super::*;
        use super::PENALTIES;

        #[test]
        fn test_same_sequence() {
            // same sequence
            let text  = "GAGATA";
            let query = "GAGATA";

            let aln = wf_align(text.as_bytes(), query.as_bytes(), PENALTIES);
            assert_eq!(aln.score, 0);
            assert_eq!(aln.cigar, String::from("6M"));
        }

        #[test]
        fn test_snp() {
            // same sequence
            let text  = "GACATA";
            let query = "GAGATA";

            let aln = wf_align(text.as_bytes(), query.as_bytes(), PENALTIES);
            assert_eq!(aln.score, 4);
            assert_eq!(aln.cigar, String::from("2M1X3M"));
        }

        #[test]
        fn test_paper_example() {
            let text  = "GATACA";
            let query = "GAGATA";

            let aln = wf_align(text.as_bytes(), query.as_bytes(), PENALTIES);
            assert_eq!(aln.score, 8);
            assert_eq!(aln.cigar, String::from("2M1X1M1X1M"));
        }

        #[test]
        fn test_long_sequences() {
            let text  = "TCTATACTGCGCGTTTGGAGAAATAAAATAGTTCTATACTGCGCGTTTGGAGAAATAAAATAGTTCTATACTGCGCGTTTGGAGAAATAAAATAGTTCTATACTGCGCGTTTGGAGAAATAAAATAGT";
            let query = "TCTTTACTCGCGCGTTGGAGAAATACAATAGTTCTTTACTCGCGCGTTGGAGAAATACAATAGTTCTTTACTCGCGCGTTGGAGAAATACAATAGTTCTTTACTCGCGCGTTGGAGAAATACAATAGT";

            let aln = wf_align(text.as_bytes(), query.as_bytes(), PENALTIES);
            assert_eq!(aln.score, 96);
            assert_eq!(aln.cigar, String::from("3M1X4M1I6M1D10M1X9M1X4M1I6M1D10M1X9M1X4M1I6M1D10M1X9M1X4M1I6M1D10M1X6M"));
            // assert_eq!(aln.cigar, String::from("3M1X4M1D7M1I9M1X9M1X4M1D7M1I9M1X9M1X4M1D7M1I9M1X9M1X4M1D7M1I9M1X6M"));
        }
    }

    mod utils {
        use super::super::*;

        #[test]
        fn test_run_length_encode() {
            let s = "DDMMMMMMMMMMMMMMMM";
            assert_eq!("2D16M", run_length_encode(s, false));

            let s = "DDMMMMMMMMMMMMMMMM";
            assert_eq!("2D16M", run_length_encode(s, false));

            let s = "DDMMMMMMXMMMMMMMMMM";
            assert_eq!("2D6M1X10M", run_length_encode(s, false));
        }
    }

    mod align {
        use super::super::*;
        use super::PENALTIES;

        #[test]
        fn test_same_sequence() {
            // same sequence
            let text  = "GAGATA";
            let query = "GAGATA";

            let aln = wf_align(text.as_bytes(), query.as_bytes(), PENALTIES);
            assert_eq!(aln.score, 0);
        }

        #[test]
        fn test_paper_test_case() {
            let text  = "GAT";
            let query = "GAG";
            let aln = wf_align(text.as_bytes(), query.as_bytes(), PENALTIES);
            assert_eq!(aln.score, 4);

            let text  = "GATACA";
            let query = "GAGATA";
            let aln = wf_align(text.as_bytes(), query.as_bytes(), PENALTIES);
            assert_eq!(aln.score, 8);
        }

        #[test]
        fn test_c_example() {
            let query = "TCTTTACTCGCGCGTTGGAGAAATACAATAGT";
            let text  = "TCTATACTGCGCGTTTGGAGAAATAAAATAGT";

            let aln = wf_align(&text.as_bytes()[..10], &query.as_bytes()[..10], PENALTIES);
            assert_eq!(aln.score, 12);

            let aln = wf_align(text.as_bytes(), query.as_bytes(), PENALTIES);
            assert_eq!(aln.score, 24);
        }

        #[test]
        fn test_same_length() {
            let text  = "TCTATACTGCGCGTTTGGAGAAATAAAATAGTTCTATACTGCGCGTTTGGAGAAATAAA\
                         ATAGTTCTATACTGCGCGTTTGGAGAAATAAAATAGTTCTATACTGCGCGTTTGGAGAA\
                         ATAAAATAGT";
            let query = "TCTTTACTCGCGCGTTGGAGAAATACAATAGTTCTTTACTCGCGCGTTGGAGAAATACA\
                         ATAGTTCTTTACTCGCGCGTTGGAGAAATACAATAGTTCTTTACTCGCGCGTTGGAGAA\
                         ATACAATAGT";

            let _p =  Penalties {
                mismatch: 7,
                matches: 0,
                gap_open: 11,
                gap_extend: 1,
            };

            let aln = wf_align(text.as_bytes(), query.as_bytes(), PENALTIES);
            assert_eq!(aln.score, 96);
        }

        #[test]
        fn test_longer_query() {
            let text  = "TCTTTACTCGCGCGTTGGAGAAATACAATAGTTCTTTACTCGCGCGTTGGAGAAATACAATAGT\
                         TCTTTACTCGCGCGTTGGAGAAATACAATAGTTCTTTACTCGCGCGTTGGAGAAATACAATAGT";
            let query = "TCTATACTGCGCGTTTATCTAGGAGAAATAAAATAGTTCTATACTGCGCGTTTGGAGAAATAAAATAGT\
                         TCTATACTGCGCGTTTGGAGAAATAACTATCAATAGTTCTATACTGCGCGTTTGGAGAAATAAAATAGT";

            let aln = wf_align(text.as_bytes(), query.as_bytes(), PENALTIES);
            assert_eq!(aln.score, 200);
        }

        #[test]
        fn test_longer_text() {
            let text  = "TCTATACTGCGCGTTTATCTAGGAGAAATAAAATAGTTCTATACTGCGCGTTTGGAGAAATAAAATAGT\
                         TCTATACTGCGCGTTTGGAGAAATAACTATCAATAGTTCTATACTGCGCGTTTGGAGAAATAAAATAGT";
            let query = "TCTTTACTCGCGCGTTGGAGAAATACAATAGTTCTTTACTCGCGCGTTGGAGAAATACAATAGT\
                         TCTTTACTCGCGCGTTGGAGAAATACAATAGTTCTTTACTCGCGCGTTGGAGAAATACAATAGT";

            let aln = wf_align(text.as_bytes(), query.as_bytes(), PENALTIES);
            assert_eq!(aln.score, 200);
        }
    }
}
