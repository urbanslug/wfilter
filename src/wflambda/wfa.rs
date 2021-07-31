#![allow(dead_code, unused_variables)]
use self::backtrace_utils::*;
use self::types::*;
use self::utils::*;
use std::cmp::max;

use super::super::types::CliArgs;

use indicatif::{ProgressBar, ProgressStyle};

const VERBOSITY_LEVEL: usize = 2;
const NULL_OFFSET: isize = -10;
const MIN_WAVEFRONT_LENGTH: isize = 10;
const MAX_DISTANCE_THRESHOLD: isize = 50;

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
                $offset + abs!($k - a_k)
            }
        }};
    }

    #[allow(unused_macros)]
    macro_rules! h {
        ($offset: expr, $k: expr, $a_k: expr) => {{
            if $k <= $a_k {
                $offset + abs!($k - a_k)
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
            offset + abs!(diagonal as isize - central_diagonal as isize) as usize
        }
    }

    pub fn h(offset: Offset, diagonal: usize, central_diagonal: usize) -> Offset {
        if diagonal <= central_diagonal {
            offset + abs!(diagonal as isize - central_diagonal as isize) as usize
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
            panic!("[wfa::utils::run_length_encode] empty cigar");
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

pub mod types {
    use super::super::super::types::Penalties;
    use std::cmp::max;

    use indicatif::{ProgressBar, ProgressStyle};
    use std::fs::OpenOptions;
    use std::io::prelude::*;

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
        MatchMismatch,
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
                offsets: vec![0; offset_len],
            }
        }
    }

    pub struct WavefrontSet {
        pub m: Option<Wavefront>,
        pub d: Option<Wavefront>,
        pub i: Option<Wavefront>,
    }

    impl WavefrontSet {
        pub fn new(offset_len: usize) -> Self {
            Self {
                m: Some(Wavefront::new(offset_len)),
                d: Some(Wavefront::new(offset_len)),
                i: Some(Wavefront::new(offset_len)),
            }
        }

        pub fn iwavefront(&self) -> Option<&Wavefront> {
            self.i.as_ref()
        }

        pub fn mwavefront(&self) -> Option<&Wavefront> {
            self.m.as_ref()
        }

        pub fn dwavefront(&self) -> Option<&Wavefront> {
            self.d.as_ref()
        }

        pub fn iwavefront_mut(&mut self) -> Option<&mut Wavefront> {
            self.i.as_mut()
        }

        pub fn dwavefront_mut(&mut self) -> Option<&mut Wavefront> {
            self.d.as_mut()
        }

        pub fn mwavefront_mut(&mut self) -> Option<&mut Wavefront> {
            self.m.as_mut()
        }

        pub fn fields(&self) -> (Option<&Wavefront>, Option<&Wavefront>, Option<&Wavefront>) {
            (self.m.as_ref(), self.i.as_ref(), self.d.as_ref())
        }

        pub fn fields_mut(
            &mut self,
        ) -> (
            Option<&mut Wavefront>,
            Option<&mut Wavefront>,
            Option<&mut Wavefront>,
        ) {
            (self.m.as_mut(), self.i.as_mut(), self.d.as_mut())
        }
    }

    pub type OptWavefrontSet = Option<WavefrontSet>;
    pub type BoxedWavefront = Box<OptWavefrontSet>;

    pub struct Wavefronts<'a> {
        pub query: &'a [u8],
        pub text: &'a [u8],

        pub wavefronts: Vec<BoxedWavefront>,

        pub diagonals: usize,

        pub central_diagonal: usize, // k_zero
        // max_diagonal: usize,
        pub min_diagonal: isize,
        pub a_offset: usize,

        pub penalties: Penalties,

        pub dp_matrix: Vec<Vec<Option<usize>>>,
    }

    impl<'a> Wavefronts<'a> {
        pub fn new(query: &'a [u8], text: &'a [u8], penalties: Penalties) -> Self {
            let qlen = query.len();
            let tlen = text.len();
            let diagonals = qlen + tlen - 1;

            let mut wavefronts = Vec::new();
            let x = Box::new(Some(WavefrontSet::new(diagonals)));
            wavefronts.push(x);
            // let mut wavefronts = Vec::with_capacity(diagonals);
            // (0..=diagonals).for_each(|_| {
            //    wavefronts.push(WavefrontSet::new(diagonals));
            // });

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
                min_diagonal: -(qlen as isize - 1),
                a_offset,
                dp_matrix,
            }
        }

        pub fn get_wavefront(&self, score: usize) -> Option<&WavefrontSet> {
            self.wavefronts.get(score).and_then(|x| (**x).as_ref())
        }

        pub fn get_wavefront_mut(&mut self, score: usize) -> Option<&mut WavefrontSet> {
            self.wavefronts.get_mut(score).and_then(|x| (**x).as_mut())
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
                        Some(o) => {
                            eprint!("{}", o)
                        }
                        None => {
                            eprint!("*")
                        }
                    }
                    eprint!("\t");
                });
                eprint!("\n");
            });

            eprint!("\n\n");
        }

        pub fn print_tsv(&self, filename: &str) {
            let dp_matrix = &self.dp_matrix;
            let query = self.query;
            let text = self.text;

            let mut file = OpenOptions::new()
                .write(true)
                .append(true)
                .create(true)
                .open(filename)
                .unwrap();

            let progress_bar = ProgressBar::new(dp_matrix.len() as u64);
            let template = "{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos:>7}/{len:7} ({eta_precise})";
            let progress_style = ProgressStyle::default_bar()
                .template(template)
                .progress_chars("=> ");
            progress_bar.set_style(progress_style);

            writeln!(file, "text\tquery\tscore\tqbase\ttbase").unwrap();

            dp_matrix.iter().enumerate().for_each(|(i, row)| {
                progress_bar.inc(1);
                row.iter().enumerate().for_each(|(j, col)| match col {
                    Some(score) => {
                        writeln!(
                            file,
                            "{}\t{}\t{}\t{}\t{}",
                            i, j, score, query[j] as char, text[i] as char
                        ).unwrap();
                    }
                    None => {}
                });
            });
        }
    }
}

mod backtrace_utils {
    use super::types::*;
    use super::utils::*;

    const NULL_OFFSET: isize = -10;

    pub fn backtrace_matches_check<T>(
        offset: &mut isize,
        cigar: &mut String,
        num_matches: usize,
        k: isize,
        central_diagonal: usize,
        backtrace_lambda: &mut T,
    ) where
        T: FnMut((i32, i32), (i32, i32)) -> (),
    {
        let k = compute_k(k, central_diagonal);

        {
            let query_start = compute_v(*offset as usize, k, central_diagonal);
            let target_start = compute_h(*offset as usize, k, central_diagonal);

            let query_stop = query_start + num_matches as usize;
            let target_stop = target_start + num_matches as usize;

            let query = (query_stop as i32, query_stop as i32);
            let target = (target_start as i32, target_stop as i32);

            backtrace_lambda(query, target);

            // println!("({}, {}) ({}, {})", query_start, query_stop, target_start, target_stop);
        }

        (0..num_matches).for_each(|_| {
            let v = compute_v(*offset as usize, k, central_diagonal);
            let h = compute_h(*offset as usize, k, central_diagonal);

            cigar.push('M');
            *offset -= 1;
        });
    }

    pub fn backtrace_deletion_extend_offset(
        wavefronts: &Wavefronts,
        score: isize,
        k: isize,
    ) -> isize {
        if score < 0 {
            return NULL_OFFSET;
        }

        let d = wavefronts
            .get_wavefront(score as usize)
            .unwrap()
            .dwavefront()
            .unwrap();
        if d.lo <= k + 1 && k + 1 <= d.hi {
            let k = compute_k(k, wavefronts.central_diagonal);
            return d.offsets[k + 1] as isize;
        } else {
            return NULL_OFFSET;
        }
    }

    pub fn backtrace_deletion_open_offset(
        wavefronts: &Wavefronts,
        score: isize,
        k: isize,
    ) -> isize {
        if score < 0 {
            return NULL_OFFSET;
        }

        let m = wavefronts
            .get_wavefront(score as usize)
            .unwrap()
            .mwavefront()
            .unwrap();

        if m.lo <= k + 1 && k + 1 <= m.hi {
            let k = compute_k(k, wavefronts.central_diagonal);
            return m.offsets[k + 1] as isize;
        } else {
            return NULL_OFFSET;
        }
    }

    pub fn backtrace_insertion_extend_offset(
        wavefronts: &Wavefronts,
        score: isize,
        k: isize,
    ) -> isize {
        if score < 0 {
            return NULL_OFFSET;
        }

        let i = wavefronts
            .get_wavefront(score as usize)
            .unwrap()
            .iwavefront()
            .unwrap();

        if i.lo <= k - 1 && k - 1 <= i.hi {
            let k = compute_k(k, wavefronts.central_diagonal);
            return i.offsets[k - 1] as isize;
        } else {
            return NULL_OFFSET;
        }
    }

    pub fn backtrace_insertion_open_offset(
        wavefronts: &Wavefronts,
        score: isize,
        k: isize,
    ) -> isize {
        if score < 0 {
            return NULL_OFFSET;
        }

        let m = wavefronts
            .get_wavefront(score as usize)
            .unwrap()
            .mwavefront()
            .unwrap();

        if m.lo <= k - 1 && k - 1 <= m.hi {
            let k = compute_k(k, wavefronts.central_diagonal);
            return m.offsets[k - 1] as isize + 1;
        } else {
            return NULL_OFFSET;
        }
    }

    pub fn backtrace_insertion_mismatch_offset(
        wavefronts: &Wavefronts,
        score: isize,
        k: isize,
    ) -> isize {
        if score < 0 {
            return NULL_OFFSET;
        }

        let m = wavefronts
            .get_wavefront(score as usize)
            .unwrap()
            .mwavefront()
            .unwrap();

        if m.lo <= k && k <= m.hi {
            let k = compute_k(k, wavefronts.central_diagonal);
            return m.offsets[k] as isize + 1;
        } else {
            return NULL_OFFSET;
        }
    }

    pub fn backtrace_mismatch_offset(wavefronts: &Wavefronts, score: isize, k: isize) -> isize {
        if score < 0 {
            return NULL_OFFSET;
        }

        let m = wavefronts
            .get_wavefront(score as usize)
            .unwrap()
            .mwavefront()
            .unwrap();

        if m.lo <= k && k <= m.hi {
            // eprintln!("\t[backtrace_mismatch_offset] score={} k={} offset={}", score, k, m.offsets[k as usize]);
            let k = compute_k(k, wavefronts.central_diagonal);
            return m.offsets[k] as isize + 1;
        } else {
            return NULL_OFFSET;
        }
    }
}

fn reduce(wavefronts: &mut Wavefronts, score: usize) {
    let qlen = wavefronts.query.len();
    let tlen = wavefronts.text.len();

    let a_k = wavefronts.central_diagonal;

    let wavefront: &mut WavefrontSet = wavefronts.get_wavefront_mut(score).unwrap();

    // fetch the m wavefront
    let m_wavefront: Option<&Wavefront> = wavefront.m.as_ref();
    if m_wavefront.is_none() {
        return;
    }
    let m_wavefront: &Wavefront = m_wavefront.unwrap();
    let hi = m_wavefront.hi;
    let lo = m_wavefront.lo;

    if (hi - lo + 1) < MIN_WAVEFRONT_LENGTH {
        return;
    }

    // Find minimum distance to (n,m)
    let mut min_distance: usize = max(qlen, tlen);

    for k in lo..hi {
        let i = k;
        let k: usize = compute_k(k, a_k);
        let s = format!("i={} k={}", i, k);

        let offset: Option<&usize> = m_wavefront.offsets.get(k);

        if offset.is_none() {
            continue
        }

        let offset: usize = *offset.unwrap();

        let left_v: usize = if offset >= k && qlen >= (offset - k) {
            qlen - (offset - k)
        } else {
            0
        };
        let left_h: usize = if tlen >= offset { tlen - offset } else { 0 };

        let distance = std::cmp::max(left_v, left_h);
        min_distance = std::cmp::max(min_distance, distance);
    }

    // Reduce WF from the bottom
    let i_wavefront: Option<&mut Wavefront> = wavefront.i.as_mut();
    if i_wavefront.is_some() {
        let i_wavefront_unwrapped: &mut Wavefront = i_wavefront.unwrap();
        if m_wavefront.lo > i_wavefront_unwrapped.lo {
            i_wavefront_unwrapped.lo = m_wavefront.lo
        }
        if m_wavefront.hi < i_wavefront_unwrapped.hi {
            i_wavefront_unwrapped.hi = m_wavefront.hi
        }
        if i_wavefront_unwrapped.lo > i_wavefront_unwrapped.hi {
            wavefront.i = None
        };
    }

    // Reduce WF from the top
    let d_wavefront: Option<&mut Wavefront> = wavefront.i.as_mut();
    if d_wavefront.is_some() {
        let d_wavefront_unwrapped: &mut Wavefront = d_wavefront.unwrap();
        if m_wavefront.lo > d_wavefront_unwrapped.lo {
            d_wavefront_unwrapped.lo = m_wavefront.lo
        }
        if m_wavefront.hi < d_wavefront_unwrapped.hi {
            d_wavefront_unwrapped.hi = m_wavefront.hi
        }
        if d_wavefront_unwrapped.lo > d_wavefront_unwrapped.hi {
            wavefront.d = None
        };
    }
}

fn wf_extend<T>(
    mwavefront: &mut Wavefront,
    match_lambda: T,
    central_diagonal: usize,
    score: usize,
    dp_matrix: &mut Vec<Vec<Option<usize>>>,
    verbosity: u8,
) where
    T: Fn(usize, usize, usize, usize, &mut Vec<Vec<Option<usize>>>) -> bool,
{
    let lo = mwavefront.lo;
    let hi = mwavefront.hi;

    if verbosity > 2 {
        eprintln!("[wfa::wf_extend] Extending wavefront with score {}", score);
        eprintln!("\tlo={}, hi={}", lo, hi);
    }

    for k in lo..=hi {
        let _lo = k;
        let k: usize = compute_k(k, central_diagonal);

        if k >= mwavefront.offsets.len() {
            if verbosity > 3 {
                eprintln!(
                    "[wfa::wf_extend] k={} is therefore out of scope. Skipping",
                    k
                );
            }
            continue;
        }
        let offset = mwavefront.offsets[k];
        let mut v: usize = v(offset, k, central_diagonal);
        let mut h: usize = h(offset, k, central_diagonal);

        if verbosity > 4 {
            eprintln!("\tk={} offset={}", k, offset);
            eprintln!("\tpre extend k={} offset={} ({},{})", k, offset, v, h);
        }

        while match_lambda(v, h, mwavefront.offsets[k], score, dp_matrix) {
            mwavefront.offsets[k] += 1;
            v += 1;
            h += 1;
        }

        if verbosity > 4 {
            eprintln!(
                "\tpost extend k={} offset={} ({},{})",
                k, mwavefront.offsets[k], v, h
            );
            eprintln!("");
        }
    }
}

fn wf_expand(wavefronts: &mut Wavefronts, score: usize) -> (isize, isize) {
    let s: isize = score as isize;

    let x: isize = wavefronts.penalties.mismatch as isize;
    let o: isize = wavefronts.penalties.gap_open as isize;
    let e: isize = wavefronts.penalties.gap_extend as isize;

    let m_hi = wavefronts
        .get_wavefront(score)
        .unwrap()
        .mwavefront()
        .unwrap()
        .hi;
    let i_hi = wavefronts
        .get_wavefront(score)
        .unwrap()
        .iwavefront()
        .unwrap()
        .hi;
    let d_hi = wavefronts
        .get_wavefront(score)
        .unwrap()
        .dwavefront()
        .unwrap()
        .hi;

    let m_lo = wavefronts
        .get_wavefront(score)
        .unwrap()
        .mwavefront()
        .unwrap()
        .lo;
    let d_lo = wavefronts
        .get_wavefront(score)
        .unwrap()
        .dwavefront()
        .unwrap()
        .lo;
    let i_lo = wavefronts
        .get_wavefront(score)
        .unwrap()
        .iwavefront()
        .unwrap()
        .lo;

    let hi = vec![
        if s - x < 0 {
            m_hi
        } else {
            wavefronts
                .get_wavefront((s - x) as usize)
                .unwrap()
                .mwavefront()
                .unwrap()
                .hi
        },
        if s - o - e < 0 {
            m_hi
        } else {
            wavefronts
                .get_wavefront((s - o - e) as usize)
                .unwrap()
                .mwavefront()
                .unwrap()
                .hi
        },
        if s - e < 0 {
            i_hi
        } else {
            wavefronts
                .get_wavefront((s - e) as usize)
                .unwrap()
                .iwavefront()
                .unwrap()
                .hi
        },
        if s - e < 0 {
            d_hi
        } else {
            wavefronts
                .get_wavefront((s - e) as usize)
                .unwrap()
                .dwavefront()
                .unwrap()
                .hi
        },
    ]
    .iter()
    .max()
    .unwrap()
        + 1;

    let lo = vec![
        if s - x < 0 {
            m_lo
        } else {
            wavefronts
                .get_wavefront((s - x) as usize)
                .unwrap()
                .mwavefront()
                .unwrap()
                .lo
        },
        if s - o - e < 0 {
            m_lo
        } else {
            wavefronts
                .get_wavefront((s - o - e) as usize)
                .unwrap()
                .mwavefront()
                .unwrap()
                .lo
        },
        if s - e < 0 {
            i_lo
        } else {
            wavefronts
                .get_wavefront((s - e) as usize)
                .unwrap()
                .iwavefront()
                .unwrap()
                .lo
        },
        if s - e < 0 {
            d_lo
        } else {
            wavefronts
                .get_wavefront((s - e) as usize)
                .unwrap()
                .dwavefront()
                .unwrap()
                .lo
        },
    ]
    .iter()
    .min()
    .unwrap()
        - 1;

    wavefronts
        .get_wavefront_mut(score)
        .unwrap()
        .mwavefront_mut()
        .unwrap()
        .hi = hi;
    wavefronts
        .get_wavefront_mut(score)
        .unwrap()
        .dwavefront_mut()
        .unwrap()
        .hi = hi;
    wavefronts
        .get_wavefront_mut(score)
        .unwrap()
        .iwavefront_mut()
        .unwrap()
        .hi = hi;

    wavefronts
        .get_wavefront_mut(score)
        .unwrap()
        .iwavefront_mut()
        .unwrap()
        .lo = lo;
    wavefronts
        .get_wavefront_mut(score)
        .unwrap()
        .mwavefront_mut()
        .unwrap()
        .lo = lo;
    wavefronts
        .get_wavefront_mut(score)
        .unwrap()
        .dwavefront_mut()
        .unwrap()
        .lo = lo;

    (lo, hi)
}

fn wf_next(wavefronts: &mut Wavefronts, score: usize, cli_args: &CliArgs) {
    let verbosity = cli_args.verbosity_level;
    if verbosity > 2 {
        eprintln!("[wf_next] Computing wavefront for score {}", score);
    }

    let s: isize = score as isize;

    let x: isize = wavefronts.penalties.mismatch as isize;
    let o: isize = wavefronts.penalties.gap_open as isize;
    let e: isize = wavefronts.penalties.gap_extend as isize;

    let num_wavefronts = wavefronts.wavefronts.len();
    if num_wavefronts <= score {
        (num_wavefronts..=score).for_each(|_| {
            let x = Box::new(Some(WavefrontSet::new(wavefronts.diagonals)));
            wavefronts.wavefronts.push(x);
        });
    }

    let (lo, hi) = wf_expand(wavefronts, score);

    if verbosity > 3 {
        eprintln!("\tk'\tk\tmmax\timax\tdmax");
    }

    for k in lo..=hi {
        if k < wavefronts.min_diagonal {
            if verbosity > 3 {
                eprintln!("[wf_next] k={} is therefore out of scope. Skipping", k);
            }
            continue;
        }

        let k: usize = compute_k(k, wavefronts.central_diagonal);

        let imax = *vec![
            if s - o - e < 0 || k == 0 || k + 1 >= wavefronts.diagonals {
                0
            } else {
                wavefronts
                    .get_wavefront((s - o - e) as usize)
                    .unwrap()
                    .mwavefront()
                    .unwrap()
                    .offsets[(k - 1) as usize] as isize
            },
            if s - e < 0 || k == 0 || k + 1 >= wavefronts.diagonals {
                0
            } else {
                wavefronts
                    .get_wavefront((s - e) as usize)
                    .unwrap()
                    .iwavefront()
                    .unwrap()
                    .offsets[(k - 1) as usize] as isize
            },
        ]
        .iter()
        .max()
        .unwrap();

        if k < wavefronts.diagonals {
            wavefronts
                .get_wavefront_mut(score)
                .unwrap()
                .iwavefront_mut()
                .unwrap()
                .offsets[k] = imax as usize;
        }

        let dmax = *vec![
            if s - o - e < 0 || k + 1 >= wavefronts.diagonals {
                0
            } else {
                wavefronts
                    .get_wavefront((s - o - e) as usize)
                    .unwrap()
                    .mwavefront()
                    .unwrap()
                    .offsets[k + 1] as isize
            },
            if s - e < 0 || k + 1 >= wavefronts.diagonals {
                0
            } else {
                wavefronts
                    .get_wavefront((s - e) as usize)
                    .unwrap()
                    .dwavefront()
                    .unwrap()
                    .offsets[k + 1] as isize
            },
        ]
        .iter()
        .max()
        .unwrap();

        if k < wavefronts.diagonals {
            wavefronts
                .get_wavefront_mut(score)
                .unwrap()
                .dwavefront_mut()
                .unwrap()
                .offsets[k] = dmax as usize;
        }

        let mmax = *vec![
            if s - x < 0 || k >= wavefronts.diagonals {
                0
            } else {
                (wavefronts
                    .get_wavefront((s - x) as usize)
                    .unwrap()
                    .mwavefront()
                    .unwrap()
                    .offsets[k]
                    + 1 as usize) as isize
            },
            if k >= wavefronts.diagonals {
                0
            } else {
                wavefronts
                    .get_wavefront(score)
                    .unwrap()
                    .iwavefront()
                    .unwrap()
                    .offsets[k] as isize
            },
            if k >= wavefronts.diagonals {
                0
            } else {
                wavefronts
                    .get_wavefront(score)
                    .unwrap()
                    .dwavefront()
                    .unwrap()
                    .offsets[k] as isize
            },
        ]
        .iter()
        .max()
        .unwrap();

        if k < wavefronts.diagonals {
            wavefronts
                .get_wavefront_mut(score)
                .unwrap()
                .mwavefront_mut()
                .unwrap()
                .offsets[k] = mmax as usize;
        }

        if verbosity > 3 {
            let k_prime = k as isize - wavefronts.central_diagonal as isize;
            eprintln!("\t{}\t{}\t{}\t{}\t{}", k_prime, k, mmax, imax, dmax);
        }
    }

    if cli_args.adapt {
        reduce(wavefronts, score);
    }
}

#[allow(unused_variables, unused_mut)]
fn backtrace<T>(
    wavefronts: &mut Wavefronts,
    score: usize,
    verbosity: u8,
    backtrace_lambda: &mut T,
) -> String
where
    T: FnMut((i32, i32), (i32, i32)),
{
    if verbosity > 1 {
        eprintln!("[wfa::backtrace]");
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
    let mut offset: isize = wf.mwavefront().unwrap().offsets[a_k] as isize;

    let mut k = compute_kk(a_k, a_k);

    let mut v: usize = compute_v(offset as usize, k as usize, a_k);
    let mut h: usize = compute_h(offset as usize, k as usize, a_k);

    let mut backtrace_type = Operation::MatchMismatch;

    let x: isize = wavefronts.penalties.mismatch as isize;
    let o: isize = wavefronts.penalties.gap_open as isize;
    let e: isize = wavefronts.penalties.gap_extend as isize;

    while v > 0 && h > 0 && score > 0 {
        if verbosity > 4 {
            eprintln!("\tbacktrace_type = {:?}", backtrace_type);
        }
        // compute scores
        let gap_open_score = score - o - e;
        let gap_extend_score = score - e;
        let mismatch_score = score - x;

        let del_ext: isize = if backtrace_type == Operation::Insertion {
            NULL_OFFSET
        } else {
            backtrace_deletion_extend_offset(&wavefronts, gap_extend_score, k)
        };
        let del_open: isize = if backtrace_type == Operation::Insertion {
            NULL_OFFSET
        } else {
            backtrace_deletion_open_offset(&wavefronts, gap_open_score, k)
        };
        let ins_ext: isize = if backtrace_type == Operation::Deletion {
            NULL_OFFSET
        } else {
            backtrace_insertion_extend_offset(&wavefronts, gap_extend_score, k)
        };
        let ins_open: isize = if backtrace_type == Operation::Deletion {
            NULL_OFFSET
        } else {
            backtrace_insertion_open_offset(&wavefronts, gap_open_score, k)
        };
        let misms: isize = if backtrace_type != Operation::MatchMismatch {
            NULL_OFFSET
        } else {
            backtrace_mismatch_offset(&wavefronts, mismatch_score, k)
        };

        // Compute maximum offset
        let max_all = *vec![del_ext, del_open, ins_ext, ins_open, misms]
            .iter()
            .max()
            .unwrap();

        if verbosity > 4 {
            eprintln!(
                "\tscore={} offset={} k={} \
                       gap_open_score={} gap_extend_score={} mismatch_score={} \
                       max_all={} del_ext={} del_open={} ins_ext={} ins_open={} mims={} \
                       backtrace_type={:?}",
                score,
                offset,
                k,
                gap_open_score,
                gap_extend_score,
                mismatch_score,
                max_all,
                del_ext,
                del_open,
                ins_ext,
                ins_open,
                misms,
                backtrace_type
            );
        }

        // Traceback Matches
        if backtrace_type == Operation::MatchMismatch && offset >= max_all {
            let num_matches = (offset - max_all) as usize;
            backtrace_matches_check(
                &mut offset,
                &mut cigar,
                num_matches,
                k,
                wavefronts.central_diagonal,
                backtrace_lambda,
            );
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
        let num_matches = offset as usize;
        backtrace_matches_check(
            &mut offset,
            &mut cigar,
            num_matches,
            k,
            wavefronts.central_diagonal,
            backtrace_lambda,
        );
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

// TODO: remove arg penalties
pub fn wf_align<T>(
    text: &[u8],
    query: &[u8],
    cli_args: &CliArgs,
    backtrace_lambda: &mut T,
) -> Alignment
where
    T: FnMut((i32, i32), (i32, i32)),
{
    let verbosity = cli_args.verbosity_level;
    let mut wavefronts = Wavefronts::new(query, text, cli_args.penalties);

    let qlen = query.len();
    let tlen = text.len();

    // central diagonal
    let a_k: usize = qlen - 1;

    // furthest offset along...
    let a_offset: usize = max(qlen, tlen);

    let mut score = 0;

    // Progress bar
    let bar = ProgressBar::new(a_offset as u64).with_prefix("\t");
    let template = "{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {bytes}/{total_bytes} ({eta_precise})";
    let progress_style = ProgressStyle::default_bar()
        .template(template)
        .progress_chars("=> ");
    bar.set_style(progress_style);
    let mut progress_value: u64 = 0;

    let mut exit_condition = |wavefronts: &Wavefronts, score: usize| {
        let current_offset = *wavefronts
            .get_wavefront(score)
            .unwrap()
            .mwavefront()
            .unwrap()
            .offsets
            .get(a_k)
            .unwrap();
        if verbosity > 1 {
            // handle progress bar
            let delta = current_offset as u64 - progress_value;
            bar.inc(delta);
            progress_value = current_offset as u64;
        }

        current_offset >= a_offset
    };

    let match_lambda =
        |v: usize, h: usize, offset: usize, score: usize, dp_matrix: &mut DpMatrix| -> bool {
            if h < tlen && v < qlen && dp_matrix[h][v] == None {
                dp_matrix[h][v] = Some(offset);
            }

            h < tlen && v < qlen && text[h] == query[v]
        };

    loop {
        let dp_matrix = &mut wavefronts.dp_matrix;
        let k = wavefronts
            .wavefronts
            .get_mut(score)
            .and_then(|x| (**x).as_mut());
        let m_s = k.unwrap().mwavefront_mut().unwrap();
        wf_extend(m_s, match_lambda, a_k, score, dp_matrix, verbosity);

        if exit_condition(&wavefronts, score) {
            if verbosity > 3 {
                eprintln!("Final state of the DP table");
                eprintln!("---------------------------");
                wavefronts.print();
            }

            if verbosity > 2 {
                eprintln!(
                    "\tscore: {}\n\
                           \tcentral diagonal (a_k): {}\n\
                           \tmaximum offset (a_offset): {}",
                    score, a_k, a_offset
                );
            }

            if cli_args.generate_alignment_tsv {
                let adapt_str: &str = if cli_args.adapt {"adapt"} else {"no_adapt"};
                let now: String = cli_args.start_time.format("%Y-%m-%d_%H-%M-%S").to_string();
                let filename = format!("wfilter-{}-{}.tsv", adapt_str, now);

                eprintln!("[wfa::align] Generating alignment tsv: {}", filename);

                wavefronts.print_tsv(&filename[..]);
            }

            let cigar = backtrace(&mut wavefronts, score, verbosity, backtrace_lambda);
            return Alignment { score, cigar };
        }

        if verbosity > 3 {
            wavefronts.print();
        }

        score += 1;

        wf_next(&mut wavefronts, score, cli_args);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::Penalties;
    use chrono::Local;

    static PENALTIES: Penalties = Penalties {
        mismatch: 4,
        matches: 0,
        gap_open: 6,
        gap_extend: 2,
    };

    static CLI: CliArgs = CliArgs {
        verbosity_level: 0,
        input_paf: String::new(),
        target_fasta: String::new(),
        query_fasta: String::new(),
        penalties: PENALTIES,
        adapt: false,
        generate_alignment_tsv: false,
        start_time: Local::now(),
    };

    static CLI_ADAPT: CliArgs = CliArgs {
        verbosity_level: 0,
        input_paf: String::new(),
        target_fasta: String::new(),
        query_fasta: String::new(),
        penalties: PENALTIES,
        adapt: true,
        generate_alignment_tsv: false,
        start_time: Local::now(),
    };

    fn mock_backtrace_lambda(_query: (i32, i32), _target: (i32, i32)) {}

    mod backtrace {
        use super::super::*;
        use super::*;

        #[test]
        fn test_same_sequence() {
            // same sequence
            let text = "GAGATA";
            let query = "GAGATA";

            let aln = wf_align(
                text.as_bytes(),
                query.as_bytes(),
                &CLI,
                &mut mock_backtrace_lambda,
            );
            assert_eq!(aln.score, 0);
            assert_eq!(aln.cigar, String::from("6M"));
        }

        #[test]
        fn test_snp() {
            // same sequence
            let text = "GACATA";
            let query = "GAGATA";

            let aln = wf_align(
                text.as_bytes(),
                query.as_bytes(),
                &CLI,
                &mut mock_backtrace_lambda,
            );
            assert_eq!(aln.score, 4);
            assert_eq!(aln.cigar, String::from("2M1X3M"));
        }

        #[test]
        fn test_paper_example() {
            let text = "GATACA";
            let query = "GAGATA";

            let aln = wf_align(
                text.as_bytes(),
                query.as_bytes(),
                &CLI,
                &mut mock_backtrace_lambda,
            );
            assert_eq!(aln.score, 8);
            assert_eq!(aln.cigar, String::from("2M1X1M1X1M"));
        }

        #[test]
        fn test_long_sequences() {
            let text = "TCTATACTGCGCGTTTGGAGAAATAAAATAGTTCTATACTGCGCGTTTGGAGAA\
                         ATAAAATAGTTCTATACTGCGCGTTTGGAGAAATAAAATAGTTCTATACTGCGC\
                         GTTTGGAGAAATAAAATAGT";
            let query = "TCTTTACTCGCGCGTTGGAGAAATACAATAGTTCTTTACTCGCGCGTTGGAGAA\
                         ATACAATAGTTCTTTACTCGCGCGTTGGAGAAATACAATAGTTCTTTACTCGCG\
                         CGTTGGAGAAATACAATAGT";

            let aln = wf_align(
                text.as_bytes(),
                query.as_bytes(),
                &CLI,
                &mut mock_backtrace_lambda,
            );
            assert_eq!(aln.score, 96);
            assert_eq!(
                aln.cigar,
                String::from(
                    "3M1X4M1I6M1D10M1X9M1X4M1I6M1D10M1X9M1X4M1I6M1D10M1X9M1X4M1I6M1D10M1X6M"
                )
            );
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
        use super::*;

        #[test]
        fn test_same_sequence() {
            // same sequence
            let text = "GAGATA";
            let query = "GAGATA";

            let aln = wf_align(
                text.as_bytes(),
                query.as_bytes(),
                &CLI,
                &mut mock_backtrace_lambda,
            );
            assert_eq!(aln.score, 0);
        }

        #[test]
        fn test_paper_test_case() {
            let text = "GAT";
            let query = "GAG";
            let aln = wf_align(
                text.as_bytes(),
                query.as_bytes(),
                &CLI,
                &mut mock_backtrace_lambda,
            );
            assert_eq!(aln.score, 4);

            let text = "GATACA";
            let query = "GAGATA";
            let aln = wf_align(
                text.as_bytes(),
                query.as_bytes(),
                &CLI,
                &mut mock_backtrace_lambda,
            );
            assert_eq!(aln.score, 8);
        }

        #[test]
        fn test_c_example() {
            let query = "TCTTTACTCGCGCGTTGGAGAAATACAATAGT";
            let text = "TCTATACTGCGCGTTTGGAGAAATAAAATAGT";

            let aln = wf_align(
                &text.as_bytes()[..10],
                &query.as_bytes()[..10],
                &CLI,
                &mut mock_backtrace_lambda,
            );
            assert_eq!(aln.score, 12);

            let aln = wf_align(
                text.as_bytes(),
                query.as_bytes(),
                &CLI_ADAPT,
                &mut mock_backtrace_lambda,
            );
            assert_eq!(aln.score, 24);
        }

        #[test]
        fn test_same_length() {
            let text = "TCTATACTGCGCGTTTGGAGAAATAAAATAGTTCTATACTGCGCGTTTGGAGAAATAAA\
                         ATAGTTCTATACTGCGCGTTTGGAGAAATAAAATAGTTCTATACTGCGCGTTTGGAGAA\
                         ATAAAATAGT";
            let query = "TCTTTACTCGCGCGTTGGAGAAATACAATAGTTCTTTACTCGCGCGTTGGAGAAATACA\
                         ATAGTTCTTTACTCGCGCGTTGGAGAAATACAATAGTTCTTTACTCGCGCGTTGGAGAA\
                         ATACAATAGT";

            let aln = wf_align(
                text.as_bytes(),
                query.as_bytes(),
                &CLI_ADAPT,
                &mut mock_backtrace_lambda,
            );
            assert_eq!(aln.score, 96);
        }

        #[test]
        fn test_longer_query() {
            let text = "TCTTTACTCGCGCGTTGGAGAAATACAATAGTTCTTTACTCGCGCGTTGGAGAAATACAATAGT\
                         TCTTTACTCGCGCGTTGGAGAAATACAATAGTTCTTTACTCGCGCGTTGGAGAAATACAATAGT";
            let query = "TCTATACTGCGCGTTTATCTAGGAGAAATAAAATAGTTCTATACTGCGCGTTTGGAGAAATAAAATAGT\
                         TCTATACTGCGCGTTTGGAGAAATAACTATCAATAGTTCTATACTGCGCGTTTGGAGAAATAAAATAGT";

            let aln = wf_align(
                text.as_bytes(),
                query.as_bytes(),
                &CLI_ADAPT,
                &mut mock_backtrace_lambda,
            );
            assert_eq!(aln.score, 200);
        }

        #[test]
        fn test_longer_text() {
            let text = "TCTATACTGCGCGTTTATCTAGGAGAAATAAAATAGTTCTATACTGCGCGTTTGGAGAAATAAAATAGT\
                         TCTATACTGCGCGTTTGGAGAAATAACTATCAATAGTTCTATACTGCGCGTTTGGAGAAATAAAATAGT";
            let query = "TCTTTACTCGCGCGTTGGAGAAATACAATAGTTCTTTACTCGCGCGTTGGAGAAATACAATAGT\
                         TCTTTACTCGCGCGTTGGAGAAATACAATAGTTCTTTACTCGCGCGTTGGAGAAATACAATAGT";

            let aln = wf_align(
                text.as_bytes(),
                query.as_bytes(),
                &CLI_ADAPT,
                &mut mock_backtrace_lambda,
            );

            assert_eq!(aln.score, 200);
        }
    }
}
