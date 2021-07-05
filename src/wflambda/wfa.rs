use std::cmp::{max};
use self::utils::*;
use self::types::*;

const VERBOSITY_LEVEL: usize = 0;

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
}

mod utils {
    use super::types::*;

    pub fn compute_k(k: isize, central_diagonal: usize) -> usize {
        (k + central_diagonal as isize) as usize
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
}

mod types {
    pub type Offset = usize;

    pub struct Alignment {
        pub score: usize,
    }

    #[derive(Copy, Clone)]
    pub struct Penalties {
        pub mismatch: usize,
        pub matches: usize,
        pub gap_open: usize,
        pub gap_extend: usize,
    }

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

            Self {
                query,
                text,
                wavefronts,
                diagonals,
                penalties,
                central_diagonal: qlen - 1,
                // max_diagonal: diagonals - 2,
                min_diagonal: -(qlen as isize -1),
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
    T: Fn(usize, usize, usize, &mut Vec<Vec<Option<usize>>>) -> bool,
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

        while match_lambda(v, h, score, dp_matrix) {
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

    let match_lambda = |v: usize, h: usize, score: usize, dp_matrix:  &mut Vec<Vec<Option<usize>>>| -> bool {
        if h < tlen && v < qlen && dp_matrix[h][v] == None {
            dp_matrix[h][v] = Some(score);
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
            return Alignment {score: score};
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

        let p =  Penalties {
            mismatch: 7,
            matches: 0,
            gap_open: 11,
            gap_extend: 1,
        };

        let aln = wf_align(text.as_bytes(), query.as_bytes(), p);
        assert_eq!(aln.score, 16);
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
