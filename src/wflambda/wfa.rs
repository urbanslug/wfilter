use std::cmp::{max};
use self::utils::*;
use self::types::*;

const VERBOSITY_LEVEL: usize = 1;

#[macro_use]
mod macros {

}

mod utils {

    use super::types::*;

    macro_rules! abs {
        ($x: expr) => {{
            if $x >= 0 {
                $x
            } else {
                -$x
            }
        }};
    }

    pub fn v(offset: usize, diagonal: usize, central_diagonal: usize) -> usize {
        offset + abs!(central_diagonal as isize - diagonal as isize) as usize
    }

    pub fn h(offset: usize) -> usize {
        offset
    }


    #[allow(dead_code)]
    fn print_matrix(wavefronts: &Wavefronts, score: usize, central_diagonal: usize) {
        let query = wavefronts.query;
        let text = wavefronts.text;
        let qlen = query.len();
        let tlen = text.len();

        let mut dp_matrix: Vec<Vec<Option<usize>>> = vec![vec![None; qlen]; tlen];

        let offsets: &Vec<Offset> = &wavefronts.get_wavefront(score).unwrap().mwavefront.offsets;

        offsets.iter().enumerate().for_each(|(diagonal, offset): (usize, &usize)| {
            let offset = *offset;
            let v: usize = v(offset, diagonal, central_diagonal);
            let h: usize = h(offset);

            if offset > 0 {
                dp_matrix[h][v] = Some(score);
            }
        });

        // print col nums
        eprint!("\t");
        dp_matrix[0].iter().enumerate().for_each(|(j, _)| {
            eprint!("{}\t", query[j] as char);
        });
        eprint!("\n");
        eprint!("\t");
        dp_matrix[0].iter().enumerate().for_each(|(j, _)| {
            eprint!("{}\t", j);
        });
        eprint!("\n");
        dp_matrix.iter().enumerate().for_each(|(i, row)| {
            eprint!("{} {}\t", text[i] as char, i);
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

        //std::io::stdout().flush().unwrap();
    }
}


mod types {
    pub struct Alignment {
        pub score: usize,
    }

    pub type Offset = usize;

    pub struct Wavefront {
        // TODO: make usize
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

        central_diagonal: usize, // k_zero
        max_diagonal: usize,


        pub dp_matrix: Vec<Vec<Option<usize>>>,
    }

    impl<'a> Wavefronts<'a> {
        pub fn new(query: &'a[u8], text: &'a[u8]) -> Self {

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
                central_diagonal: qlen - 1,
                max_diagonal: diagonals - 2,
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
    }
}


fn wf_extend<T>(mwavefront: &mut Wavefront,
                match_lambda: T,
                central_diagonal: usize,
                score: usize,
                dp_matrix: &mut Vec<Vec<Option<usize>>>)
where
    T: Fn(usize, usize) -> bool,
{
    for k in mwavefront.lo..=mwavefront.hi {
        let k = (k + central_diagonal as isize) as usize;
        if  k >= mwavefront.offsets.len() {
            if VERBOSITY_LEVEL > 3 {
                eprintln!("[wf_extend] k={} is therefore out of scope. Skipping", k);
            }
            continue;
        }
        let mut v: usize = v(mwavefront.offsets[k], k, central_diagonal);
        let mut h: usize = h(mwavefront.offsets[k]);

        while match_lambda(v,h) {
            // println!("k={} ({}, {})", k , v, h);
            dp_matrix[h][v] = Some(score);
            mwavefront.offsets[k] += 1;
            v += 1;
            h += 1;
        }
    }
}


fn wf_next(wavefronts: &mut Wavefronts, score: usize) {

    let x: isize = 4; let o: isize = 6; let e: isize = 1; let s: isize = score as isize;

    let num_wavefronts = wavefronts.wavefronts.len();
    if num_wavefronts <= score {
        (num_wavefronts..=score).for_each(|_| {
            wavefronts.wavefronts.push(WavefrontSet::new(wavefronts.diagonals));
        });
    }

    let m_hi =  wavefronts.get_wavefront(score).unwrap().mwavefront.hi;
    let m_lo =  wavefronts.get_wavefront(score).unwrap().mwavefront.lo;
    let d_hi =  wavefronts.get_wavefront(score).unwrap().dwavefront.hi;
    let d_lo =  wavefronts.get_wavefront(score).unwrap().dwavefront.lo;
    let i_hi =  wavefronts.get_wavefront(score).unwrap().iwavefront.hi;
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
    ].iter().max().unwrap() - 1;

    // set new lo
    wavefronts.get_wavefront_mut(score).unwrap().mwavefront.hi = hi;
    wavefronts.get_wavefront_mut(score).unwrap().mwavefront.lo = lo;
    wavefronts.get_wavefront_mut(score).unwrap().iwavefront.hi = hi;
    wavefronts.get_wavefront_mut(score).unwrap().iwavefront.lo = lo;
    wavefronts.get_wavefront_mut(score).unwrap().dwavefront.hi = hi;
    wavefronts.get_wavefront_mut(score).unwrap().dwavefront.lo = lo;

    for k in lo..=hi {
        if k < 0 {
            if VERBOSITY_LEVEL > 3 {
                eprintln!("[wf_next] k={} is therefore out of scope. Skipping", k);
            }
            continue;
        }

        let k = k as usize;

        let imax = *vec![
            if s-o-e < 0 || k <= 0  || k+1 >= wavefronts.diagonals { 0 } else { wavefronts.get_wavefront((s-o-e) as usize).unwrap().mwavefront.offsets[(k-1) as usize] },
            if s-e < 0   || k <= 0 || k+1 >= wavefronts.diagonals { 0 } else { wavefronts.get_wavefront((s-e) as usize).unwrap().iwavefront.offsets[(k-1) as usize] },
        ].iter().max().unwrap() + 1;

        let dmax = *vec![
            if s-o-e < 0 || k+1 >= wavefronts.diagonals { 0 } else { wavefronts.get_wavefront((s-o-e) as usize).unwrap().mwavefront.offsets[k+1] },
            if s-e < 0 || k+1 >= wavefronts.diagonals { 0 } else { wavefronts.get_wavefront((s-e) as usize).unwrap().dwavefront.offsets[k+1] },
        ].iter().max().unwrap();

        let mmax = *vec![
            if s-x < 0 || k >= wavefronts.diagonals { 0 } else { wavefronts.get_wavefront((s-x) as usize).unwrap().mwavefront.offsets[k] + 1 },
            if  k >= wavefronts.diagonals { 0 } else { wavefronts.get_wavefront(score).unwrap().iwavefront.offsets[k] },
            if  k >= wavefronts.diagonals { 0 } else { wavefronts.get_wavefront(score).unwrap().dwavefront.offsets[k] },
        ].iter().max().unwrap();

        if k >= wavefronts.diagonals { continue; }
        wavefronts.get_wavefront_mut(score).unwrap().iwavefront.offsets[k] = imax;
        wavefronts.get_wavefront_mut(score).unwrap().dwavefront.offsets[k] = dmax;
        wavefronts.get_wavefront_mut(score).unwrap().mwavefront.offsets[k] = mmax;
    }

}

fn wf_align(text: &[u8], query: &[u8]) -> Alignment {
    let mut wavefronts = Wavefronts::new(query, text);

    let qlen = query.len();
    let tlen = text.len();

    // central diagonal
    let a_k: usize = qlen - 1;

    // furthest offset along...
    let a_offset: usize = max(qlen, tlen);

    let mut score = 0;

    // let exit_condition = || { *wavefronts.get_wavefront(score).unwrap().mwavefront.offsets.get(a_k).unwrap() >= a_offset };

    let match_lambda = |v: usize, h: usize| {
        h < tlen && v < qlen &&  text[h] == query[v]
    };

    loop {

        let dp_matrix = &mut wavefronts.dp_matrix;
        let m_s = &mut wavefronts.wavefronts.get_mut(score).unwrap().mwavefront;
        wf_extend(m_s, match_lambda, a_k, score, dp_matrix);

        if wavefronts.get_wavefront(score).unwrap().mwavefront.offsets[a_k] >= a_offset {
            let score = if score > 0 {score - 1} else { score };

            if VERBOSITY_LEVEL > 0 {
                eprintln!("score: {}, central_diagonal(a_k): {}, max_offset(a_offset): {} ",
                          score,
                          a_k,
                          a_offset);
            }

            return Alignment {score: score};
        }

        if VERBOSITY_LEVEL > 3 {
            wavefronts.print();
        }

        score += 1;

        wf_next(&mut wavefronts, score);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_same_sequence() {
        // same sequence
        let text  = "GAGATA";
        let query = "GAGATA";

        let aln = wf_align(text.as_bytes(), query.as_bytes());
        assert_eq!(aln.score, 0);
    }


    #[test]
    fn test_paper_test_case() {
        let text  = "GATACA";
        let query = "GAGATA";
        let aln = wf_align(text.as_bytes(), query.as_bytes());
        assert_eq!(aln.score, 8);
    }

    #[test]
    fn test_same_length() {
        let text  = "TCTTTACTCGCGCGTTGGAGAAATACAATAGTTCTTTACTCGCGCGTTGGAGAAATACAATAGT\
                     TCTTTACTCGCGCGTTGGAGAAATACAATAGTTCTTTACTCGCGCGTTGGAGAAATACAATAGT";
        let query = "TCTATACTGCGCGTTTGGAGAAATAAAATAGTTCTATACTGCGCGTTTGGAGAAATAAAATAGT\
                     TCTATACTGCGCGTTTGGAGAAATAAAATAGTTCTATACTGCGCGTTTGGAGAAATAAAATAGT";
        let aln = wf_align(text.as_bytes(), query.as_bytes());
        assert_eq!(aln.score, 250);
    }

    #[test]
    fn test_longer_query() {
        let text  = "TCTTTACTCGCGCGTTGGAGAAATACAATAGTTCTTTACTCGCGCGTTGGAGAAATACAATAGT\
                     TCTTTACTCGCGCGTTGGAGAAATACAATAGTTCTTTACTCGCGCGTTGGAGAAATACAATAGT";
        let query = "TCTATACTGCGCGTTTATCTAGGAGAAATAAAATAGTTCTATACTGCGCGTTTGGAGAAATAAAATAGT\
                     TCTATACTGCGCGTTTGGAGAAATAACTATCAATAGTTCTATACTGCGCGTTTGGAGAAATAAAATAGT";

        let aln = wf_align(text.as_bytes(), query.as_bytes());
        assert_eq!(aln.score, 536);
    }

    #[test]
    fn test_longer_text() {
        let text  = "TCTATACTGCGCGTTTATCTAGGAGAAATAAAATAGTTCTATACTGCGCGTTTGGAGAAATAAAATAGT\
                     TCTATACTGCGCGTTTGGAGAAATAACTATCAATAGTTCTATACTGCGCGTTTGGAGAAATAAAATAGT";
        let query = "TCTTTACTCGCGCGTTGGAGAAATACAATAGTTCTTTACTCGCGCGTTGGAGAAATACAATAGT\
                     TCTTTACTCGCGCGTTGGAGAAATACAATAGTTCTTTACTCGCGCGTTGGAGAAATACAATAGT";

        let aln = wf_align(text.as_bytes(), query.as_bytes());
        assert_eq!(aln.score, 526);
    }
}
