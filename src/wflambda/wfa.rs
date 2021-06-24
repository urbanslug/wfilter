#![allow(dead_code, unused_macros)]
/*
A diagonal is a number around the 0/center diagnal. Could be greater or less than 0
Diagonal number: h - v
v = offset - diagonal
h = offset
k is the diagnal
 */

use std::cmp;


const MATCH_SCORE: isize = 0;
const MIS_MATCH_SCORE: isize = 4;
const GAP_OPENING_SCORE: isize = 6;
const GAP_EXTENSION_SCORE: isize = 2;


/*
k is the name of the diagonal
We name diagonals based on the row - column
Example: a diagonal passing through 12, 14 is diagonal -2
*/
macro_rules! k {
    ($h: expr, $v: expr) => {{
        $h - $v
    }};
}

macro_rules! v {
    ($k: expr, $offset: expr) => {{
        $offset - $k
    }}
}

macro_rules! h {
    ($k: expr, $offset: expr) => {{
        $offset
    }}
}

// will implicitly change type
macro_rules! negate {
    ($x: expr) => {{
        -($x as isize)
    }};
}

macro_rules! abs {
    ($x: expr) => {{
        if $x >= 0 {
            $x
        } else {
            -$x
        }
    }};
}

#[derive(Debug)]
struct Alignment {
    t_len: usize,
    q_len: usize,
    alignment_score: i16,
    num_wavefronts: usize,
    a_k: usize, // starting diagonal
}



type Offset = i16;
#[derive(Debug)]
struct WaveFront {
    low: isize,  // lowest diagonal that this wavfront touches
    high: isize, // top diagonal that this wavfront touches

    /*
    offsets is a vector of diagonals with the index being the diagonal and the
    value being the offset i.e how far we've moved along that diagonal

    Below we have offsets for 7 diagonals along which we've moved 6, 4, 2 ... 0 steps

      0   1   2   3   4   5   6
     ----------------------------
    | 6 | 4 | 2 | 2 | 0 | 0 | 0 |
     ----------------------------
     */
    offsets: Vec<Offset>, // matches & mismatches
}

impl WaveFront {
    // create a new set of wavefronts
    fn new(max_offset: usize) -> Self {
        Self {
            low: 0,
            high: 0,
            // TODO: make this more efficient.
            // It is wasteful as not all reach here but helps the exit condition
            offsets: vec![0; max_offset],
        }
    }

    // check in the offsets
    // never fail. Could be dangerous
    fn at_offset(&self, o: usize) -> Offset {
        if o < self.offsets.len() {
            *self.offsets.get(o).unwrap()
        } else {
            0
        }
    }
}

// Wavefronts with score s (WF_s)
#[derive(Debug)]
struct WFs {
    iwavefront: WaveFront,
    dwavefront: WaveFront,
    mwavefront: WaveFront,
}

impl WFs {
    // Create the initial WFs set at the origin
    fn new(max_offset: usize) -> Self {
        Self {
            iwavefront: WaveFront::new(max_offset),
            dwavefront: WaveFront::new(max_offset),
            mwavefront: WaveFront::new(max_offset),
        }
    }
}

// wavefront at index 0 is the wavefront for score 0 i.e WF_0
struct WaveFronts {
    wavefronts: Vec<WFs>,
    diagonals: usize,
    max_offset: usize,
    lowest_diagonal: isize,
}

impl WaveFronts {
    // create a new set of wavefronts
    fn new(_max_offset: usize, text_len: usize, query_len: usize) -> Self {

        let m = (text_len+query_len)*4;
        let mut wfs: Vec<WFs> = Vec::with_capacity(query_len+text_len);
        (0..=text_len+query_len).for_each(|_| {
            wfs.push(WFs::new(m))
        });

        Self {
            // create it with at least the zero score wavefront
            wavefronts:  wfs,
            diagonals: text_len+query_len,
            max_offset: m, // just trying to set a large number
            lowest_diagonal: -(query_len as isize * 2),
        }
    }

    fn diagnoal_index(&self, diagonal: isize) -> usize {
        if diagonal >= self.lowest_diagonal {
            let index = diagonal + abs!(self.lowest_diagonal);
            /*
            println!("offset of diagonal {}, with lowest diagonal {}, is {}",
                     diagonal, self.lowest_diagonal, index);
             */
            index as usize
        } else {
            panic!("Diagonal out of range")
        }
    }

    fn get_wavefront_mut<'a>(&'a mut self, score: usize) -> Option<&'a mut WFs> {
        self.wavefronts.get_mut(score)
    }

    fn get_wavefront(&self, score: usize) -> Option<&WFs> {
        self.wavefronts.get(score)
    }

    #[allow(dead_code)]
    fn print_matrix(&self) {
        ()
    }
}



/*
Take m-wavefront for a score s
the score s is determined in wf_align
 */
fn wf_extend<T>(mwavefront: &mut WaveFront, match_lambda: T)
where
    T: Fn(usize, usize) -> bool,
{
    // k is the diagonal
    // we are going from the lowest to highest diagonal
    mwavefront.offsets.iter_mut().enumerate().for_each(|(diagonal, offset): (usize, &mut i16)| {
        // we expect diagonal to run from 0...max_diagonals
        let mut v: usize = v!(diagonal as i16, *offset) as usize;
        let mut h: usize = h!(diagonal as i16, *offset) as usize;

        while match_lambda(v, h) {
            // extend the wavefront along this diagonal
            // TODO: too expensive?
            /*
            match mwavefront.offsets.get_mut(diagonal) {
                Some(o) => { *o += 1 }
                _ => panic!("[wfa::wf_extend] ...")
            }
             */
            *offset += 1;
            //mwavefront.offsets[diagonal] = mwavefront.offsets[diagonal] + 1;
            v += 1;
            h += 1
        }
    });
}

/*
Compute the next wavefront for the score
A wavfront will touch a set of diagonals
from hi to low and be extended along them
 */
fn wf_next(
    wavefronts: &mut WaveFronts,
    _query: &str,
    _text: &str,
    score: usize,
    _max_offset: usize,
) {
    eprintln!("[wfa::wf_next]");
    if score <= 0 {
        panic!("[wfa::wf_next] score={}, is too low.", score);
    }

    eprint!("\tCreating wavefront for score {}.\n", score);
    eprint!("\t\tPrevious mlo {}, mhi {}",
            wavefronts.get_wavefront(score-1).unwrap().mwavefront.low,
            wavefronts.get_wavefront(score-1).unwrap().mwavefront.high);

    // If a wavefront for that score doesn't exit, create it.
    match wavefronts.get_wavefront_mut(score) {
        None => {
            let c = wavefronts.wavefronts.len();
            (c..=score).for_each(|_| { wavefronts.wavefronts.push(WFs::new(wavefronts.max_offset)) })
        },
        _ => {},
    }

    let x = MIS_MATCH_SCORE as usize;
    let o = GAP_OPENING_SCORE as usize;
    let e = GAP_EXTENSION_SCORE as usize;
    let s = score;

    let ix = x as isize;
    let iscore = score as isize;
    let io = o as isize;
    let ie = e as isize;

    let hi: isize = vec![
        if iscore - ix < 0 { 0 } else { wavefronts.get_wavefront(score - x).unwrap().mwavefront.high },
        if iscore - io - ie < 0 { 0 } else {wavefronts.get_wavefront(score - o -e).unwrap().mwavefront.high},
        if iscore - ie < 0 { 0 } else {wavefronts.get_wavefront(score - e).unwrap().iwavefront.high},
        if iscore - ie < 0 { 0 } else {wavefronts.get_wavefront(score - e).unwrap().dwavefront.high},
    ].iter().max().unwrap() + 1;

    if hi <= 0 { eprintln!("[wfa::wf_next] Computed a hi of {}. This indicates serious a problem", hi)}

    let lo: isize = vec![
        if iscore - ix < 0 { 0 } else { wavefronts.get_wavefront(score - x).unwrap().mwavefront.low },
        if iscore - io - ie < 0 { 0 } else {wavefronts.get_wavefront(score - o -e).unwrap().mwavefront.low},
        if iscore - ie < 0 { 0 } else {wavefronts.get_wavefront(score - e).unwrap().iwavefront.low},
        if iscore - ie < 0 { 0 } else {wavefronts.get_wavefront(score - e).unwrap().dwavefront.low},
    ].iter().min().unwrap() - 1;

    if lo >= 0 { eprintln!("[wfa::wf_next] Computed a lo of {}. This indicates serious a problem", lo)}

    wavefronts.wavefronts[score].iwavefront.low = lo;
    wavefronts.wavefronts[score].iwavefront.high = hi;

    wavefronts.wavefronts[score].dwavefront.low = lo;
    wavefronts.wavefronts[score].dwavefront.high = hi;

    wavefronts.wavefronts[score].mwavefront.low = lo;
    wavefronts.wavefronts[score].mwavefront.high = hi;

    eprint!("\t\tCurrent lo {}, hi {}\n", lo, hi);
    eprint!("\t\tCurrent wavefronts.max_offset {},  dlen {}, mlen {}, ilen {}\n",
            wavefronts.max_offset,
            wavefronts.get_wavefront(score).unwrap().mwavefront.offsets.len(),
            wavefronts.get_wavefront(score).unwrap().dwavefront.offsets.len(),
            wavefronts.get_wavefront(score).unwrap().iwavefront.offsets.len()
    );

    // k is a diagnal
    (lo..=hi).for_each(|diagonal: isize| {
        let k: usize = wavefronts.diagnoal_index(diagonal);

        eprint!("\t\tk={}, diagonal={}\n", k , diagonal);

        let imax: i16 = vec![
            if iscore - io - ie < 0 || k == 0 { 0 } else { wavefronts.get_wavefront(s-o-e).unwrap().mwavefront.offsets[k-1] },
            if iscore - ie < 0 || k == 0 { 0 } else { wavefronts.get_wavefront(s-e).unwrap().iwavefront.offsets[k-1] }
        ].iter().max().unwrap() + 1;

        let dmax: i16 = *vec![
            if iscore - io - ie < 0 { 0 } else { wavefronts.get_wavefront(s-o-e).unwrap().mwavefront.offsets[k+1] },
            if iscore - ie < 0 { 0 } else { wavefronts.get_wavefront(s-e).unwrap().dwavefront.offsets[k+1] }
        ].iter().max().unwrap();

        let mmax: i16 = *vec![
            if iscore - ix < 0 { 0 } else { wavefronts.get_wavefront(s-x).unwrap().mwavefront.offsets[k] + 1 },
            wavefronts.get_wavefront(s).unwrap().iwavefront.offsets[k],
            wavefronts.get_wavefront(s).unwrap().dwavefront.offsets[k]
        ].iter().max().unwrap();

        wavefronts.wavefronts[score].iwavefront.offsets[k] = imax;
        wavefronts.wavefronts[score].dwavefront.offsets[k] = dmax;
        wavefronts.wavefronts[score].mwavefront.offsets[k] = mmax;
    });
}

fn wf_align(query: &str, text: &str) -> Alignment {
    // TODO: is this the best idea? it lets us easily access the wavefronts
    // assume m >= n
    let query_len = query.len(); // n
    let text_len = text.len(); // m

    // starting diagonal
    // this is the diagonal that reaches point n,m
    let a_k: usize = abs!(text_len as i32 - query_len as i32) as usize;

    // max offset on diagonal a_k
    // it is also the furthest point (offset) from the starting diagonal (a_k) in the DP matrix
    // this offset reaches point n,m
    let a_offset = text_len;

    let offset_count  = query_len+text_len+2;

    // start conditions
    // not actually a set
    // set the score at the origin to 0
    // we start with the match wavefront
    let start_score: usize = 0;

    // this is the origin and assumes query_len == text_len
    // at offset 0
    let start_offset = 0;

    // initialize wavefronts
    // pass max offsets to set aside capacity for all possible offsets
    let mut wavefronts = WaveFronts::new(offset_count, text_len, query_len);
    // println!("{:?}", wavefronts.wavefronts); TODO: remove
    // initial wavefront set
    let initial_wavefront_set: &mut WFs = match wavefronts.get_wavefront_mut(start_score) {
        Some(wf) => wf,
        None => panic!("Initial wavefront not set"),
    };
    initial_wavefront_set
        .mwavefront
        .offsets
        .insert(start_offset, 0);

    let mut score = start_score;

    let match_lambda =
        |v: usize, h: usize| {
            v < query_len && h < text_len && (query.chars().nth(v).unwrap() == text.chars().nth(h).unwrap())
        };

    eprint!("[wfa::wfalign]\n");
    eprint!("\tStats:\n");
    eprint!("\t------\n");
    eprint!("\ttext_len(m) = {}, query_len(n) = {}, diagonals={}, max_offset={}, lowest_diagonal={}\n\n",
            query_len, text_len,
            wavefronts.diagonals, wavefronts.max_offset, wavefronts.lowest_diagonal);

    loop {
        let mwavefront = &mut wavefronts.get_wavefront_mut(score).unwrap().mwavefront;
        wf_extend(
            mwavefront,
            match_lambda,
        );

        if score > 50 { break; }

        // if wavefront of the current score is at offset
        // if we have reaches a_k
        if wavefronts
            .get_wavefront(score)
            .unwrap()
            .mwavefront
            .at_offset(a_k)
            >= a_offset as i16
        {
            break;
        }

        // compute the wavefront for the next score
        score += 1;
        wf_next(&mut wavefronts, query, text, score, offset_count);
    }


    let aln = Alignment {
        num_wavefronts: score+1,
        alignment_score: wavefronts.get_wavefront(score).unwrap().mwavefront.at_offset(a_k),
        q_len: query_len,
        t_len: text_len,
        a_k: a_k
    };

    eprintln!("{:#?}", aln);
    eprintln!("\n\n.................................................................................................................\n\n\n");

    aln
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_align_equal_length() {
        // same sequence
        let text  = "GAGATA";
        let query = "GAGATA";
        let aln = wf_align(query, text);
        assert_eq!(aln.alignment_score, 6);
        assert_eq!(aln.num_wavefronts, 1);

        // different sequences
        let text  = "GAGAAT";
        let query = "GAGATA";
        let aln = wf_align(query, text);
        assert_eq!(aln.alignment_score, 6);

    }

    #[test]
    fn test_align_longer_text() {

        // different sequences
        let text  = "GAGGTACAAT";
        let query = "GAGATA";
        let aln = wf_align(query, text);
        assert_eq!(aln.alignment_score, 6);
    }

    #[test]
    fn test_align_longer_query() {

        // different sequences
        let text  = "GAGAAT";
        let query = "GAGATAGATA";
        let aln = wf_align(query, text);
        assert_eq!(aln.alignment_score, 6);
    }

    #[test]
    fn test_wf_diagonal_index() {
        let text = "GATACA";
        let query = "GAGATA";

        let wf = WaveFronts::new(text.len(), text.len(), query.len());

        assert_eq!(wf.lowest_diagonal, -(query.len() as isize));
        assert_eq!(wf.diagnoal_index(text.len() as isize), text.len()+query.len());

        let negative_query_len = negate!(query.len());
        assert_eq!(wf.diagnoal_index(negative_query_len), 0);
    }

    #[test]
    fn test_marcos() {
        assert_eq!(k!(14, 12), 2);
        assert_eq!(abs!(-14), 14);
        assert_eq!(abs!(5), 5);
    }
}
