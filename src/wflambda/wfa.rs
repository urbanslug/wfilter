#![allow(dead_code, unused_macros)]

const MATCH_SCORE: isize = 0;
const MIS_MATCH_SCORE: isize = 4;
const GAP_OPENING_SCORE: isize = 6;
const GAP_EXTENSION_SCORE: isize = 2;

const VERBOSITY_LEVEL: usize = 4;

/*
A diagonal is a named/numbered based on the row - column
All diagonals are located relative to a central diagonal m - n.
where m is the text length and n is the query length
if less the central diagonal is less than zero, we flip its sign

Diagonal number: h - v
v = offset - diagonal
h = offset

k = h - v (k is the diagonal)

Example: a diagonal passing through (v, h) (12, 14) is diagonal 2
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

type Offset = i16;

#[derive(Debug)]
struct Alignment {
    t_len: usize,
    q_len: usize,
    alignment_score: Offset,
    num_wavefronts: usize,
    central_diagonal: usize, // starting diagonal
}


#[derive(Debug)]
struct WaveFront {
    low: isize,  // lowest diagonal that this wavfront touches
    high: isize, // top diagonal that this wavfront touches

    /*
    offsets is a vector of diagonals with the index being the diagonal and the
    value being the offset i.e how far we've moved along that diagonal

    Below we have offsets for 7 diagonals along which we've moved 6, 4, 2 ... 0 steps

      0   1   2   3   4   5   6
     ---------------------------
    | 6 | 4 | 2 | 2 | 0 | 0 | 0 |
     ---------------------------
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
    }}

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

    max_offset: usize, // the furthest that the longest diagonal can go

    highest_diagonal: isize, // diagonal at the top right corner of the matrix
    lowest_diagonal: isize,  // diagonal at the bottom left corner of the matrix

    current_high: isize,
    current_low: isize,
}

impl WaveFronts {
    // create a new set of wavefronts
    fn new(text_len: usize, query_len: usize) -> Self {


        // The furthest from the origin that we expect the wavefront that travels the furthest to go.
        // This will be along the longest diagonal, the diagonal at the origin.
        // Should we add two because we start along 1 on the row and col of the matrix?
        let max_offset = query_len+text_len+2;


        let expected_wavefront_count = query_len+text_len;
        let mut wfs: Vec<WFs> = Vec::with_capacity(expected_wavefront_count);

        (0..=text_len+query_len).for_each(|_| {
            wfs.push(WFs::new(max_offset))
        });

        Self {
            // create it with at least the zero score wavefront
            wavefronts:  wfs, // initialize with the expected number of wavefronts
            diagonals: text_len+query_len,
            max_offset,
            lowest_diagonal: -(query_len as isize),
            highest_diagonal: text_len as isize,
            current_low: 0,
            current_high: 0,
        }
    }

    fn diagnoal_index(&self, diagonal: isize) -> Option<usize> {

        if diagonal < self.lowest_diagonal || diagonal > self.highest_diagonal {
            return None;
        }

        let index = diagonal + abs!(self.lowest_diagonal);
        Some(index as usize)
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

        if VERBOSITY_LEVEL > 2 {
            eprintln!("[wfa::wf_extend]");
            eprintln!("\tdiagonal={} offset={}", diagonal, offset);
        }


        while match_lambda(v, h) {
            if VERBOSITY_LEVEL > 2 {
                eprintln!("\t\tv={} h={}", v, h);
            }

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
) {
    if VERBOSITY_LEVEL > 2 {
        eprintln!("\n[wfa::wf_next]");
    }

    if score <= 0 {
        panic!("[wfa::wf_next] score={}, is too low.", score);
    }

    if VERBOSITY_LEVEL > 2 {
        eprint!("\tCreating wavefront for score {}.\n", score);
        eprint!("\t\tPrevious mlo {}, mhi {}",
                wavefronts.get_wavefront(score-1).unwrap().mwavefront.low,
                wavefronts.get_wavefront(score-1).unwrap().mwavefront.high);
    }

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

    if VERBOSITY_LEVEL > 2 {
        eprint!("\t\tCurrent lo {}, hi {}\n", lo, hi);
        eprint!("\t\tCurrent wavefronts.max_offset {}\n", wavefronts.max_offset);

        /*
        eprint!("\t\tCurrent wavefronts.max_offset {},  dlen {}, mlen {}, ilen {}\n",
                wavefronts.max_offset,
                wavefronts.get_wavefront(score).unwrap().mwavefront.offsets.len(),
                wavefronts.get_wavefront(score).unwrap().dwavefront.offsets.len(),
                wavefronts.get_wavefront(score).unwrap().iwavefront.offsets.len()
        );
         */
    }

    wavefronts.current_high = hi;
    wavefronts.current_low = lo;

    for diagonal in lo..=hi {

        let k: usize = match wavefronts.diagnoal_index(diagonal) {
            Some(k) => { k },
            None => {
                if VERBOSITY_LEVEL > 2 {
                    eprint!("\t\tdiagonal={} is out of bounds\n", diagonal);
                }

                continue
            }
        };

        let imax: i16 = vec![
            if iscore - io - ie < 0 || k == 0 { iscore as i16 } else { wavefronts.get_wavefront(s-o-e).unwrap().mwavefront.offsets[k-1] },
            if iscore - ie < 0 || k == 0 { iscore as i16 } else { wavefronts.get_wavefront(s-e).unwrap().iwavefront.offsets[k-1] }
        ].iter().max().unwrap() + 1;

        let dmax: i16 = *vec![
            if iscore - io - ie < 0 { iscore as i16 } else { wavefronts.get_wavefront(s-o-e).unwrap().mwavefront.offsets[k+1] },
            if iscore - ie < 0 { iscore as i16 } else { wavefronts.get_wavefront(s-e).unwrap().dwavefront.offsets[k+1] }
        ].iter().max().unwrap();

        let mmax: i16 = *vec![
            if iscore - ix < 0 { iscore as i16 } else { wavefronts.get_wavefront(s-x).unwrap().mwavefront.offsets[k] + 1 },
            wavefronts.get_wavefront(s).unwrap().iwavefront.offsets[k],
            wavefronts.get_wavefront(s).unwrap().dwavefront.offsets[k]
        ].iter().max().unwrap();

        if VERBOSITY_LEVEL > 3 {
            eprint!("\t\tk={}, diagonal={}, imax={}, dmax={}, mmax={}\n",
                    k, diagonal, imax, dmax, mmax);
        }

        wavefronts.wavefronts[score].iwavefront.offsets[k] = imax;
        wavefronts.wavefronts[score].dwavefront.offsets[k] = dmax;
        wavefronts.wavefronts[score].mwavefront.offsets[k] = mmax;
    };
}

fn wf_align(query: &str, text: &str) -> Alignment {
    // TODO: is this the best idea? it lets us easily access the wavefronts
    // assume m >= n
    let query_len = query.len(); // n
    let text_len = text.len();   // m

    // starting diagonal
    // this is the diagonal that reaches point n,m
    let a_k: usize = abs!(text_len as i32 - query_len as i32) as usize;

    // max offset on diagonal a_k
    // it is also the furthest point (offset) from the starting diagonal (a_k) in the DP matrix
    // this offset reaches point n,m
    let a_offset = text_len;

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
    let mut wavefronts = WaveFronts::new(text_len, query_len);

    // initial wavefront set
    let initial_wavefront_set: &mut WFs = match wavefronts.get_wavefront_mut(start_score) {
        Some(wf) => wf,
        None => panic!("Initial wavefront not set"),
    };
    initial_wavefront_set
        .mwavefront
        .offsets
        .insert(start_offset, 0);

    let match_lambda = |v: usize, h: usize| {
            v < query_len && h < text_len && (query.chars().nth(v).unwrap() == text.chars().nth(h).unwrap())
        };

    let end_reached = |score: usize, wavefronts: &WaveFronts| -> bool {
        (wavefronts.current_high > wavefronts.highest_diagonal  &&
         wavefronts.current_low < wavefronts.lowest_diagonal) ||
        wavefronts.get_wavefront(score).unwrap().mwavefront.at_offset(a_k) >= a_offset as i16
    };

    if VERBOSITY_LEVEL > 0 {
        eprintln!("[wfa::wfalign]");
        eprintln!("\tquery= {}", if query_len < 10 { &query[..] } else { &query[..10] });
        eprintln!("\ttext=  {}", if text_len < 10 { &text[..] } else { &text[..10] });
        eprintln!("\ttext_len(m)={}, query_len(n)={}, a_k={}, a_offset={}, diagonals={}, max_offset={}, lowest_diagonal={}, highest diagonal={}\n",
                query_len, text_len,
                a_k, a_offset,
                wavefronts.diagonals, wavefronts.max_offset,
                wavefronts.lowest_diagonal, wavefronts.highest_diagonal);
    }

    let mut score = start_score;
    loop {
        let mwavefront = &mut wavefronts.get_wavefront_mut(score).unwrap().mwavefront;
        wf_extend(
            mwavefront,
            match_lambda,
        );

        // TODO: remove
        // prevents an infinite loop
        let fake_max = 3*std::cmp::max(query_len, text_len);
        if score > fake_max {
            eprintln!("\n\n");
            eprintln!("[wfa::wfalign]");
            eprintln!("\tAbort: current score {} is beyond 3*{} = {}.",
                    score, std::cmp::max(query_len, text_len), fake_max);
            break;
        }

        if VERBOSITY_LEVEL > 2 {
            println!("a_k={1}, score[{}].offset[{}]={} a_offset={} {:?}",
                     score, a_k,
                     wavefronts.get_wavefront(score).unwrap().mwavefront.at_offset(a_k),
                     a_offset,
                     wavefronts.get_wavefront(score).unwrap().mwavefront
            );
        }

        // if we have computed the score for the max cell (n,m)
        if end_reached(score, &wavefronts) { break; }

        score += 1;

        // compute the wavefront for the newly incremented score
        wf_next(&mut wavefronts, query, text, score);
    }

    let aln = Alignment {
        num_wavefronts: score+1,
        alignment_score: wavefronts.get_wavefront(score).unwrap().mwavefront.at_offset(a_k),
        q_len: query_len,
        t_len: text_len,
        central_diagonal: a_k
    };

    if VERBOSITY_LEVEL > 0 {
        eprintln!("\n");
        eprintln!("{:#?}", aln);
        eprintln!("\n\n\n\
                   ................................................................\
                   ................................................................\
                   \n\n\n");
    }

    aln
}

#[cfg(test)]
mod tests {
    use super::*;
    use super::super::needleman_wunsch as nw;

    #[test]
    fn test_align_equal_length_same_sequence() {
        // same sequence
        let text  = "GAGATA";
        let query = "GAGATA";
        let aln_wfa = wf_align(query, text);

        assert_eq!(aln_wfa.num_wavefronts, 1);

        let aln_nw = nw::compute_alginment_score(query, text);
        assert_eq!(aln_wfa.alignment_score as i32, aln_nw);
    }

    #[test]
    fn test_align_equal_length_different_sequences() {
        let text  = "GAGAAT";
        let query = "GAGATA";
        let aln_wfa = wf_align(query, text).alignment_score as i32;
        let aln_nw = nw::compute_alginment_score(query, text);
        //assert_eq!(aln_wfa, aln_nw);

        let text  = "GATACA";
        let query = "GAGATA";
        let aln_wfa = wf_align(query, text).alignment_score as i32;
        assert_eq!(aln_wfa, 8);

        let text  = "GAGGTACAAT";
        let query = "GAGGTACAAG";
        let aln_wfa = wf_align(query, text).alignment_score as i32;
        let aln_nw = nw::compute_alginment_score(query, text);
        assert_eq!(aln_wfa, aln_nw);

        let text  = "TCTTTACTCGCGCGTTGGAGAAATACAATAGTTCTTTACTCGCGCGTTGGAGAAATACAATAGT\
                     TCTTTACTCGCGCGTTGGAGAAATACAATAGTTCTTTACTCGCGCGTTGGAGAAATACAATAGT";
        let query = "TCTATACTGCGCGTTTGGAGAAATAAAATAGTTCTATACTGCGCGTTTGGAGAAATAAAAT\
                     AGTTCTATACTGCGCGTTTGGAGAAATAAAATAGTTCTATACTGCGCGTTTGGAGAAATAA\
                     AATAGT";
        let aln = wf_align(query, text);
        assert_eq!(aln.alignment_score as i32, nw::compute_alginment_score(query, text));
    }

    #[test]
    fn test_align_longer_text() {
        let text  = "GAGGTACAAT";
        let query = "GAGATA";
        let aln = wf_align(query, text);
        assert_eq!(aln.alignment_score, 0);
    }

    #[test]
    fn test_align_longer_query() {
        let text  = "GAGA";
        let query = "GAGGT";
        let aln = wf_align(query, text);
        let aln_nw = nw::compute_alginment_score(query, text);
        assert_eq!(aln.alignment_score, 0);
    }

    #[test]
    fn test_wf_diagonal_index() {
        let text  = "GATACA";
        let query = "GAGATA";

        let wf = WaveFronts::new(text.len(), query.len());

        assert_eq!(wf.lowest_diagonal, negate!(query.len()));
        assert_eq!(wf.diagnoal_index(text.len() as isize), Some(text.len()+query.len()));

        let negative_query_len = negate!(query.len());
        assert_eq!(wf.diagnoal_index(negative_query_len), Some(0));
    }

    #[test]
    fn test_marcos() {
        assert_eq!(k!(14, 12), 2);
        assert_eq!(abs!(-14), 14);
        assert_eq!(abs!(5), 5);
    }
}
