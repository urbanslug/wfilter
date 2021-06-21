use std::cmp;

// const MATCH_SCORE: usize = 0;
const MIS_MATCH_SCORE: usize = 4;
const GAP_OPENING_SCORE: usize = 6;
const GAP_EXTENSION_SCORE: usize = 2;

/*
k is the name of the diagonal
We name diagonals based on the row - column
Example: a diagonal passing through 14, 12 is diagonal 2
*/
macro_rules! k {
    ($h: expr, $v: expr) => {{
        $h - $v
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
struct WaveFront {
    low: usize,  // lowest diagonal that this wavfront touches
    high: usize, // top diagonal that this wavfront touches

    /*
    offsets is a vector of diagonals with the index being the diagonal and the
    value being the offset i.e how far we've moved along that diagonal

    Below we have offsets for 7 diagonals along which we've moved 6, 4, 2 ... 0 steps

      0   1   2   3   4   5   6
     ----------------------------
    | 6 | 4 | 2 | 2 | 0 | 0 | 0 |
     ----------------------------
     */
    offsets: Vec<usize>, // matches & mismatches
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
    fn at_offset(&self, o: usize) -> usize {
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
}

impl WaveFronts {
    // create a new set of wavefronts
    fn new(max_offset: usize) -> Self {
        Self {
            // create it with at least the zero score wavefront
            wavefronts: vec![WFs::new(max_offset)],
        }
    }

    fn get_wavefront_mut<'a>(&'a mut self, score: usize) -> Option<&'a mut WFs> {
        self.wavefronts.get_mut(score)
    }

    fn get_wavefront(&self, score: usize) -> Option<&WFs> {
        self.wavefronts.get(score)
    }
}

/*
Take m-wavefront for a score s
the score s is determined in wf_align
 */
fn wf_extend<T>(mwavefront: &mut WaveFront, _query: &str, _text: &str, match_lambda: T)
where
    T: Fn(usize, usize) -> bool,
{
    // k is the diagonal
    // we are going from the lowest to highest diagonal
    for k in mwavefront.low..=mwavefront.high {
        //let mut diagonal_k_offset = *mwavefront.offsets.get(k).unwrap();
        // mwavefront.offsets[k] is M_s_k
        let mut v: usize = mwavefront.offsets[k] - k; // columns
        let mut h: usize = mwavefront.offsets[k];

        while match_lambda(v, h) {
            // extend the wavefront along this diagonal
            mwavefront.offsets[k] = mwavefront.offsets[k] + 1; // why this?
            v += 1;
            h += 1
        }
    }
}

// Compute the next wavefront for the score
// A wavfront will touch a set of diagonals
// from hi to low and be extended along them
fn wf_next(
    wavefronts: &mut WaveFronts,
    _query: &str,
    _text: &str,
    score: usize,
    max_offset: usize,
) {
    let wavefront_set: &WFs = match wavefronts.get_wavefront(score) {
        Some(wf) => wf,
        // create the wavefront
        None => {
            // hopefully this doesn't break
            // make sure we have enough wavefronts
            while wavefronts.wavefronts.len() <= score {
                wavefronts.wavefronts.push(WFs::new(max_offset));
            }
            match wavefronts.get_wavefront(score) {
                Some (w) => w,
                None => panic!("[wfa::wf_next] could not find wavefront for score {} despite pushing enough wavefronts", score),
            }
        }
    };

    let mwavefront = &wavefront_set.mwavefront;
    let iwavefront = &wavefront_set.iwavefront;
    let dwavefront = &wavefront_set.dwavefront;

    let high = cmp::max(
        cmp::max(
            mwavefront.high - MIS_MATCH_SCORE,
            mwavefront.high - GAP_OPENING_SCORE - GAP_EXTENSION_SCORE,
        ),
        cmp::max(
            iwavefront.high - GAP_EXTENSION_SCORE,
            dwavefront.high - GAP_EXTENSION_SCORE,
        ),
    ) + 1;

    let low = cmp::max(
        cmp::max(
            mwavefront.low - MIS_MATCH_SCORE,
            mwavefront.low - GAP_OPENING_SCORE - GAP_EXTENSION_SCORE,
        ),
        cmp::max(
            iwavefront.low - GAP_EXTENSION_SCORE,
            dwavefront.low - GAP_EXTENSION_SCORE,
        ),
    ) + 1;

    for k in low..=high {
        let imax = cmp::max(
            wavefronts
                .get_wavefront(score - GAP_OPENING_SCORE - GAP_EXTENSION_SCORE)
                .unwrap()
                .mwavefront
                .offsets[k - 1],
            wavefronts
                .get_wavefront(score - GAP_EXTENSION_SCORE)
                .unwrap()
                .iwavefront
                .offsets[k - 1],
        ) + 1;

        let dmax = cmp::max(
            wavefronts
                .get_wavefront(score - GAP_OPENING_SCORE - GAP_EXTENSION_SCORE)
                .unwrap()
                .mwavefront
                .offsets[k + 1],
            wavefronts
                .get_wavefront(score - GAP_EXTENSION_SCORE)
                .unwrap()
                .dwavefront
                .offsets[k + 1],
        );

        let mmax = cmp::max(
            wavefronts
                .get_wavefront(score - MIS_MATCH_SCORE)
                .unwrap()
                .mwavefront
                .offsets[k]
                + 1,
            cmp::max(
                wavefronts.get_wavefront(score).unwrap().iwavefront.offsets[k],
                wavefronts.get_wavefront(score).unwrap().dwavefront.offsets[k],
            ),
        );

        // We have to borrow it as mutable at the end
        // unwrap should just work because we got the wavefront at that score above
        let wavefront_set = wavefronts.get_wavefront_mut(score).unwrap();
        wavefront_set.iwavefront.offsets[k] = imax;
        wavefront_set.dwavefront.offsets[k] = dmax;
        wavefront_set.mwavefront.offsets[k] = mmax;
    }
}

fn wf_align(query: &str, text: &str) {
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
    let mut wavefronts = WaveFronts::new(a_offset);
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
        |v: usize, h: usize| query.chars().nth(v).unwrap() == text.chars().nth(h).unwrap();

    loop {
        wf_extend(
            &mut wavefronts.get_wavefront_mut(score).unwrap().mwavefront,
            query,
            text,
            match_lambda,
        );

        // if wavefront of the current score is at offset
        // if we have reaches a_k
        if wavefronts
            .get_wavefront(score)
            .unwrap()
            .mwavefront
            .at_offset(a_k)
            >= a_k
        {
            break;
        }

        // compute the wavefront for the next score
        score += 1;
        wf_next(&mut wavefronts, query, text, score, a_offset);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_align_equal_length() {
        let text = "GATACA";
        let query = "GAGATA";
        assert_eq!(wf_align(query, text), ());
    }

    #[test]
    fn test_marcos() {
        assert_eq!(k!(14, 12), 2);
    }
}
