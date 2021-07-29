#![allow(dead_code, unused_variables)]
use std::cmp::max;
use std::io::{self, Write};

const MISMATCH_SCORE: i32 = -1;
const GAP_SCORE: i32 = -2;
const MATCH_SCORE: i32 = 1;

#[allow(dead_code)]
struct TraceBack {
    cells: Vec<Cell>,
    pattern: String,
    text: String,
}
type Cell = (usize, usize);
#[derive(Debug)]
struct Score {
    score: i32,
    cell: Cell,
}

type AlignmentMatrix = Vec<Vec<i32>>;
struct AlignmentMatrixWrapper<'a> {
    pattern: &'a [u8],
    text: &'a [u8],
    m: AlignmentMatrix,
    score: Option<Score>,
}

impl<'a> AlignmentMatrixWrapper<'a> {
    fn init(text: &'a str, pattern: &'a str) -> Self {
        // initialize a text_length by pattern_length matrix
        // row are text and columns are pattern
        // v is query, h is target. M(h,v)

        /*

             Matrix

                      0 <--------- pattern (j) -------------> max(text_length+1, text_length+1)

                 --------------------------------------------
                 | 0 |  |  |  |  |  |  |  |  |  |  |  |  |  |
             ^   --------------------------------------------
             |   |   |  |  |  |  |  |  |  |  |  |  |  |  |  |
                 --------------------------------------------
             t   |   |  |  |  |  |  |  |  |  |  |  |  |  |  |
             e   --------------------------------------------
             x   |   |  |  |  |  |  |  |  |  |  |  |  |  |  |
             t   --------------------------------------------
                 |   |  |  |  |  |  |  |  |  |  |  |  |  |  |
            (i)  --------------------------------------------
             |   |   |  |  |  |  |  |  |  |  |  |  |  |  |  |
             |   --------------------------------------------
             |   |   |  |  |  |  |  |  |  |  |  |  |  |  |  |
             |   --------------------------------------------

        max(text_length+1, text_length+1)

              */

        // add extra col and row at the top and left for init
        let text_length = text.len() + 1;
        let pattern_length = pattern.len() + 1;

        let v = Vec::with_capacity(pattern_length);
        let mut m = vec![v; max(text_length, pattern_length)]; // h is implicit here

        // init the M(0,0) to 0
        m[0].insert(0, 0);

        // init top row and left most column to save time
        // init top row
        for i in 1..max(pattern_length, text_length) {
            let left_cell = m[0][i - 1];
            m.get_mut(0).unwrap().insert(i, left_cell + GAP_SCORE);
        }

        // init left most column
        for j in 1..max(pattern_length, text_length) {
            let cell_above = m[j - 1][0];
            m.get_mut(j).unwrap().insert(0, cell_above + GAP_SCORE);
        }

        Self {
            pattern: &pattern.as_bytes()[..],
            text: &text.as_bytes()[..],
            m,
            score: None,
        }
    }

    fn get_score(&self) -> i32 {
        self.score.as_ref().unwrap().score
    }
}

/*
Align
-----

Max(
  M[i-1, j-1] + match_or_mismatch(M[i,j]),
  M[i-1, j] + gap_penalty,                  cell to the left
  M[i, j-1] + gap_penalty,                  cell above
)
 */
fn align(matrix_wrapper: &mut AlignmentMatrixWrapper) {
    let mut max_score = Score {
        score: 0,
        cell: (0, 0),
    };

    let text: &[u8] = matrix_wrapper.text;
    let pattern: &[u8] = matrix_wrapper.pattern;
    let aln: &mut AlignmentMatrix = &mut matrix_wrapper.m;

    // Start from M[1,1] because top row and left most columns are full
    let match_or_mismatch = |i: usize, j: usize| -> i32 {
        let i = match text.get(i) {
            Some(v) => v,
            _ => return MISMATCH_SCORE,
        };
        let j = match pattern.get(j) {
            Some(h) => h,
            _ => return MISMATCH_SCORE,
        };

        if i == j {
            MATCH_SCORE
        } else {
            MISMATCH_SCORE
        }
    };

    // does this fill the entire matrix?
    let h = aln.len();
    let v = aln[0].len();
    for i in 1..h {
        for j in 1..v {
            let cell_above = aln[i - 1][j];
            let left_cell = aln[i][j - 1];
            let prev_diagonal = aln[i - 1][j - 1];

            let best_score = max(
                prev_diagonal + match_or_mismatch(i - 1, j - 1),
                max(left_cell + GAP_SCORE, cell_above + GAP_SCORE),
            );
            aln.get_mut(i).unwrap().insert(j, best_score);

            if best_score > max_score.score {
                max_score = Score {
                    score: best_score,
                    cell: (i, j),
                };
            }
        }
    }

    matrix_wrapper.score = Some(max_score);
}

fn compute_traceback(matrix_wrapper: &AlignmentMatrixWrapper) -> TraceBack {
    let matrix = &matrix_wrapper.m;
    let pattern = matrix_wrapper.pattern;
    let text = matrix_wrapper.text;
    let mut h = matrix_wrapper.text.len();
    let mut v = matrix_wrapper.pattern.len();
    let mut traceback: Vec<Cell> = Vec::with_capacity(max(v, h));

    let mut h_s = String::new();
    let mut v_s = String::new();

    loop {
        traceback.insert(0, (h, v));

        if h == 0 && v == 0 {
            break;
        }

        let cell_above = if h == 0 { i32::MIN } else { matrix[h - 1][v] };
        let left_cell = if v == 0 { i32::MIN } else { matrix[h][v - 1] };
        let prev_diagonal = matrix[h - 1][v - 1];

        let lowest = max(prev_diagonal, max(left_cell, cell_above));

        if lowest == cell_above {
            v_s.push(pattern[v - 1] as char);
            h_s.push('-');

            h = h - 1;
        } else if lowest == left_cell {
            v_s.push('-');
            h_s.push(text[h - 1] as char);

            v = v - 1;
        } else {
            h_s.push(text[h - 1] as char);
            v_s.push(pattern[v - 1] as char);
            h = h - 1;
            v = v - 1;
        };
    }

    h_s = h_s.chars().rev().collect();
    v_s = v_s.chars().rev().collect();

    traceback.shrink_to_fit(); // TODO: Is this necessary?
    TraceBack {
        pattern: v_s,
        text: h_s,
        cells: traceback,
    }
}

fn print_alignment(traceback: &TraceBack) {
    println!("{}", traceback.text);
    traceback
        .text
        .chars()
        .zip(traceback.pattern.chars())
        .for_each(|(t, p)| if t == p { print!("|") } else { print!(" ") });
    print!("\n");
    io::stdout().flush().unwrap();

    println!("{}", traceback.pattern);
    println!("");
}

fn needleman_wunsch(text: &str, pattern: &str) {
    let mut aln: &mut AlignmentMatrixWrapper = &mut AlignmentMatrixWrapper::init(text, pattern);
    align(&mut aln);
    let tr = compute_traceback(&aln);
    print_alignment(&tr);
}

pub fn compute_alginment_score(text: &str, pattern: &str) -> i32 {
    let mut aln: &mut AlignmentMatrixWrapper = &mut AlignmentMatrixWrapper::init(text, pattern);
    align(&mut aln);

    aln.get_score()
}

#[cfg(test)]
mod tests {
    use super::*;

    const TEXT: &'static str = "TCTTTACTCGCGCGTTGGAGAAATACAATAGTTCTTTACTCGCGCGTTGGAGAAATACAATAGT\
                                TCTTTACTCGCGCGTTGGAGAAATACAATAGTTCTTTACTCGCGCGTTGGAGAAATACAATAGT";

    const PATTERN: &'static str = "TCTATACTGCGCGTTTGGAGAAATAAAATAGTTCTATACTGCGCGTTTGGAGAAATAAAAT\
                                   AGTTCTATACTGCGCGTTTGGAGAAATAAAATAGTTCTATACTGCGCGTTTGGAGAAATAA\
                                   AATAGT";

    #[ignore]
    #[test]
    fn test_align_equal_length() {
        let text_window = &TEXT[..12];
        let pattern_window = &PATTERN[..12];
        assert_eq!(needleman_wunsch(text_window, pattern_window), ());

        let text_window = &TEXT[..12];
        let pattern_window = &PATTERN[12..24];
        assert_eq!(needleman_wunsch(text_window, pattern_window), ());


        let text = &TEXT[..];
        let pattern = &PATTERN[..];
        let mut aln: &mut AlignmentMatrixWrapper = &mut AlignmentMatrixWrapper::init(text, pattern);
        align(&mut aln);

        assert_eq!(aln.get_score(), 92);
    }

    #[ignore]
    #[test]
    fn test_align_varying_length() {
        let text_window = &TEXT[..14];
        let pattern_window = &PATTERN[..18];
        assert_eq!(needleman_wunsch(text_window, pattern_window), ());

        let text_window = &TEXT[..18];
        let pattern_window = &PATTERN[..14];
        assert_eq!(needleman_wunsch(text_window, pattern_window), ());

        assert_eq!(needleman_wunsch(TEXT, PATTERN), ());
    }
}
