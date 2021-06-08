use crate::paf::PafAlignment;
use crate::types::{Interval, SequenceType};
use std::str::FromStr;

// compute matc
fn compute_match_intervals(seq_type: SequenceType, start: u32, stop: u32, cigar: &str) -> Vec<Interval> {
    let mut intervals: Vec<Interval> = Vec::new();

    // ordered
    // let mut matches = Vec::<u32>::new();
    // let mut matches = Vec::<u32>::new();
    let mut buffer = String::new();
    let mut cursor = start;
    // let mut high = start;

    cigar.chars().for_each(|c: char| {
        match c {
            'M' => {
                let m: u32 = u32::from_str(&buffer[..]).unwrap();

                intervals.push(Interval(cursor, cursor+m));
                cursor+=m;
                buffer.clear();
            },
            'X' =>  {
                let x: u32 = u32::from_str(&buffer[..]).unwrap();
                cursor += x;
                buffer.clear();
            },
            'I' =>  {
                let i: u32 = u32::from_str(&buffer[..]).unwrap();
                if seq_type == SequenceType::Target { cursor -= i } else {cursor += i};
                buffer.clear();
            },
            'D' =>  {
                let d: u32 = u32::from_str(&buffer[..]).unwrap();
                if seq_type == SequenceType::Target { cursor += d } else {cursor -= d};
                buffer.clear();
            },
            _ => {
                //TODO: check that it's a number in a char
                buffer.push(c);
            }
        }
    });

    intervals
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_compute_match_intervals() {
        let intervals_computed: Vec<Interval> = compute_match_intervals(SequenceType::Query, 0, 330243, "330243M");
        let intervals: Vec<Interval> = vec![Interval(0, 330243)];
        assert_eq!(intervals, intervals_computed);


        let cg = "15M1I158M1I24M1I169M1I1147M1I24M1I851M1I13M1I3900M1D25M1I874M4I10847M3D4400M1I1494M1D4041M1I8577M14I1340M2D21138M2I7776M6D3563M2I83120M10D5541M2D27729M1I2M13I49698M1I5030M2I17541M1D22531M1I187M1D458M1D80M1I75M1I266M1I48M1I269M1I460M1D240M";
        let intervals_computed: Vec<Interval> = compute_match_intervals(SequenceType::Query, 41052,  324759, cg);
        let intervals: Vec<Interval> = vec![];
        assert_eq!(intervals, intervals_computed);

        let intervals_computed: Vec<Interval> = compute_match_intervals(SequenceType::Target, 0, 283680, cg);
        let intervals: Vec<Interval> = vec![];
        assert_eq!(intervals, intervals_computed);
    }
}
