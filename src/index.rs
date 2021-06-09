use coitrees;
use std::convert::TryFrom;
use std::str::FromStr;

use super::paf;
use super::types;

// TODO: account for strand & stop
fn compute_match_intervals(
    seq_type: types::SequenceType,
    _strand: types::Strand,
    start: u32,
    _stop: u32,
    cigar: &str,
) -> Vec<types::Interval> {
    let mut intervals: Vec<types::Interval> = Vec::new();
    let mut buffer = String::new();
    let mut cursor = start;

    cigar.chars().for_each(|c: char| {
        match c {
            'M' | '=' => {
                // TODO: consider the ambiguity of M being match/mismatch
                let m: u32 = u32::from_str(&buffer[..]).unwrap();
                intervals.push(types::Interval(cursor, cursor + m));
                cursor += m;
                buffer.clear();
            }
            'X' => {
                let x: u32 = u32::from_str(&buffer[..]).unwrap();
                cursor += x;
                buffer.clear();
            }
            'I' => {
                let i: u32 = u32::from_str(&buffer[..]).unwrap();
                if seq_type == types::SequenceType::Target {
                    cursor -= i
                } else {
                    cursor += i
                };
                buffer.clear();
            }
            'D' => {
                let d: u32 = u32::from_str(&buffer[..]).unwrap();
                if seq_type == types::SequenceType::Target {
                    cursor += d
                } else {
                    cursor -= d
                };
                buffer.clear();
            }
            _ => {
                // At this point we expect the char to be a base 10 digit i.e '0', '1' .. '9'
                match c {
                    b if b.is_digit(10) => buffer.push(b),
                    b if b.is_ascii_alphabetic() => panic!(
                        "[index::coompute_match_intervals] Unexpected char {} in CIGAR string",
                        b
                    ),
                    _ => panic!(
                        "[index::coompute_match_intervals] Unknown char {} in CIGAR string",
                        c
                    ),
                }
            }
        }
    });

    intervals
}

pub fn index_paf_matches(
    p: &paf::PAF,
) -> (
    coitrees::COITree<types::AlignmentMetadata, u32>,
    coitrees::COITree<types::AlignmentMetadata, u32>,
) {
    let alignments: &Vec<paf::PafAlignment> = p.get_alignments();
    let mut query_intervals: Vec<types::Interval> = Vec::new();
    let mut target_intervals: Vec<types::Interval> = Vec::new();

    alignments.iter().for_each(|a: &paf::PafAlignment| {
        let mut t = compute_match_intervals(
            types::SequenceType::Target,
            a.strand,
            a.target_start,
            a.target_end,
            &a.cigar[..],
        );
        let mut q = compute_match_intervals(
            types::SequenceType::Query,
            a.strand,
            a.query_start,
            a.query_end,
            &a.cigar[..],
        );

        query_intervals.append(&mut q);
        target_intervals.append(&mut t);
    });

    let gen_coitree =
        |intervals: Vec<types::Interval>| -> coitrees::COITree<types::AlignmentMetadata, u32> {
            // Generate coitrees::IntervalNodes
            let interval_nodes: Vec<coitrees::IntervalNode<types::AlignmentMetadata, u32>> =
                intervals
                    .iter()
                    .map(|types::Interval(start, stop): &types::Interval| {
                        let start = i32::try_from(*start)
                            .expect("[index::index_paf] Could not convert start u32 to i32");
                        let end = i32::try_from(*stop)
                            .expect("[index::index_paf] Could not convert end u32 to i32");
                        coitrees::IntervalNode::<types::AlignmentMetadata, u32>::new(start, end, ())
                    })
                    .collect();

            coitrees::COITree::new(interval_nodes)
        };

    let query_coitree = gen_coitree(query_intervals);
    let target_coitree = gen_coitree(target_intervals);

    (query_coitree, target_coitree)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[ignore]
    #[test]
    fn test_compute_match_intervals() {
        // Forward
        let intervals_computed: Vec<types::Interval> = compute_match_intervals(
            types::SequenceType::Query,
            types::Strand::Forward,
            0,
            330243,
            "330243M",
        );
        let intervals: Vec<types::Interval> = vec![types::Interval(0, 330243)];
        assert_eq!(intervals, intervals_computed);

        let cg = "15M1I158M1I24M1I169M1I1147M1I24M1I851M1I13M1I3900M1D25M1I874M4I10847M3D4400M1I1494M1D4041M1I8577M14I1340M2D21138M2I7776M6D3563M2I83120M10D5541M2D27729M1I2M13I49698M1I5030M2I17541M1D22531M1I187M1D458M1D80M1I75M1I266M1I48M1I269M1I460M1D240M";
        let intervals_computed: Vec<types::Interval> = compute_match_intervals(
            types::SequenceType::Query,
            types::Strand::Forward,
            41052,
            324759,
            cg,
        );
        let intervals: Vec<types::Interval> = vec![];
        assert_eq!(intervals, intervals_computed);

        let intervals_computed: Vec<types::Interval> = compute_match_intervals(
            types::SequenceType::Target,
            types::Strand::Forward,
            0,
            283680,
            cg,
        );
        let intervals: Vec<types::Interval> = vec![];
        assert_eq!(intervals, intervals_computed);

        // Reverse strand
    }

    #[test]
    fn test_index_paf() {
        static TEST_PAF_STRING: &str = "\
        qry\t330243\t0\t330243\t+\ttgt\t330243\t0\t330243\t330243\t330243\t60\tNM:i:0\tms:i:660486\tAS:i:660486\tnn:i:0\ttp:A:P\tcm:i:62290\ts1:i:329202\ts2:i:262341\tde:f:0\trl:i:2730\tcg:Z:330243M\
        \n\
        qry\t329347\t41052\t324759\t+\ttgt\t283680\t0\t283680\t283613\t283736\t0\tNM:i:123\tms:i:566760\tAS:i:566760\tnn:i:0\ttp:A:S\tcm:i:53397\ts1:i:282348\tde:f:0.0003\trl:i:2765\tcg:Z:15M1I158M1I24M1I169M1I1147M1I24M1I851M1I13M1I3900M1D25M1I874M4I10847M3D4400M1I1494M1D4041M1I8577M14I1340M2D21138M2I7776M6D3563M2I83120M10D5541M2D27729M1I2M13I49698M1I5030M2I17541M1D22531M1I187M1D458M1D80M1I75M1I266M1I48M1I269M1I460M1D240M
";
        let alignments: paf::PAF = paf::PAF::from_str(TEST_PAF_STRING);
        let (query_index, target_index): (
            coitrees::COITree<types::AlignmentMetadata, u32>,
            coitrees::COITree<types::AlignmentMetadata, u32>,
        ) = index_paf_matches(&alignments);

        // should apply to all of them
        assert_eq!(38, query_index.query_count(0, 330_243));
        // the first match in the second alignment plus the first alignment which covers everything
        assert_eq!(2, query_index.query_count(41_052, 41_067));
        assert_eq!(query_index.len(), target_index.len());
    }
}
