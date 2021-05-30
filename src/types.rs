/*
PAF: a Pairwise mApping Format
Spec https://github.com/lh3/miniasm/blob/master/PAF.md

PAF is a text format describing the approximate mapping positions between two
set of sequences. PAF is TAB-delimited with each line consisting of the
following predefined fields:

|Col|Type  |Description                               |
|--:|:----:|:-----------------------------------------|
|1  |string|Query sequence name                       |
|2  |int   |Query sequence length                     |
|3  |int   |Query start (0-based; BED-like; closed)   |
|4  |int   |Query end (0-based; BED-like; open)       |
|5  |char  |Relative strand: "+" or "-"               |
|6  |string|Target sequence name                      |
|7  |int   |Target sequence length                    |
|8  |int   |Target start on original strand (0-based) |
|9  |int   |Target end on original strand (0-based)   |
|10 |int   |Number of residue matches                 |
|11 |int   |Alignment block length                    |
|12 |int   |Mapping quality (0-255; 255 for missing)  |

If PAF is generated from an alignment, column 10 equals the number of sequence
matches, and column 11 equals the total number of sequence matches, mismatches,
insertions and deletions in the alignment. If alignment is not available,
column 10 and 11 are still required but may be highly inaccurate.

A PAF file may optionally contain SAM-like typed key-value pairs at the end of
each line.
 */
use std::str;
use std::str::FromStr;

use std::convert::TryFrom;

use crate::io;

// -------
// Structs
// --------
// A struct over a single alignment, that is, a single line of a PAF file
#[derive(Debug, Clone)]
pub struct PafAlignment {
    query:           String, // Query sequence name
    query_length:    u32,    // Query sequence length
    query_start:     u32,    // Query start (0-based; BED-like; closed)
    query_end:       u32,    // Query end (0-based; BED-like; open)
    strand:          char,   // Relative strand: "+" or "-"
    target:          String, // target sequence name
    target_length:   u32,    // Target sequence length
    target_start:    u32,    // Target start on original strand (0-based)
    target_end:      u32,    // Target end on original strand (0-based)
    // residue_matches: u32,    // Number of residue matches
    // block_len:       u32,    // Alignment block length
    // quality:         String, // Mapping quality (0-255; 255 for missing)
    // cigar:           String, // SAM style CIGAR string TODO: specify CIGAR version
}

// A struct over the entire PAF file
#[derive(Debug)]
pub struct PAF {
    alignments: Vec<PafAlignment>
}

// -------
// Traits
// -------

impl PafAlignment {
    pub fn from_lines(lines: Vec<String>) -> Vec<PafAlignment> {
        lines
            .iter()
            .map(|line| PafAlignment::from_line(&line[..]))
            .collect::<Vec<PafAlignment>>()
    }

    pub fn from_line(line: &str) -> Self {
        let it: Vec<&str> = line.split_whitespace().collect();

        //need a more robust way to index into the vector
        PafAlignment{
            query:         it[0].to_string(),
            query_length:  u32::from_str(it[1]).unwrap(),
            query_start:   u32::from_str(it[2]).unwrap(),
            query_end:     u32::from_str(it[3]).unwrap(),
            strand:        char::from_str(it[4]).unwrap(),
            target:        it[5].to_string(),
            target_length: u32::from_str(it[6]).unwrap(),
            target_start:  u32::from_str(it[7]).unwrap(),
            target_end:    u32::from_str(it[8]).unwrap(),
        }
    }

    pub fn first(&self) -> i32 {
        i32::try_from(self.query_start).expect("Couldn't get first")
    }

    pub fn last(&self) -> i32 {
        i32::try_from(self.target_end).expect("Couldn't get last")
    }
}


impl PAF {
    pub fn from_file(file_name: &str) -> PAF {
        // read PAF file & return a vector of Strings for each line
        let lines: Vec<String> = io::read_file(&file_name[..]);

        // Parse each line into a PafAlignment
        let alignments: Vec<PafAlignment> = PafAlignment::from_lines(lines);

        PAF{ alignments }
    }

    pub fn get_alignments(&self) -> &Vec<PafAlignment> {
        &self.alignments
    }
}
