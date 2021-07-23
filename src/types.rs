use coitrees;
use std::fmt;

#[derive(Copy, Clone)]
pub struct Penalties {
    pub mismatch: u64,
    pub matches: u64,
    pub gap_open: u64,
    pub gap_extend: u64,
}

pub struct CliArgs {
    pub verbosity_level: u8,
    pub input_paf: String,
    pub target_fasta: String,
    pub query_fasta: String,
    pub penalties: Penalties,
}

impl CliArgs {
    #[allow(dead_code)] // TODO: exists for testing purposes
    pub fn new(verbosity_level: u8,
               paf_filepath: &str,
               target_filepath: &str,
               query_filepath: &str,
               penalties: Option<Penalties>) -> Self {
        let penalties = match penalties {
            Some(p) => p,
            _ => Penalties {
                mismatch: 4,
                matches: 0,
                gap_open: 6,
                gap_extend: 2,
            }
        };

        CliArgs {
            verbosity_level,
            input_paf: String::from(paf_filepath),
            target_fasta: String::from(target_filepath),
            query_fasta: String::from(query_filepath),
            penalties,
        }
    }
}

#[derive(PartialEq, Debug)]
pub struct Interval(pub u32, pub u32, pub usize);

// dummy type to hold metadata
pub type AlignmentMetadata = usize;

pub type Index = coitrees::COITree<AlignmentMetadata, u32>;

#[derive(PartialEq, Clone, Copy)]
pub enum SequenceType {
    Target,
    Query,
}

#[derive(PartialEq, Clone, Copy)]
pub enum Strand {
    Forward,
    Reverse,
}

impl fmt::Debug for Strand {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let x: char = if *self == Strand::Forward { '+' } else { '-' };
        f.write_str(&x.to_string())
    }
}
