use std::fmt;

pub struct CliArgs {
    pub verbosity_level: u8,
    pub input_paf: String,
    pub target_fasta: String,
    pub query_fasta: String,
}

#[derive(PartialEq, Debug)]
pub struct Interval(pub u32, pub u32);

// dummy type to hold metadata
pub type AlignmentMetadata = ();

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

impl CliArgs {
    #[allow(dead_code)] // TODO: remove dead code
    pub fn new(verbosity_level: u8,
               paf_filepath: &str,
               target_filepath: &str,
               query_filepath: &str) -> Self {
        CliArgs {
            verbosity_level,
            input_paf: String::from(paf_filepath),
            target_fasta: String::from(target_filepath),
            query_fasta: String::from(query_filepath),
        }
    }
}
