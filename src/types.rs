use std::fmt;

pub struct CliArgs {
    pub debug: bool,
    pub input_paf: String,
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
    pub fn new(debug: bool, paf_filepath: &str) -> Self {
        CliArgs {
            debug,
            input_paf: String::from(paf_filepath),
        }
    }
}
