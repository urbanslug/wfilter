pub struct CliArgs {
    pub debug: bool,
    pub input_paf: String,
}

#[derive(PartialEq, Debug)]
pub struct Interval(pub u32, pub u32);

#[derive(PartialEq)]
pub enum SequenceType {
    Target,
    Query
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
