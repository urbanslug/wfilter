mod cli;
mod types;
mod io;

use types::PAF;
use types::PafAlignment;

fn main() {
    // Parse CLI args
    let paf_file_path: String = cli::start();

    // Parse the PAF input file
    eprintln!("Parsing input PAF: {}", paf_file_path);
    let paf = PAF::from_file(&paf_file_path[..]);
    let alignments: &Vec<PafAlignment> = paf.get_alignments();

    // Read the PAF into interval trees
}
