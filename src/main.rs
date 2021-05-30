mod cli;
mod types;
mod io;

use types::PAF;
use types::PafAlignment;

fn main() {
    let paf_file_path: String = cli::start();

    let paf = PAF::from_file(&paf_file_path[..]);

    let y: &PafAlignment = &paf.get_alignments()[0];
    println!("{:#?}", y);
}
