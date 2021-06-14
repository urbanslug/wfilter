use coitrees;
use std::time::Instant;

// local
mod wflambda;

mod cli;
mod index;
mod io;
mod paf;
mod types;

fn main() {
    // Parse CLI args
    let args: types::CliArgs = cli::start();
    let paf_file_path: String = args.input_paf;

    // Parse the PAF input file
    let now = Instant::now();
    if args.debug {
        eprintln!("Parsing input PAF: {}", paf_file_path)
    }
    let paf = paf::PAF::from_file(&paf_file_path[..]);
    if args.debug {
        eprintln!(
            "Done parsing. Time taken {}s",
            now.elapsed().as_millis() as f64 / 1000.0
        )
    }

    let (_query_index, _target_index): (
        coitrees::COITree<types::AlignmentMetadata, u32>,
        coitrees::COITree<types::AlignmentMetadata, u32>,
    ) = index::index_paf_matches(&paf);
}
