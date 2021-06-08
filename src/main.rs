use std::time::Instant;
use coitrees::COITree;

mod cli;
mod io;
mod types;
// mod wflambda;
mod paf;
mod index;

// use wflambda::types;
use types::CliArgs;
use paf::PafAlignment;
use paf::PAF;

fn main() {
    // Parse CLI args
    let args: CliArgs = cli::start();
    let paf_file_path: String = args.input_paf;

    // Parse the PAF input file
    let now = Instant::now();
    if args.debug {
        eprintln!("Parsing input PAF: {}", paf_file_path)
    }
    let paf = PAF::from_file(&paf_file_path[..]);
    if args.debug {
        eprintln!(
            "Done parsing. Time taken {}s",
            now.elapsed().as_millis() as f64 / 1000.0
        )
    }

    // Read the query and target into a COITrees
    if args.debug {
        eprintln!("Generating query interval tree")
    }
    let now = Instant::now();
    let _query_coitree: COITree<PafAlignment, u32> = paf.gen_query_coitree();
    if args.debug {
        eprintln!(
            "Done generating query interval tree. Time taken {}s",
            now.elapsed().as_millis() as f64 / 1000.0
        )
    }

    if args.debug {
        eprintln!("Generating target interval tree")
    }
    let now = Instant::now();
    let _target_coitree: COITree<PafAlignment, u32> = paf.gen_target_coitree();
    if args.debug {
        eprintln!(
            "Done generating target interval tree. Time taken {}s",
            now.elapsed().as_millis() as f64 / 1000.0
        )
    }

    // TODO: remove
    // eprintln!("{}\t{}", _target_coitree.len(), _query_coitree.len());
}
