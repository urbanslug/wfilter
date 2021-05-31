use coitrees::IntervalNode;
use coitrees::COITree;

use std::time::Instant;

mod cli;
mod types;
mod io;

use types::PAF;
use types::PafAlignment;
use types::CliArgs;

fn main() {
    // Parse CLI args
    let args: CliArgs = cli::start();
    let paf_file_path: String =  args.input_paf;

    // Parse the PAF input file
    let now = Instant::now();
    if args.debug { eprintln!("Parsing input PAF: {}", paf_file_path) }
    let paf = PAF::from_file(&paf_file_path[..]);
    let alignments: &Vec<PafAlignment> = paf.get_alignments();
    if args.debug { eprintln!("Done parsing. Time taken {}s", now.elapsed().as_millis() as f64 / 1000.0) }

    // Read the PAF into a COITree
    if args.debug { eprintln!("Building interval tree") }
    let now = Instant::now();
    let interval_nodes:  Vec<IntervalNode<PafAlignment, u32>> = alignments
        .into_iter()
        .map(|a| IntervalNode::<PafAlignment, u32>::new(a.first(), a.last(), a.clone()))
        .collect();
    let _interval_tree: COITree<PafAlignment, u32> = COITree::new(interval_nodes);
    if args.debug { eprintln!("Done generating interval tree. Time taken {}s", now.elapsed().as_millis() as f64 / 1000.0) }
}
