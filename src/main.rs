use coitrees::IntervalNode;
use coitrees::COITree;

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
    let interval_nodes:  Vec<IntervalNode<PafAlignment, u32>> = alignments
        .into_iter()
        .map(|a| IntervalNode::<PafAlignment, u32>::new(a.first(), a.last(), a.clone()))
        .collect();

    let interval_tree: COITree<PafAlignment, u32> = COITree::new(interval_nodes);

    println!("{}", interval_tree.len());
}
