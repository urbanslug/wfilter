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
    let paf_file_path: &str = &args.input_paf[..];
    let verbosity = args.verbosity_level;

    // Parse target fasta
    let now = Instant::now();
    if verbosity > 0 {
        eprintln!("[wfilter::main] parsing target: {}", args.target_fasta);
    }
    let target = io::read_fasta(&args.target_fasta[..]);
    if verbosity  > 0 {
        eprintln!(
            "[wfilter::main] done parsing target. Time taken {} seconds",
            now.elapsed().as_millis() as f64 / 1000.0
        );
        eprintln!("target length {}", target.len());
    }

    // Parse query fasta
    let now = Instant::now();
    if verbosity > 0 {
        eprintln!("[wfilter::main] parsing query: {}", args.query_fasta);
    }
    let query = io::read_fasta(&args.query_fasta[..]);
    if verbosity > 0 {
        eprintln!(
            "[wfilter::main] done parsing query. Time taken {} seconds",
            now.elapsed().as_millis() as f64 / 1000.0
        );
        eprintln!("query length {}", query.len());
    }

    // TODO: put these into default args
    if verbosity > 0 {
        eprintln!("[wfilter::main] align");
    }
    let p =  wflambda::wfa::types::Penalties {
        mismatch: 4,
        matches: 0,
        gap_open: 6,
        gap_extend: 2,
    };
    // TODO: pass &p
    let aln = wflambda::wfa::wf_align(target.as_bytes(), query.as_bytes(), p, &args);
    if verbosity > 0 {
        eprintln!(
            "[wfilter::main] done aligning. Time taken {} seconds",
            now.elapsed().as_millis() as f64 / 1000.0
        )
    };

    if verbosity > 1 {
        eprintln!("score {}", aln.score);
        eprintln!("{}", aln.cigar);
    }


    // Parse the PAF input file
    let now = Instant::now();
    if verbosity > 0 {
        eprintln!("Parsing PAF: {}", paf_file_path)
    }
    let paf = paf::PAF::from_file(paf_file_path);
    if verbosity > 0 {
        eprintln!(
            "[wfilter::main] done parsing PAF. Time taken {} seconds",
            now.elapsed().as_millis() as f64 / 1000.0
        )
    }

    let (_query_index, _target_index): (
        coitrees::COITree<types::AlignmentMetadata, u32>,
        coitrees::COITree<types::AlignmentMetadata, u32>,
    ) = index::index_paf_matches(&paf);
}
