use coitrees;
use std::collections::HashSet;
use std::time::Instant;

// local
mod wflambda;

mod cli;
mod fasta;
mod index;
mod io;
mod paf;
mod types;
mod utils;

fn filter(
    target: &fasta::FastaFile,
    target_index: &types::Index,
    query: &fasta::FastaFile,
    query_index: &types::Index,
    cli_args: &types::CliArgs,
) -> HashSet<usize> {
    let verbosity = cli_args.verbosity_level;

    let mut query_lines: HashSet<usize> = HashSet::new();
    let mut target_lines: HashSet<usize> = HashSet::new();

    let mut backtrace_lambda = |query: (i32, i32), target: (i32, i32)| {
        let y = |i: &coitrees::IntervalNode<types::AlignmentMetadata, u32>| {
            target_lines.insert(i.metadata);
        };
        let z = |i: &coitrees::IntervalNode<types::AlignmentMetadata, u32>| {
            query_lines.insert(i.metadata);
        };

        target_index.query(target.0, target.1, y);
        query_index.query(query.0, query.1, z);
    };

    let now = Instant::now();
    for t in target.iter() {
        for q in query.iter() {
            if verbosity > 1 {
                eprintln!(
                    "[wfilter::main::filter] Aligning {} length: {} bases and {} length: {} bases",
                    std::str::from_utf8(&t.header[..]).unwrap(),
                    utils::pretty_print_int(t.seq.len() as isize),
                    std::str::from_utf8(&q.header[..]).unwrap(),
                    utils::pretty_print_int(q.seq.len() as isize)
                );
            };

            let aln =
                wflambda::wfa::wf_align(&t.seq[..], &q.seq[..], cli_args, &mut backtrace_lambda);

            if verbosity > 3 {
                eprintln!("score {}", aln.score);
                eprintln!("{}", aln.cigar);
            }
        }
    }

    let mut filtered = target_lines;
    filtered.extend(&query_lines);

    if verbosity > 0 {
        eprintln!(
            "\t[wfilter::main::align] finished all alignments. Time taken {} seconds",
            now.elapsed().as_millis() as f64 / 1000.0
        )
    };

    filtered
}

fn main() {
    // ------------
    //    CLI
    // ------------

    // Parse CLI args
    let args: types::CliArgs = cli::start();
    let paf_file_path: &str = &args.input_paf[..];
    let verbosity = args.verbosity_level;

    // ------------
    //     PAF
    // ------------

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

    // index
    let (query_index, target_index): (types::Index, types::Index) = index::index_paf_matches(&paf);

    // ------------
    //     FASTA
    // ------------

    // Parse target fasta
    let now = Instant::now();
    if verbosity > 0 {
        eprintln!("[wfilter::main] parsing target: {}", args.target_fasta);
    }

    let target: fasta::FastaFile = fasta::Fasta::from_path(&args.target_fasta[..]);
    if verbosity > 0 {
        eprintln!(
            "[wfilter::main] done parsing target. Time taken {} seconds",
            now.elapsed().as_millis() as f64 / 1000.0
        );
        eprintln!(
            "[wfilter::main] Number of sequences in the target {}",
            target.len()
        );
    }

    // Parse query fasta
    let now = Instant::now();
    if verbosity > 0 {
        eprintln!("[wfilter::main] parsing query: {}", args.query_fasta);
    }
    let query: fasta::FastaFile = fasta::Fasta::from_path(&args.query_fasta[..]);
    if verbosity > 0 {
        eprintln!(
            "[wfilter::main] done parsing query. Time taken {} seconds",
            now.elapsed().as_millis() as f64 / 1000.0
        );
        eprintln!(
            "[wfilter::main] Number of sequences in the query {}",
            query.len()
        );
    }

    // ------------
    //     Filter
    // ------------

    let now = Instant::now();
    if verbosity > 0 {
        eprintln!("[wfilter::main] Filtering");
    }

    let lines: HashSet<usize> = filter(&target, &target_index, &query, &query_index, &args);
    let mut lines = lines.into_iter().collect::<Vec<usize>>();
    lines.sort();

    io::copy_filtered(paf_file_path, &lines);
    if verbosity > 0 {
        eprintln!(
            "[wfilter::main] done filtering. Time taken {} seconds",
            now.elapsed().as_millis() as f64 / 1000.0
        );
    }
}

#[cfg(test)]
mod tests {
    use super::fasta::Fasta;
    use super::index;
    use super::paf;
    use super::types::{CliArgs, Index};
    use super::*;

    static PAF_STRING: &str = "\
    qry\t330243\t0\t330243\t+\ttgt\t330243\t0\t330243\t330243\t330243\t60\tNM:i:0\tms:i:660486\
    \tAS:i:660486\tnn:i:0\ttp:A:P\tcm:i:62290\ts1:i:329202\ts2:i:262341\tde:f:0\trl:i:2730\
    \tcg:Z:330243M";

    static TEXT: &str = ">species_y\n\
                         TCTATACTGCGCGTTTATCTAGGAGAAATAAAATAGTTCTATACTGCGCGTTTGGAGAAATAAAATAGT\
                         TCTATACTGCGCGTTTGGAGAAATAACTATCAATAGTTCTATACTGCGCGTTTGGAGAAATAAAATAGT";
    static QUERY: &str = ">species_x\n\
                          TCTTTACTCGCGCGTTGGAGAAATACAATAGTTCTTTACTCGCGCGTTGGAGAAATACAATAGT\
                          TCTTTACTCGCGCGTTGGAGAAATACAATAGTTCTTTACTCGCGCGTTGGAGAAATACAATAGT";

    #[test]
    fn test_filter() {
        let alignments: paf::PAF = paf::PAF::from_str(PAF_STRING);
        let (query_index, target_index): (Index, Index) = index::index_paf_matches(&alignments);

        let penalties = types::Penalties {
            mismatch: 4,
            matches: 0,
            gap_open: 6,
            gap_extend: 2,
        };

        let args = CliArgs {
            verbosity_level: 0,
            input_paf: String::new(),
            target_fasta: String::new(),
            query_fasta: String::new(),
            penalties: penalties,
            adapt: false,
        };

        let text = Fasta::from_str(TEXT);
        let query = Fasta::from_str(QUERY);

        filter(&text, &target_index, &query, &query_index, &args);
    }
}
