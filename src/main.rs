use coitrees;
use std::time::Instant;

// local
mod wflambda;

mod cli;
mod index;
mod io;
mod paf;
mod fasta;
mod types;

fn filter(target: &fasta::FastaFile,
          target_index: &types::Index,
          query: &fasta::FastaFile,
          query_index: &types::Index,
          cli_args: &types::CliArgs) {
    let verbosity = cli_args.verbosity_level;

    let backtrace_lambda = |query: (i32, i32), target: (i32, i32)| -> bool {
        let mut oq: Vec<usize> = Vec::new();
        let y = |i: &coitrees::IntervalNode<types::AlignmentMetadata, u32>| { oq.push(i.metadata) };
        target_index.query(target.0, target.1, y);
        let res = target_index.query_count(target.0, target.1) > 0 || query_index.query_count(query.0, query.1) > 0;

        if res {
            println!("{} {:?}", res, oq);
        }

        res
    };

    let now = Instant::now();
    for t in target.iter() {
        for q in query.iter() {
            if verbosity > 1 {
                eprintln!("\t[wfilter::main::align] Aligning {} and {}",
                          std::str::from_utf8(&t.header[..]).unwrap(),
                          std::str::from_utf8(&q.header[..]).unwrap()
                );
            };

            let aln = wflambda::wfa::wf_align(&t.seq[..], &q.seq[..], cli_args, &backtrace_lambda);

            if verbosity > 3 {
                eprintln!("score {}", aln.score);
                eprintln!("{}", aln.cigar);
            }
        }
    }

    if verbosity > 0 {
        eprintln!(
            "\t[wfilter::main::align] finished all alignments. Time taken {} seconds",
            now.elapsed().as_millis() as f64 / 1000.0
        )
    };
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
    if verbosity  > 0 {
        eprintln!(
            "[wfilter::main] done parsing target. Time taken {} seconds",
            now.elapsed().as_millis() as f64 / 1000.0
        );
        eprintln!("Number of sequences in target {}", target.len());
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
        eprintln!("Number of sequences in query {}", query.len());
    }

    if verbosity > 0 {
        eprintln!("[wfilter::main] filter");
    }

    filter(&target, &target_index, &query, &query_index, &args);

}

#[cfg(test)]
mod tests {
    use super::*;
    use super::paf;
    use super::index;
    use super::types::{Index, CliArgs};
    use super::fasta::Fasta;

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

        let penalties =  types::Penalties {
            mismatch: 4,
            matches: 0,
            gap_open: 6,
            gap_extend: 2,
        };

        let args = CliArgs::new(2, "", "", "", Some(penalties));
        let text = Fasta::from_str(TEXT);
        let query = Fasta::from_str(QUERY);

        filter(&text, &target_index, &query, &query_index, &args);
    }
}
