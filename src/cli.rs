use clap::{App, Arg};
use std::env;

use crate::types;

// Env vars
const NAME: &str = env!("CARGO_PKG_NAME");
const VERSION: &str = env!("CARGO_PKG_VERSION");
const DESCRIPTION: &str = env!("CARGO_PKG_DESCRIPTION");
const AUTHORS: &str = env!("CARGO_PKG_AUTHORS");

pub fn start() -> types::CliArgs {
    let matches = App::new(NAME)
        .version(VERSION)
        .author(AUTHORS)
        .about(DESCRIPTION)
        .arg(
            Arg::with_name("input_paf")
                .short("p")
                .long("paf")
                .value_name("FILE")
                .help("Path to input PAF file")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("target_fasta")
                .short("t")
                .long("target")
                .value_name("FILE")
                .help("Path to input target fasta file")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("query_fasta")
                .short("q")
                .long("query")
                .value_name("FILE")
                .help("Path to input target fasta file")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("adapt")
                .short("a")
                .long("adapt")
                .multiple(false)
                .help("To apply adaptive wavefront alignment [Default: false]"),
        )
        .arg(
            Arg::with_name("alignment_matrix")
                .short("m")
                .long("alignment-matrix")
                .multiple(false)
                .help("Generate alignment matrix TSV file `wfilter-<no_adapt/adapt>-<now>.tsv` [Default: false]"),
        )
        .arg(
            Arg::with_name("v")
                .short("v")
                .multiple(true)
                .help("Sets the level of verbosity [Default: 0]"),
        )
        .get_matches();

    // Gets a value for config if supplied by user, or defaults to "default.conf"
    let paf_file_path: &str = matches.value_of("input_paf").unwrap();
    let target_file_path: &str = matches.value_of("target_fasta").unwrap();
    let query_file_path: &str = matches.value_of("query_fasta").unwrap();
    let adapt: bool = matches.is_present("adapt");
    let generate_alignment_tsv: bool = matches.is_present("alignment_matrix");
    let verbosity_level: u8 = matches.occurrences_of("v") as u8;

    types::CliArgs::new(
        verbosity_level,
        paf_file_path,
        target_file_path,
        query_file_path,
        None, // TODO: implement penalties
        adapt,
        generate_alignment_tsv
    )
}
