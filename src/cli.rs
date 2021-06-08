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
                .short("f")
                .long("input-paf")
                .value_name("FILE")
                .help("Path to input PAF file")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("debug")
                .short("d")
                .long("debug")
                .help("Show debug statements [Default: false]"),
        )
        .get_matches();

    // Gets a value for config if supplied by user, or defaults to "default.conf"
    let paf_file_path: &str = matches.value_of("input_paf").unwrap();
    let debug: bool = matches.is_present("debug");

    types::CliArgs::new(debug, paf_file_path)
}
