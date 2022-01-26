mod process;

use clap::{App, Arg};
use std::error::Error;
use std::str::FromStr;

fn main() -> Result<(), Box<dyn Error>> {
    let matches = App::new("fastqc-rs")
        .version("0.3.1-alpha.1")
        .author("Felix W. <fxwiegand@wgdnet.de>")
        .about("A quality control tool for FASTQ files written in rust")
        .arg(
            Arg::new("fastq")
                .help("The input FASTQ file to use.")
                .takes_value(true)
                .required(true)
                .short('q')
                .long("fastq"),
        )
        .arg(
            Arg::new("k")
                .help("The length k of k-mers for k-mer counting.")
                .takes_value(true)
                .required(false)
                .default_value("5")
                .short('k'),
        )
        .arg(
            Arg::new("summary")
                .help("Creates an output file for usage with MultiQC under the given path.")
                .takes_value(true)
                .required(false)
                .short('s'),
        )
        .get_matches();

    let k = u8::from_str(matches.value_of("k").unwrap())?;
    let summary = matches.value_of("summary");
    crate::process::process(matches.value_of("fastq").unwrap(), k, summary)
}
