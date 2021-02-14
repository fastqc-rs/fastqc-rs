mod process;

use clap::{App, Arg};
use std::error::Error;
use std::str::FromStr;

fn main() -> Result<(), Box<dyn Error>> {
    let matches = App::new("fastqc-rs")
        .version("0.1")
        .author("Felix W. <fxwiegand@wgdnet.de>")
        .about("A quality control tool for FASTQ files written in rust")
        .arg(
            Arg::new("fastq")
                .about("The input FASTQ file to use.")
                .takes_value(true)
                .required(true)
                .short('q')
                .long("fastq"),
        )
        .arg(
            Arg::new("k")
                .about("The length k of k-mers for k-mer counting.")
                .takes_value(true)
                .required(false)
                .default_value("5")
                .short('k'),
        )
        .arg(
            Arg::new("summary")
                .about("Creates an output file for usage with MultiQC under the given path. The default path is the current working directory.")
                .takes_value(true)
                .required(false)
                .default_value(".")
                .short('s'),
        )
        .get_matches();

    let k = u8::from_str(matches.value_of("k").unwrap())?;
    let summary = matches.value_of("summary");
    crate::process::process(matches.value_of("fastq").unwrap(), k, summary)
}
