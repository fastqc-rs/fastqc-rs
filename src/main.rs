mod process;

use clap::{Arg, Command};
use std::error::Error;
use env_logger::Builder;

pub fn init_log() -> u64 {
    Builder::from_default_env().init();
    println!("\n ************** initializing logger *****************\n");
    1
}

fn main() -> Result<(), Box<dyn Error>> {
    let _ = init_log();

    let matches = Command::new("fastqc-rs")
        .about("A FASTQ quality control tool inspired by fastQC")
        .version("0.3.3")
        .author("Felix W. <fxwiegand@wgdnet.de>")
        .arg(
            Arg::new("fastq")
                .short('q')
                .long("fastq")
                .value_name("FILE")
                .help("The input FASTQ file to use.")
                .required(true)
                .value_parser(clap::value_parser!(String))
        )
        .arg(
            Arg::new("k")
                .short('k')
                .long("kmer")
                .value_name("K")
                .help("The length k of k-mers for k-mer counting.")
                .default_value("5")
                .value_parser(clap::value_parser!(String))
        )
        .arg(
            Arg::new("summary")
                .short('s')
                .long("summary")
                .value_name("FILE")
                .help("Creates an output file for usage with MultiQC under the given path.")
                .value_parser(clap::value_parser!(String))
        )
        .get_matches();

    let fastq_file = matches.get_one::<String>("fastq").unwrap();
    let k: u8 = matches.get_one::<String>("k")
                       .unwrap()
                       .parse()
                       .expect("K should be a valid number");
    let summary = matches.get_one::<String>("summary");

    crate::process::process(fastq_file, k, summary)
}
