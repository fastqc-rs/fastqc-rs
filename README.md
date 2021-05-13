# <img src="./img/fastqc-rs-ferris.svg" width=100em alt="fastqc-rs logo" /> fastqc-rs

![Rust](https://github.com/fxwiegand/fastqc-rs/workflows/Rust/badge.svg)
[![Crates.io](https://img.shields.io/crates/d/fastqc-rs.svg?label=crates.io%20downloads)](https://crates.io/crates/fastqc-rs)
[![Crates.io](https://img.shields.io/crates/v/fastqc-rs.svg)](https://crates.io/crates/fastqc-rs)
[![Crates.io](https://img.shields.io/crates/l/fastqc-rs.svg)](https://crates.io/crates/fastqc-rs)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/fastqc-rs/README.html)
[![Bioconda downloads](https://img.shields.io/conda/dn/bioconda/fastqc-rs?label=bioconda%20downloads)](https://anaconda.org/bioconda/fastqc-rs)

A fast quality control tool for FASTQ files written in rust inspired by [fastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/). Results are written to `stdout` as a self containing html report with visualizations for all statistics. Summary files for usage with [MultiQC](https://multiqc.info) can also be generated.

Available statistics are:
- Read length
- Sequence quality score
- Sequence quality per base
- Sequence content per base
- k-mer content
- GC content

For a detailed list of changes, take a look at the [CHANGELOG](CHANGELOG.md).

## Installation

There are multiple ways to install fastqc-rs:

#### Bioconda

fastqc-rs is available via [Bioconda](https://bioconda.github.io).
With Bioconda set up, installation is as easy as

    conda install fastqc-rs

#### Cargo

If the [Rust](https://www.rust-lang.org/tools/install) compiler and associated [Cargo](https://github.com/rust-lang/cargo/) are installed, fastqc-rs can be installed via

    cargo install fastqc-rs

#### Source

Download the source code and within the root directory of source run

    cargo install

## Usage

```
fqc -q path/to/my_sequence.fastq > report.html
```

Arguments: 

| Parameter                 | Default       | Description   |	
| :------------------------ |:-------------:| :-------------|
| -q --fastq 	       |	-           |The path to the FASTQ file to use
| -k          | 5           |The length k of k-mers for k-mer counting
| -s --summary          | -           |Creates an output file for usage with [MultiQC](https://multiqc.info) under the given path
