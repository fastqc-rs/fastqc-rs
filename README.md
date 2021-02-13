# fastqc-rs

![Rust](https://github.com/fxwiegand/fastqc-rs/workflows/Rust/badge.svg)
[![Crates.io](https://img.shields.io/crates/d/fastqc-rs.svg)](https://crates.io/crates/fastqc-rs)
[![Crates.io](https://img.shields.io/crates/v/fastqc-rs.svg)](https://crates.io/crates/fastqc-rs)
[![Crates.io](https://img.shields.io/crates/l/fastqc-rs.svg)](https://crates.io/crates/fastqc-rs)

A fast quality control tool for FASTQ files written in rust inspired by [fastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).

Available statistics are:
- Read length
- Sequence quality score
- Sequence quality per base
- Sequence content per base
- k-mer content
- GC content

## Usage

```
fqc -q path/to/my_sequence.fastq > report.html
```

Arguments: 

| Parameter                 | Default       | Description   |	
| :------------------------ |:-------------:| :-------------|
| -q --fastq 	       |	-           |The path to the FASTQ file to use
| -k          | 5           |The length k of k-mers for k-mer counting
| -s --summary          | .           |Creates an output file for usage with MultiQC under the given path
