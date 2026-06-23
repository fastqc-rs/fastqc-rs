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

### Basic Usage

Generate an HTML report (output to stdout):

```bash
fqc -q path/to/my_sequence.fastq > report.html
```

### Command-Line Arguments

| Parameter        | Default | Description |
| :--------------- | :-----: | :---------- |
| `-q --fastq`     | -       | **(Required)** The path to the FASTQ file to analyze (supports `.fastq`, `.fastq.gz`, `.fq`, `.fq.gz`) |
| `-k --kmer`      | 5       | The length k of k-mers for k-mer counting (detects over-represented sequences) |
| `-s --summary`   | -       | Output directory for MultiQC summary file (`fastqc_data.txt` will be created). If the directory doesn't exist, it will be created automatically. |
| `--no-html`      | false   | Skip HTML report output to stdout (useful when network is unavailable or only MultiQC summary is needed) |

### Examples

#### 1. Basic HTML report

```bash
fqc -q sample.fastq > report.html
```

#### 2. Compressed FASTQ input

```bash
# Automatically handles .gz compressed files
fqc -q sample.fastq.gz > report.html
```

#### 3. Custom k-mer length

```bash
# Use k=7 for more specific sequence detection
fqc -q sample.fastq -k 7 > report.html
```

#### 4. Generate MultiQC summary

```bash
# Creates fastqc_data.txt in the specified directory
# If the directory doesn't exist, it will be created automatically
fqc -q sample.fastq -s output_dir/
```

#### 5. Generate MultiQC summary only (skip HTML)

```bash
# Skip HTML output when network is unavailable
# Useful for batch processing without needing external JavaScript
fqc -q sample.fastq -s output_dir/ --no-html
```

#### 6. Batch processing multiple files

```bash
# Process all FASTQ files in a directory
for file in *.fastq.gz; do
    name=$(basename "$file" .fastq.gz)
    fqc -q "$file" -s "output/$name/" > "reports/$name.html"
done
```

#### 7. Parallel batch processing

```bash
# Using GNU parallel for faster processing
ls *.fastq.gz | parallel -j 4 'fqc -q {} -s output/{/.}/ > reports/{/.}.html'
```

### Output

- **HTML Report**: Written to stdout, contains interactive visualizations for all statistics
- **MultiQC Summary** (`-s`): Creates `fastqc_data.txt` compatible with [MultiQC](https://multiqc.info) for aggregating multiple samples
