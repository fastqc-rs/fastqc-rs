# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

fastqc-rs is a Rust-based FASTQ quality control tool inspired by fastQC. It analyzes sequencing files and generates HTML reports with visualizations and optional MultiQC-compatible summary files.

## Build and Development Commands

```bash
# Build the project
cargo build

# Build in release mode (with LTO enabled)
cargo build --release

# Run the binary
cargo run -- -q path/to/file.fastq > report.html

# Run tests
cargo test

# Run a specific test
cargo test test_report
cargo test test_quartiles1

# Install locally
cargo install
```

## Running the Tool

```bash
# Basic usage (generates HTML report to stdout)
fqc -q path/to/my_sequence.fastq > report.html

# With custom k-mer length and MultiQC summary
fqc -q input.fastq -k 7 -s summary_dir/
```

Arguments:
- `-q/--fastq`: Input FASTQ file (required)
- `-k/--kmer`: K-mer length for k-mer counting (default: 5)
- `-s/--summary`: Output directory for MultiQC summary file

## Code Architecture

### Entry Point
- **src/main.rs**: CLI argument parsing using clap, forwards to `process::process()`

### Core Processing
- **src/process.rs**: Single-file module containing all analysis logic
  - `process()`: Main function that orchestrates the entire pipeline
  - Uses `needletail` for FASTQ parsing
  - Collects statistics in HashMaps: read lengths, base counts, quality scores, k-mers, GC content
  - Computes quartiles using custom `quartiles()` function for box plots

### Report Generation
The project uses embedded files (via `include_str!()`) in src/report/:
- **Vega-Lite JSON specs**: Define visualization charts (base_per_pos_specs.json, counter_specs.json, etc.)
  - These are loaded, mutated with data, and serialized back to JSON
- **Tera templates**: report.html.tera and fastqc_summary.txt.tera
  - Custom filter `embed_source` fetches external JavaScript dependencies

### Analysis Pipeline Flow
1. Parse FASTQ file and collect raw statistics
2. Compute derived metrics (percentages, averages, quartiles)
3. Determine warning levels (pass/warn/fail) for each metric
4. Populate Vega-Lite specs with collected data
5. Render HTML report via Tera template to stdout
6. Optionally generate MultiQC summary file

### Data Structures
- `rustc_hash::FxHashMap` used throughout for performance (aka `HashMap` in code)
- Position-based data keyed by `usize` (base position in read)
- Quality histograms: `Vec<usize>` of length 94 (Phred scores 0-93)

## Test Structure

- **tests/lib.rs**: Integration test that runs the binary and compares HTML output
- **tests/resources/example.fastq**: Test input file
- **tests/expected/report.html**: Expected output (with certain lines excluded from comparison)
- Unit tests in process.rs for `quartiles()` function

## Dependencies

- **needletail**: FASTQ parsing
- **clap**: CLI argument parsing
- **tera**: HTML template rendering
- **serde_json**: JSON manipulation for Vega-Lite specs
- **rustc-hash**: Fast hasher for HashMaps
- **itertools**: Iterator utilities
