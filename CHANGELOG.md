# Change Log
All notable changes to this project will be documented in this file.
This project adheres to [Semantic Versioning](http://semver.org/).

## [0.3.2] - 2022-06-07
### Changed
- Update dependencies

## [0.3.1] - 2022-01-26
### Changed
- Compiling problems with recent clap versions ([#8](https://github.com/fastqc-rs/fastqc-rs/pull/8)).
- Removed `Cargo.lock` from `.gitignore`.

## [0.3.0] - 2021-11-19
### Changed
- Major performance upgrade by using a more performant hash function ([#6](https://github.com/fastqc-rs/fastqc-rs/pull/6)). 

## [0.2.2] - 2021-05-14
### Changed
- Minor cosmetic changes and transferring fastqc-rs to its own organization.

## [0.2.1] - 2021-03-18
### Changed
- Fix for the [MultiQC](https://multiqc.info) summary.

## [0.2.0] - 2021-02-15
### Changed
- Various fixes for the [MultiQC](https://multiqc.info) summary.
- Removed default value for `-s`.

## [0.1.1] - 2021-02-14
### Added
- New parameter `-s` that allows the user to generate summary files for usage with [MultiQC](https://multiqc.info).

## [0.1.0] - 2021-02-12
Initial release of fastqc-rs including processing FASTQ files and generating a self contained html report with visualizations.
