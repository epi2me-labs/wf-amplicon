# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [v0.6.1]
### Changed
- More informative warnings for failed samples in the report intro section.

## [v0.6.0]
### Added
- Assembly step to de-novo consensus mode to handle longer amplicons. README was updated accordingly.

### Removed
- Default local executor CPU and RAM limits

## [v0.5.0]
### Changed
- Names of barcoded directories in the sample sheet now need to be of format `barcodeXY`.

### Added
- Support for specifying reference sequences for individual samples with an extra `"ref"` column in the sample sheet.

## [v0.4.1]
### Changed
- The workflow now also emits VCF index files as well as BAM / VCF index files for combined outputs when running with `--combine_results`.

### Fixed
- Emitting an empty consensus FASTA file when consensus generation failed in certain situations. Instead, no file is emitted.
- Now uses the downsampled BAM for `medaka annotate`.

### Removed
- Misleading allelic balance statistic in report.

## [v0.4.0]
### Added
- Running the workflow without a reference will switch the workflow to "no-reference mode" and will use SPOA to construct a consensus sequence _de novo_.

## [v0.3.5]
### Changed
- The workflow now also outputs BAM index files.

## [v0.3.4]
### Changed
- The workflow now downsamples reads for each amplicon to be in the suitable depth range for Medaka.

## [v0.3.3]
### Added
- MacOS ARM64 support.

## [v0.3.2]
### Changed
- Updated Medaka to 1.9.1.

## [v0.3.1]
### Changed
- No longer publishes empty result files (BAM, VCF, consensus FASTA) for samples which do not have any reads left after pre-processing and filtering.
- Now uses Medaka v1.8.2. Options for `basecaller_cfg` were updated accordingly. The default now is `dna_r10.4.1_e8.2_400bps_sup@v4.2.0`.
- The per-sample summary table in the report no longer shows sample metadata columns unless metadata was provided by the user via a sample sheet.

## [v0.3.0]
### Changed
- VCF files now use the sample alias as sample name instead of `SAMPLE`.
- The reference FASTA file with sanitized sequence headers (with `:`, `*`, and whitespace replaced with `_`) which is used by the workflow internally due to some tools not tolerating these symbols in the sequence IDs is now also published alongside the other results.

### Added
- Parameter `--combine_results` to also output merged BAM and VCF files.

## [v0.2.3]
### Fixed
- The workflow now prints a warning when there are no reads after filtering / preprocessing and produces a truncated report showing only the pre-processing stats table.

## [v0.2.2]
### Fixed
- Bug where the `mosdepth` process would fail if the length of a reference sequence was smaller than the number of depth windows requested.

## [v0.2.1]
### Changed
- GitHub issue templates
- Example command to use demo data.

## [v0.2.0]
### Changed
- Instead of dropping variants with `DP < min_coverage`, set their `FILTER` column to `LOW_DEPTH` in the results VCF.
- Bumped minimum required Nextflow version to 22.10.8.
- Enum choices are enumerated in the `--help` output.
- Enum choices are enumerated as part of the error message when a user has selected an invalid choice.

### Fixed
- Replaced `--threads` option in `fastqingress.nf` with hardcoded values to remove warning about undefined `param.threads`.

## [v0.1.3]
### Added
- Consensus sequences and BAM files as output files.

## [v0.1.2]
### Added
- Configuration for running demo data in AWS.
- Tags for `epi2melabs` desktop app.

## [v0.1.1]
### Fixed
- Reference not being required by the schema.
- Bug in documentation that prevented blog post from building.

## [v0.1.0]
* First release.

