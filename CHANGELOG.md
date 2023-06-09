# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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

