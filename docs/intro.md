## Introduction

This [Nextflow](https://www.nextflow.io/) workflow provides a simple way to
analyse Oxford Nanopore reads generated from amplicons.<br>

It requires the raw reads and a reference FASTA file containing one sequence per
amplicon. After filtering (based on read length and quality) and trimming,
reads are aligned to the reference using
[minimap2](https://github.com/lh3/minimap2). Variants are then called with
[Medaka](https://github.com/nanoporetech/medaka). Results include an interactive
HTML report and VCF files containing the called variants.<br>

As mentioned above, the reference FASTA file needs to contain one sequence per
amplicon for now. An option to provide a whole-genome reference file and pairs
of primers might be added in the future if requested by users.
