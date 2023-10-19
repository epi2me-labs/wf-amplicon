## Introduction

This [Nextflow](https://www.nextflow.io/) workflow provides a simple way to
analyse Oxford Nanopore reads generated from amplicons.

The workflow requires raw reads in FASTQ format and can be run in two modes:
* With a reference: After initial filtering (based on read length and quality)
  and adapter trimming, [minimap2](https://github.com/lh3/minimap2) is used to
  align the reads to a reference FASTA file (please note that the reference
  should only contain the expected sequences of the individual amplicons).
  Variants are then called with
  [Medaka](https://github.com/nanoporetech/medaka).
* Without a reference: Like for the "reference mode", reads are first filtered
  and trimmed. Then, `medaka smolecule` is used to generate the consensus of the
  reads of each sample _de novo_. Reads are re-aligned against the consensus to
  produce coverage plots for the report.

The results of the workflow include an interactive HTML report, FASTA files with
the consensus sequences of the amplicons, BAM files with the alignments, and VCF
files containing the variants (if run in "reference mode").
