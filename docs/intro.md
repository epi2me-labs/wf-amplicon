## Introduction

This [Nextflow](https://www.nextflow.io/) workflow provides a simple way to
analyse Oxford Nanopore reads generated from amplicons.

The workflow requires raw reads in FASTQ format and can be run in two modes:
* Variant calling mode: Trigger this mode by passing a reference FASTA file.
  After initial filtering (based on read length and quality) and adapter
  trimming, [minimap2](https://github.com/lh3/minimap2) is used to align the
  reads to the reference (please note that the reference should only contain the
  expected sequences of the individual amplicons). Variants are then called with
  [Medaka](https://github.com/nanoporetech/medaka). This mode allows for
  multiple amplicons per barcode (for details on how to map specific target
  amplicons to individual samples / barcodes, see below).
* De-novo consensus mode: This mode is run when no reference file is passed.
  Like for the "variant calling mode", reads are first filtered and trimmed.
  Then, a consensus sequence is generated _de novo_ from the reads of each
  sample. Reads are re-aligned against the draft consensus for polishing with
  [Medaka](https://github.com/nanoporetech/medaka). Please note that only one
  amplicon per barcode is supported in de-novo consensus mode.

The results of the workflow include an interactive HTML report, FASTA files with
the consensus sequences of the amplicons, BAM files with the alignments, and VCF
files containing the variants (if run in variant calling mode).
