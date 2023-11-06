# wf-amplicon






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




## Quickstart

The workflow relies on the following dependencies:

- [Nextflow](https://www.nextflow.io/) for managing compute and software
  resources.
- Either [Docker](https://www.docker.com/products/docker-desktop) or
  [Singularity](https://docs.sylabs.io/guides/latest/user-guide/) to provide
  isolation of the required software.

It is not necessary to clone or download the git repository in order to run the
workflow. For more information on running EPI2ME Labs workflows, [visit our
website](https://labs.epi2me.io/wfindex).

**Workflow options**

If you have [Nextflow](https://www.nextflow.io/) installed, you can run the
following to obtain the workflow:

```
nextflow run epi2me-labs/wf-amplicon --help
```

This will show the workflow's command line options.

### Usage

This section covers how to use the workflow when providing a reference file with
the expected sequences of the amplicons. For details on how to run the workflow
without a reference (i.e. for _de novo_ construction of the amplicon consensus
sequence), see the relevant section below.

Basic usage is as follows

```bash
nextflow run epi2me-labs/wf-amplicon \
    --fastq $input \
    --reference references.fa \
    --threads 4
```

`$input` can be a single FASTQ file, a directory containing FASTQ files, or a
directory containing barcoded sub-directories which in turn contain FASTQ files.
A sample sheet can be included with `--sample_sheet` and a sample name for an
individual sample with `--sample`. If a sample sheet is provided, it needs to
contain at least these three columns:

- `barcode`: names of the sub-directories with the FASTQ files of the different
  samples
- `alias`: sample names
- `type`: sample types (needs to be one of `test_sample`, `positive_control`,
  `negative_control`, `no_template_control`)

Additionally, it can also contain a column named `ref`, which can be used to
specify one or more sequences in the reference FASTA file for the individual
samples. If a sample should be aligned against multiple reference sequences,
their IDs can be included as a space-delimited list. If a specific sample should
be aligned against all references, you can leave the corresponding cell in the
`ref` column blank. As an example, the following sample sheet tells the workflow
that the first two barcodes (`sample1` and `sample2`) are expected to contain
reads from both amplicons in the reference, whereas the other two samples only
contain one amplicon each.

```
barcode,alias,type,ref
barcode01,sample1,test_sample,katG::NC_000962.3:2154725-2155670 rpoB::NC_000962.3:760285-761376
barcode02,sample2,test_sample,
barcode03,sample3,positive_control,katG::NC_000962.3:2154725-2155670
barcode04,sample4,test_sample,rpoB::NC_000962.3:760285-761376
```

Relevant options for filtering of raw reads are

- `--min_read_length`
- `--max_read_length`
- `--min_read_qual`

After filtering and trimming with
[Porechop](https://github.com/rrwick/Porechop), reads can optionally be
downsampled. You can control the number of reads to keep per sample with
`--reads_downsampling_size`. Samples with fewer than `--min_n_reads` are
ignored.

After alignment, haploid variants are called with
[Medaka](https://github.com/nanoporetech/medaka). You can set the minimum
coverage a variant needs to exceed in order to be included in the results with
`--min_coverage`. Variants with lower coverage will still be included in the
resulting VCF files, but with `LOW_DEPTH` instead of `PASS` in the `FILTER`
column.

The workflow selects the appropriate
[Medaka models](https://github.com/nanoporetech/medaka#models) based on the basecaller
configuration that was used to process the signal data. You can use the
parameter `--basecaller_cfg` to provide this information (e.g.
`dna_r10.4.1_e8.2_400bps_hac`). Alternatively, you can choose the
[Medaka](https://github.com/nanoporetech/medaka) model directly with
`--medaka_model`.

If you want to use
[Singularity](https://docs.sylabs.io/guides/latest/user-guide/) instead of
[Docker](https://www.docker.com/products/docker-desktop), add `-profile
singularity` to the nextflow invocation.

### Key outputs

- Interactive HTML report detailing the results.
- BAM files with reads aligned against the reference.
- VCF files with variants called by
  [Medaka](https://github.com/nanoporetech/medaka).
- FASTA files with consensus sequences generated by applying the variants to the
  provided reference.

A sub-directory with these output files will be created for each sample.
Additionally, a single BAM and VCF file containing the results of all samples
will be produced if the option `--combine_results` is provided.

### Known issues and limitations

The consensus sequence is generated by applying the variants found by
[Medaka](https://github.com/nanoporetech/medaka) to the provided reference. This
comes with the caveat that deletions at either end of the amplicon are not going
to be reflected in the consensus (i.e. it will still contain the deleted
regions). For this reason we suggest to always have a look at the coverage plots
in the generated HTML report.

### Running without a reference

The workflow can also be run without a reference. It will then use
`medaka smolecule` to construct a consensus sequence for each barcode.
`smolecule` relies on [SPOA](https://github.com/rvaser/spoa) and the workflow
takes advantage of the recently added option in SPOA to prune low-coverage nodes
from the partial order graph. The threshold for this is relative to the number
of reads per sample (after filtering) and can be modified with
`--spoa_minimum_relative_coverage` (e.g. set to `0.2` to require a coverage of at least
20% of reads along the whole consensus).

[SPOA](https://github.com/rvaser/spoa) works best for LSK data. The performance
of the consensus generation may be poor for RBK data in some cases, especially
if no long reads spanning most of the amplicon were found. We therefore
recommend to only use "no reference mode" for short amplicons (4 kb or shorter).
To improve robustness on RBK, an assembly-based approach might be implemented in
the future.





## Useful links

- [Nextflow](https://www.nextflow.io/)
- [Docker](https://www.docker.com/products/docker-desktop)
- [Singularity](https://docs.sylabs.io/guides/latest/user-guide/)
- [minimap2](https://github.com/lh3/minimap2)
- [Medaka](https://github.com/nanoporetech/medaka)
- [Porechop](https://github.com/rrwick/Porechop)
