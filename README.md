# wf-amplicon






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
the expected sequences of the amplicons ("variant calling mode"). For details on
how to run the workflow without a reference ("de-novo consensus mode"), see the
relevant section below.

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

After initial filtering, reads can optionally be downsampled
(`--reads_downsampling_size` controls the number of reads to keep per sample).
The workflow supports random downsampling as well as selecting reads based on
read length. Subsetting based on read length tends to work better for de-novo
consensus generation and is thus set as default (see the section for the de-novo
consensus mode below for details and for how to disable this). However, it
should not detriment the performance of the variant-calling mode. The selected
reads are then trimmed with [Porechop](https://github.com/rrwick/Porechop)
before being aligned against the reference.

> Note: Samples with fewer reads than `--min_n_reads` after preprocessing are
> ignored.

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
[miniasm](https://github.com/lh3/miniasm) or
[SPOA](https://github.com/rvaser/spoa) to create a draft consensus. `miniasm`
tends to work better for longer amplicons and `spoa` for shorter amplicons. The
workflow attempts to create an assembly with `miniasm` first. If this succeeds,
the draft assembly is polished with [Racon](https://github.com/lbcb-sci/racon) and
the input reads are re-aligned against it. This is followed by some post-processing
(trimming low-coverage regions from the ends) and quality control (QC). During QC, the following
sanity-checks are performed:

- The mean coverage of the consensus is above `--minimum_mean_depth` (default: 30X).
- The fraction of primary alignments for the re-aligned reads is more than
  `--primary_alignment_threshold` (default: 70%). Large numbers of secondary /
  supplementary alignments indicate that something might have gone wrong during
  assembly. If your amplicons contain long repetitive regions, you can lower
  this number.

If more than one contig from the draft assembly passes these filters, the one
with the highest coverage is selected.

If assembly generation fails or none of the contigs produced by `miniasm` pass
QC, `spoa` is tried. Regardless of how the draft consensus was generated
(`miniasm` or `spoa`), a final polishing step with
[Medaka](https://github.com/nanoporetech/medaka) is performed.

Note: Since `spoa` tends to produce a better consensus than `miniasm` for short
amplicons, users can force the workflow to run `spoa` if the assembly created by
`miniasm` is shorter than `--force_spoa_length_threshold` (even if it passes
QC).

**Random downsampling vs selecting reads by length:**

Despite the workflow dropping reads containing adapter sequences in the middle
(with `porechop --discard_middle`), some reads in the input data stemming from
concatemers or other PCR artifacts might still make it through preprocessing and
could thus interfere with consensus generation. A simple way to avoid this is to
drop a small fraction (e.g. 5%) of longest reads. Additionally, `spoa` depends
on the order of reads in the input and benefits from "seeing" longer reads
first. Therefore, the following options are used per default
`--drop_frac_longest_reads 0.05 --take_longest_remaining_reads`. To disable this
and enforce random downsampling of input reads, use `--drop_frac_longest_reads 0
--take_longest_remaining_reads false`.




## Useful links

- [Nextflow](https://www.nextflow.io/)
- [Docker](https://www.docker.com/products/docker-desktop)
- [Singularity](https://docs.sylabs.io/guides/latest/user-guide/)
- [minimap2](https://github.com/lh3/minimap2)
- [miniasm](https://github.com/lh3/miniasm)
- [SPOA](https://github.com/rvaser/spoa)
- [Medaka](https://github.com/nanoporetech/medaka)
- [Porechop](https://github.com/rrwick/Porechop)
