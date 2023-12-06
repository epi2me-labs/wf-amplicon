> Note: The initial preprocessing steps are the same for the variant calling and the de-novo consensus mode.

#### 1. Read preprocessing

The input data can be a single FASTQ file, a directory containing FASTQ files, or a directory containing barcoded sub-directories which in turn contain FASTQ files.
A sample sheet can be included with `--sample_sheet` and a sample name for an individual sample with `--sample`.
If a sample sheet is provided, it needs to contain at least these three columns:

- `barcode`: names of the sub-directories with the FASTQ files of the different
  samples
- `alias`: sample names
- `type`: sample types (needs to be one of `test_sample`, `positive_control`,
  `negative_control`, `no_template_control`)

After parsing the sample sheet, the raw reads are filtered. Relevant options for filtering are

- `--min_read_length`
- `--max_read_length`
- `--min_read_qual`

Reads can optionally also be downsampled (`--reads_downsampling_size` controls the number of reads to keep per sample).
The workflow supports random downsampling as well as selecting reads based on read length.
Subsetting based on read length tends to work better for de-novo consensus generation and usually does not decrease performance of the variant calling mode.
It is thus set as default (see the FAQ section below for details and for how to disable this).
The selected reads are then trimmed with [Porechop](https://github.com/rrwick/Porechop) prior to downstream analysis.

> Note: Samples with fewer reads than `--min_n_reads` after preprocessing are ignored.


### Variant calling mode

#### 2. Align reads

Preprocessed reads are aligned against the provided reference with [Minimap2](https://github.com/lh3/minimap2) (please note that the reference should only contain the expected sequences of the individual amplicons).
[Bamstats](https://github.com/epi2me-labs/fastcat#bamstats) is used to collate alignment statistics and [mosdepth](https://github.com/brentp/mosdepth) to calculate the depth of coverage along the amplicon.

#### 3. Call variants

After alignment, haploid variants are called with [Medaka](https://github.com/nanoporetech/medaka).
You can set the minimum coverage a variant needs to exceed in order to be included in the results with `--min_coverage`.
Variants with lower coverage will still be listed in the resulting VCF files, but with `LOW_DEPTH` instead of `PASS` in the `FILTER` column.

The workflow selects the appropriate [Medaka models](https://github.com/nanoporetech/medaka#models) based on the basecaller configuration that was used to process the signal data.
You can use the parameter `--basecaller_cfg` to provide this information (e.g. `dna_r10.4.1_e8.2_400bps_hac`).
Alternatively, you can choose the [Medaka](https://github.com/nanoporetech/medaka) model directly with `--medaka_model`.

#### 4. Use the variants to generate a consensus

The variants passing the depth filter are then incorporated in the reference to create the consensus sequence.

### De-novo consensus mode

#### 2. Create a draft consensus

When running in de-novo mode, [miniasm](https://github.com/lh3/miniasm) or [SPOA](https://github.com/rvaser/spoa) are used to create a draft consensus after the reads were selected and trimmed.
`miniasm` tends to work better for longer amplicons and `spoa` for shorter amplicons.
The workflow attempts to create an assembly with `miniasm` first.
If this succeeds, the draft assembly is polished with [Racon](https://github.com/lbcb-sci/racon) and the input reads are re-aligned against it.
This is followed by some post-processing (trimming low-coverage regions from the ends) and quality control (QC).
During QC, the following sanity-checks are performed:

- The mean coverage of the consensus needs to be above `--minimum_mean_depth` (default: 30X).
- The fraction of primary alignments for the re-aligned reads needs to be more than
  `--primary_alignment_threshold` (default: 70%). Large numbers of secondary /
  supplementary alignments indicate that something might have gone wrong during
  assembly. If your amplicons contain long repetitive regions, you can lower
  this number.

If more than one contig from the draft assembly passes these filters, the one with the highest coverage is selected.

If assembly generation fails or none of the contigs produced by `miniasm` pass QC, `spoa` is tried.
Regardless of how the draft consensus was generated (`miniasm` or `spoa`), a final polishing step with
[Medaka](https://github.com/nanoporetech/medaka) is performed.

> Note: Since `spoa` tends to produce a better consensus than `miniasm` for short amplicons, users can force the workflow to run `spoa` if the assembly created by `miniasm` is shorter than `--force_spoa_length_threshold` (even if it passes QC).
