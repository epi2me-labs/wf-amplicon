# Amplicon workflow

Nextflow workflow for analysing Oxford Nanopore reads created by amplicon sequencing.



## Introduction

This [Nextflow](https://www.nextflow.io/) workflow provides a simple way to
analyse Oxford Nanopore reads generated from amplicons.

The workflow requires raw reads in FASTQ format and can be run in two modes:
* Variant calling mode: Trigger this mode by passing a reference FASTA file.
  After initial filtering (based on read length and quality) and adapter
  trimming, [minimap2](https://github.com/lh3/minimap2) is used to align the
  reads to the reference. Variants are then called with
  [Medaka](https://github.com/nanoporetech/medaka). This mode allows for
  multiple amplicons per barcode (before running with multiple amplicons per
  barcode, please read the relevant sections in the [Pipeline Overview](#pipeline-overview) and [FAQs](#FAQs) below).
* De-novo consensus mode: This mode is run when no reference file is passed.
  Like for the "variant calling mode", reads are first filtered and trimmed.
  Then, a consensus sequence is generated _de novo_ from the reads of each
  sample. Reads are re-aligned against the draft consensus for polishing with
  [Medaka](https://github.com/nanoporetech/medaka). Please note that only one
  amplicon per barcode is supported in de-novo consensus mode.

> Note: This workflow is *not* intended for marker gene sequencing of mixtures / communities of different organisms (e.g. 16S sequencing).
> In de-novo consensus mode it expects a single amplicon per barcode.
> When running in variant calling mode, multiple amplicons per barcode can be processed, but their sequences need to be sufficiently different from each other so that most reads only align to one of the provided references.




## Compute requirements

Recommended requirements:

+ CPUs = 12
+ Memory = 32GB

Minimum requirements:

+ CPUs = 6
+ Memory = 16GB

Approximate run time: 0.5-5 minutes per sample (depending on number of reads, length of amplicons, and available compute).

ARM processor support: True




## Install and run

These are instructions to install and run the workflow on command line. You can also access the workflow via the [EPI2ME application](https://labs.epi2me.io/downloads/).

The workflow uses [Nextflow](https://www.nextflow.io/) to manage compute and software resources, therefore Nextflow will need to be installed before attempting to run the workflow.

The workflow can currently be run using either [Docker](https://www.docker.com/products/docker-desktop) or
[Singularity](https://docs.sylabs.io/guides/3.0/user-guide/index.html) to provide isolation of
the required software. Both methods are automated out-of-the-box provided
either Socker or Singularity is installed. This is controlled by the [`-profile`](https://www.nextflow.io/docs/latest/config.html#config-profiles) parameter as exemplified below.

It is not required to clone or download the git repository in order to run the workflow.
More information on running EPI2ME workflows can be found on our [website](https://labs.epi2me.io/wfindex).

The following command can be used to obtain the workflow. This will pull the repository in to the assets folder of Nextflow and provide a list of all parameters available for the workflow as well as an example command:

```
nextflow run epi2me-labs/wf-amplicon –-help
```

A demo dataset is provided for testing of the workflow. It can be downloaded using:

```
wget https://ont-exd-int-s3-euwst1-epi2me-labs.s3.amazonaws.com/wf-amplicon/wf-amplicon-demo.tar.gz
tar -xzvf wf-amplicon-demo.tar.gz
```

The workflow can be run with the demo data using:

```
nextflow run epi2me-labs/wf-amplicon \
    --fastq wf-amplicon-demo/fastq \
    --reference wf-amplicon-demo/reference.fa \
    -profile standard
```

For further information about running a workflow on the command line see https://labs.epi2me.io/wfquickstart/




## Related protocols

This workflow is designed to take input sequences that have been produced from [Oxford Nanopore Technologies](https://nanoporetech.com/) devices.

Find related protocols in the [Nanopore community](https://community.nanoporetech.com/docs/).




## Input example

<!---Example of input directory structure, delete and edit as appropriate per workflow.--->
This workflow accepts FASTQ files as input.

The FASTQ input parameters for this workflow accept one of three cases: (i) the path to a single FASTQ file; (ii) the path to a top-level directory containing FASTQ files; (iii) the path to a directory containing one level of sub-directories which in turn contain FASTQ files. In the first and second cases (i and ii), a sample name can be supplied with `--sample`. In the last case (iii), the data is assumed to be multiplexed with the names of the sub-directories as barcodes. In this case, a sample sheet can be provided with `--sample_sheet`.

```
(i)                     (ii)                 (iii)
input_reads.fastq   ─── input_directory  ─── input_directory
                        ├── reads0.fastq     ├── barcode01
                        └── reads1.fastq     │   ├── reads0.fastq
                                             │   └── reads1.fastq
                                             ├── barcode02
                                             │   ├── reads0.fastq
                                             │   ├── reads1.fastq
                                             │   └── reads2.fastq
                                             └── barcode03
                                              └── reads0.fastq
```



## Input parameters

### Input Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| fastq | string | FASTQ files to use in the analysis. | This accepts one of three cases: (i) the path to a single FASTQ file; (ii) the path to a top-level directory containing FASTQ files; (iii) the path to a directory containing one level of sub-directories which in turn contain FASTQ files. In the first and second case, a sample name can be supplied with `--sample`. In the last case, the data is assumed to be multiplexed with the names of the sub-directories as barcodes. In this case, a sample sheet can be provided with `--sample_sheet`. |  |
| analyse_unclassified | boolean | Analyse unclassified reads from input directory. By default the workflow will not process reads in the unclassified directory. | If selected and if the input is a multiplex directory the workflow will also process the unclassified directory. | False |
| reads_downsampling_size | integer | Downsample to this number of reads per sample. | Downsampling is performed after filtering. Downsampling is performed by default to improve performance. Set to 0 to disable downsampling and use all data for analysis. Increasing the downsampling size (or setting it to zero) may lead to excessive memory consumption. | 1500 |
| reference | string | Path to a reference FASTA file. | The reference file should contain one sequence per amplicon. |  |


### Sample Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| sample_sheet | string | A CSV file used to map barcodes to sample aliases. The sample sheet can be provided when the input data is a directory containing sub-directories with FASTQ files. | The sample sheet is a CSV file with, minimally, columns named `barcode` and `alias`. Extra columns are allowed. In variant calling mode, a `ref` column can be added to tell the workflow which reference sequences should be used for which samples (please see the FAQs section in the documentation for details). |  |
| sample | string | A single sample name for non-multiplexed data. Permissible if passing a single .fastq(.gz) file or directory of .fastq(.gz) files. |  |  |


### Pre-processing Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| min_read_length | integer | Shorter reads will be removed. |  | 300 |
| max_read_length | integer | Longer reads will be removed. |  |  |
| min_read_qual | number | Reads with a lower mean quality will be removed. |  | 10 |
| drop_frac_longest_reads | number | Drop fraction of longest reads. | The very longest reads might be concatemers or contain other artifacts. In many cases removing them simplifies de novo consensus generation. When running variant calling mode with multiple amplicons per sample, it is recommended to set this to 0.0. | 0.05 |
| take_longest_remaining_reads | boolean | Whether to use the longest (remaining) reads. | With this option enabled, reads are not randomly selected during downsampling. Instead, after dropping the longest reads (unless the `drop_frac_longest_reads` parameter is set to 0) the longest remaining reads are kept. This is recommended for de-novo mode when working with long amplicons. When running variant calling mode with multiple amplicons per sample, it is recommended to disable this option. | True |
| min_n_reads | number | Samples / barcodes with fewer reads will not be processed. |  | 40 |


### Variant Calling Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| min_coverage | integer | Minimum coverage for variants to keep. | Only variants covered by more than this number of reads are reported in the resulting VCF file. | 20 |
| basecaller_cfg | string | Name of the basecaller model that processed the signal data; used to select an appropriate Medaka model. | The basecaller configuration is used to automatically select the appropriate Medaka model. The automatic selection can be overridden with the 'medaka_model' parameters. Available models are: 'dna_r10.4.1_e8.2_400bps_hac@v3.5.2', 'dna_r10.4.1_e8.2_400bps_sup@v3.5.2', 'dna_r9.4.1_e8_fast@v3.4', 'dna_r9.4.1_e8_hac@v3.3', 'dna_r9.4.1_e8_sup@v3.3', 'dna_r10.4.1_e8.2_400bps_hac_prom', 'dna_r9.4.1_450bps_hac_prom', 'dna_r10.3_450bps_hac', 'dna_r10.3_450bps_hac_prom', 'dna_r10.4.1_e8.2_260bps_hac', 'dna_r10.4.1_e8.2_260bps_hac_prom', 'dna_r10.4.1_e8.2_400bps_hac', 'dna_r9.4.1_450bps_hac', 'dna_r9.4.1_e8.1_hac', 'dna_r9.4.1_e8.1_hac_prom'. | dna_r10.4.1_e8.2_400bps_sup@v4.2.0 |
| medaka_model | string | The name of the Medaka model to use. This will override the model automatically chosen based on the provided basecaller configuration. | The workflow will attempt to map the basecaller model (provided with 'basecaller_cfg') used to a suitable Medaka model. You can override this by providing a model with this option instead. |  |


### De-novo Consensus Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| force_spoa_length_threshold | integer | Consensus length below which to force SPOA consensus generation. | If the consensus generated by `miniasm` is shorter than this value, force consensus generation with SPOA (regardless of whether the sequence produced by `miniasm` passed QC or not). The rationale for this parameter is that `miniasm` sometimes gives slightly truncated assemblies for short amplicons from RBK data, whereas SPOA tends to be more robust in this regard. | 2000 |
| spoa_minimum_relative_coverage | number | Minimum coverage (relative to the number of reads per sample after filtering) when constructing the consensus with SPOA. | Needs to be a number between 0.0 and 1.0. The result of multiplying this number with the number of reads for the corresponding sample (after filtering) is passed on to SPOA's `--min-coverage` option. | 0.15 |
| minimum_mean_depth | integer | Mean depth threshold to pass consensus quality control. | Draft consensus sequences with a lower average depth of coverage after re-aligning the input reads will fail QC. | 30 |
| primary_alignments_threshold | number | Fraction of primary alignments to pass quality control. | Draft consensus sequences with a lower fraction of primary alignments after re-aligning the input reads will fail QC. | 0.7 |


### Output Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| out_dir | string | Directory for output of all workflow results. |  | output |
| combine_results | boolean | Whether to merge per-sample results into a single BAM / VCF file. | Per default, results are grouped per sample. With this option, an additional BAM and VCF file are produced which contain the alignments / variants for all samples and amplicons. | False |


### Advanced Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| number_depth_windows | integer | Number of windows used during depth of coverage calculations. | Depth of coverage is calculated for each sample across each amplicon split into this number of windows. A higher number will produce more fine-grained plots at the expense of run time. | 100 |
| medaka_target_depth_per_strand | integer | Downsample each amplicon to this per-strand depth before running Medaka. | Medaka performs best with even strand coverage and depths between 80X and 400X. To avoid too high coverage, the workflow downsamples the reads for each amplicon to this per-strand depth before running Medaka. Changing this value is discouraged as it might cause decreased performance. | 150 |


### Miscellaneous Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| threads | integer | Maximum number of CPU threads to use per workflow task. | Several tasks in this workflow benefit from using multiple CPU threads. This option sets the maximum number of CPU threads for such processes. The total CPU resources used by the workflow are constrained by the Nextflow executor configuration. | 4 |
| disable_ping | boolean | Enable to prevent sending a workflow ping. |  | False |






## Outputs

Output files may be aggregated including information for all samples or provided per sample. Per-sample files will be prefixed with respective aliases and represented below as {{ alias }}.

| Title | File path | Description | Per sample or aggregated |
|-------|-----------|-------------|--------------------------|
| workflow report | ./wf-amplicon-report.html | Report for all samples. | aggregated |
| Sanitized reference file | ./reference_sanitized_seqIDs.fasta | Some programs used by the workflow don't like special characters (like colons) in the sequence IDs in the reference FASTA file. The reference is thus "sanitized" by replacing these characters with underscores. This file is only generated when the workflow is run in variant calling mode. | aggregated |
| Alignments BAM file | ./{{ alias }}/alignments/aligned.sorted.bam | BAM file with alignments of input reads against the references (in variant calling mode) or the created consensus (in de-novo consensus mode). | per-sample |
| Alignments index file | ./{{ alias }}/alignments/aligned.sorted.bam.bai | Index for alignments BAM file. | per-sample |
| De-novo consensus FASTQ file | ./{{ alias }}/consensus/consensus.fastq | Consensus file generated by de-novo consensus pipeline. | per-sample |
| Consensus FASTA file | ./{{ alias }}/consensus/consensus.fasta | Consensus file generated variant calling pipeline. | per-sample |
| Variants VCF file | ./{{ alias }}/variants/medaka.annotated.vcf.gz | VCF file of variants detected against the provided reference. | per-sample |
| Variants index file | ./{{ alias }}/variants/medaka.annotated.vcf.gz.csi | Index for variants VCF file. | per-sample |




## Pipeline overview

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

By default, reads are also downsampled.
The workflow supports random downsampling as well as selecting reads based on read length.
Subsetting based on read length tends to work better for de-novo consensus generation and does not decrease performance of variant calling mode when running with only one amplicon per sample.

- `--reads_downsampling_size`: This controls the number of reads to keep per sample
- `--drop_frac_longest_reads`: Drop a fraction (e.g. 0.05 for 5%) of the longest reads.
  As the very longest reads might stem from concatemers or other sequencing artifacts, they might interfere with consensus generation and it is generally better to remove them in de-novo consensus mode.
- `--take_longest_remaining_reads`: Use the longest reads instead of random downsampling; note that this is done after dropping the longest reads according to `--drop_frac_longest_reads`.
  This can be beneficial when assembling long amplicons in de-novo consensus mode in some cases.

The default values for these parameters are listed in the [Inputs section](#pre-processing-options).
The defaults work well for

- De-novo consensus mode.
- Variant calling mode with one amplicon per sample.

When running variant calling mode with more than one amplicon per sample we recommend changing the parameters to `--drop_frac_longest_reads 0 --take_longest_remaining_reads false` in order to perform random downsampling.
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




## Troubleshooting

+ If the workflow fails please run it with the demo data set to ensure the workflow itself is working. This will help us determine if the issue is related to the environment, input parameters or a bug.
+ Please see [here](https://labs.epi2me.io/trouble-shooting/) for how to resolve some common Nextflow issues and [here](https://labs.epi2me.io/how-to-exits/) for how to interpret command exit codes.




## FAQ's

*Do I need to change any parameters or can I stick to the defaults?* - The default parameters have been set for optimal results when running de-novo consensus mode or variant calling mode with a single amplicon per sample with amplicon sizes between 500 and 5,000 bp.
When running variant calling mode with several amplicons per sample, we recommend setting `--drop_frac_longest_reads 0 --take_longest_remaining_reads false`.
Also, when your amplicons are very long or short, it might be a good idea to change the minimum read length filter cutoff.
For more details, please refer to the [Pipeline Overview](#pipeline-overview).

*How can I select which reference is paired with which sample when running variant calling mode with multiple amplicons?* - If you want to map specific reference sequences to specific amplicons, you can add a column named `ref` containing the respective target reference IDs to the sample sheet.
If a barcode should be aligned against multiple reference sequences, their IDs can be included as a space-delimited list.
If a specific sample should be aligned against all references, you can leave the corresponding cell in the `ref` column blank.
Consider the following example with 4 barcodes and 2 sequences in the reference file. The sample sheet tells the workflow that the first two barcodes (`sample1` and `sample2`) are expected to contain
reads from both amplicons in the reference, whereas the other two samples only contain one amplicon each.

```
barcode,alias,type,ref
barcode01,sample1,test_sample,katG::NC_000962.3:2154725-2155670 rpoB::NC_000962.3:760285-761376
barcode02,sample2,test_sample,
barcode03,sample3,positive_control,katG::NC_000962.3:2154725-2155670
barcode04,sample4,test_sample,rpoB::NC_000962.3:760285-761376
```

*How much coverage do I need?* - We recommend >150X average coverage across the individual amplicons.
1500 reads per amplicon should thus be enough in the vast majority of cases.
You can speed up execution time by setting `--reads_downsampling_size 1500` (or to a smaller number if your amplicons are not longer than a few kb).

*Why does the workflow select reads based on length rather than random downsampling per default?* - Despite the workflow dropping reads containing adapter sequences in the middle (using `porechop --discard_middle`), some reads in the input data stemming from concatemers or other PCR artifacts might still make it through preprocessing and could thus interfere with consensus generation.
A simple way to avoid this is to drop a small fraction (e.g. 5%) of longest reads.
Furthermore, `spoa` depends on the order of reads in the input and benefits from "seeing" longer reads first.
Therefore, the following workflow options are used per default `--drop_frac_longest_reads 0.05 --take_longest_remaining_reads`.
This essentially drops the longest 5% of reads and then takes selects the next longest reads (either all of them or the number specified with `--reads_downsampling_size`).
To disable the default behaviour and enforce random downsampling of input reads, use `--drop_frac_longest_reads 0 --take_longest_remaining_reads false`.
Disabling the default behaviour in favour of random downsampling is generally recommended when running variant calling mode with multiple amplicons per sample.

*Is there a difference between the consensus sequences created in variant calling and de-novo consensus mode?* - In variant calling mode, the consensus sequence is generated by applying the variants found by [Medaka](https://github.com/nanoporetech/medaka) to the provided reference.
This comes with the caveat that structural variants too large to be detected by Medaka are not going to be represented in the consensus.
For this reason we suggest to always have a look at the coverage plots in the generated HTML report.
In de-novo consensus mode, however, the consensus is generated either via assembly or [POA](https://doi.org/10.1093/bioinformatics/btg109) and no prior information besides the reads is taken into account.

If your question is not answered here, please report any issues or suggestions on the [github issues](https://github.com/epi2me-labs/wf-alignment/issues) page or start a discussion on the [community](https://community.nanoporetech.com/).




## Related blog posts

## Related blog posts

- [How to align your data](https://labs.epi2me.io/how-to-align/)

See the [EPI2ME website](https://labs.epi2me.io/) for lots of other resources and blog posts.




