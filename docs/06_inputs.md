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
| sample_sheet | string | A CSV file used to map barcodes to sample aliases. The sample sheet can be provided when the input data is a directory containing sub-directories with FASTQ files. | The sample sheet is a CSV file with, minimally, columns named `barcode` and `alias`. Extra columns are allowed. A `type` column is required for certain workflows and should have the following values; `test_sample`, `positive_control`, `negative_control`, `no_template_control`. |  |
| sample | string | A single sample name for non-multiplexed data. Permissible if passing a single .fastq(.gz) file or directory of .fastq(.gz) files. |  |  |


### Pre-processing Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| min_read_length | integer | Shorter reads will be removed. |  | 300 |
| max_read_length | integer | Longer reads will be removed. |  |  |
| min_read_qual | number | Reads with a lower mean quality will be removed. |  | 10 |
| drop_frac_longest_reads | number | Drop fraction of longest reads. | The very longest reads might be concatemers or contain other artifacts. In many cases removing them simplifies de novo consensus generation. | 0.05 |
| take_longest_remaining_reads | boolean | Whether to use the longest (remaining) reads. | With this option, reads are not randomly selected during downsampling (potentially after the longest reads have been removed), but instead the longest remaining reads are taken. This generally improves performance on long amplicons. | True |
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
| medaka_target_depth_per_strand | integer | Downsample each amplicon to this per-strand depth before running Medaka. | Medaka performs best with even strand coverage and depths between 80X and 400X. To avoid too high coverage, the workflow downsamples the reads for each amplicon to this per-strand depth before running Medaka. Changing this value is discouraged as it might cause decreased performance. | 75 |


### Miscellaneous Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| threads | integer | Maximum number of CPU threads to use per workflow task. | Several tasks in this workflow benefit from using multiple CPU threads. This option sets the maximum number of CPU threads for such processes. The total CPU resources used by the workflow are constrained by the executor configuration in `nextflow.config`. | 4 |
| disable_ping | boolean | Enable to prevent sending a workflow ping. |  | False |


