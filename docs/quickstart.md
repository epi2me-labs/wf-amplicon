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
individual sample with `--sample`.

Relevant options for filtering of raw reads are

- `--min_read_length`
- `--max_read_length`
- `--min_read_qual`

After filtering and trimming with
[Porechop](https://github.com/rrwick/Porechop), reads can optionally be
downsampled. You can control the number of reads to keep per sample with
`--reads_downsampling_size`.

Haploid variants are then called with
[Medaka](https://github.com/nanoporetech/medaka). You can set the minimum
coverage a variant needs to exceed in order to be included in the results with
`--min_coverage`. The workflow selects the appropriate
[Medaka](https://github.com/nanoporetech/medaka) model based on the basecaller
configuration that was used to process the signal data. You can use the
parameter `--basecaller_cfg` to provide this information (e.g.
`dna_r10.4.1_e8.2_400bps_hac`) to the workflow. Alternatively, you can choose
the [Medaka](https://github.com/nanoporetech/medaka) model directly with
`--medaka_model`.

If you want to use
[Singularity](https://docs.sylabs.io/guides/latest/user-guide/) instead of
[Docker](https://www.docker.com/products/docker-desktop), add `-profile
singularity`.

### Key outputs

- Interactive HTML report detailing the results.
- VCF files with variants called by
  [Medaka](https://github.com/nanoporetech/medaka).
