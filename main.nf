#!/usr/bin/env nextflow

import groovy.json.JsonBuilder
nextflow.enable.dsl = 2

include { fastq_ingress } from "./lib/ingress"
include { pipeline as variantCallingPipeline } from "./modules/local/reference-based"


OPTIONAL_FILE = file("$projectDir/data/OPTIONAL_FILE")


process getVersions {
    label "wfamplicon"
    cpus 1
    output: path "versions.txt"
    script:
    """
    python --version | tr ' ' ',' | sed 's/P/p/' > versions.txt
    fastcat --version | sed 's/^/fastcat,/' >> versions.txt
    python -c "import ezcharts; print(f'ezcharts,{ezcharts.__version__}')" >> versions.txt
    python -c "import dominate; print(f'dominate,{dominate.__version__}')" >> versions.txt
    python -c "import numpy; print(f'numpy,{numpy.__version__}')" >> versions.txt
    python -c "import pandas; print(f'pandas,{pandas.__version__}')" >> versions.txt
    python -c "import pysam; print(f'pysam,{pysam.__version__}')" >> versions.txt
    python -c "import si_prefix; print(f'si-prefix,{si_prefix.__version__}')" >> versions.txt
    seqkit version | tr ' ' ',' >> versions.txt
    printf "porechop,%s\\n" \$(porechop --version) >> versions.txt
    samtools --version | head -n1 | tr ' ' ',' >> versions.txt
    printf "minimap2,%s\\n" \$(minimap2 --version) >> versions.txt
    """
}

process addMedakaToVersionsFile {
    label "medaka"
    cpus 1
    input: path "old_versions.txt"
    output: path "versions.txt"
    script:
    """
    medaka --version | tr ' ' ',' >> old_versions.txt
    mv old_versions.txt versions.txt
    """
}

process getParams {
    label "wfamplicon"
    cpus 1
    output:
        path "params.json"
    script:
    String paramsJSON = new JsonBuilder(params).toPrettyString()
    """
    # Output nextflow params object to JSON
    echo '$paramsJSON' > params.json
    """
}

process downsampleReads {
    label "wfamplicon"
    cpus Math.min(params.threads, 3)
    input:
        tuple val(meta), path("reads.fastq.gz")
        val n_reads
    output: tuple val(meta), path("downsampled.fastq.gz")
    script:
    int bgzip_threads = task.cpus == 1 ? 1 : task.cpus - 1
    """
    seqkit sample -2 reads.fastq.gz -n $n_reads \
    | bgzip -@ $bgzip_threads > downsampled.fastq.gz
    """
}

process porechop {
    label "wfamplicon"
    cpus Math.min(params.threads, 5)
    input: tuple val(meta), path("reads.fastq.gz")
    output:
        tuple val(meta), path("porechopped.fastq.gz"), emit: seqs
        tuple val(meta), path("porechopped-per-file-stats.tsv"), emit: stats
    script:
    int porechop_threads = Math.max(task.cpus - 2, 1)
    //  run fastcat on porechopped reads so that we can include the post-trimming stats
    //  in the report
    """
    fastcat <(porechop -i reads.fastq.gz -t $porechop_threads --discard_middle) \
        -s $meta.alias \
        -f porechopped-per-file-stats.tsv \
    | bgzip > porechopped.fastq.gz
    """
}

process collectFilesInDir {
    label "wfamplicon"
    cpus 1
    input: tuple path("staging_dir/*"), val(dirname)
    output: path(dirname)
    script:
    """
    mv staging_dir $dirname
    """
}

process makeReport {
    label "wfamplicon"
    cpus 1
    input:
        path "data/*"
        path "reference.fasta"
        path sample_sheet, stageAs: "sample_sheet/*"
        path "versions.txt"
        path "params.json"
    output:
        path "wf-amplicon-report.html"
    script:
    // we need to check `.fileName.name` here because `TaskPath.name` includes the
    // staging directory (https://github.com/nextflow-io/nextflow/issues/3574)
    String sample_sheet_arg = sample_sheet.fileName.name == OPTIONAL_FILE.name ? "" : \
        "--sample-sheet \"$sample_sheet\""
    """
    workflow-glue report \
        --report-fname wf-amplicon-report.html \
        --data data \
        --reference reference.fasta \
        $sample_sheet_arg \
        --versions versions.txt \
        --params params.json
    """
}

// See https://github.com/nextflow-io/nextflow/issues/1636. This is the only way to
// publish files from a workflow whilst decoupling the publish from the process steps.
// The process takes a tuple containing the filename and the name of a sub-directory to
// put the file into. If the latter is `null`, puts it into the top-level directory.
process output {
    // publish inputs to output directory
    label "wfamplicon"
    cpus 1
    publishDir (
        params.out_dir,
        mode: "copy",
        saveAs: { dirname ? "$dirname/$fname" : fname }
    )
    input:
        tuple path(fname), val(dirname)
    output:
        path fname
    """
    """
}


// workflow module
workflow pipeline {
    take:
        ch_reads
    main:
        def software_versions = getVersions() | addMedakaToVersionsFile
        def workflow_params = getParams()

        // Make channel to hold the results required for the report. It will have the
        // shape `[meta, file(s)]` where the second item can be a list of files or an
        // individual file.
        ch_results_for_report = Channel.empty()
        // Make a channel for the files to be published. It will look like `[file, name
        // of sub-dir to be published in | null]` (when `null`, the file will be
        // published in the top level output directory).
        ch_to_publish = Channel.empty()
        | mix(
            software_versions | map { [it, null] },
            workflow_params | map { [it, null] },
        )

        // get reference
        Path ref = file(params.reference, checkIfExists: true)

        // drop elements in `ch_reads` that only contain a metamap (this occurs when
        // there was an entry in a sample sheet but no corresponding barcode dir)
        ch_reads = ch_reads.filter { meta, reads, stats_dir -> reads }

        // add fastcat stats of raw reads to channel with results for report
        ch_results_for_report = ch_reads
        | map { meta, reads, stats_dir -> [meta, *file(stats_dir.resolve("*"))] }

        // remove fastcat stats from reads channel
        ch_reads = ch_reads.map { meta, reads, stats_dir -> [meta, reads] }

        // downsample
        if (params.reads_downsampling_size) {
            ch_reads = downsampleReads(ch_reads, params.reads_downsampling_size)
        }

        // trim adapters
        ch_reads = porechop(ch_reads).seqs

        // add the post-trim stats to the report channel
        ch_results_for_report = ch_results_for_report.join(porechop.out.stats)

        // parse the post-porechop fastcat per-file stats to get the number of reads
        // after filtering for each sample (`splitCsv` actually works on lists as well;
        // it will emit `[meta, row]` for each row in the CSV, but in this case there is
        // only one row in each CSV file anyway)
        ch_post_filter_n_reads = porechop.out.stats
        .splitCsv(sep: "\t", header: true)
        .map { meta, stats -> [meta, stats["n_seqs"] as int] }

        // make sure that there are reads left after pre-processing and print a warning
        // otherwise
        ch_post_filter_n_reads
        .collect { meta, n_seqs -> n_seqs }
        // sum up the `n_seqs`  column and throw an error if the total is `0`
        .map { if (it.sum() == 0) {
                log.warn "No reads left after pre-processing; report will only show " +
                    "limited results. Perhaps consider relaxing the filtering " +
                    "criteria (`--min_read_qual` etc.)."
            }
        }

        // drop samples for which all reads have been removed during filtering
        ch_reads = ch_reads
        | join(ch_post_filter_n_reads)
        | filter { meta, reads, n_seqs -> n_seqs > 0 }
        | map { meta, reads, n_seqs -> [meta, reads] }

        variantCallingPipeline(ch_reads, ref)

        // add variant calling results to channels for the report + to publish
        ch_results_for_report = ch_results_for_report
        | join(variantCallingPipeline.out.mapping_stats, remainder: true)
        | join(variantCallingPipeline.out.depth, remainder: true)
        | join(variantCallingPipeline.out.variants, remainder: true)
        // drop any `null`s from the joined list
        | map { it -> it.findAll { it } }

        ch_to_publish = ch_to_publish
        | mix(
            variantCallingPipeline.out.sanitized_ref | map { [it, null] },
            variantCallingPipeline.out.variants
            | map { meta, vcf -> [vcf, "$meta.alias/variants"] },
            variantCallingPipeline.out.mapped
            | map { meta, bam, bai -> [[bam, bai], "$meta.alias/alignments"] }
            | transpose,
            variantCallingPipeline.out.consensus
            | map { meta, cons -> [cons, "$meta.alias/consensus"] },
        )
        // combine VCF and BAM files if requested
        if (params.combine_results) {
            ch_to_publish = ch_to_publish
            | mix(
                variantCallingPipeline.out.combined_vcfs | map { [it, null] },
                variantCallingPipeline.out.combined_bams | map { [it, null] },
            )
        }

        // collect files for report in directories (one dir per sample)
        ch_results_for_report = ch_results_for_report
        | map {
            meta = it[0]
            files = it[1..-1]
            [files, meta.alias]
        }
        | collectFilesInDir

        // create the report and add it to the publishing channel
        makeReport(
            ch_results_for_report | collect,
            ref,
            file(params.sample_sheet ?: OPTIONAL_FILE),
            software_versions,
            workflow_params,
        )
        ch_to_publish = ch_to_publish
        | mix(makeReport.out | map { [it, null] } )

    emit:
        combined_results_to_publish = ch_to_publish
}


// entrypoint workflow
WorkflowMain.initialise(workflow, params, log)
workflow {

    Pinguscript.ping_start(nextflow, workflow, params)

    ArrayList fastcat_extra_args = []
    if (params.min_read_length) { fastcat_extra_args << "-a $params.min_read_length" }
    if (params.max_read_length) { fastcat_extra_args << "-b $params.max_read_length" }
    if (params.min_read_qual) { fastcat_extra_args << "-q $params.min_read_qual" }

    samples = fastq_ingress([
        "input":params.fastq,
        "sample":params.sample,
        "sample_sheet":params.sample_sheet,
        "analyse_unclassified":params.analyse_unclassified,
        "stats": true,
        "fastcat_extra_args": fastcat_extra_args.join(" ")])

    // run workflow
    pipeline(samples)

    // publish results
    pipeline.out.combined_results_to_publish
    | toList
    | flatMap
    | output
}

workflow.onComplete {
    Pinguscript.ping_complete(nextflow, workflow, params)
}
workflow.onError {
    Pinguscript.ping_error(nextflow, workflow, params)
}
