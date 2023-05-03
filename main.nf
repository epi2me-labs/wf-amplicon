#!/usr/bin/env nextflow

import groovy.json.JsonBuilder
nextflow.enable.dsl = 2

include { fastq_ingress } from "./lib/fastqingress"
include { pipeline as variantCallingPipeline } from "./modules/local/pipeline"


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
    //  run fastcat on porechopped reads so that we can include the post-trimming stats
    //  in the report
    """
    fastcat <(porechop -i reads.fastq.gz -t $task.cpus --discard_middle) \
        -s $meta.alias \
        -f porechopped-per-file-stats.tsv \
    | bgzip > porechopped.fastq.gz
    """
}

process collectFilesInDir {
    label "wfamplicon"
    cpus 1
    input: tuple val(meta), path("staging_dir/*"), val(dirname)
    output: tuple val(meta), path(dirname)
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
        path "metadata.json"
        path "versions.txt"
        path "params.json"
    output:
        path "wf-amplicon-report.html"
    script:
    """
    workflow-glue report \
        --report-fname wf-amplicon-report.html \
        --data data \
        --reference reference.fasta \
        --meta-json metadata.json \
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

        // get reference
        Path ref = file(params.reference, checkIfExists: true)

        // put fastcat stats into results channel
        def ch_per_sample_results = ch_reads
        | map { meta, reads, stats_dir -> [meta, *file(stats_dir.resolve("*"))] }

        // remove fastcat stats from reads channel
        ch_reads = ch_reads.map { meta, reads, stats_dir -> [meta, reads] }

        // downsample
        if (params.reads_downsampling_size) {
            ch_reads = downsampleReads(ch_reads, params.reads_downsampling_size)
        }

        // trim adapters
        ch_reads = porechop(ch_reads).seqs

        // add the post-trim stats to the results channel
        ch_per_sample_results = ch_per_sample_results.join(porechop.out.stats)

        // run variant calling pipeline
        variantCallingPipeline(ch_reads, ref)

        // add variant calling results to main results channel
        ch_per_sample_results = ch_per_sample_results
        | join(variantCallingPipeline.out.mapping_stats)
        | join(variantCallingPipeline.out.depth)
        | join(variantCallingPipeline.out.variants)

        // collect files for report in directories (one per sample)
        def ch_results_for_report = ch_per_sample_results
        | map {
            meta = it[0]
            rest = it[1..-1]
            [meta, rest, meta.alias]
        }
        | collectFilesInDir
        | map { meta, dirname -> dirname }

        // get all metadata and combine into single JSON file
        def metadata_json = ch_reads
        | map { meta, reads -> meta }
        | collect
        | map { new JsonBuilder(it).toPrettyString() }
        | collectFile(name: "metadata.json", newLine: true)

        // create channel with files to publish; the channel will have the shape `[file,
        // name of sub-dir to be published in]`.
        def ch_to_publish = Channel.empty()
        | mix(
            software_versions | map { [it, null] },
            workflow_params | map { [it, null] },
            variantCallingPipeline.out.variants
            | map { meta, vcf -> [vcf, "$meta.alias/variants"] }
        )

        makeReport(
            ch_results_for_report | collect,
            ref,
            metadata_json,
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

    if (params.disable_ping == false) {
        Pinguscript.ping_post(workflow, "start", "none", params.out_dir, params)
    }

    ArrayList fastcat_extra_args = []
    if (params.min_read_length) { fastcat_extra_args << "-a $params.min_read_length" }
    if (params.max_read_length) { fastcat_extra_args << "-b $params.max_read_length" }
    if (params.min_read_qual) { fastcat_extra_args << "-q $params.min_read_qual" }

    samples = fastq_ingress([
        "input":params.fastq,
        "sample":params.sample,
        "sample_sheet":params.sample_sheet,
        "analyse_unclassified":params.analyse_unclassified,
        "fastcat_stats": true,
        "fastcat_extra_args": fastcat_extra_args.join(" ")])

    // run workflow
    pipeline(samples)

    // publish results
    pipeline.out.combined_results_to_publish
    | toList
    | flatMap
    | output
}

if (params.disable_ping == false) {
    workflow.onComplete {
        Pinguscript.ping_post(workflow, "end", "none", params.out_dir, params)
    }

    workflow.onError {
        Pinguscript.ping_post(workflow, "error", "$workflow.errorMessage", params.out_dir, params)
    }
}
