#!/usr/bin/env nextflow

import groovy.json.JsonBuilder
nextflow.enable.dsl = 2

include { fastq_ingress } from './lib/fastqingress'
include { clusterReads } from './subworkflows/clustering/vsearch'
include { draftAssembly } from './subworkflows/assembly/flye'

OPTIONAL_FILE = file("$projectDir/data/OPTIONAL_FILE")

process getVersions {
    label "wfamplicon"
    cpus 1
    output:
        path "versions.txt"
    script:
    """
    python -c "import pysam; print(f'pysam,{pysam.__version__}')" >> versions.txt
    fastcat --version | sed 's/^/fastcat,/' >> versions.txt
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



// See https://github.com/nextflow-io/nextflow/issues/1636. This is the only way to
// publish files from a workflow whilst decoupling the publish from the process steps.
// The process takes a tuple containing the filename and the name of a sub-directory to
// put the file into. If the latter is `null`, puts it into the top-level directory.
process output {
    // publish inputs to output directory
    label "wfamplicon"
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
        reads
    main:
        software_versions = getVersions()
        workflow_params = getParams()
    
        // the reads have already been filtered by `fastcat` --> cluster next
        clustering = clusterReads(
            reads.map {it[0..1]},
            params.min_cluster_size,
        )

        
    emit:
        workflow_params
}


params.min_read_length = 300
params.max_read_length = 3600
params.min_read_qual = 8
params.min_cluster_size = 0.2

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
        "fastcat_stats": false,
        "fastcat_extra_args": fastcat_extra_args.join(" ")])
    
    // looks like this is the most robust way to check if a `param` coming from the
    // command line is a number
    if (
        params.min_cluster_size instanceof String ||
        !params.min_cluster_size.toString().isNumber()
    ) {
        error "`--min_cluster_size` must be a float or integer."
    }

    pipeline(samples)
}

if (params.disable_ping == false) {
    workflow.onComplete {
        Pinguscript.ping_post(workflow, "end", "none", params.out_dir, params)
    }

    workflow.onError {
        Pinguscript.ping_post(workflow, "error", "$workflow.errorMessage", params.out_dir, params)
    }
}
