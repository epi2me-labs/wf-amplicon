#!/usr/bin/env nextflow

import groovy.json.JsonBuilder
nextflow.enable.dsl = 2

include { fastq_ingress } from "./lib/ingress"
include { pipeline as variantCallingPipeline } from "./modules/local/variant-calling"
include {
    pipeline as deNovoPipeline_asm; pipeline as deNovoPipeline_spoa;
} from "./modules/local/de-novo"


OPTIONAL_FILE = file("$projectDir/data/OPTIONAL_FILE")


process getVersions {
    label "wfamplicon"
    cpus 1
    memory "2 GB"
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
    mosdepth --version | tr ' ' ',' >> versions.txt
    printf "miniasm,%s\\n" \$(miniasm -V) >> versions.txt
    printf "racon,%s\\n" \$(racon --version) >> versions.txt
    printf "csvtk,%s\\n" \$(csvtk version | grep -oP '\\d+\\.\\d+.\\d+') >> versions.txt
    """
}

process addMedakaToVersionsFile {
    label "medaka"
    cpus 1
    memory "2 GB"
    input: path "old_versions.txt"
    output: path "versions.txt"
    script:
    """
    cp old_versions.txt versions.txt
    medaka --version | tr ' ' ',' >> versions.txt
    """
}

process getParams {
    label "wfamplicon"
    cpus 1
    memory "2 GB"
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
    memory "2 GB"
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

/*
Instead of downsampling randomly, this process subsets input reads based on length.
Depending on parameters, it drops a fraction (e.g. `0.05` for 5%) of longest reads and
then selects a number of the next-longest reads (if `take_longest_remaining_reads` is
`true`; if `false`, it will randomly downsample the reads after dropping the longest).
Selected reads will always be sorted by length.
*/
process subsetReads {
    label "wfamplicon"
    cpus 1
    memory "2 GB"
    input:
        tuple val(meta), path("reads.fastq.gz"), path("fastcat-stats")
        val drop_longest_frac
        val take_longest_remaining_reads
        val n_reads
    output: tuple val(meta), path("subset.fastq.gz")
    script:
    String take_longest_remaining_arg = ""
    if (take_longest_remaining_reads) {
        take_longest_remaining_arg = "--take-longest-remaining"
    }
    """
    echo $meta.alias  # sometimes useful for debugging

    workflow-glue subset_reads \
        fastcat-stats/per-read-stats.tsv.gz \
        $drop_longest_frac \
        $n_reads \
        $take_longest_remaining_arg \
        > target_ids.txt

    # only attempt to subset reads if neither the FASTQ or the file with target IDs is
    # empty (`samtools fqidx` would fail otherwise)
    if [[ ( -n \$(zcat reads.fastq.gz | head -c1) ) && ( -s target_ids.txt ) ]]; then
        samtools fqidx reads.fastq.gz -r target_ids.txt | bgzip > subset.fastq.gz
    else
        echo -n | bgzip > subset.fastq.gz
    fi
    """
}

process porechop {
    label "wfamplicon"
    cpus Math.min(params.threads, 5)
    memory {
        // `porechop` quite annoyingly loads the whole FASTQ into memory
        // (https://github.com/rrwick/Porechop/issues/77) and uses up to 4x the size of
        // the `.fastq.gz` file (depending on compression ratio); we give another factor
        // of 2 for extra margin to be on the safe side
        def fastq_size = fastq.size()
        fastq_size > 2e9 ? "32 GB" : (fastq_size > 2.5e8 ? "16 GB" : "2 GB")
    }
    input: tuple val(meta), path(fastq, stageAs: "reads.fastq.gz")
    output:
        tuple val(meta), path("porechopped.fastq.gz"), emit: seqs
        tuple val(meta), path("porechopped-per-file-stats.tsv"), emit: stats
    script:
    int porechop_threads = Math.max(task.cpus - 2, 1)
    //  run fastcat on porechopped reads so that we can include the post-trimming stats
    //  in the report
    """
    (
        # only run porechop if the FASTQ isn't empty
        if [[ -n \$(zcat reads.fastq.gz | head -c1) ]]; then
            porechop -i reads.fastq.gz -t $porechop_threads --discard_middle
        else
            echo -n
        fi
    ) \
    | fastcat /dev/stdin -s $meta.alias -f porechopped-per-file-stats.tsv \
    | bgzip > porechopped.fastq.gz
    """
}

process collectFilesInDir {
    label "wfamplicon"
    cpus 1
    memory "2 GB"
    input: tuple path("staging_dir/*"), val(dirname)
    output: path(dirname)
    script:
    """
    mv staging_dir $dirname
    """
}

process concatTSVs {
    label "wfamplicon"
    cpus 1
    memory "2 GB"
    input:
        tuple val(meta), path("input/f*.tsv")
        val fname
    output: tuple val(meta), path(fname)
    script:
    """
    csvtk concat -t input/* > $fname
    """

}

process makeReport {
    label "wfamplicon"
    cpus 1
    memory "8 GB"
    input:
        path "data/*"
        path ref, stageAs: "ref/*"
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
    String ref_arg = ref.fileName.name == OPTIONAL_FILE.name ? "" : \
        "--reference \"$ref\""
    """
    workflow-glue report \
        --report-fname wf-amplicon-report.html \
        --data data \
        $ref_arg \
        $sample_sheet_arg \
        --downsampling-size $params.reads_downsampling_size \
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
    memory "2 GB"
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

        // drop elements in `ch_reads` that only contain a metamap (this occurs when
        // there was an entry in a sample sheet but no corresponding barcode dir)
        ch_reads = ch_reads.filter { meta, reads, stats_dir -> reads }

        // add fastcat stats of raw reads to channel with results for report
        ch_results_for_report = ch_reads
        | map { meta, reads, stats_dir -> [meta, *file(stats_dir.resolve("*"))] }

        // either subset reads (i.e. drop the longest and then take the next longest) or
        // downsample naively or take all reads
        if (params.drop_frac_longest_reads || params.take_longest_remaining_reads) {
            // subset reads
            ch_reads = subsetReads(
                ch_reads,
                params.drop_frac_longest_reads,
                params.take_longest_remaining_reads,
                params.reads_downsampling_size,
            )
        } else if (params.reads_downsampling_size) {
            // downsample reads
            ch_reads = downsampleReads(
                ch_reads | map { meta, reads, stats -> [meta, reads] },
                params.reads_downsampling_size
            )
        } else {
            // take all reads
            ch_reads = ch_reads.map { meta, reads, stats_dir -> [meta, reads] }
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
        // sum up the `n_seqs` column and print warning if the total is `0`
        .map { if (it.sum() == 0) {
                log.warn "No reads left after pre-processing; report will only show " +
                    "limited results. Perhaps consider relaxing the filtering " +
                    "criteria (`--min_read_qual` etc.)."
            }
        }

        // drop samples with fewer than `params.min_n_reads` reads after filtering
        ch_reads = ch_reads
        | join(ch_post_filter_n_reads)
        | filter { meta, reads, n_seqs -> n_seqs >= params.min_n_reads }
        | map { meta, reads, n_seqs -> [meta, reads] }

        Path ref = OPTIONAL_FILE
        if (params.reference) {
            // variant calling mode; make sure the ref file exists and run the variant
            // calling pipeline
            ref = file(params.reference, checkIfExists: true)

            variantCallingPipeline(ch_reads, ref)

            // add variant calling results to channels for the report + to publish
            ch_results_for_report = ch_results_for_report
            | join(variantCallingPipeline.out.mapping_stats, remainder: true)
            | join(variantCallingPipeline.out.depth, remainder: true)
            | join(variantCallingPipeline.out.variants, remainder: true)
            // drop any `null`s from the joined list
            | map { it - null }

            ch_to_publish = ch_to_publish
            | mix(
                variantCallingPipeline.out.sanitized_ref | map { [it, null] },
                variantCallingPipeline.out.variants
                | map { meta, vcf, idx -> [[vcf, idx], "$meta.alias/variants"] }
                | transpose,
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
                    variantCallingPipeline.out.combined_vcfs | map { [it, null] }
                    | transpose,
                    variantCallingPipeline.out.combined_bams | map { [it, null] }
                    | transpose,
                )
            }
        } else {
            // de-novo consensus mode: run pipeline for de-novo consensus; we try
            // `miniasm` first and then run `spoa` on all samples that failed
            deNovoPipeline_asm(ch_reads, "miniasm")

            deNovoPipeline_spoa(
                deNovoPipeline_asm.out.metas_failed
                // use `combine` instead of `join` here because `join` with `remainder:
                // true` does not work as expected if one channel contains only the key
                // (as is the case with the `metas_failed` channel here)
                | combine(ch_reads, by: 0),
                "spoa"
            )

            // combine the results for `miniasm` and `spoa`
            ch_no_ref_results = deNovoPipeline_asm.out.passed
            | mix(deNovoPipeline_spoa.out.passed)
            | multiMap { meta, cons, bam, bai, bamstats, flagstat, depth ->
                consensus: [meta, cons]
                mapped: [meta, bam, bai]
                for_report: [meta, bamstats, flagstat, depth]
            }

            // get the QC summary TSVs and concat them for each sample (if the `miniasm`
            // and `spoa` pipelines were run for a sample, it will have two QC summary
            // TSVs)
            ch_concat_qc_summaries = concatTSVs(
                deNovoPipeline_asm.out.qc_summaries
                | mix(deNovoPipeline_spoa.out.qc_summaries)
                | groupTuple,
                "qc-summary.tsv"
            )

            // add de-novo results to main results channel
            ch_results_for_report = ch_results_for_report
            | join(ch_no_ref_results.for_report, remainder: true)
            | join(ch_concat_qc_summaries, remainder: true)
            // drop any `null`s from the joined list
            | map { it -> it.findAll { it } }

            // add the consensus and alignments (re-aligned against the consensus) to
            // the output channel for publishing
            ch_to_publish = ch_to_publish
            | mix(
                ch_no_ref_results.consensus
                | map { meta, cons -> [cons, "$meta.alias/consensus"] },
                ch_no_ref_results.mapped
                | map { meta, bam, bai -> [[bam, bai], "$meta.alias/alignments"] }
                | transpose,
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
