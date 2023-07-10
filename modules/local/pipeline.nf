process alignReads {
    label "wfamplicon"
    cpus params.threads
    input:
        tuple val(meta), path("reads.fastq.gz")
        path "reference.fasta"
    output:
        tuple val(meta), path("*.bam"), path("*.bai")
    script:
    """
    minimap2 -t $task.cpus -ax map-ont reference.fasta reads.fastq.gz \
    | samtools sort -@ $task.cpus -o aligned.sorted.bam -

    samtools index aligned.sorted.bam
    """
}

process bamstats {
    label "wfamplicon"
    cpus Math.min(params.threads, 2)
    input: tuple val(meta), path("input.bam"), path("input.bam.bai")
    output: tuple val(meta), path("bamstats.tsv"), path("bamstats-flagstat.tsv")
    script:
    """
    bamstats -u input.bam -s $meta.alias -f bamstats-flagstat.tsv -t $task.cpus \
    > bamstats.tsv
    """
}

process mosdepth {
    label "wfamplicon"
    cpus Math.min(params.threads, 3)
    input:
        tuple val(meta), path("input.bam"), path("input.bam.bai"), val(ref_id)
        val n_windows
    output: tuple val(meta), path("depth.regions.bed.gz")
    script:
    int mosdepth_extra_threads = task.cpus - 1
    """
    # get the length of the reference
    ref_length=\$(samtools idxstats input.bam | awk '\$1 == "$ref_id" {print \$2}')

    # calculate the corresponding window length (check `ref_length` first because
    # `expr a / b` returns non-zero exit code when `a < b`)
    if [ "\$ref_length" -lt "$n_windows" ]; then
        window_length=1
    else
        window_length=\$(expr \$ref_length / $n_windows)
    fi

    # get the depths (we could add `-x`, but it loses a lot of detail from the depth
    # curves)
    mosdepth -b \$window_length -n -c "$ref_id" depth input.bam
    """
}

process concatMosdepthResultFiles {
    label "wfamplicon"
    cpus 1
    input: tuple val(meta), path("depth.*.bed.gz")
    output: tuple val(meta), path("per-window-depth.tsv.gz")
    script:
    """
    # add a header line and concatenate the depth .bed files (note that gzipped data
    # can be concatenated just like regular data)
    cat <(printf "ref\\tstart\\tend\\tdepth\\n" | gzip) depth.*.bed.gz \
    > per-window-depth.tsv.gz
    """
}

process sanitizeRefFile {
    // some tools can't deal with `:` or `*` in ref FASTA ID lines --> replace them (and
    // whitespace) with underscores
    label "wfamplicon"
    cpus 1
    input: path "reference.fasta"
    output: path "reference_sanitized_seqIDs.fasta"
    script:
    """
    sed '/^>/s/:\\|\\*\\| /_/g' reference.fasta > reference_sanitized_seqIDs.fasta
    """
}

process lookupMedakaVariantModel {
    label "wfamplicon"
    input:
        path("lookup_table")
        val basecall_model
    output:
        stdout
    shell:
    '''
    medaka_model=$(workflow-glue resolve_medaka_model \
        lookup_table '!{basecall_model}' "medaka_variant")
    echo -n $medaka_model
    '''
}

process medakaConsensus {
    label "medaka"
    cpus Math.min(params.threads, 2)
    input:
        tuple val(meta), path("input.bam"), path("input.bam.bai"), val(reg)
        val medaka_model
    output: tuple val(meta), path("consensus_probs.hdf")
    script:
    """
    medaka consensus input.bam consensus_probs.hdf \
        --threads $task.cpus --regions "$reg" --model $medaka_model
    """
}

process medakaVariant {
    label "medaka"
    cpus 1
    input:
        tuple val(meta),
            path("consensus_probs*.hdf"),
            path("input.bam"),
            path("input.bam.bai")
        path "reference.fasta"
        val min_coverage
    output:
        tuple val(meta), path("medaka.annotated.vcf.gz"), emit: filtered
        tuple val(meta), path("medaka.annotated.unfiltered.vcf"), emit: unfiltered
        tuple val(meta), path("medaka.consensus.fasta"), emit: consensus
    """
    medaka variant reference.fasta consensus_probs*.hdf medaka.vcf
    medaka tools annotate --dpsp medaka.vcf reference.fasta input.bam \
        medaka.annotated.unfiltered.vcf

    # filter variants
    bcftools filter medaka.annotated.unfiltered.vcf \
        -e 'INFO/DP < $min_coverage' \
        -s LOW_DEPTH \
        -Oz -o medaka.annotated.vcf.gz

    # make consensus seqs
    bcftools index medaka.annotated.vcf.gz
    bcftools consensus -f reference.fasta medaka.annotated.vcf.gz \
        -i 'FILTER="PASS"' \
        -o medaka.consensus.fasta
    """
}

// workflow module
workflow pipeline {
    take:
        // reads have already been downsampled and trimmed
        ch_reads
        ref
    main:
        // sanitize seq IDs in ref FASTA
        def ref = sanitizeRefFile(ref)

        // get medaka model (look up or override)
        def medaka_model
        if (params.medaka_model) {
            log.warn "Overriding Medaka model with ${params.medaka_model}."
            medaka_model = params.medaka_model
        } else {
            Path lookup_table = file(
                "${projectDir}/data/medaka_models.tsv", checkIfExists: true)
            medaka_model = lookupMedakaVariantModel(lookup_table, params.basecaller_cfg)
        }

        // align to reference
        alignReads(ch_reads, ref)

        // get the seq IDs of the amplicons from the ref file
        def ch_amplicon_seqIDs = ref | splitFasta(record: [id: true]) | map { it.id }

        // run medaka consensus (the process will run once for each sample--amplicon
        // combination)
        def ch_medaka_consensus_probs = medakaConsensus(
            alignReads.out | combine(ch_amplicon_seqIDs),
            medaka_model
        ) | groupTuple

        // get the variants
        medakaVariant(
            ch_medaka_consensus_probs | join(alignReads.out),
            ref,
            params.min_coverage
        )

        // get mapping stats for report
        bamstats(alignReads.out)

        // run mosdepth on each sample--amplicon combination
        mosdepth(
            alignReads.out | combine(ch_amplicon_seqIDs),
            params.number_depth_windows
        )
        // concat the depth files for each sample
        mosdepth.out
        | groupTuple
        | concatMosdepthResultFiles
    emit:
        mapped = alignReads.out | map { meta, bam, bai -> [meta, bam] }
        mapping_stats = bamstats.out
        depth = concatMosdepthResultFiles.out
        variants = medakaVariant.out.filtered
        consensus = medakaVariant.out.consensus
}
