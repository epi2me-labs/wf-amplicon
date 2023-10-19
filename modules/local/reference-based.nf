include {
    alignReads;
    bamstats;
    mosdepth;
    concatMosdepthResultFiles;
    lookupMedakaModel
} from "./common"


process downsampleBAMforMedaka {
    label "wfamplicon"
    cpus 1
    input:
        tuple val(meta), path("input.bam"), path("input.bam.bai"), path("bamstats.tsv")
    output: tuple val(meta), path("downsampled.bam"), path("downsampled.bam.bai")
    script:
    """
    workflow-glue get_reads_with_longest_aligned_ref_len \
        bamstats.tsv \
        $params.medaka_target_depth_per_strand\
    > downsampled.read_IDs

    samtools view input.bam -N downsampled.read_IDs -o downsampled.bam
    samtools index downsampled.bam
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

    # use the sample alias as sample name in the VCF and filter variants
    bcftools reheader medaka.annotated.unfiltered.vcf -s <(echo '$meta.alias') \
    | bcftools filter \
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

process mergeVCFs {
    label "medaka"
    cpus 1
    input: path "VCFs/file*.vcf.gz"
    output: path "combined.vcf.gz"
    script:
    """
    (
        cd VCFs
        ls | xargs -n1 bcftools index
    )
    bcftools merge VCFs/file*.vcf.gz -Oz -o combined.vcf.gz
    """
}

process mergeBAMs {
    label "wfamplicon"
    cpus 1
    input:
        path "BAMs/file*.bam"
        path "indices/file*.bam.bai"
    output: path "combined.bam"
    script:
    """
    samtools merge BAMs/* indices/* -pXo combined.bam
    """
}


// workflow module
workflow pipeline {
    take:
        // reads have already been downsampled and trimmed
        // expects shape: `[meta, reads]`
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
            medaka_model = lookupMedakaModel(
                lookup_table, params.basecaller_cfg, "medaka_variant"
            )
        }

        // align to reference
        ch_reads
        | combine(ref)
        | alignReads

        // get mapping stats for report and pre-Medaka downsampling
        bamstats(alignReads.out)

        // downsample to target depth before running Medaka
        alignReads.out
        | join(bamstats.out)
        | map { meta, bam, bai, bamstats, flagstat ->
            [meta, bam, bai, bamstats]
        }
        | downsampleBAMforMedaka

        // get the seq IDs of the amplicons from the ref file
        def ch_amplicon_seqIDs = ref | splitFasta(record: [id: true]) | map { it.id }

        // run medaka consensus (the process will run once for each sample--amplicon
        // combination)
        def ch_medaka_consensus_probs = medakaConsensus(
            downsampleBAMforMedaka.out | combine(ch_amplicon_seqIDs),
            medaka_model
        ) | groupTuple

        // get the variants
        medakaVariant(
            ch_medaka_consensus_probs | join(alignReads.out),
            ref,
            params.min_coverage
        )

        // run mosdepth on each sample--amplicon combination
        mosdepth(
            alignReads.out | combine(ch_amplicon_seqIDs),
            params.number_depth_windows
        )
        // concat the depth files for each sample
        mosdepth.out
        | groupTuple
        | concatMosdepthResultFiles

        // combine per-sample BAM and VCF files if requested
        combined_vcfs = null
        combined_bams = null
        if (params.combine_results) {
            combined_vcfs = mergeVCFs(
                medakaVariant.out.filtered.collect { meta, vcf -> vcf }
            )
            combined_bams = mergeBAMs(
                alignReads.out.collect { meta, bam, bai -> bam },
                alignReads.out.collect { meta, bam, bai -> bai }
            )
        }
    emit:
        sanitized_ref = ref
        mapped = alignReads.out
        mapping_stats = bamstats.out
        depth = concatMosdepthResultFiles.out
        variants = medakaVariant.out.filtered
        consensus = medakaVariant.out.consensus
        combined_vcfs
        combined_bams
}
