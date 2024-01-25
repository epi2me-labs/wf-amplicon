include {
    alignReads;
    bamstats;
    mosdepth;
    concatMosdepthResultFiles;
    lookupMedakaModel
} from "./common"


process sanitizeRefFile {
    // some tools don't like `:` or `*` in FASTA header lines lines (e.g. medaka will
    // try to parse the sequence ID as region string if it contains `:`) --> replace
    // them (and whitespace) with underscores
    label "wfamplicon"
    cpus 1
    memory "2 GB"
    input: path "reference.fasta"
    output: path "reference_sanitized_seqIDs.fasta"
    script:
    """
    sed '/^>/s/[:\\*\\t ]/_/g' reference.fasta > reference_sanitized_seqIDs.fasta
    """
}

process subsetRefFile {
    label "wfamplicon"
    cpus 1
    memory "2 GB"
    input: tuple val(metas), path("reference.fasta"), val(target_seqs)
    output: tuple val(metas), path("reference_subset.fasta"), val(target_seqs)
    script:
    // need to manually join with " " here as this will be an `ArrayList` and not a
    // `BlankSeparatedList`
    String refs_str = target_seqs.join(" ")
    """
    samtools faidx reference.fasta $refs_str > reference_subset.fasta
    """
}

process downsampleBAMforMedaka {
    label "wfamplicon"
    cpus 1
    memory "8 GB"
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

process medakaConsensus {
    label "medaka"
    cpus Math.min(params.threads, 2)
    memory "8 GB"
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
    cpus 2
    memory "8 GB"
    input:
        tuple val(meta),
            path("consensus_probs*.hdf"),
            path("input.bam"),
            path("input.bam.bai"),
            path("reference.fasta")
        val min_coverage
    output:
        tuple(
            val(meta),
            path("medaka.annotated.vcf.gz"),
            path("medaka.annotated.vcf.gz.csi"),
            emit: filtered
        )
        tuple val(meta), path("medaka.consensus.fasta"), emit: consensus
    script:
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
    memory "4 GB"
    input: path "VCFs/file*.vcf.gz"
    output: tuple path("combined.vcf.gz"), path("combined.vcf.gz.csi")
    script:
    """
    (
        cd VCFs
        ls | xargs -n1 bcftools index
    )
    bcftools merge VCFs/file*.vcf.gz -Oz -o combined.vcf.gz
    bcftools index combined.vcf.gz
    """
}

process mergeBAMs {
    label "wfamplicon"
    cpus 1
    memory "16 GB"
    input:
        path "BAMs/file*.bam"
        path "indices/file*.bam.bai"
    output: tuple path("combined.bam"), path("combined.bam.bai")
    script:
    """
    samtools merge BAMs/* indices/* -pXo combined.bam
    samtools index combined.bam
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

        // some tools don't like `:` or `*` in FASTA header lines lines (e.g. medaka
        // will try to parse the sequence ID as region string if it contains `:`) and we
        // need to sanitize the ref IDs
        san_ref = sanitizeRefFile(ref)

        // create a map from refID to sanitized refID; we'll use this to translate
        // target ref IDs from the sample sheet to sanitized ref IDs so that we can
        // subset the sanitized ref file
        ref_id_map = Channel.empty()
        | concat(
            Channel.of(ref).splitFasta(record: [id: true]).map{ it.id }.collect(),
            san_ref.splitFasta(record: [id: true]).map{ it.id }.collect()
        )
        | toList
        | map { it.transpose().collectEntries() as LinkedHashMap }

        // We don't want to write out sanitized FASTA files more often than we have to.
        // We thus group the samples by target ref IDs and subset only once per group
        ch_branched = ch_reads
        | map { meta, reads -> [meta["ref"]?.split() as Set, meta] }
        | groupTuple
        | combine(san_ref)
        | combine(ref_id_map)
        | branch { refs, metas, san_ref, ref_id_map ->
            with_target_refs: refs as boolean
            no_target_refs: true
        }

        // For samples with target ref IDs, subset the sanitized ref file and add the
        // subset of sanitized ref IDs to the channel. For samples without target ref
        // IDs, add the full sanitized ref file and all sanitized ref IDs.
        ch_sanitized_refs_and_ids = ch_branched.with_target_refs
        | map { refs, metas, san_ref, ref_id_map ->
            ArrayList ordered_target_ref_ids = ref_id_map.keySet().findAll { it in refs }
            ArrayList sanitized_target_ref_ids = \
                ref_id_map.subMap(ordered_target_ref_ids).values()
            [metas, san_ref, sanitized_target_ref_ids]
        }
        | subsetRefFile
        | mix(
            ch_branched.no_target_refs | map { _, metas, san_ref, ref_id_map ->
                [metas, san_ref, ref_id_map.values() as ArrayList]
            }
        )
        // transpose the channel for shape `[meta, sanitized ref, sanitized ref IDs]`
        | transpose(by: 0)

        // put the sanitized + subsetted refs and the ref IDs into separate channels
        ch_sanitized_refs = ch_sanitized_refs_and_ids
        | map { meta, ref, ids -> [meta, ref] }
        ch_sanitized_ids = ch_sanitized_refs_and_ids
        | map { meta, ref, ids -> [meta, ids] }

        // align the reads
        ch_reads
        | join(ch_sanitized_refs)
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

        // run medaka consensus (the process will run once for each sample--amplicon
        // combination)
        ch_medaka_consensus_probs = medakaConsensus(
            // join with and transpose on the list of sanitized IDs for each sample
            downsampleBAMforMedaka.out | join(ch_sanitized_ids) | transpose(by: 3),
            medaka_model
        ) | groupTuple

        // get the variants
        medakaVariant(
            ch_medaka_consensus_probs
            | join(downsampleBAMforMedaka.out)
            | join(ch_sanitized_refs),
            params.min_coverage
        )

        // run mosdepth on each sample--amplicon combination
        mosdepth(
            alignReads.out | join(ch_sanitized_ids) | transpose(by: 3),
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
                medakaVariant.out.filtered.collect { meta, vcf, idx -> vcf }
            )
            combined_bams = mergeBAMs(
                alignReads.out.collect { meta, bam, bai -> bam },
                alignReads.out.collect { meta, bam, bai -> bai }
            )
        }
    emit:
        sanitized_ref = san_ref
        mapped = alignReads.out
        mapping_stats = bamstats.out
        depth = concatMosdepthResultFiles.out
        variants = medakaVariant.out.filtered
        consensus = medakaVariant.out.consensus
        combined_vcfs
        combined_bams
}
