include {
    lookupMedakaModel;
    alignReads;
    bamstats;
    mosdepth;
    concatMosdepthResultFiles;
} from "./common"


process medakaSmolecule {
    label = "medaka"
    cpus Math.min(params.threads, 2)
    errorStrategy { task.exitStatus == 65 ? "retry" : "finish" }
    maxRetries 3
    input:
        tuple val(meta), path("reads.fastq.gz"), val(n_reads)
        val medaka_model
    output: tuple val(meta),
        path("reads.fastq.gz"),
        path("consensus.fasta"),
        optional: true
    script:
    // POA is dependent on the order of input reads. Therefore, in case of failure due
    // to producing multiple output sequences, shuffling the input reads might lead to a
    // single consensus. If no single consensus is generated on the third try, no output
    // will be produced.
    String cat_cmd = task.attempt == 1 ? "zcat" : "seqkit shuffle -s $task.attempt"
    int spoa_min_cov = Math.round(n_reads * params.spoa_minimum_relative_coverage)
    //  When passing a single FASTx file to `medaka smolecule`, it expects multiple
    //  amplicons in that file (with the amplicon to which a read belongs being denoted
    //  by the part of the id line before the first underscore; e.g. amp1_r1, amp1_r2,
    //  ..., ampX_rY). We need to format the ID lines accordingly.
    """
    medaka smolecule out \
        <($cat_cmd reads.fastq.gz | awk '{
            if (NR % 4 == 1) {printf ">read_%d\\n", ++c} else if (NR % 4==2) print;
        }') \
        --length 0 \
        --model $medaka_model \
        --threads $task.cpus \
        --spoa_min_coverage $spoa_min_cov

    # If `smolecule` produces more than one sequence, fail the process with exit code
    # `65` (it will shuffle the reads upon retry). If already in the third attempt, exit
    # with code `0` (no retry and no output).
    if [[ \$(grep -c '^>' out/consensus.fasta) -gt 1 ]]; then
        echo "QUITTING: More than one sequence in 'out/consensus.fasta'."
        [[ $task.attempt -ge 3 ]] && exit_code=0 || exit_code=65
        exit \$exit_code
    fi

    # give consensus sequence a more meaningful ID before emitting
    sed -E '1s/^>.*/>${meta["alias"]}/' out/consensus.fasta > consensus.fasta
    """
}


// workflow module
workflow pipeline {
    take:
        // reads have already been downsampled and trimmed
        // expects shape: `[meta, reads, n_reads]`
        ch_reads
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
                lookup_table, params.basecaller_cfg, "medaka_consensus"
            )
        }

        // run `medaka smolecule` to get the consensus
        medakaSmolecule(ch_reads, medaka_model)
        // realign the reads against the consensus (for coverage stats etc.)
        | alignReads

        // get mapping stats for report
        bamstats(alignReads.out)

        // run mosdepth on each sample
        mosdepth(
            alignReads.out | combine(Channel.of(null)),
            params.number_depth_windows
        )
        // concat the depth files for each sample
        mosdepth.out
        | groupTuple
        | concatMosdepthResultFiles
    emit:
        consensus = medakaSmolecule.out.map { meta, reads, cons -> [meta, cons] }
        mapped = alignReads.out
        mapping_stats = bamstats.out
        depth = concatMosdepthResultFiles.out
}
