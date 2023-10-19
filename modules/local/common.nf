process alignReads {
    label "wfamplicon"
    cpus params.threads
    input: tuple val(meta), path("reads.fastq.gz"), path("reference.fasta")
    output: tuple val(meta), path("*.bam"), path("*.bai")
    script:
    """
    minimap2 -t $task.cpus -ax map-ont reference.fasta reads.fastq.gz \
        -R '@RG\\tID:$meta.alias\\tSM:$meta.alias' \
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
    output: tuple val(meta), path("depth.regions.bed.gz"), optional: true
    script:
    int mosdepth_extra_threads = task.cpus - 1
    ref_id = ref_id ?: ""
    """
    REF_ID="$ref_id"
    # get ref IDs and lengths with `samtools idxstats`
    samtools idxstats input.bam > idxstats

    # if `ref_id` was `null`, assume that there was only one reference and look up its
    # ID from the idxstats
    if [ -z "\$REF_ID" ]; then
        REF_ID=\$(head -n1 idxstats | cut -f1)
        if [[ \$REF_ID = '*' ]]; then
            echo "QUITTING: Only unmapped reads in 'input.bam'."
            exit 0
        fi
    fi

    # get the length of the reference
    REF_LENGTH=\$(grep -w "\$REF_ID" idxstats | cut -f2)

    # calculate the corresponding window length (check `REF_LENGTH` first because
    # `expr a / b` returns non-zero exit code when `a < b`)
    window_length=1
    if [ "\$REF_LENGTH" -gt "$n_windows" ]; then
        window_length=\$(expr \$REF_LENGTH / $n_windows)
    fi

    # get the depths (we could add `-x`, but this loses a lot of detail from the depth
    # curves)
    mosdepth -t $mosdepth_extra_threads -b \$window_length -n -c "\$REF_ID" depth input.bam
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

process lookupMedakaModel {
    label "wfamplicon"
    input:
        path("lookup_table")
        val basecall_model
        val model_type
    output:
        stdout
    shell:
    '''
    medaka_model=$(workflow-glue resolve_medaka_model \
        lookup_table '!{basecall_model}' !{model_type})
    echo -n $medaka_model
    '''
}