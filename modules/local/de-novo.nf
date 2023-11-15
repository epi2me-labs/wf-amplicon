include {
    lookupMedakaModel;
    alignReads as alignDraft;
    alignReads as alignPolished;
    bamstats as bamstatsDraft;
    bamstats as bamstatsPolished;
} from "./common"


/*
Get draft consensus with SPOA. Like `medaka smolecule`, the python script interleaves
reads so that forward and reverse reads are uniformly distributed before being passed to
SPOA.
*/
process spoa {
    label = "medaka"
    cpus 1
    input: tuple val(meta), path("reads.fastq.gz")
    output: tuple val(meta), path("reads.fastq.gz"), path("asm.fasta")
    script:
    String min_cov_args = ""
    if (params.spoa_minimum_relative_coverage) {
        min_cov_args = "--relative-min-coverage $params.spoa_minimum_relative_coverage"
    }
    """
    echo $meta.alias  # makes some debugging easier

    workflow-glue run_spoa reads.fastq.gz $min_cov_args > asm.fasta
    """
}

/*
Assemble draft consensus with `miniasm`. `miniasm` params are taken from pomoxis
`mini_assemble`. Emits the env variable `STATUS` (can be 'failed' or 'passed') alongside
the draft asssembly. If none of the assembled contigs is longer than
`params.force_spoa_length_threshold`, `STATUS` is set to 'failed'.
*/
process miniasm {
    label = "wfamplicon"
    cpus params.threads
    input: tuple val(meta), path("reads.fastq.gz")
    output: tuple val(meta), path("reads.fastq.gz"), path("asm.fasta"), env(STATUS)
    script:
    int mapping_threads = Math.max(1, task.cpus - 1)
    """
    STATUS=failed
    echo $meta.alias  # makes some debugging easier

    # overlap reads and assemble draft consensus (note that `miniasm` outputs the
    # assembly in GFA format)
    minimap2 -L -x ava-ont -t $mapping_threads reads.fastq.gz reads.fastq.gz \
    | miniasm -s 100 -e 3 -f reads.fastq.gz - \
    | awk '/^S/{print ">"\$2"\\n"\$3}' > asm.fasta  # extract header + seq from GFA

    if [[ -s asm.fasta ]]; then
        # check if at least one of the contigs is longer than the force-SPOA-threshold
        samtools faidx asm.fasta
        longest=\$(sort -k2rn asm.fasta.fai | head -n1 | cut -f2)
        if [[ \$longest -gt $params.force_spoa_length_threshold ]]; then
            STATUS=passed
        fi
    fi
    """
}

/*
Polish draft consensus with `racon` (one round only). `racon` params are taken from
pomoxis `mini_assemble`. `--no-trimming` is added as it sometimes trims the consensus
too aggressively. We trim the sequences downstream instead.
*/
process racon {
    label = "wfamplicon"
    cpus params.threads
    input: tuple val(meta), path("reads.fastq.gz"), path("draft.fasta")
    output: tuple val(meta), path("reads.fastq.gz"), path("polished.fasta")
    script:
    int mapping_threads = Math.max(1, task.cpus - 1)
    """
    echo $meta.alias  # makes some debugging easier

    # align against draft
    minimap2 -L -x ava-ont -t $mapping_threads draft.fasta reads.fastq.gz \
    | bgzip > pre-racon.paf.gz

    # run racon
    racon -m 8 -x -6 -g -8 -w 500 -t $task.cpus -q -1 --no-trimming \
        reads.fastq.gz pre-racon.paf.gz draft.fasta \
        > polished.fasta
    """
}

/*
Polish draft consensus with `medaka`. Also, emit the polished sequences as FASTQ
(with consensus quality scores calculated by medaka) and not FASTA.
*/
process medakaConsensus {
    label "medaka"
    cpus Math.min(params.threads, 3)
    input:
        tuple val(meta), path("input.bam"), path("input.bam.bai"), path("draft.fasta")
        val medaka_model
    output: tuple val(meta), path("consensus.fastq")
    script:
    """
    medaka consensus input.bam consensus_probs.hdf \
        --threads $task.cpus --model $medaka_model

    medaka stitch consensus_probs.hdf draft.fasta consensus.fastq \
        --threads $task.cpus --qualities
    """
}

/*
Use `mosdepth` to calculate per-base depths. This is used for trimming the consensus and
not for the plots in the report.
*/
process mosdepthPerBase {
    label "wfamplicon"
    cpus Math.min(params.threads, 3)
    input: tuple val(meta), path("input.bam"), path("input.bam.bai")
    output: tuple val(meta), path("depth.per-base.bed.gz")
    script:
    int mosdepth_extra_threads = task.cpus - 1
    """
    echo $meta.alias  # makes some debugging easier

    mosdepth -t $mosdepth_extra_threads depth input.bam
    """
}

/*
Get `mosdepth` depths for a given number of windows. Expects only a single reference to
be present in the input BAM. Emits a TSV with header instead of the original BED file as
this is expected by the report code.
*/
process mosdepthWindows {
    label "wfamplicon"
    cpus Math.min(params.threads, 3)
    input:
        tuple val(meta), path("input.bam"), path("input.bam.bai")
        val n_windows
    output: tuple val(meta), path("per-window-depth.tsv.gz")
    script:
    int mosdepth_extra_threads = task.cpus - 1
    """
    echo $meta.alias  # makes some debugging easier

    # get ref IDs and lengths with `samtools idxstats`
    samtools idxstats input.bam | grep -v '^*' > idxstats

    # we expect only a single ref seq
    if [[ \$(wc -l < idxstats) -ne 1 ]]; then
        echo "Found unexpected number of references in input BAM."
        exit 1
    fi

    # get the length of the reference
    REF_LENGTH=\$(cut -f2 idxstats)

    # calculate the corresponding window length (check `REF_LENGTH` first because
    # `expr a / b` returns non-zero exit code when `a < b`)
    window_length=1
    if [ "\$REF_LENGTH" -gt "$n_windows" ]; then
        window_length=\$(expr \$REF_LENGTH / $n_windows)
    fi

    # get the depths (we could add `-x`, but this loses a lot of detail from the depth
    # curves)
    mosdepth -t $mosdepth_extra_threads -b \$window_length -n depth input.bam

    # add a header to the BED file
    cat <(printf "ref\\tstart\\tend\\tdepth\\n" | gzip) depth.regions.bed.gz \
    > per-window-depth.tsv.gz

    # clean up BED file to save some space
    rm depth.regions.bed.gz
    """
}

/*
Performs QC on the draft consensus. The script filters out contigs that have either too
little mean depth or a too large fraction of secondary + supplementary alignments. If
more than one contig passes both filters, the one with higher mean depth is chosen. It
is then trimmed based on `mosdepth` per-base depths. The process outputs two optional
channels (`passed` and `failed`), depending on whether QC was passed. It makes sure that
it emits in only one of the channels.

`qc-summary.tsv` contains info like mean depth, the number of primary alignments, etc.
for each contig. It also contains the assembly method so that files from the SPOA and
`miniasm` pipelines can be concatenated downstream.
*/
process trimAndQC {
    label "wfamplicon"
    cpus params.threads
    cpus 1
    input:
        tuple val(meta),
            path("consensus.fastq"),
            path("flagstat.tsv"),
            path("depth.per-base.bed.gz")
        val asm_method
    output:
        tuple val(meta), path("passed/consensus.fastq"), path("passed/qc-summary.tsv"),
            optional: true, emit: passed
        tuple val(meta), path("failed/qc-summary.tsv"), optional: true, emit: failed
    script:
    """
    echo $meta.alias  # makes some debugging easier

    workflow-glue trim_and_qc \
        --alias $meta.alias \
        --asm-method $asm_method \
        --depth depth.per-base.bed.gz \
        --flagstat flagstat.tsv \
        --consensus consensus.fastq \
        --outdir-pass passed \
        --outdir-fail failed \
        --minimum-depth $params.minimum_mean_depth \
        --primary-threshold $params.primary_alignments_threshold \
        --relative-depth-trim-threshold $params.spoa_minimum_relative_coverage \
        --qc-summary-tsv qc-summary.tsv

    # sanity check: make sure that `qc-summary.tsv` was only written for either passed
    # or failed (and not both)
    if [[ ( -f passed/qc-summary.tsv ) && ( -f failed/qc-summary.tsv ) ]]; then
        echo "Found 'qc-summary.tsv' in both 'passed' and 'failed' directories."
        exit 1
    elif [[ ( ! -f passed/qc-summary.tsv ) && ( ! -f failed/qc-summary.tsv ) ]]; then
        echo "Found 'qc-summary.tsv' in neither 'passed' nor 'failed' directories."
        exit 1
    fi
    """
}


// workflow module
workflow pipeline {
    take:
        // expected shape: `[meta, reads]` (reads have already been downsampled and
        // trimmed)
        ch_reads
        method
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

        if (method !in ["miniasm", "spoa"]) {
            error "Invalid de novo method '$method' (needs to be either 'miniasm' " +
                "or 'spoa')."
        }

        // run miniasm -> racon or SPOA depending on `method`
        if (method == "miniasm") {
            ch_branched = miniasm(ch_reads)
            // branch off samples for which the assembly has failed (or was shorter than
            // `--force_spoa_length_threshold`)
            | branch { meta, reads, asm, status ->
                passed: status == "passed"
                failed: status == "failed"
                err: error "Post-assembly status is neither 'passed' nor 'failed' " +
                    "for sample '$meta.alias'."
            }

            ch_draft = ch_branched.passed
            | map { meta, reads, asm, status -> [meta, reads, asm] }
            | racon
            // failed samples will be re-tried with SPOA in `main.nf`
            ch_failed = ch_branched.failed
        } else {
            ch_draft = spoa(ch_reads)
            ch_failed = Channel.empty()
        }

        // re-align reads against the draft consensus
        alignDraft(ch_draft)

        // polish with medaka
        medakaConsensus(
            alignDraft.out
            | join(ch_draft)
            | map { meta, bam, bai, reads, cons -> [meta, bam, bai, cons] },
            medaka_model
        )

        // get mapping stats and mosdepth per-base results for trimming the consensus
        bamstatsDraft(alignDraft.out)
        mosdepthPerBase(alignDraft.out)

        // select and trim the consensus sequence with the largest mean coverage
        trimAndQC(
            medakaConsensus.out
            | join(
                bamstatsDraft.out | map { meta, bamstats, flagstat -> [meta, flagstat] }
            )
            | join(mosdepthPerBase.out),
            method
        )

        // if QC passed, re-align again against the selected polished contig and get
        // stats + depths for report
        ch_reads
        | join(trimAndQC.out.passed, remainder: true)
        | filter { null !in it }
        | map { meta, reads, cons, qc_summary ->
            [meta, reads, cons]
        }
        | alignPolished

        bamstatsPolished(alignPolished.out)
        mosdepthWindows(alignPolished.out, params.number_depth_windows)
    emit:
        // emit three channels:
        // * passed: consensus, stats, and depths for successful samples
        // * metas_failed: meta maps of samples that failed (either the `miniasm` step
        //   or QC)
        // * qc_summaries: all QC summaries (regardless of fail or pass)
        passed = trimAndQC.out.passed
        | map { meta, cons, qc_summary -> [meta, cons] }
        | join(alignPolished.out)
        | join(bamstatsPolished.out)
        | join(mosdepthWindows.out)

        metas_failed = ch_failed
        | map { meta, reads, asm, status -> meta }
        | mix(trimAndQC.out.failed | map { meta, qc_summary -> meta } )

        qc_summaries = trimAndQC.out.passed
        | map { meta, cons, qc_summary -> [meta, qc_summary] }
        | mix(trimAndQC.out.failed)
}
