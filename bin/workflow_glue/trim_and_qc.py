"""Select a contig if multiple were created, do some QC, modify depth / stats files."""
from pathlib import Path
import sys

import pandas as pd
import pysam

from .util import get_named_logger, wf_parser  # noqa: ABS101


def main(args):
    """Run the entry point."""
    logger = get_named_logger("mostCovRef")

    logger.info("Read input files.")

    flagstat = pd.read_csv(args.flagstat, sep="\t", index_col=0)

    depths = pd.read_csv(
        args.depth, sep="\t", header=None, names=["ref", "start", "end", "depth"]
    )

    with pysam.FastxFile(args.consensus) as f:
        fastq = {entry.name: entry for entry in f}

    # make sure we got the same seq IDs in flagstat and the FASTQ file
    seq_ids = [x for x in flagstat.index if x != "*"]
    if set(fastq.keys()) != set(seq_ids):
        ValueError("Sequence IDs in FASTQ and in flagstat file don't agree.")

    # if all reads in the BAM are unmapped, `mosdepth` annoyingly outputs an empty
    # per-base depth file; we need to account for that
    if depths.empty:
        mean_depths = pd.Series(0, index=seq_ids)
    else:
        # we got a valid depth file; calculate the mean depths of the candidate
        # consensus contigs
        mean_depths = depths.groupby("ref").apply(
            lambda df: df.eval("(end - start) * depth").sum() / df["end"].iloc[-1]
        )

    # perform QC checks on each contig
    qc_stats = pd.DataFrame(
        columns=[
            "method",
            "length",
            "mean_depth",
            "primary",
            "secondary",
            "supplementary",
            "primary_ratio",
            "status",
            "fail_reason",
        ]
    )
    qc_stats.index.name = "contig"
    for contig, mean_depth in mean_depths.items():
        status = "passed"
        contig_len = len(fastq[contig].sequence)
        fail_reasons = []
        if mean_depth < args.minimum_depth:
            status = "failed"
            fail_reasons.append("low depth")
        primary_ratio = flagstat.loc[contig, "primary"] / flagstat.loc[contig, "total"]
        if primary_ratio < args.primary_threshold:
            status = "failed"
            fail_reasons.append("low primary ratio")

        qc_stats.loc[contig] = [
            args.asm_method,
            contig_len,
            round(mean_depth, 3),
            *flagstat.loc[contig, ["primary", "secondary", "supplementary"]],
            round(primary_ratio, 3),
            status,
            ", ".join(fail_reasons),
        ]

    qc_passed = qc_stats.query("status == 'passed'")

    failed = qc_passed.empty

    # create the output directory
    outdir = Path(args.outdir_fail if failed else args.outdir_pass)
    outdir.mkdir()

    qc_stats.to_csv(outdir / args.qc_summary_tsv, sep="\t")

    if failed:
        logger.warning(
            "Consensus generation failed for this sample. No consensus candidate "
            "passed the quality checks."
        )
        sys.exit()

    # select the contig with the largest mean depth
    selected = qc_passed["mean_depth"].sort_values(ascending=False).index[0]
    logger.info(
        f"Selected contig '{selected}'. "
        f"Write modified depth and stats files to {outdir}."
    )

    # only keep the depths of the selected contig
    depths.query("ref == @selected", inplace=True)

    # determine the limits for trimming off low-coverage regions at the ends
    first_and_last_idx = depths.query(
        "depth > depth.max() * @args.relative_depth_trim_threshold"
    ).index[[0, -1]]
    depths = depths.loc[slice(*first_and_last_idx)]
    trim_start = depths.iloc[0]["start"]
    trim_end = depths.iloc[-1]["end"]

    # extract the selected consensus sequence and trim low-coverage regions from the
    # ends
    logger.info("Extract and write selected consensus sequence.")
    with pysam.FastxFile(args.consensus) as f:
        for entry in f:
            if entry.name == selected:
                break
        else:
            raise ValueError(
                f"Selected consensus contig {selected} not found in"
                f"consensus FASTx file '{args.consensus}'"
            )
    entry.name = args.alias
    entry.sequence = entry.sequence[trim_start:trim_end]
    entry.quality = entry.quality[trim_start:trim_end]
    with open(outdir / args.consensus.name, "w") as f:
        f.write(str(entry) + "\n")

    logger.info("Done")


def argparser():
    """Argument parser for entrypoint."""
    parser = wf_parser("get_ref_with_largest_mean_depth")
    parser.add_argument(
        "--alias",
        type=str,
        help="Sample alias.",
        required=True,
    )
    parser.add_argument(
        "--asm-method",
        type=str,
        help="Assembly method (is added to QC summary TSV).",
        required=True,
    )
    parser.add_argument(
        "--depth",
        type=Path,
        help="Path to mosdepth per-base BED file.",
        required=True,
    )
    parser.add_argument(
        "--flagstat",
        type=Path,
        help="Path to bamstats flagstat TSV file.",
        required=True,
    )
    parser.add_argument(
        "--consensus",
        type=Path,
        help="Path to consensus FASTQ file.",
        required=True,
    )
    parser.add_argument(
        "--outdir-pass",
        type=Path,
        help="Name of output directory if passing QC.",
        required=True,
    )
    parser.add_argument(
        "--outdir-fail",
        type=Path,
        help="Name of output directory if failing QC.",
        required=True,
    )
    parser.add_argument(
        "--relative-depth-trim-threshold",
        type=float,
        help=(
            "Trim ends of consensus sequence if the relative coverage drops below "
            "this value."
        ),
        required=True,
    )
    parser.add_argument(
        "--minimum-depth",
        type=int,
        help="Minimum required depth to trust a consensus sequence.",
        required=True,
    )
    parser.add_argument(
        "--primary-threshold",
        type=float,
        help=(
            "Threshold for proportion of primary alignments (will write a "
            "warning for samples with lower fraction)."
        ),
        required=True,
    )
    parser.add_argument(
        "--qc-summary-tsv",
        type=str,
        help="Name of TSV file into which to write QC summary.",
        required=True,
    )
    return parser
