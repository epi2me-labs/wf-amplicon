"""Get IDs of reads with largest `aligned_ref_len` from `bamstats` results TSV file."""
import sys

import pandas as pd

from .util import get_named_logger, wf_parser  # noqa: ABS101


def argparser():
    """Argument parser for entrypoint."""
    parser = wf_parser("get_reads_with_longest_aligned_ref")
    parser.add_argument(
        "bamstats_results",
        help="`bamstats` results TSV file.",
    )
    parser.add_argument(
        "num_reads_per_strand",
        help="Number of reads to select.",
        type=int,
    )
    return parser


def main(args):
    """Run the entry point."""
    logger = get_named_logger("longAlnRds")
    logger.info(
        f"Extract IDs of {args.num_reads_per_strand} reads with "
        "longest `aligned_ref_len`."
    )

    # we only load the relevant columns and use categorical dtype where possible to
    # reduce the memory footprint
    dtypes = {
        "name": str,
        "ref": "category",
        "direction": "category",
        "aligned_ref_len": int,
    }
    df = pd.read_csv(
        args.bamstats_results, sep="\t", usecols=dtypes.keys(), dtype=dtypes
    ).query('ref != "*"')

    # for each ref--direction combo, write out the read IDs of the longest
    # `args.num_reads_per_strand` reads
    df.groupby(["ref", "direction"]).apply(
        lambda grp_df: sys.stdout.write(
            "\n".join(
                grp_df.sort_values("aligned_ref_len", ascending=False)["name"].head(
                    args.num_reads_per_strand
                )
            )
            + "\n"
        )
    )

    logger.info("Finished writing read IDs with longest `aligned_ref_len` to STDOUT.")
