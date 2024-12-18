"""Subset reads based on length (from fastcat per-read stats)."""
from pathlib import Path
import sys

import pandas as pd

from .util import get_named_logger, wf_parser  # noqa: ABS101


def main(args):
    """Run the entry point."""
    logger = get_named_logger("subsetReads")

    logger.info("Read per-read stats and sort lengths.")
    sorted_lengths = pd.read_csv(
        args.per_read_stats,
        sep="\t",
        index_col="read_id",
        usecols=["read_id", "read_length"],
    )["read_length"].sort_values(ascending=False)

    drop_longest_n = 0
    if args.drop_longest_frac:
        drop_longest_n = int(round(len(sorted_lengths) * args.drop_longest_frac))
        logger.info(f"Drop {drop_longest_n} longest reads.")
        sorted_lengths = sorted_lengths.iloc[drop_longest_n:]

    if args.take_longest_remaining:
        if args.n_reads == 0:
            logger.info("Take all remaining reads (sorted by length).")
            target_ids = sorted_lengths.index
        else:
            logger.info(f"Take {args.n_reads} longest remaining reads.")
            target_ids = sorted_lengths.iloc[: args.n_reads].index
    else:
        logger.info("Randomly select reads and write out IDs.")
        # select a subset and sort by length again (but only if we still got enough
        # reads)
        if args.n_reads > len(sorted_lengths) or args.n_reads == 0:
            target_ids = sorted_lengths.index
        else:
            target_ids = (
                sorted_lengths.sample(args.n_reads).sort_values(ascending=False).index
            )

    sys.stdout.write("\n".join(target_ids) + "\n")

    logger.info(
        f"Finished extracting {len(target_ids)} reads from '{args.per_read_stats}'."
    )


def argparser():
    """Argument parser for entrypoint."""
    parser = wf_parser("subset_reads")
    parser.add_argument(
        "per_read_stats", type=Path, help="Path to per-read stats TSV file"
    )
    parser.add_argument(
        "drop_longest_frac", type=float, help="Fraction of longest reads to drop"
    )
    parser.add_argument(
        "n_reads", type=int, help="Number of reads to select. If 0, take all."
    )
    parser.add_argument(
        "--take-longest-remaining",
        action="store_true",
        help=(
            "Whether to select the longest remaining reads or just randomly "
            "sample after dropping the longest reads"
        ),
    )
    return parser
