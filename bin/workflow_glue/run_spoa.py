"""Interleave forward and reverse reads to improve SPOA performance."""

from pathlib import Path
import sys

from medaka.common import reverse_complement
import parasail
import pysam
import spoa

from .util import get_named_logger, wf_parser  # noqa: ABS101


def get_seqs_from_fastx_and_check_lengths(file, logger, max_len=None):
    """Get sequences from FASTx file while making sure they are not too long.

    We don't want to run SPOA on very long reads as it was not intended for this and
    will take an exorbitant amount of memory.

    :param file: FASTx file
    :param logger: logger for error message
    :param max_len: maximum allowed read length, defaults to None
    :yield: sequences from FASTx file
    """
    with pysam.FastxFile(file, "r") as f:
        for entry in f:
            seq = entry.sequence
            if max_len is not None and len(seq) > max_len:
                # we shouldn't run SPOA on long reads
                logger.error(
                    "Tried to run SPOA with reads longer than "
                    f"{max_len}. This is what assemblers are for. "
                    "Aborting..."
                )
                sys.exit(0)
            yield seq


def align(s1, s2):
    """Use parasail to align two sequences.

    :param s1: string with first DNA sequence
    :param s2: string with second DNA sequence
    :return: parasail alignment Result object
    """
    return parasail.sw_trace_striped_16(
        s1, s2, open=8, extend=4, matrix=parasail.dnafull
    )


def interleave_lists(l1, l2):
    """Uniformly interleave two lists.

    For example:
    ```
    >>> interleave_lists([1, 2, 3, 4], list('ABCDEF'))
    [1, 'A', 'B', 2, 'C', 3, 'D', 'E', 4, 'F']
    ```

    :param l1: first list
    :param l2: second list
    :return: list containing items of `l1` and `l2` uniformly interleaved
    """
    interleaved = []
    rate_1 = 1 / len(l1)
    rate_2 = 1 / len(l2)

    c_1, c_2 = 0, 0

    itr_1 = iter(l1)
    itr_2 = iter(l2)

    while True:
        try:
            if c_1 > c_2:
                interleaved.append(next(itr_2))
                c_2 += rate_2
            else:
                interleaved.append(next(itr_1))
                c_1 += rate_1
        except StopIteration:
            break

    return interleaved


def main(args):
    """Run the entry point."""
    logger = get_named_logger("runSPOA")

    logger.info("Get read orientations...")
    # read input file and determine the orientation of the reads
    fwd = []
    rev = []
    seqs = get_seqs_from_fastx_and_check_lengths(
        args.fastq, logger, args.max_allowed_read_length
    )
    try:
        first_seq = next(seqs)
    except StopIteration:
        # looks like the input file was empty; this shouldn't happen, but we might as
        # well catch it and exit gracefully
        logger.error(f"Input file '{args.fastq}' appears to be empty. Aborting...")
        sys.exit(0)
    for seq in seqs:
        rc = reverse_complement(seq)
        fwd_score = align(seq, first_seq).score
        rev_score = align(rc, first_seq).score
        if fwd_score > rev_score:
            fwd.append(seq)
        else:
            rev.append(rc)

    # uniformly interleave the reads (as long as we got fwd and rev reads)
    if fwd and rev:
        interleaved_reads = interleave_lists(fwd, rev)
        logger.info("Finished interleaving reads.")
    else:
        # we know that we got either fwd or rev reads as otherwise the file must have
        # been empty (which we checked above)
        logger.info("Only got reads from one strand; not interleaving...")
        interleaved_reads = fwd or rev

    # determine minimum coverage if param was provided
    min_cov = None
    if args.relative_min_coverage is not None:
        min_cov = int(round(len(interleaved_reads) * args.relative_min_coverage))
        logger.info(f"SPOA min coverage: {min_cov}.")

    # run spoa
    cons, _ = spoa.poa(interleaved_reads, genmsa=False, min_coverage=min_cov)

    # run spoa a second time with the previous result as first read
    cons, _ = spoa.poa([cons, *interleaved_reads], genmsa=False, min_coverage=min_cov)

    # write out the result
    with open(args.output, "w") as outfile:
        outfile.write(f">consensus\n{cons}\n")

    logger.info("Wrote consensus to STDOUT.")


def argparser():
    """Argument parser for entrypoint."""
    parser = wf_parser("run_spoa")
    parser.add_argument("fastq", type=Path, help="Path to input FASTQ file")
    parser.add_argument(
        "-o", required=True, dest="output", help="Output FASTA file name"
    )
    parser.add_argument(
        "--relative-min-coverage",
        type=float,
        help="Minimum relative coverage of POA graph",
    )
    parser.add_argument(
        "--max-allowed-read-length",
        type=int,
        help=(
            "SPOA can't deal with very long reads; "
            "don't run it when encountering reads longer than this"
        ),
    )
    return parser
