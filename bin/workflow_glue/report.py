"""Create workflow report."""
from pathlib import Path

import dominate.tags as html_tags
import dominate.util as dom_util
import ezcharts as ezc
from ezcharts.components import fastcat
from ezcharts.components.ezchart import EZChart
from ezcharts.components.reports import labs
from ezcharts.layout.snippets import Tabs
from ezcharts.layout.snippets.stats import Stats
from ezcharts.layout.snippets.table import DataTable
import pandas as pd
import pysam

from . import report_util as util  # noqa: ABS101

from .util import get_named_logger, wf_parser  # noqa: ABS101

# number of points in the depth-across-amplicon line plots
N_DEPTH_WINDOWS = 100


def argparser():
    """Argument parser for entrypoint."""
    parser = wf_parser("report")
    parser.add_argument(
        "--report-fname",
        required=True,
        help="Name of report HTML file.",
    )
    parser.add_argument(
        # all relevant files will be in sub-directories (one per sample) in this
        # directory
        "--data",
        required=True,
        type=Path,
        help=(
            "Directory with one sub-directory per sample containing files "
            "required for the report."
        ),
    )
    parser.add_argument(
        "--reference",
        type=Path,
        help="FASTA file with reference sequences for the individual amplicon.",
    )
    parser.add_argument(
        "--sample-sheet",
        type=Path,
        help="Sample sheet with metadata.",
    )
    parser.add_argument(
        "--versions",
        help="CSV with program versions.",
        required=True,
    )
    parser.add_argument(
        "--params",
        required=True,
        help="JSON file containing the workflow parameters.",
    )
    parser.add_argument(
        "--revision",
        default="unknown",
        help="git branch/tag of the executed workflow",
    )
    parser.add_argument(
        "--commit",
        default="unknown",
        help="git commit of the executed workflow",
    )
    return parser


def main(args):
    """Run the entry point."""
    logger = get_named_logger("Report")

    # in case there was a sample sheet, read it so that we can show the metadata in the
    # per-sample summary table
    metadata = None
    if args.sample_sheet:
        metadata = pd.read_csv(args.sample_sheet, index_col="alias")

    # read data for report
    datasets = sorted(
        [util.ReportDataSet(d) for d in args.data.glob("*")],
        key=lambda x: x.sample_alias,
    )

    # if there was a sample sheet, make the samples in `datasets` have the same order as
    # in the sample sheet
    if metadata is not None:
        datasets_dict = {d.sample_alias: d for d in datasets}
        datasets = [
            datasets_dict[alias] for alias in metadata.index if alias in datasets_dict
        ]

    # create, fill, and write out the report
    report = labs.LabsReport(
        "Workflow Amplicon Sequencing report", "wf-amplicon", args.params, args.versions
    )
    populate_report(report, metadata, datasets, args.reference)

    report.write(args.report_fname)
    logger.info(f"Report written to '{args.report_fname}'.")


def populate_report(report, metadata, datasets, ref_fasta):
    """Fill the report with content."""
    # put whether we got a ref file into a variable
    ref_mode = ref_fasta is not None
    analysis_type_str = "Variant calling" if ref_mode else "Consensus generation"
    samples = [d.sample_alias for d in datasets]
    # check if any samples are missing any of the required inputs and filter them out
    # (we will mention them in the report below); this could only happen if there were
    # no reads for that particular sample
    bad_datasets = [d for d in datasets if not d.all_inputs_valid]
    datasets = [d for d in datasets if d.all_inputs_valid]
    ref_seqs = []
    if ref_mode:
        with pysam.FastxFile(ref_fasta) as f:
            for entry in f:
                ref_seqs.append(entry.name)

    # in case there are no "good" datasets, print a short notice to the report and exit
    # early
    if not datasets:
        with report.add_section("Warning: Analysis incomplete!", "Warning!"):
            if not bad_datasets:
                # there were no datasets at all
                html_tags.p(
                    """
                    No analysis was run as no valid input data were found!
                    """
                )
                return
            # we got invalid data for all samples --> warn user and create truncated
            # report
            html_tags.p(
                f"""{analysis_type_str} failed for all samples.
                Therefore, only a limited report is available. Consider relaxing the
                filtering criteria (e.g. using a lower value for
                """,
                html_tags.kbd("--min_read_qual"),
                """) as this is sometimes caused by filtering out all reads from the
                input data.
                """,
            )
        # show the pre-processing stats and omit all other sections
        preprocessing_section(report, bad_datasets)
        return

    # intro and "at a glance" stats
    with report.add_section("Introduction", "Intro."):
        intro_str = (
            "This report contains tables and plots to help interpret the results "
            "of wf-amplicon. The workflow was run in"
        )
        if not ref_mode:
            # we're in no-ref mode
            html_tags.p(
                intro_str,
                html_tags.b("no-reference mode"),
                """. The individual sections of the report summarize the outcomes of the
                different steps of the workflow (read filtering, consensus generation,
                re-mapping against the consensus for depth of coverage calculation).
                """,
            )
        else:
            # we're in ref mode
            html_tags.p(
                intro_str,
                html_tags.b("reference mode"),
                """. The individual sections of the report summarize the
                outcomes of the different steps of the workflow (read filtering, mapping
                against the reference file containing the amplicon sequences, variant
                calling).
                """,
            )
        html_tags.p("The input data contained:")
        # brief summary of how many samples / refs there were
        for name, items in zip(["sample", "amplicon"], [samples, ref_seqs]):
            if name == "amplicon" and not ref_mode:
                continue
            if len(items) > 1:
                name += "s"
            items_str = ", ".join(items[:7]) + (", ..." if len(items) > 7 else "")
            dom_util.raw(
                f"""
                <b>{len(items)} {name}:</b><br>
                {items_str}<br><br>"""
            )
        # mention if there were any datasets with invalid input files
        if bad_datasets:
            dom_util.raw(
                f"""
                <b>Warning: {analysis_type_str} failed for the following sample(s) (and
                these will be omitted from the rest of the report):</b><br>
                {", ".join([d.sample_alias for d in bad_datasets])}<br><br>
                """
            )
        html_tags.h5("At a glance")
        html_tags.p(
            """
            Key results for the individual samples are shown below. You can use the
            dropdown menu to visualize the results for a different sample.
            """
        )
        # "at a glance" stats in cards (with one dropdown tab per sample)
        tabs = Tabs()
        with tabs.add_dropdown_menu():
            for d in datasets:
                basic_summary = d.get_basic_summary_stats(ref_mode)
                with tabs.add_dropdown_tab(d.sample_alias):
                    # add values + titles of cards for stats reported in either case
                    # (no-ref mode and with reference) first
                    stats_card_items = [
                        (f'{basic_summary["reads"]:,g}', "Reads"),
                        (f'{basic_summary["bases"]:,g}', "Bases"),
                        (basic_summary["mean_length"].round(1), "Mean length"),
                        (basic_summary["mean_quality"].round(1), "Mean quality"),
                    ]
                    if not ref_mode:
                        # add the no-ref-only cards
                        stats_card_items += [
                            (f'{basic_summary["unmapped_reads"]:,g}', "Unmapped reads"),
                            (
                                f'{basic_summary["consensus_length"]:,g}',
                                "Consensus length",
                            ),
                        ]
                    else:
                        n_amplicons_stats_str = (
                            f'{basic_summary["amplicons"]:g} / {len(ref_seqs)}'
                        )
                        n_snps_stats_str = f'{basic_summary["snps"]:g}'
                        n_indels_stats_str = f'{basic_summary["indels"]:g}'
                        stats_card_items += [
                            (n_amplicons_stats_str, "Amplicons detected"),
                            (
                                basic_summary["overall_mean_depth"].round(1),
                                "Mean coverage across all amplicons",
                            ),
                            (
                                basic_summary["min_mean_depth"].round(1),
                                "Smallest mean coverage for any amplicon",
                            ),
                            (n_snps_stats_str, "SNVs"),
                            (n_indels_stats_str, "Indels"),
                        ]
                    Stats(
                        columns=3,
                        items=stats_card_items,
                    )

    # show preprocessing stats
    preprocessing_section(report, datasets)

    # summarize bamstats results of all samples for the following report sections
    bamstats_summary = util.summarize_bamstats(datasets)

    # section for alignment stats and (if in reference mode) detected variants summary
    with report.add_section("Summary", "Summary"):
        if not ref_mode:
            # only show re-alignment stats: make a single table with two rows per
            # sample; one with stats of successfully re-aligned reads and one with stats
            # of unaligned reads
            html_tags.p(
                """
                The table below summarizes the main results of re-mapping the reads of
                each barcode against the generated consensus.
                """
            )
            with dom_util.container() as c:
                no_ref_summary_table = util.format_no_ref_summary_table(
                    bamstats_summary, datasets
                )
            c.clear()
            DataTable.from_pandas(no_ref_summary_table)
        else:
            with dom_util.container() as c:
                # create summary table with one row per sample
                per_sample_summary_table = util.format_per_sample_summary_table(
                    bamstats_summary, datasets
                )
                # add sample meta data
                per_sample_summary_table = (
                    pd.concat((metadata, per_sample_summary_table), axis=1)
                    .fillna("-")
                    .rename(columns=str.capitalize)
                )
                per_sample_summary_table.index.name = "Sample alias"

                # summary table with one row per amplicon
                per_amplicon_summary_table = util.format_per_amplicon_summary_table(
                    bamstats_summary, datasets
                )
            c.clear()

            html_tags.p(
                """
                The two tables below (one per tab) briefly summarize the main results of
                mapping the reads to the provided amplicon references and subsequent
                variant calling.
                """
            )
            tabs = Tabs()
            with tabs.add_tab("Per-sample summary"):
                DataTable.from_pandas(per_sample_summary_table)
            with tabs.add_tab("Per-amplicon summary"):
                DataTable.from_pandas(per_amplicon_summary_table)

            html_tags.p(
                dom_util.raw(
                    """
                    The following table breaks the results down further (one
                    sample&ndash;amplicon combination per row).
                    """
                )
            )
            with dom_util.container() as c:
                summary_table = DataTable.from_pandas(
                    util.format_stats_table(bamstats_summary.fillna(0)).reset_index(),
                    use_index=False,
                )
            c.clear()
            dom_util.raw(str(summary_table))
    # depth coverage plots
    with report.add_section("Depth of coverage", "Depth"):
        html_tags.p(
            """
            Coverage along the individual amplicon, (use the dropdown menu to view the
            plots for the individual amplicons).
            """
        )
        tabs = Tabs()
        with tabs.add_dropdown_menu():
            # get a table in long format with the amplicon ID, centre position of the
            # depth windows, depth values, and sample alias
            per_amplicon_depths = pd.concat(
                [
                    d.depth.assign(
                        pos=lambda df: df[["start", "end"]].mean(axis=1),
                        sample=d.sample_alias,
                    )[["ref", "pos", "depth", "sample"]]
                    for d in datasets
                ]
            )
            for amplicon, depth_df in per_amplicon_depths.groupby("ref"):
                # drop samples with zero depth along the whole amplicon
                total_depths = depth_df.groupby('sample')['depth'].sum()
                samples_with_nonzero_depth = list(total_depths.index[total_depths > 0])
                if not samples_with_nonzero_depth:
                    continue
                depth_df = depth_df.query('sample.isin(@samples_with_nonzero_depth)')
                with tabs.add_dropdown_tab(amplicon):
                    plt = ezc.lineplot(
                        data=depth_df.round(2),
                        x="pos",
                        y="depth",
                        hue="sample",
                    )
                    plt.title = {"text": "Coverage along amplicon"}
                    plt.xAxis.max = depth_df["pos"].max()
                    plt.xAxis.name = "Position along amplicon"
                    plt.yAxis.name = "Sequencing depth"
                    for s in plt.series:
                        # only show the line and no circles
                        s.showSymbol = False
                    EZChart(plt, "epi2melabs")

    # add variant tables (skip if there were no VCFs)
    if not ref_mode:
        return
    with report.add_section("Variants", "Variants"):
        html_tags.p(
            "Haploid variant calling was performed with Medaka. Variants with low ",
            "depth (i.e. smaller than ",
            html_tags.kbd("--min_coverage"),
            ') are shown under the "Low depth" tab.',
        )
        # combine all variants into a single dataframe with `[sample, amplicon]` as
        # multi-level index
        comb_variants = (
            (
                pd.concat(
                    (
                        d.variants.reset_index().assign(sample=d.sample_alias)
                        for d in datasets
                    )
                )
                .astype({"pos": int, "DP": int})
                .rename(
                    columns={
                        "amp": "amplicon",
                        "pos": "position",
                        "ref": "ref. allele",
                        "alt": "alt. allele",
                        "DP": "Depth",
                    }
                )
                .rename(columns=lambda x: x.capitalize())
            )
            .sort_values(["Sample", "Amplicon", "Position"])
            .set_index(["Sample", "Amplicon"])
        )
        # one tab + table for the variants that passed the depth filter and one for the
        # others
        tabs = Tabs()
        with tabs.add_tab("Filters passed"):
            DataTable.from_pandas(
                comb_variants.query("Filter == 'PASS'")
                .drop(columns="Filter")
                .reset_index(),
                use_index=False,
            )
        with tabs.add_tab("Low depth"):
            DataTable.from_pandas(
                comb_variants.query("Filter == 'LOW_DEPTH'")
                .drop(columns="Filter")
                .reset_index(),
                use_index=False,
            )


def preprocessing_section(report, datasets):
    """Add a section with a table of pre-processing stats and the SeqSummary plots."""
    # `fastcat` omits the reads removed during filtering from the per-read stats
    # file, but includes them in the per-file summary. We'll get some basic summary
    # stats (number of reads, number of bases, mean qual, etc.) for both sets of
    # files and put them into a table
    preprocessing_stats = pd.DataFrame(
        index=["Raw", "Filtered", "Downsampled, trimmed"],
        columns=[
            "reads",
            "bases",
            "min_read_length",
            "max_read_length",
            "mean_quality",
        ],
        dtype=float,
    )
    fastcat_per_file_stats = pd.concat((d.fastcat_per_file_stats for d in datasets))
    preprocessing_stats.loc["Raw"] = [
        (total_reads := fastcat_per_file_stats["n_seqs"].sum()),
        fastcat_per_file_stats["n_bases"].sum(),
        fastcat_per_file_stats["min_length"].min(),
        fastcat_per_file_stats["max_length"].max(),
        (fastcat_per_file_stats.eval("n_seqs * mean_quality").sum() / total_reads)
        if total_reads != 0
        else 0,
    ]
    # now get the same stats for after filtering
    fastcat_per_read_stats = pd.concat((d.fastcat_per_read_stats for d in datasets))
    preprocessing_stats.loc["Filtered"] = [
        fastcat_per_read_stats.shape[0],
        fastcat_per_read_stats["read_length"].sum(),
        fastcat_per_read_stats["read_length"].min(),
        fastcat_per_read_stats["read_length"].max(),
        fastcat_per_read_stats["mean_quality"].mean(),
    ]
    # and for after trimming
    post_trim_per_file_stats = pd.concat((d.post_trim_per_file_stats for d in datasets))
    preprocessing_stats.loc["Downsampled, trimmed"] = [
        (total_reads := post_trim_per_file_stats["n_seqs"].sum()),
        post_trim_per_file_stats["n_bases"].sum(),
        post_trim_per_file_stats["min_length"].min(),
        post_trim_per_file_stats["max_length"].max(),
        (post_trim_per_file_stats.eval("n_seqs * mean_quality").sum() / total_reads)
        if total_reads != 0
        else 0,
    ]
    # there could be NaNs in case some datasets didn't have any reads
    preprocessing_stats.fillna(0, inplace=True)
    # some formatting
    for col in ("reads", "bases"):
        preprocessing_stats[col] = (
            preprocessing_stats[col].astype(int).apply(util.si_format)
        )
    for col in ("min_read_length", "max_read_length"):
        preprocessing_stats[col] = (
            preprocessing_stats[col].astype(int).apply(lambda x: f"{x:,d}")
        )
    preprocessing_stats["mean_quality"] = preprocessing_stats["mean_quality"].round(1)
    preprocessing_stats.index.name = "Condition"
    preprocessing_stats.rename(
        columns=lambda col: col.capitalize().replace("_", " "), inplace=True
    )
    with report.add_section("Preprocessing", "Preproc."):
        # show the stats in a DataTable
        html_tags.p(
            """
            Some basic stats covering the raw reads and the reads remaining after the
            initial filtering step (based on length and mean quality) as well as after
            downsampling and trimming are illustrated in the table below.
            """
        )
        DataTable.from_pandas(preprocessing_stats)

        # combine all post-filter per-read stats for the SeqSummary
        comb_per_read_stats = pd.concat(
            data.fastcat_per_read_stats for data in datasets
        )
        # return early if there were no per-read stats (i.e. there were no reads left
        # after filtering)
        if comb_per_read_stats.empty:
            return
        # `SeqSummary` plots
        html_tags.p(
            """
            The following plots show the read quality and length distributions as well
            as the base yield after filtering (but before downsampling / trimming) for
            each sample (use the dropdown menu to view the plots for the
            individual samples).
            """
        )
        fastcat.SeqSummary(comb_per_read_stats)
