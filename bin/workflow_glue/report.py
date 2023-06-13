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

# set order to sort samples by type (as defined in the sample sheet)
SAMPLE_TYPE_ORDER = {
    "test_sample": 0,
    "positive_control": 1,
    "negative_control": 2,
    "no_template_control": 3,
}


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
        "--meta-json",
        type=Path,
        help="JSON file containing the meta data for all samples.",
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

    # read the meta data for all files
    meta_df = (
        pd.read_json(args.meta_json)
        .set_index("alias")
        .rename(columns=lambda col: col.capitalize())
        .fillna("-")
    )

    report = labs.LabsReport(
        "Workflow Amplicon Sequencing report", "wf-amplicon", args.params, args.versions
    )

    # read data for report (sort by sample type and alias)
    datasets = sorted(
        [util.ReportDataSet(d) for d in args.data.glob("*")],
        key=lambda x: (
            SAMPLE_TYPE_ORDER[meta_df.loc[x.sample_alias, "Type"]],
            x.sample_alias,
        ),
    )
    samples = [d.sample_alias for d in datasets]
    # check if any samples are missing any of the required inputs and filter them out
    # (we will mention them in the report below); this could only happen if there were
    # no reads for that particular sample
    bad_datasets = [d for d in datasets if not d.all_inputs_valid]
    datasets = [d for d in datasets if d.all_inputs_valid]
    ref_seqs = []
    with pysam.FastxFile(args.reference) as ref_file:
        for entry in ref_file:
            ref_seqs.append(entry.name)

    # intro and "at a glance" stats
    with report.add_section("Introduction", "Intro."):
        html_tags.p(
            """
            This report contains tables and plots to help interpret the results of the
            wf-amplicon workflow. The individual sections of the report summarize the
            outcomes of the different steps of the workflow (read filtering, mapping
            against the reference file containing the amplicon sequences, variant
            calling).
            """
        )
        html_tags.p("The input data contained:")
        # brief summary of how many samples / refs there were
        for name, items in zip(["sample", "amplicon"], [samples, ref_seqs]):
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
                <b>No valid results were found for the following sample(s) (and these
                will be omitted from the rest of the report):</b><br>
                {", ".join([d.sample_alias for d in bad_datasets])}<br><br>
                """
            )
        html_tags.h5("At a glance")
        html_tags.p(
            """
            Key results for the individual samples are shown below. You can use the
            dropdown menu to select individual samples.
            """
        )
        # "at a glance" stats in cards (with one dropdown tab per sample)
        tabs = Tabs()
        with tabs.add_dropdown_menu():
            for d in datasets:
                basic_summary = d.get_basic_summary_stats()
                with tabs.add_dropdown_tab(d.sample_alias):
                    Stats(
                        columns=3,
                        items=[
                            (f'{int(basic_summary["reads"]):,d}', "Reads"),
                            (f'{int(basic_summary["bases"]):,d}', "Bases"),
                            (basic_summary["mean_length"].round(1), "Mean length"),
                            (basic_summary["mean_quality"].round(1), "Mean quality"),
                            (
                                f'{basic_summary["amplicons"]:g} / {len(ref_seqs)}',
                                "Amplicons detected",
                            ),
                            (
                                basic_summary["overall_mean_depth"].round(1),
                                "Mean coverage across all amplicons",
                            ),
                            (
                                basic_summary["min_mean_depth"].round(1),
                                "Smallest mean coverage for any amplicon",
                            ),
                            (f'{int(basic_summary["snps"])}', "SNVs"),
                            (f'{int(basic_summary["indels"])}', "Indels"),
                        ],
                    )

    # raw reads stats
    with report.add_section("Preprocessing", "Preproc."):
        # `fastcat` omits the reads removed during filtering from the per-read stats
        # file, but includes them in the per-file summary. We'll get some basic summary
        # stats (number of reads, number of bases, mean qual, etc.) for both sets of
        # files and put them into a table
        basic_summary_stats = pd.DataFrame(
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
        basic_summary_stats.loc["Raw"] = [
            (total_reads := fastcat_per_file_stats["n_seqs"].sum()),
            fastcat_per_file_stats["n_bases"].sum(),
            fastcat_per_file_stats["min_length"].min(),
            fastcat_per_file_stats["max_length"].max(),
            fastcat_per_file_stats.eval("n_seqs * mean_quality").sum() / total_reads,
        ]
        # now get the same stats for after filtering
        fastcat_per_read_stats = pd.concat((d.fastcat_per_read_stats for d in datasets))
        basic_summary_stats.loc["Filtered"] = [
            fastcat_per_read_stats.shape[0],
            fastcat_per_read_stats["read_length"].sum(),
            fastcat_per_read_stats["read_length"].min(),
            fastcat_per_read_stats["read_length"].max(),
            fastcat_per_read_stats["mean_quality"].mean(),
        ]
        # and for after trimming
        post_trim_per_file_stats = pd.concat(
            (d.post_trim_per_file_stats for d in datasets)
        )
        basic_summary_stats.loc["Downsampled, trimmed"] = [
            (total_reads := post_trim_per_file_stats["n_seqs"].sum()),
            post_trim_per_file_stats["n_bases"].sum(),
            post_trim_per_file_stats["min_length"].min(),
            post_trim_per_file_stats["max_length"].max(),
            post_trim_per_file_stats.eval("n_seqs * mean_quality").sum() / total_reads,
        ]
        # some formatting
        for col in ("reads", "bases"):
            basic_summary_stats[col] = (
                basic_summary_stats[col].astype(int).apply(util.si_format)
            )
        for col in ("min_read_length", "max_read_length"):
            basic_summary_stats[col] = (
                basic_summary_stats[col].astype(int).apply(lambda x: f"{x:,d}")
            )
        basic_summary_stats["mean_quality"] = basic_summary_stats["mean_quality"].round(
            1
        )
        basic_summary_stats.index.name = "Condition"
        basic_summary_stats.rename(
            columns=lambda col: col.capitalize().replace("_", " "), inplace=True
        )
        # show the stats in a DataTable
        html_tags.p(
            """
            Some basic stats covering the raw reads and the reads remaining after the
            initial filtering step (based on length and mean quality) as well as after
            downsampling and trimming are illustrated in the table below.
            """
        )
        DataTable.from_pandas(basic_summary_stats)

        # `SeqSummary` plots
        html_tags.p(
            """
            The following plots show the read quality and length distributions as well
            as the base yield after filtering (but before downsampling / trimming) for
            each sample (use the dropdown menu to view the plots for the
            individual samples).
            """
        )
        fastcat.SeqSummary(
            pd.concat((data.fastcat_per_read_stats for data in datasets))
        )

    # summarize bamstats results of all samples for the following report sections
    bamstats_summary = util.summarize_bamstats(datasets)

    # create summary table with one row per sample
    per_sample_summary_table = util.format_per_sample_summary_table(
        bamstats_summary, datasets
    )
    # add sample meta data
    per_sample_summary_table = pd.concat((meta_df, per_sample_summary_table), axis=1)
    per_sample_summary_table.index.name = "Sample alias"

    # summary table with one row per amplicon
    per_amplicon_summary_table = util.format_per_amplicon_summary_table(
        bamstats_summary, datasets
    )

    # section for alignment stats and detected variants summary
    with report.add_section("Summary", "Summary"):
        html_tags.p(
            """
            The two tables below (one per tab) briefly summarize the main results of
            mapping the reads to the provided amplicon references and subsequent variant
            calling.
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
                The following table breaks the results down further (one row per
                sample&ndash;amplicon combination).
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

    # variant tables
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
                .astype({"pos": int, "DP": int, "AB": float})
                .rename(
                    columns={
                        "amp": "amplicon",
                        "pos": "position",
                        "ref": "ref. allele",
                        "alt": "alt. allele",
                        "DP": "Depth",
                        "AB": "Allelic balance",
                    }
                )
                .rename(columns=lambda x: x.capitalize())
            )
            .sort_values(["Sample", "Amplicon", "Position"])
            .set_index(["Sample", "Amplicon"])
        )
        comb_variants["Allelic balance"] = (
            comb_variants["Allelic balance"].astype(str) + "%"
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

    report.write(args.report_fname)
    logger.info(f"Report written to '{args.report_fname}'.")
