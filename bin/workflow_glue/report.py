"""Create workflow report."""
import itertools
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
from ezcharts.plots.util import choose_palette
import pandas as pd
import pysam

from . import report_util as util  # noqa: ABS101

from .util import get_named_logger, wf_parser  # noqa: ABS101


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
        help="FASTA file with reference sequences for the individual amplicons.",
    )
    parser.add_argument(
        "--sample-sheet",
        type=Path,
        help="Sample sheet with metadata.",
    )
    parser.add_argument(
        "--downsampling-size",
        type=int,
        default=0,
        help="Downsampled number of reads (per sample).",
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
        "--wf-version",
        required=True,
        help="version of the executed workflow",
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
        "Workflow Amplicon Sequencing report",
        "wf-amplicon",
        args.params,
        args.versions,
        args.wf_version,
    )
    populate_report(report, metadata, datasets, args.reference, args.downsampling_size)

    report.write(args.report_fname)
    logger.info(f"Report written to '{args.report_fname}'.")


def populate_report(report, metadata, all_datasets, ref_fasta, downsampling_size):
    """Fill the report with content."""
    # put whether we got a ref file into a variable
    de_novo = ref_fasta is None
    analysis_type = "de-novo consensus" if de_novo else "variant calling"
    samples = [d.sample_alias for d in all_datasets]
    # check if any samples are missing any of the required inputs and filter them out
    # (we will mention them in the report below); this could only happen if there were
    # no reads for that particular sample
    bad_datasets = [d for d in all_datasets if not d.all_inputs_valid]
    datasets = [d for d in all_datasets if d.all_inputs_valid]
    ref_seqs = []
    if de_novo:
        # make sure that there was only one amplicon per sample in de-novo mode
        for d in datasets:
            if len(d.detected_amplicons) != 1:
                raise ValueError(
                    "Found unexpected number of amplicons "
                    f"in de-novo consensus mode for sample {d.sample_alias}."
                )
        # some samples might have failed because of consensus QC and some for
        # other reasons
        bad_datasets_failed_qc = [d for d in bad_datasets if d.qc_summary is not None]
        bad_datasets_other = [d for d in bad_datasets if d.qc_summary is None]

    else:
        # in variant calling mode we can have many amplicons
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
                f"""The {analysis_type} pipeline failed for all samples.
                Therefore, only a limited report is available. Consider relaxing the
                filtering criteria (e.g. using a lower value for
                """,
                html_tags.code("--min_read_qual"),
                """) as this is sometimes caused by filtering out all reads from the
                input data.
                """,
            )
        # show the pre-processing stats and consensus QC summaries (if in de-novo mode);
        # omit all other sections
        preprocessing_section(report, bad_datasets)

        if de_novo:
            de_novo_qc_section(report, bad_datasets_failed_qc)

        return

    # intro and "at a glance" stats
    with report.add_section("Introduction", "Intro."):
        intro_str = (
            "This report contains tables and plots to help interpret the results "
            "of wf-amplicon. The workflow was run in"
        )
        if de_novo:
            # we're in de-novo mode
            html_tags.p(
                intro_str,
                html_tags.b("de-novo consensus mode"),
                """. The individual sections of the report summarize the outcomes of the
                different steps of the workflow (read filtering, consensus generation,
                re-mapping against the consensus for depth of coverage calculation).
                """,
            )
        else:
            # we're in ref mode
            html_tags.p(
                intro_str,
                html_tags.b("variant calling mode"),
                """. The individual sections of the report summarize the
                outcomes of the different steps of the workflow (read filtering, mapping
                against the reference file containing the amplicon sequences, variant
                calling).
                """,
            )
            with html_tags.div(cls="alert alert-info"):
                html_tags.p(
                    html_tags.b("Note: "),
                    "If the sequence IDs in the reference file contained special "
                    "characters, they were replaced with underscores."
                )
        html_tags.p("The input data contained:")
        # brief summary of how many samples / refs there were
        for name, items in zip(["sample", "amplicon"], [samples, ref_seqs]):
            if name == "amplicon" and de_novo:
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
            if de_novo:
                # mention datasets with failed QC and those that failed for other
                # reasons separately
                failed_qc_aliases = [d.sample_alias for d in bad_datasets_failed_qc]
                failed_assembly_aliases = [d.sample_alias for d in bad_datasets_other]
                if failed_assembly_aliases:
                    html_tags.b("Warning:")
                    html_tags.p(
                        f"The {analysis_type} pipeline failed for the following ",
                        (
                            "sample (and it"
                            if (n_failed := len(failed_assembly_aliases)) == 1
                            else f"{n_failed} samples (and they"
                        ),
                        " will be omitted in some sections of the report):",
                        html_tags.br(),
                        ", ".join(failed_assembly_aliases),
                    )
                    html_tags.p(
                        "A potential reason for this might be that not enough reads "
                        "remained for analysis after filtering and preprocessing."
                    )
                if failed_qc_aliases:
                    html_tags.b("Warning:")
                    html_tags.p(
                        " Consensus quality control failed for the following ",
                        util.format_number_and_plural(len(failed_qc_aliases), "sample"),
                        ":",
                        html_tags.br(),
                        ", ".join(failed_qc_aliases),
                    )
            else:
                html_tags.b("Warning:")
                html_tags.p(
                    f"The {analysis_type} pipeline failed for the following ",
                    util.format_number_and_plural(len(bad_datasets), "sample"),
                    " (and these will be omitted in some sections of the report):",
                    html_tags.br(),
                    ", ".join([d.sample_alias for d in bad_datasets]),
                )
                html_tags.p(
                    "A potential reason for this might be that not enough reads "
                    "remained for analysis after filtering and preprocessing."
                )
            html_tags.br()
        if downsampling_size:
            with html_tags.div(cls="alert alert-warning"):
                html_tags.p(
                    html_tags.b("Note: "),
                    "The data was downsampled to "
                    f"{downsampling_size} reads per sample.",
                )
        html_tags.h5("At a glance")
        html_tags.p(
            """
            Key results for the individual samples are shown below. You can use the
            dropdown menu to view the results for a different sample.
            """
        )
        # "at a glance" stats in cards (with one dropdown tab per sample)
        tabs = Tabs()
        with tabs.add_dropdown_menu():
            for d in all_datasets:
                basic_summary = d.get_basic_summary_stats(de_novo)
                with tabs.add_dropdown_tab(d.sample_alias):
                    if d in bad_datasets:
                        if de_novo and d.sample_alias in failed_qc_aliases:
                            html_tags.p(
                                html_tags.b("Note: "),
                                "Consensus QC failed for this sample.",
                            )
                        else:
                            html_tags.p(
                                html_tags.b("Note: "),
                                "The analysis pipeline failed for this sample "
                                "(e.g. potentially due to too few reads).",
                            )
                    # add values + titles of cards for stats reported in both modes
                    # (de-novo and variant calling) first
                    stats_card_items = [
                        (f'{basic_summary["reads"]:,g}', "Reads"),
                        (f'{basic_summary["bases"]:,g}', "Bases"),
                        (basic_summary["mean_length"].round(1), "Mean length"),
                        (basic_summary["mean_quality"].round(1), "Mean quality"),
                    ]
                    if de_novo:
                        # add the de-novo-only cards
                        stats_card_items += [
                            (f'{basic_summary["unmapped_reads"]:,g}', "Unmapped reads"),
                            (
                                f'{basic_summary["consensus_length"]:,g}',
                                "Consensus length",
                            ),
                        ]
                    else:
                        # Define the number of amplicons we expect to find. If a `ref`
                        # column was in the sample sheet, it is the number of ref IDs
                        # listed there. Otherwise, it is the total number of references.
                        # Note that the ref column can also contain an empty cell.
                        try:
                            n_exp_amps = len(
                                metadata.loc[d.sample_alias, "ref"].split()
                            )
                        except (AttributeError, KeyError):
                            n_exp_amps = len(ref_seqs)
                        n_amplicons_stats_str = (
                            f'{basic_summary["amplicons"]:g} / {n_exp_amps}'
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
    preprocessing_section(report, all_datasets)

    # brief section for QC summaries (de-novo mode only)
    if de_novo:
        de_novo_qc_section(
            report,
            sorted(datasets + bad_datasets_failed_qc, key=lambda d: d.sample_alias),
        )

    # summarize bamstats results of all samples for the following report sections
    bamstats_summary = util.summarize_bamstats(datasets)

    # only keep the sample amplicon combinations that had alignments
    bamstats_summary.query("reads > 0", inplace=True)

    # section for alignment stats (per sample and per-amplicon if variant calling mode;
    # otherwise only per-sample)
    if de_novo:
        with report.add_section("Re-alignment summary", "Re-align"):
            # only show re-alignment stats: make a single table with two rows per
            # sample; one with stats of successfully re-aligned reads and one with stats
            # of unaligned reads
            html_tags.p(
                """
                The table below summarizes the main results of re-mapping the reads of
                each barcode against the generated consensus. Percentages of unmapped
                reads are relative to the number of reads for that particular sample.
                Other percentages are relative to the total number of reads / bases
                including all samples.
                """
            )
            with dom_util.container() as c:
                de_novo_summary_table = util.format_de_novo_summary_table(
                    bamstats_summary, datasets
                )
            c.clear()
            DataTable.from_pandas(de_novo_summary_table)
    else:
        with report.add_section("Summary", "Summary"):
            with dom_util.container() as c:
                # create summary table with one row per sample
                per_sample_summary_table = util.format_per_sample_summary_table(
                    bamstats_summary, datasets
                )
                # add sample meta data
                per_sample_summary_table = (
                    # for this table, drop `ref` column from sample sheet if present
                    pd.concat(
                        (
                            metadata.drop(columns="ref", errors="ignore")
                            if metadata is not None
                            else None,
                            per_sample_summary_table,
                        ),
                        axis=1,
                    )
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
                variant calling. Percentages of unmapped reads are relative to the
                number of reads for that particular sample. Other percentages are
                relative to the total number of reads / bases including all samples.
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
            depth_dfs = []
            for d in datasets:
                depth_df = d.depth.assign(
                    pos=lambda df: df[["start", "end"]].mean(axis=1),
                    sample=d.sample_alias,
                )[["ref", "pos", "depth", "sample"]]
                # if an amplicon has 0 depth along its whole length, drop it
                nonzero_depth_amplicons = [
                    k
                    for k, v in depth_df.groupby("ref")["depth"].sum().items()
                    if v > 0
                ]
                if not nonzero_depth_amplicons:
                    continue
                depth_df = depth_df.query("ref.isin(@nonzero_depth_amplicons)")
                depth_dfs.append(depth_df)
            # in de-novo mode, we re-aligned against the consensus and use the sample
            # name as amplicon name
            if de_novo:
                depth_dfs = [df.eval("ref = sample") for df in depth_dfs]
            # combine the data and group by amplicon
            per_amplicon_depths = pd.concat(depth_dfs).groupby("ref")
            # if there are amplicons with more than one sample (i.e. plots with more
            # than one line), define colours for each sample; otherwise, just go with
            # the default
            palette = None
            if not de_novo:
                for amplicon, depth_df in per_amplicon_depths:
                    if len(depth_df["sample"].unique()) > 1:
                        palette = dict(zip(samples, itertools.cycle(choose_palette())))
                        break
            # draw the plots
            for amplicon, depth_df in per_amplicon_depths:
                with tabs.add_dropdown_tab(amplicon):
                    plt = ezc.lineplot(
                        data=depth_df.round(2),
                        x="pos",
                        y="depth",
                        hue="sample",
                        palette=palette,
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
    if de_novo:
        return
    with report.add_section("Variants", "Variants"):
        html_tags.p(
            "Haploid variant calling was performed with Medaka. Variants with low ",
            "depth (i.e. smaller than ",
            html_tags.code("--min_coverage"),
            ') are shown under the "Low depth" tab. The numbers in the ',
            '"depth" column relate to the sequencing depth used to perform '
            "variant calling.",
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
    fastcat_comb_length_hists = pd.concat((d.fastcat_length_hist for d in datasets))
    fastcat_comb_qual_hists = pd.concat((d.fastcat_qual_hist for d in datasets))
    preprocessing_stats.loc["Filtered"] = [
        n_reads := fastcat_comb_length_hists['count'].sum(),
        fastcat_comb_length_hists.eval('start * count').sum(),
        fastcat_comb_length_hists['start'].min(),
        fastcat_comb_length_hists['start'].max(),
        fastcat_comb_qual_hists.eval('(start + end) / 2 * count').sum() / n_reads
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

        # return early if there were no per-read stats (i.e. there were no reads left
        # after filtering)
        if preprocessing_stats.loc["Filtered", "Reads"] == 0:
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
        # give all histograms to the fastcat SeqSummary
        fastcat.SeqSummary(
            tuple(d.data_dir for d in datasets),
            sample_names=tuple(d.sample_alias for d in datasets),
        )


def de_novo_qc_section(report, datasets):
    """Add a section with tables of post-assembly consensus QC stats."""
    # make sure datasets are in the right order
    datasets = sorted(datasets, key=lambda d: d.sample_alias)
    with report.add_section("Quality Control", "QC"):
        html_tags.p(
            "After creating a draft assembly (either with ",
            html_tags.a("Miniasm", href="https://github.com/lh3/miniasm"),
            " or with ",
            html_tags.a("SPOA", href="https://github.com/rvaser/spoa"),
            """) basic quality control is performed on the consensus candidates (i.e.
            the contigs in the assembly). The QC stats for the contigs produced by the
            assembly pipeline are listed in the table below. If there are no contigs for
            the "miniasm" method in the table, the Miniasm step either
            failed or produced contigs that were shorter than the """,
            html_tags.code("--force_spoa_length_threshold"),
            """ parameter (and thus SPOA was run regardless of assembly quality). If
            there are no contigs for "spoa", SPOA was not run as a contig produced by
            Miniasm passed all filters and was selected. If the sample is missing
            altogether, both assembly methods failed.
            """,
        )
        tabs = Tabs()
        with tabs.add_dropdown_menu():
            for d in datasets:
                qc_summary = d.qc_summary
                qc_summary["fail_reason"] = qc_summary["fail_reason"].fillna("-")
                qc_summary["mean_depth"] = qc_summary["mean_depth"].round(1)
                qc_summary["length"] = [f"{x:,d}" for x in qc_summary["length"]]
                # give the contigs more meaningful names
                method_counts = {"miniasm": 0, "spoa": 0}
                for idx, row in qc_summary.iterrows():
                    method = row["method"]
                    method_counts[method] += 1
                    qc_summary.loc[idx, "contig"] = f"{method}_{method_counts[method]}"

                qc_summary.set_index("contig", inplace=True)
                qc_summary.index.name = "Contig"
                qc_summary.rename(
                    columns=lambda col: col.replace("_", " ").capitalize(), inplace=True
                )
                with tabs.add_dropdown_tab(d.sample_alias):
                    DataTable.from_pandas(qc_summary)
