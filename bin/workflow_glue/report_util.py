"""Common functions and baseclass for holding report data."""
import dominate.tags as html_tags
from ezcharts.layout.snippets import Progress
import pandas as pd
import pysam
import si_prefix


class ReportDataSet:
    """Holds the data required for the report for an individual sample."""

    def __init__(self, data_dir):
        """Parse data for sample found in `data_dir`."""
        self.data_dir = data_dir
        self.sample_alias = data_dir.name
        # the below will be set to `False` when one or more of the input files are empty
        # (i.e. they only contain the CSV header)
        self.all_inputs_valid = True

        # read the fastcat per-file stats (these include all reads pre-filtering)
        self.fastcat_per_file_stats = self.read_csv(
            data_dir / "per-file-stats.tsv", sep="\t"
        )
        # read the fastcat per-read stats (these include all reads present after
        # filtering and before downsampling or trimming)
        self.fastcat_per_read_stats = self.read_csv(
            data_dir / "per-read-stats.tsv.gz", sep="\t", index_col="read_id"
        )
        # read the post-trimming fastcat per-file stats
        self.post_trim_per_file_stats = self.read_csv(
            data_dir / "porechopped-per-file-stats.tsv", sep="\t"
        )

        # read bamstats files
        self.bamstats = self.read_csv(data_dir / "bamstats.tsv", sep="\t", index_col=0)
        self.bamstats_flagstat = self.read_csv(
            data_dir / "bamstats-flagstat.tsv", sep="\t", index_col=0
        )
        # get the list of amplicons with at least one primary alignment
        self.detected_amplicons = (
            list(self.bamstats_flagstat.query("primary > 0").index)
            if self.bamstats_flagstat is not None
            else []
        )
        if not self.detected_amplicons:
            # no amplicons were deteced (i.e. not a single read mapped to any of the
            # sequences in the reference file)
            self.all_inputs_valid = False

        # read the depth data
        self.depth = self.read_csv(data_dir / "per-window-depth.tsv.gz", sep="\t")

        # if there is a VCF, read the variants
        vcf_path = data_dir / "medaka.annotated.vcf.gz"
        self.variants = pd.DataFrame(
            columns=["amp", "pos", "ref", "alt", "type", "filter", "DP"]
        ).set_index(["amp", "pos"])
        if vcf_path.exists():
            # check that the VCF only contains a single sample
            vcf_file = pysam.VariantFile(vcf_path)
            if len(vcf_file.header.samples) != 1:
                raise ValueError(f"Found multiple samples in {vcf_path}.")
            # process variants
            for entry in vcf_file.fetch():
                # we should only have one alt allele
                if len(entry.alts) != 1:
                    raise ValueError(
                        f"Found entry with multiple ALT alleles in {vcf_path}."
                    )
                self.variants.loc[(entry.chrom, entry.pos), :] = [
                    entry.ref,
                    entry.alts[0],
                    # the field `entry.alleles_variant_types` contains a tuple `('REF',
                    # alt-type)`
                    entry.alleles_variant_types[1],
                    entry.filter[0].name,
                    entry.info["DP"],
                ]

        # read the QC summaries if they exist
        qc_summary = data_dir / "qc-summary.tsv"
        if qc_summary.exists():
            self.qc_summary = pd.read_csv(qc_summary, sep="\t")
        else:
            self.qc_summary = None

    def __repr__(self):
        """Return human-readable string (simply the sample alias)."""
        return f"ReportDataSet('{self.sample_alias}')"

    def read_csv(self, file, **kwargs):
        """Read a CSV and set `self.all_inputs_valid` to `False` if it's empty.

        This is a wrapper around `pd.read_csv` that sets `self.all_inputs_valid = False`
        if the CSV contains only a header line and no data.
        """
        try:
            df = pd.read_csv(file, **kwargs)
        except FileNotFoundError:
            self.all_inputs_valid = False
            return None
        if df.empty:
            self.all_inputs_valid = False
        return df

    def get_basic_summary_stats(self, de_novo):
        """Collect basic summary stats.

        This is used in the "at a glance" section at the beginning of the report.
        """
        # create series to hold the summary stats and add stats that are needed in
        # de-novo and variant calling mode
        extra_fields = (
            ["unmapped_reads", "consensus_length"]
            if de_novo
            else ["amplicons", "overall_mean_depth", "min_mean_depth", "snps", "indels"]
        )
        basic_summary = pd.Series(
            0,
            index=["reads", "bases", "mean_length", "mean_quality"] + extra_fields,
        )
        basic_summary["reads"] = self.post_trim_per_file_stats["n_seqs"].sum()
        basic_summary["bases"] = self.post_trim_per_file_stats["n_bases"].sum()
        basic_summary["mean_quality"] = (
            self.post_trim_per_file_stats.eval("n_seqs * mean_quality")
            / basic_summary["reads"]
        )
        basic_summary["mean_length"] = self.post_trim_per_file_stats.eval(
            "n_bases / n_seqs"
        )
        # return early if this sample didn't have all valid inputs (the post-trim
        # per-file stats should always be there)
        if not self.all_inputs_valid:
            return basic_summary
        # if in de-novo mode, add the de-novo-only stats and return
        if de_novo:
            basic_summary["unmapped_reads"] = self.bamstats.eval('ref == "*"').sum()
            basic_summary["consensus_length"] = self.depth["end"].iloc[-1]
            return basic_summary
        basic_summary["amplicons"] = len(self.detected_amplicons)
        # now get the mean depths. we got different window sizes for the different
        # amplicons. thus, we'll
        # 1. calculate the "summed depths" first (i.e. the sum of window size * mean
        #    depth for each window) for each detected amplicons
        # 2. get the ref lengths and divide the "summed depths" by them
        grpd = self.depth.query("ref.isin(@self.detected_amplicons)").groupby("ref")
        summed_depths = grpd.apply(lambda df: df.eval("(end - start) * depth").sum())
        ref_lengths = grpd.last()["end"]
        mean_depths = summed_depths / ref_lengths
        basic_summary["overall_mean_depth"] = summed_depths.sum() / ref_lengths.sum()
        basic_summary["min_mean_depth"] = mean_depths.min()
        basic_summary["snps"] = self.variants.query("type == 'SNP'").shape[0]
        basic_summary["indels"] = self.variants.query("type == 'INDEL'").shape[0]
        return basic_summary


def custom_progress_bar(value, display_val=None):
    """Make a bootstrap progress bar with different styling to be used in tables.

    One thing to keep in mind is that the constructor for `Progress` uses `dominate`
    context managers and therefore the progress bar will be added to the enclosing
    context regardless of whether this desired or not.
    """
    p = Progress(
        value_min=0,
        value_max=100,
        value_now=round(value, 1),
        display_val=display_val,
        height=30,
    )
    # set some style attributes to get a bordered light-grey bar with black font
    p.children[0].children[0].children[0]["style"] += ";color:black"
    p.children[0].children[0]["style"] += ";background-color:rgb(210, 210, 210)"
    p.children[0]["style"] += ";width:120px;border:1px solid black"
    return p


def summarize_bamstats(datasets):
    """Summarize the bamstats data from all the samples.

    The resulting DataFrame contains summed / averaged values for each sample--amplicon
    combination.
    """
    summary_stats = pd.DataFrame(
        columns=[
            "sample",
            "amplicon",
            "reads",
            "bases",
            "median_read_length",
            "mean_cov",
            "mean_acc",
            "variants",
            "indels",
        ],
    ).set_index(["sample", "amplicon"])

    for d in datasets:
        # order refs so that the row for unmapped is at the end for each sample
        refs = d.bamstats["ref"].unique()
        if "*" in refs:
            refs = [ref for ref in refs if ref != "*"] + ["*"]
        for amplicon in refs:
            df = d.bamstats.query("ref == @amplicon")
            n_reads = df.shape[0]
            n_bases = df["read_length"].sum()
            med_read_length = df["read_length"].median()
            mean_cov = df["ref_coverage"].mean()
            mean_acc = df["acc"].mean()
            try:
                variants = d.variants.loc[amplicon].shape[0]
                indels = d.variants.loc[amplicon].query('type == "INDEL"').shape[0]
            except KeyError:
                variants = indels = 0
            summary_stats.loc[(d.sample_alias, amplicon), :] = [
                n_reads,
                n_bases,
                med_read_length,
                mean_cov,
                mean_acc,
                variants,
                indels,
            ]
    return summary_stats


def si_format(n):
    """Use `si-prefix`, but don't add a decimal when `n` is smaller than 1000.

    `si_prefix.si_format()` always adds the decimals specified by `precision`, even when
    calling it on a small integer. This wrapper prevents that.
    """
    return si_prefix.si_format(n, precision=0 if n < 1000 else 1)


def format_stats_table(df):
    """Pretty-format a table with bamstats data.

    Individual columns will be formatted differently depending on their column name
    (e.g. values in columns named `reads` will be replaced by the HTML code for a
    progress bar showing the number of reads as well as the proportion of the total).
    The resulting DataFrame is expected to be passed to `DataTable.from_pandas()` to be
    added to the report.

    Note that progress bars will implicitly added to the enclosing `dominate` context,
    which is something that needs to be taken care of.
    """
    total_reads = df["reads"].sum()
    total_bases = df["bases"].sum()

    # used later for making sure the table columns containing HTML code for progress
    # bars are sortable
    num_leading_zeros = len(str(df.shape[0])) + 1

    # initialise resulting DataFrame
    df_for_table = df.copy()

    # check if columns with a certain type of data are present (e.g. read length, mean
    # coverage, etc.); formatting of the individual columns depends on column name (e.g.
    # for 'reads' there will be a progress bar and the column 'median_read_length' will
    # be cast to `int`)
    if "reads" in df_for_table.columns:
        # progress bars for number of reads column
        for n, (idx, n_reads) in enumerate(
            df["reads"].astype(int).sort_values().items()
        ):
            # get percentage of total number of reads
            perc = 100 * n_reads / total_reads if total_reads != 0 else 0
            # some number formatting
            display_val = f"{si_format(n_reads)} ({perc:.0f}%)"
            # add the progress bar to the table cell
            df_for_table.loc[idx, "reads"] = str(
                # We want the column to be sortable by the actual values even though it
                # contains the HTML progress bar elements. A hacky way to achieve this
                # is by prepending a hidden span containing the zero-filled index from
                # the `enumerate()` above before the progress bar to maintain order. The
                # table will be sorted based on this number and ignore whatever is in
                # the progress bar.
                html_tags.span(str(n).zfill(num_leading_zeros), hidden=True)
            ) + str(custom_progress_bar(perc, display_val))
    if "bases" in df_for_table.columns:
        # progress bars for number of bases column
        for n, (idx, n_bases) in enumerate(df["bases"].sort_values().items()):
            perc = 100 * n_bases / total_bases if total_bases != 0 else 0
            display_val = f"{si_format(n_bases)} ({perc:.0f}%)"
            df_for_table.loc[idx, "bases"] = str(
                # same trick to keep the column sortable as above
                html_tags.span(str(n).zfill(num_leading_zeros), hidden=True)
            ) + str(custom_progress_bar(perc, display_val))
    if "variants" in df_for_table.columns:
        # make sure the input DataFrame got a column with indels as well
        if "indels" not in df.columns:
            raise ValueError(
                "`format_stats_table()` was asked to format a 'variants' column, "
                "but there was no accompanying 'indels' column."
            )
        # variants summary column (simply showing total number of variants and indels)
        df_for_table["variants"] = df.apply(
            lambda row: f"{int(row['variants'])} ({int(row['indels'])})", axis=1
        )
        # we don't need the `indels` column
        df_for_table.drop(columns="indels", inplace=True)

    # round mean coverage and alignment accuracy
    for col in ("mean_cov", "mean_acc"):
        if col in df_for_table.columns:
            df_for_table[col] = df[col].round(1)

    # median read length, number of amplicons, and number of samples should be ints
    for col in ("median_read_length", "amplicons", "samples"):
        if col in df_for_table.columns:
            df_for_table[col] = df[col].astype(int)

    if "unmapped" in df_for_table.columns:
        # progress bars for number of unmapped reads column; same as for "reads" and
        # "bases" above
        for idx, (sample, (n_reads, n_unmapped)) in enumerate(
            df[["reads", "unmapped"]].astype(int).sort_values("unmapped").iterrows()
        ):
            perc = 100 * n_unmapped / n_reads
            display_val = f"{si_format(n_unmapped)} ({perc:.0f}%)"
            df_for_table.loc[sample, "unmapped"] = str(
                html_tags.span(str(idx).zfill(num_leading_zeros), hidden=True)
            ) + str(custom_progress_bar(perc, display_val))

    # Format column names, replace '*' in the index with 'Unmapped', and return
    if df_for_table.index.name is not None:
        df_for_table.index.name = df_for_table.index.name.capitalize()
    return (
        df_for_table.rename(index={"*": "Unmapped"})
        .rename_axis(index=str.capitalize)
        .rename(
            columns=lambda col: col.replace("variants", "variants (indels)")
            .replace("mean_acc", "mean_acc.")
            .replace("mean_cov", "mean_cov.")
            .replace("_", " ")
            .capitalize()
        )
    )


def format_per_sample_summary_table(summary_stats, datasets):
    """Summarize and format bamstats results on a per-sample level.

    `summary_stats` contains summary stats for each sample--amplicon combination (i.e.
    it's got a multi-index of `['sample', 'amplicon']`)
    """
    # initialise results dataframe with zeros
    per_sample_summary_stats = pd.DataFrame(
        0,
        index=summary_stats.index.unique("sample"),
        columns=[
            "reads",
            "bases",
            "median_read_length",
            "amplicons",
            "unmapped",
            "variants",
            "indels",
        ],
    )
    # group by sample and aggregate the stats
    for sample, df in summary_stats.groupby("sample"):
        per_sample_summary_stats.loc[sample, ["reads", "bases"]] = df[
            ["reads", "bases"]
        ].sum()
        # for median read length we have a look at the whole bamstats DataFrame
        # belonging to that sample
        (dataset,) = [x for x in datasets if x.sample_alias == sample]
        per_sample_summary_stats.loc[sample, "median_read_length"] = dataset.bamstats[
            "read_length"
        ].median()
        per_sample_summary_stats.loc[sample, "amplicons"] = (
            df.query("reads > 0 and amplicon != '*'")
        ).shape[0]
        per_sample_summary_stats.loc[sample, "unmapped"] = (
            df.loc[(sample, "*"), "reads"]
            if "*" in df.index.get_level_values("amplicon")
            else 0
        )
        per_sample_summary_stats.loc[sample, ["variants", "indels"]] = df[
            ["variants", "indels"]
        ].sum()
    per_sample_summary_stats.fillna(0, inplace=True)
    return format_stats_table(per_sample_summary_stats)


def format_per_amplicon_summary_table(summary_stats, datasets):
    """Summarize and format bamstats results on a per-amplicon level.

    `summary_stats` contains summary stats for each sample--amplicon combination (i.e.
    it's got a multi-index of `['sample', 'amplicon']`)
    """
    per_amplicon_summary_stats = pd.DataFrame(
        0,
        index=summary_stats.index.unique("amplicon"),
        columns=[
            "reads",
            "bases",
            "median_read_length",
            "samples",
            "mean_cov",
            "mean_acc",
        ],
    )
    # group by amplicon and aggregate the stats
    for amplicon, df in summary_stats.groupby("amplicon"):
        per_amplicon_summary_stats.loc[amplicon, ["reads", "bases"]] = df[
            ["reads", "bases"]
        ].sum()
        # for median read length we have to combine all reads belonging to this amplicon
        # from the individual per-sample datasets
        per_amplicon_summary_stats.loc[amplicon, "median_read_length"] = pd.concat(
            (d.bamstats.query(f"ref == '{amplicon}'")["read_length"] for d in datasets)
        ).median()
        per_amplicon_summary_stats.loc[amplicon, "samples"] = (df["reads"] > 0).sum()
        # multiply the mean cov + acc. per sample--amplicon combination by the number of
        # reads for that combo and divide by total number of reads for this amplicon to
        # get the overall per-amplicon average
        per_amplicon_summary_stats.loc[amplicon, ["mean_cov", "mean_acc"]] = (
            df[["mean_cov", "mean_acc"]].multiply(df["reads"], axis="index").sum()
            / df["reads"].sum()
        )
        per_amplicon_summary_stats.loc[amplicon, ["variants", "indels"]] = df[
            ["variants", "indels"]
        ].sum()
    per_amplicon_summary_stats.fillna(0, inplace=True)
    return format_stats_table(per_amplicon_summary_stats)


def format_de_novo_summary_table(summary_stats, datasets):
    """Summarize and format per-sample bamstats results in de-novo mode.

    `summary_stats` contains summary stats for each sample--amplicon combination. In the
    de-novo case this means there are two rows for each sample. One for the reads that
    were successfully re-aligned against the consensus (`summary_stats.loc[sample_name,
    sample_name]`) and one for the reads that did not re-align against the consensus
    (`summary_stats[sample_name, "*"]`). This function combines both entries into one
    row containing all the relevant information for the summary table.
    """
    # initialise results dataframe with zeros
    de_novo_summary_stats = pd.DataFrame(
        0,
        index=summary_stats.index.unique("sample"),
        columns=[
            "reads",
            "unmapped",
            "bases",
            "median_read_length",
            "mean_cov",
            "mean_acc",
        ],
    )
    # group by sample and aggregate the stats
    for sample, df in summary_stats.groupby("sample"):
        de_novo_summary_stats.loc[sample, ["reads", "bases"]] = df[
            ["reads", "bases"]
        ].sum()
        try:
            de_novo_summary_stats.loc[sample, "unmapped"] = df.loc[
                (sample, "*"), "reads"
            ]
        except KeyError:
            de_novo_summary_stats.loc[sample, "unmapped"] = 0
        # for median read length we have a look at the whole bamstats DataFrame
        # belonging to that sample
        (dataset,) = [x for x in datasets if x.sample_alias == sample]
        de_novo_summary_stats.loc[sample, "median_read_length"] = dataset.bamstats[
            "read_length"
        ].median()
        # there should only be one amplicon per sample
        (amplicon,) = set(df.reset_index()["amplicon"]) - set("*")
        de_novo_summary_stats.loc[sample, ["mean_cov", "mean_acc"]] = df.loc[
            (sample, amplicon), ["mean_cov", "mean_acc"]
        ]
    de_novo_summary_stats.fillna(0, inplace=True)
    return format_stats_table(de_novo_summary_stats)


def format_number_and_plural(num, singular):
    """Format a word as plural unless `num` is `1`."""
    word = singular if num == 1 else singular + "s"
    return f"{num} {word}"
