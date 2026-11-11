"""
===================
DRAM Visualizations
===================

Script that generates a product visualization from the DRAM output.
"""

from __future__ import annotations

import logging
from pathlib import Path

import click
import panel as pn
import polars as pl

from dram_rules.apps.heatmap import Dashboard
from dram_rules.definitions import (
    BACKUP_ALIAS_COLUMNS,
    BACKUP_FASTA_COLUMNS,
    BACKUP_GROUPBY_COLUMNS,
    BACKUP_LABEL_COLUMNS,
    DEFAULT_ALIAS_COLUMN,
    DEFAULT_FASTA_COLUMN,
    DEFAULT_GROUPBY_COLUMN,
    DEFAULT_LABEL_COLUMN,
)
from dram_rules.processing.process_annotations import (
    build_tax_edge_df_pl,
    build_taxonomy_df_pl,
    build_tree_pl,
    load_mapping,
    scan_tsv,
)
from dram_rules.rule_parser.rules import (
    ID_PARSERS,
    CompiledRules,
    build_present_map,
    evaluate_cycles,
    prepare_present_map_df,
    get_val_cols,
)

logger = logging.getLogger("dram.viz")

pn.config.reuse_sessions = True
pn.config.global_loading_spinner = True

RULES_SYSTEMS = {
    "legacy": Path(__file__).parent / "data/rules.tsv",
    "ag": Path(__file__).parent / "data/rules_ag.tsv",
    "bgc": Path(__file__).parent / "data/rules_bgc.tsv",
    "aquatic": Path(__file__).parent / "data/rules_aquatic.tsv",
    "eng_sys": Path(__file__).parent / "data/rules_eng_sys.tsv",
    "gut": Path(__file__).parent / "data/rules_gut.tsv",
    "default": Path(__file__).parent / "data/rules_gut.tsv",
}


def get_column_name(
    name: str,
    backup_names: list[str],
    name_of_data: str,
    df: pl.LazyFrame | pl.DataFrame = None,
    cols: list[str] = None,
    raise_error: bool = True,
) -> str:
    assert df is not None or cols is not None, (
        "Must provide either a dataframe or a list of columns to search for the column name"
    )
    if cols is None:
        if isinstance(df, pl.LazyFrame):
            cols = df.collect_schema().names()
        else:
            cols = df.columns
    if name in cols:
        return name
    for backup_name in backup_names:
        if backup_name in cols:
            return backup_name
    if raise_error:
        raise ValueError(
            f"None of the following column names were found in the {name_of_data}: {[name] + backup_names}"
        )
    else:
        return name


@click.command()
@click.option("--annotations", "-a", type=Path, help="Path to the annotations tsv file")
@click.option(
    "--fasta_column",
    "-f",
    type=str,
    default=DEFAULT_FASTA_COLUMN,
    help="Column from annotations file to use as fasta names",
)
@click.option(
    "--output_dir",
    "-o",
    type=Path,
    help="Path to the output directory",
    default=Path.cwd().resolve(),
)
@click.option(
    "--dashboard",
    "-d",
    is_flag=True,
    show_default=True,
    default=False,
)
@click.option(
    "--rules_system",
    "-rs",
    type=click.Choice(list(RULES_SYSTEMS.keys()), case_sensitive=True),
    help="Choose alternative in-built rules for the visualization. This option cannot be specified if `--rule_tsv` is specified. Will default to default if neither this or rule_tsv is provided",
)
@click.option(
    "--rules_tsv",
    "-r",
    type=click.Path(exists=True),
    help="This is an optional path to a custom rules file with strict formatting. It will over write the original rules file that is stored with the script. This option cannot be specified if `--rule_tsv` is specified.",
)
@click.option(
    "--common_rules_tsv",
    "-c",
    type=click.Path(exists=True),
    help="Path to rules sheet with common rules to be added to all rules sheet. Useful for common aliases across multiple rules. If none is provided, will use default common rules system defined in Rule Parser.",
)
@click.option(
    "--functions",
    "-fn",
    "function_files",
    type=click.Path(exists=True, dir_okay=False),
    multiple=True,
    help="Python file of extra rule functions (decorated with @rule_function from dram_rules.rules_api). Can be given more than once.",
)
@click.option(
    "--mapping",
    "-m",
    type=Path,
    help="Path to mapping file of gene abundances in samples. Must have a column named 'Geneid' that matches the query_id column in the annotations file.",
)
@click.option(
    "--label_column",
    "-l",
    type=str,
    default=DEFAULT_LABEL_COLUMN,
    help="Column in rules file to use as labels in visualization",
)
@click.option(
    "--alias_column",
    "-ac",
    type=str,
    default=DEFAULT_ALIAS_COLUMN,
    help="Column from rules file to use as alias in rules. This allows you to alias a rule on one line and then use that alias in a different line.",
)
@click.option(
    "--group_column",
    "-g",
    type=str,
    help="Column in rules file to group by in visualization (creates separate heatmap for each unique value in this column)",
    default=DEFAULT_GROUPBY_COLUMN,
)
@click.option(
    "--port",
    default=5006,
    type=click.INT,
    help="Port to launch dashboard on if launching dashboard",
)
@click.option(
    "--save_dataframes",
    "-sd",
    is_flag=True,
    show_default=True,
    default=False,
)
@click.option(
    "--save_needed_features",
    "-sf",
    is_flag=True,
    show_default=True,
    default=False,
    help="Save txt file of the gene_ids in your rules needed. Potentially useful for debugging.",
)
def main(
    annotations,
    fasta_column,
    output_dir,
    dashboard,
    rules_system,
    rules_tsv,
    common_rules_tsv,
    function_files,
    mapping,
    label_column,
    alias_column,
    group_column,
    port,
    save_dataframes,
    save_needed_features,
):
    """
    Make a product heatmap visualization from the DRAM output.
    """
    import time

    s = time.time()
    if common_rules_tsv is None:
        common_rules_tsv = Path(__file__).parent / "data/rules_common.tsv"

    if rules_system and rules_tsv:
        raise click.BadArgumentUsage(
            "You may only supply either a rules system or a custom rules tsv, not both"
        )
    if rules_system:
        rules_tsv = RULES_SYSTEMS[rules_system]
    if not rules_system and not rules_tsv:
        rules_tsv = RULES_SYSTEMS["default"]

    rules_lf = scan_tsv(rules_tsv, infer_schema_length=None, truncate_ragged_lines=True)

    rules_cols = rules_lf.collect_schema().names()
    group_column = get_column_name(
        group_column, BACKUP_GROUPBY_COLUMNS, name_of_data="rules", cols=rules_cols
    )
    label_column = get_column_name(
        label_column, BACKUP_LABEL_COLUMNS, name_of_data="rules", cols=rules_cols
    )
    # Downstream (evaluation, heatmaps, common rules) reads rule names from "name".
    if label_column != DEFAULT_LABEL_COLUMN:
        rules_lf = rules_lf.drop(DEFAULT_LABEL_COLUMN, strict=False).rename(
            {label_column: DEFAULT_LABEL_COLUMN}
        )
        label_column = DEFAULT_LABEL_COLUMN
    alias_column = get_column_name(
        alias_column,
        BACKUP_ALIAS_COLUMNS,
        name_of_data="rules",
        cols=rules_cols,
        raise_error=False,
    )
    if "long_name" not in rules_cols:
        rules_lf = rules_lf.with_columns(long_name=pl.col(label_column))

    rules_lf = rules_lf.with_columns(
        pl.col("long_name").fill_null(pl.col(label_column)).alias("long_name")
    )

    # Lazy: each use below reads only the columns it needs, and the full rows
    # load only if a rule reads them (Rows.df).
    raw_anno = scan_tsv(annotations, infer_schema_length=10_000)
    fasta_column = get_column_name(
        fasta_column, BACKUP_FASTA_COLUMNS, name_of_data="annotations", df=raw_anno
    )
    # If fasta names are all ints (the samples from users were ints themselves),
    # we need to cast to strings for later
    raw_anno = raw_anno.rename({fasta_column: "genome"}).with_columns(
        pl.col("genome").cast(pl.String)
    )
    anno_columns = raw_anno.collect_schema().names()

    kw = dict(
        rules=rules_lf,
        label_col=label_column,
        alias_col=alias_column,
        rules_col="rule",
        common_rules_path=common_rules_tsv,
        function_files=function_files,
    )
    # kw = dict(rules_path=rules_path, label_col="module", parent_col=alias_column, rules_col="rule")
    compiled = CompiledRules.from_rules(**kw)
    logger.info(f"Compiled rules in {time.time() - s} seconds")

    genes = raw_anno.select("query_id", "genome").collect()
    sample_names = genes["genome"].unique().sort().to_list()
    anno_df = prepare_present_map_df(
        df=raw_anno,
        count_col="query_id",
        group_col="genome",
        besthit_cols=list(ID_PARSERS),
        needed_features=compiled.needed_features,
    )
    dfs = {}
    eval_cycles_kw = {}
    if mapping:
        means, abundances, sample_names_mapping = load_mapping(mapping, genes)
        key = means.columns[0]  # query_id or genome
        # query_id, genome, hit, sample, sample_abundance: each sample's hits
        sample_hits_df = anno_df.join(abundances, on=key)
        anno_df = anno_df.join(means, on=key)
        # Annotation rows per sample, for rules that read rows (count_rows,
        # tax, ...): a gene belongs to every sample where its abundance isn't 0.
        sample_annotations = raw_anno.join(
            abundances.lazy().select(key, "sample"), on=key
        )

        present_map = build_present_map(
            df=sample_hits_df,
            sample_col="sample",
            needed_features=compiled.needed_features,
            sample_names=sample_names_mapping,
        )
        dfs["sample"] = evaluate_cycles(
            compiled=compiled,
            samples=sample_names_mapping,
            present_map=present_map,
            annotations=sample_annotations,
            sample_col="sample",
            group_col=group_column,
            anno_df=sample_hits_df,
            value_col="sample_abundance",
        )
        eval_cycles_kw["anno_df"] = anno_df
        eval_cycles_kw["value_col"] = "mean_sample_abundance"

    present_map = build_present_map(
        df=anno_df,
        sample_col="genome",
        needed_features=compiled.needed_features,
        sample_names=sample_names,
    )

    logger.info(f"Built present map in {time.time() - s} seconds")

    dfs["genome"] = evaluate_cycles(
        compiled=compiled,
        samples=sample_names,
        present_map=present_map,
        annotations=raw_anno,
        sample_col="genome",
        additional_cols=["long_name"],
        group_col=group_column,
        **eval_cycles_kw,
    )
    logger.info("Evaluated all rules in:")
    logger.info(time.time() - s)

    extra_cols = [
        col for col in ["Completeness", "Contamination"] if col in anno_columns
    ]
    if extra_cols:
        df = (
            raw_anno.unpivot(index="genome", on=extra_cols, variable_name=label_column)
            .unique()
            .collect()
        )
        # df = df.rename({groupby_column: "genome"})
        # We do this reorder the dfs dictionary to ensure that the metadata dataframe is the first one
        d = {"Meta": df}
        d.update(dfs["genome"])
        dfs["genome"] = d

    tax_tree_data = None
    if "taxonomy" in anno_columns:
        tax_df = build_taxonomy_df_pl(raw_anno, "genome")

        tax_edge_df, tax_df = build_tax_edge_df_pl(tax_df)

        tax_tree_data = build_tree_pl(
            tax_edge_df,
            id_cb=lambda source, child, parent_id: f"{parent_id};{child}",
        )

        for key, df in dfs["genome"].items():
            dfs["genome"][key] = tax_df.join(df, on="genome", how="left")

    output_dir.mkdir(parents=True, exist_ok=True)

    if save_dataframes:
        for df_type, dfs_dict in dfs.items():
            y_col = df_type  # the table's row id column: "genome" or "sample"
            df = pl.concat(
                [d for key, d in dfs_dict.items() if key != "Meta"],
                how="diagonal_relaxed",
            )

            val_cols = get_val_cols(df, allow_multiple=True)
            df.select([y_col, pl.all().exclude(y_col)]).write_csv(
                output_dir / f"df_vertical_{df_type}.tsv", separator="\t"
            )
            df.pivot(
                on=[y_col],
                values=val_cols,
                # drop taxonomy because taxonomy and genome have a relationship, so either we would
                # need to concatenate columns awkwardly, or drop it, or make a taxonomy column which would
                # blow up the dataframe
                index=[
                    col
                    for col in df.columns
                    if col
                    not in [
                        *val_cols,
                        y_col,
                        "taxonomy",
                        "gene_ids_present",
                        "steps",
                        "steps_present",
                        "subunits",
                        "subunits_present",
                    ]
                ],
            ).write_csv(output_dir / f"df_{df_type}.tsv", separator="\t")
            logger.info(f"Saved dataframe to {output_dir / f'df_{df_type}.tsv'}")

    if save_needed_features:
        (output_dir / "needed_features.txt").write_text(
            "\n".join(compiled.needed_features_raw)
        )
    rule_type = (
        (Path(rules_tsv).stem)
        .lower()
        .removeprefix("rule_")
        .removeprefix("rules_")
        .removeprefix("distill_")
        .removesuffix("_rules")
        .removesuffix("_rule")
    ).upper()
    kw = dict(
        dfs=dfs,
        taxanomy_tree_data=tax_tree_data,
        output_dir=output_dir,
        title=f"DRAM {rule_type} VISUALIZATION",
        # mapping=bool(mapping),
    )
    logger.info(
        f"Finished all processing in: {time.time() - s} seconds, starting visualization"
    )
    if dashboard:
        pn.serve(
            lambda: Dashboard(**kw),
            port=port,
        )
    else:
        d = Dashboard(**kw)
        column_options = d.column_options.options
        current_y = d.y_col
        current_color_col = d.color_col
        for y_col, c_cols in column_options.items():
            for c_col in c_cols:
                if y_col == current_y and c_col == current_color_col:
                    continue
                d.y_col = y_col
                d.color_col = c_col
                d.update_plot()
                d.download_heatmap()


if __name__ == "__main__":
    main()
