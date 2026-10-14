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
import polars.selectors as cs

from dram_viz.apps.heatmap import Dashboard
from dram_viz.definitions import (
    BACKUP_ALIAS_COLUMNS,
    BACKUP_FASTA_COLUMNS,
    BACKUP_GROUPBY_COLUMNS,
    BACKUP_LABEL_COLUMNS,
    DEFAULT_ALIAS_COLUMN,
    DEFAULT_FASTA_COLUMN,
    DEFAULT_GROUPBY_COLUMN,
    DEFAULT_LABEL_COLUMN,
)
from dram_viz.processing.process_annotations import (
    build_tax_edge_df_pl,
    build_tax_tree_selected_recurse,
    build_taxonomy_df_pl,
    build_tree_pl,
)
from dram_viz.rule_parser.src.rules import (
    ID_EXPR_DICT,
    CompiledRules,
    build_present_map,
    evaluate_cycles,
    prepare_present_map_df,
    get_val_col,
)

logger = logging.getLogger("dram.viz")

pn.config.reuse_sessions = True
pn.config.global_loading_spinner = True

RULES_SYSTEMS = {
    "default": Path(__file__).parent / "data/rules.tsv",
    "ag": Path(__file__).parent / "data/rules_ag.tsv",
    "bgc": Path(__file__).parent / "data/rules_bgc.tsv",
    # "aquatic": Path(__file__).parent / "data/rules_aquatic.tsv",
    # "eng_sys": Path(__file__).parent / "data/rules_eng_sys.tsv",
    "gut": Path(__file__).parent / "data/rules_gut.tsv",
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


def join_present_map_df_to_mapping_df(
    df: pl.DataFrame,
    group_col: str,
    count_col: str,
    besthit_cols: list[str],
    needed_features: set[str],
    additional_cols: list[str] = None,
    mapping_df: pl.DataFrame = None,
):

    sample_names = [
        col for col in mapping_df.columns if col not in [count_col, group_col]
    ]
    hit_col = "hit"
    df = (  # Here we join mapping_df to df, to get the mean value accross all samples
        df.join(mapping_df, on=count_col)
        .with_columns(mean_sample_abundance=pl.mean_horizontal(sample_names))
        .select([count_col, group_col, hit_col, "mean_sample_abundance"])
    )

    mapping_df = (
        # (inner) join anno df to unpivotted mapping df to get `query_id, genome, hit, sample`
        # to link sample and query_id in mapping file to specific hits in anno file.
        # anno file should already be filtered to anno needed features, so this should be pretty small
        df.select([count_col, group_col, hit_col]).join(
            (
                # unpivot mapping df to get `count_col, sample, sample_abundance` as a long table of all nonnull nonzero abundances
                # we consider the sample_abundance the abundance of a specific sample for a specific count_col, which
                # we will match later by hit from the anno df
                mapping_df.unpivot(
                    index=count_col,
                    on=sample_names,
                    variable_name="sample",
                    value_name="sample_abundance",
                ).filter(
                    pl.col("sample_abundance").is_not_null()
                    & (pl.col("sample_abundance") != 0)
                )
            ),
            on=count_col,
            how="inner",
            # validate="1:m"
        )
    )
    return df, mapping_df, sample_names


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
    "--group_colunm",
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
    mapping,
    label_column,
    alias_column,
    group_colunm,
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
        common_rules_tsv = Path(__file__).parent / "rule_parser/src/rules_common.tsv"

    if rules_system and rules_tsv:
        raise click.BadArgumentUsage(
            "You may only supply either a rules system or a custom rules tsv, not both"
        )
    if rules_system:
        rules_tsv = RULES_SYSTEMS[rules_system]
    if not rules_system and not rules_tsv:
        rules_tsv = RULES_SYSTEMS["default"]

    rules_lf = pl.scan_csv(
        rules_tsv, separator="\t", infer_schema_length=None, truncate_ragged_lines=True
    ).fill_null("")

    rules_cols = rules_lf.collect_schema().names()
    group_colunm = get_column_name(
        group_colunm, BACKUP_GROUPBY_COLUMNS, name_of_data="rules", cols=rules_cols
    )
    label_column = get_column_name(
        label_column, BACKUP_LABEL_COLUMNS, name_of_data="rules", cols=rules_cols
    )
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

    raw_anno = pl.read_csv(
        annotations,
        separator="\t",
        infer_schema_length=10_000,
        # columns=list(ID_EXPR_DICT.keys()) + [groupby_column]
    )
    fasta_column = get_column_name(
        fasta_column, BACKUP_FASTA_COLUMNS, name_of_data="annotations", df=raw_anno
    )
    raw_anno = raw_anno.rename({fasta_column: "genome"})
    # If fasta names are all ints (the samples from users were ints themselves),
    # we need to cast to strings for later
    raw_anno = raw_anno.with_columns(pl.col("genome").cast(pl.String))

    kw = dict(
        rules=rules_lf,
        label_col=label_column,
        alias_col=alias_column,
        rules_col="rule",
        allow_visualize_functions=True,
        common_rules_path=common_rules_tsv,
    )
    # kw = dict(rules_path=rules_path, label_col="module", parent_col=alias_column, rules_col="rule")
    compiled = CompiledRules.from_rules(**kw)
    logger.info(f"Compiled rules in {time.time() - s} seconds")

    sample_names = (
        raw_anno.select("genome").unique().sort("genome").to_series().to_list()
    )
    besthit_cols = list(ID_EXPR_DICT.keys())
    dfs = {}
    eval_cycles_kw = {}
    if mapping:
        mapping_df = pl.read_csv(mapping, separator="\t", ignore_errors=True).fill_null(
            0
        )
        # This renames the first column to query_id regardless of initial form (Genome, GeneId, reference, etc.)
        # It then drops known additional meta columns that could be there depending on the format
        # finally does a catch all of keeping only the label column and all floats (sample columns) as a backup
        mapping_df = (
            mapping_df.rename({mapping_df.columns[0]: "query_id"})
            .drop(
                ["Chr", "Start", "End", "Strand", "Length", "KO", "Description"],
                strict=False,
            )
            .select(pl.col("query_id"), cs.numeric())
        )
        # if mapping_df col 0 maps to anno query_id col
        if mapping_df.select(
            pl.col("query_id")
            .is_in(raw_anno.select(pl.col("query_id")).to_series().implode())
            .all()
        ).item():
            mapping_df = mapping_df.join(
                raw_anno.select(["query_id", "genome"]).unique(),
                on="query_id",
                validate="1:1",
            )
        # if mapping_df col 0 maps to anno input_fasta/genome col
        elif mapping_df.select(
            pl.col("query_id")
            .is_in(raw_anno.select(pl.col("genome")).to_series().implode())
            .all()
        ).item():
            mapping_df = mapping_df.rename({"query_id": "genome"}).join(
                raw_anno.select(["query_id", "genome"]).unique(),
                on="genome",
                validate="1:m",
            )
        else:
            raise ValueError(
                "First Column in Mapping file can't be mapped to raw annotation file either to query_id column (gene level) or input_fasta (fasta file name, stand in for genome)"
            )

        anno_df = prepare_present_map_df(
            df=raw_anno,
            count_col="query_id",
            group_col="genome",
            besthit_cols=besthit_cols,
            needed_features=compiled.needed_features,
        )
        anno_df, mapping_df, sample_names_mapping = join_present_map_df_to_mapping_df(
            df=anno_df,
            group_col="genome",
            count_col="query_id",
            besthit_cols=besthit_cols,
            needed_features=compiled.needed_features,
            mapping_df=mapping_df,
        )

        present_map = build_present_map(
            df=mapping_df,
            sample_col="sample",
            needed_features=compiled.needed_features,
            sample_names=sample_names_mapping,
        )
        dfs["sample"] = evaluate_cycles(
            compiled=compiled,
            samples=sample_names_mapping,
            present_map=present_map,
            annotations=raw_anno,
            sample_col="sample",
            # additional_cols=["long_name"],
            group_col=group_colunm,
            anno_df=mapping_df,
            value_col="sample_abundance",
        )
        eval_cycles_kw["anno_df"] = anno_df
        eval_cycles_kw["value_col"] = "mean_sample_abundance"
    else:
        anno_df = prepare_present_map_df(
            df=raw_anno,
            count_col="query_id",
            group_col="genome",
            besthit_cols=besthit_cols,
            needed_features=compiled.needed_features,
        )

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
        group_col=group_colunm,
        **eval_cycles_kw,
    )
    logger.info("Evaluated all rules in:")
    logger.info(time.time() - s)

    extra_cols = [
        col for col in ["Completeness", "Contamination"] if col in raw_anno.columns
    ]
    if extra_cols:
        df = raw_anno.unpivot(
            index="genome", on=extra_cols, variable_name=label_column
        ).unique()
        # df = df.rename({groupby_column: "genome"})
        # We do this reorder the dfs dictionary to ensure that the metadata dataframe is the first one
        d = {"Meta": df}
        d.update(dfs["genome"])
        dfs["genome"] = d

    tax_tree_data = None
    selected_tax_tree = None
    if "taxonomy" in raw_anno:
        tax_df = build_taxonomy_df_pl(raw_anno, "genome")

        tax_edge_df, tax_df = build_tax_edge_df_pl(tax_df)

        tax_tree_data = build_tree_pl(
            tax_edge_df,
            id_cb=lambda source, child, parent_id: f"{parent_id};{child}",
        )
        selected_tax_tree = build_tax_tree_selected_recurse(tax_tree_data)

        for key, df in dfs["genome"].items():
            dfs["genome"][key] = tax_df.join(df, on="genome", how="left")

    del raw_anno

    if save_dataframes:
        for df_type, dfs_dict in dfs.items():
            df = pl.concat(
                [d for key, d in dfs_dict.items() if key != "Meta"],
                how="diagonal_relaxed",
            )

            val_col = get_val_col(df)
            df.pivot(
                on=["genome"],
                values=val_col,
                index=[col for col in df.columns if col not in [val_col, "genome"]],
            ).write_csv(output_dir / f"df_{df_type}.tsv", separator="\t")
            logger.info(f"Saved dataframe to {output_dir / f'df_{df_type}.tsv'}")

    if save_needed_features:
        (output_dir / "needed_features.txt").write_text(
            "\n".join(compiled.needed_features_raw)
        )

    kw = dict(
        dfs=dfs,
        taxanomy_tree_data=tax_tree_data,
        selected_tax_tree=selected_tax_tree,
        output_dir=output_dir,
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
