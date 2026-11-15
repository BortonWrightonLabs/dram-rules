"""The visualize step: heatmaps and tables of each rule's Presence, Coverage
and abundance per genome (and per sample, with a mapping file)."""

from __future__ import annotations

import logging
import time
from pathlib import Path

import panel as pn
import polars as pl

from dram_rules.apps.heatmap import Dashboard
from dram_rules.definitions import (
    BACKUP_ALIAS_COLUMNS,
    BACKUP_GROUPBY_COLUMNS,
    DEFAULT_ALIAS_COLUMN,
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
    CompiledRules,
    RuleError,
    evaluate,
    find_label_column,
    get_val_cols,
)

logger = logging.getLogger("dram.viz")

# Significant figures of abundances in the wide tables' row key.
ABUNDANCE_SIG_FIGS = 10

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


def load_sheet(
    path: Path,
    common_rules_tsv: Path,
    function_files: tuple[str, ...],
    label_col: str | None = None,
    alias_col: str | None = None,
    group_col: str | None = None,
) -> CompiledRules:
    """Compile a visualize or traits sheet. Its label column becomes `name`,
    which the heatmaps read, and `long_name` falls back to the label. The
    label, alias and Group columns are looked up unless given."""
    rules_lf = scan_tsv(path, infer_schema_length=None, truncate_ragged_lines=True)
    rules_cols = rules_lf.collect_schema().names()
    label_column = label_col or find_label_column(rules_cols)
    # Downstream (evaluation, heatmaps, common rules) reads rule names from "name".
    if label_column != DEFAULT_LABEL_COLUMN:
        rules_lf = rules_lf.drop(DEFAULT_LABEL_COLUMN, strict=False).rename(
            {label_column: DEFAULT_LABEL_COLUMN}
        )
    alias_column = alias_col or next(
        (c for c in [DEFAULT_ALIAS_COLUMN, *BACKUP_ALIAS_COLUMNS] if c in rules_cols),
        DEFAULT_ALIAS_COLUMN,
    )
    if "long_name" not in rules_cols:
        rules_lf = rules_lf.with_columns(long_name=pl.col(DEFAULT_LABEL_COLUMN))
    rules_lf = rules_lf.with_columns(
        pl.col("long_name").fill_null(pl.col(DEFAULT_LABEL_COLUMN))
    )
    compiled = CompiledRules.from_rules(
        rules=rules_lf,
        label_col=DEFAULT_LABEL_COLUMN,
        alias_col=alias_column,
        rules_col="rule",
        group_col=group_col,
        common_rules_path=common_rules_tsv,
        function_files=function_files,
    )
    if compiled.group_col is None:
        raise RuleError(
            "No Group column: expected one of"
            f" {[DEFAULT_GROUPBY_COLUMN, *BACKUP_GROUPBY_COLUMNS]}"
        )
    return compiled


def visualize(
    sheets: dict[Path, tuple[Path, CompiledRules]],
    raw_anno: pl.LazyFrame,
    hits: pl.DataFrame,
    mapping: Path | None,
) -> dict[Path, dict]:
    """Evaluate each sheet, `{output folder: (sheet path, compiled sheet)}`,
    for its heatmaps and tables. `hits` are the annotations' Hits for every
    sheet; the mapping file and the taxonomy are loaded once and shared.
    Returns each folder's `Dashboard` arguments; nothing is written."""
    s = time.time()
    anno_columns = raw_anno.collect_schema().names()
    genes = raw_anno.select("query_id", "genome").collect()
    sample_names = genes["genome"].unique().sort().to_list()
    anno_df = hits
    abundance_col = None
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
        abundance_col = "mean_sample_abundance"

    extra_cols = [
        col for col in ["Completeness", "Contamination"] if col in anno_columns
    ]
    meta = None
    if extra_cols:
        meta = (
            raw_anno.unpivot(
                index="genome", on=extra_cols, variable_name=DEFAULT_LABEL_COLUMN
            )
            .unique()
            .collect()
        )

    tax_df = tax_tree_data = None
    if "taxonomy" in anno_columns:
        tax_df = build_taxonomy_df_pl(raw_anno, "genome")
        tax_edge_df, tax_df = build_tax_edge_df_pl(tax_df)
        tax_tree_data = build_tree_pl(
            tax_edge_df,
            id_cb=lambda source, child, parent_id: f"{parent_id};{child}",
        )
    logger.info(
        f"Loaded the shared data for {len(sheets)} sheets in {time.time() - s:.1f} s"
    )

    results = {}
    for output_dir, (rules_tsv, compiled) in sheets.items():
        dfs = {}
        if mapping:
            dfs["sample"] = evaluate(
                compiled,
                sample_hits_df,
                sample_annotations,
                "sample",
                sample_names_mapping,
                abundance_col="sample_abundance",
            )
        dfs["genome"] = evaluate(
            compiled,
            anno_df,
            raw_anno,
            "genome",
            sample_names,
            additional_cols=["long_name"],
            abundance_col=abundance_col,
        )
        if meta is not None:
            # The metadata frame goes first.
            dfs["genome"] = {"Meta": meta, **dfs["genome"]}
        if tax_df is not None:
            for key, df in dfs["genome"].items():
                dfs["genome"][key] = tax_df.join(df, on="genome", how="left")
        results[output_dir] = {
            "dfs": dfs,
            "taxanomy_tree_data": tax_tree_data,
            "title": title(rules_tsv),
        }
    return results


def save_visualized(
    sheets: dict[Path, tuple[Path, CompiledRules]],
    results: dict[Path, dict],
    save_needed_features: bool,
) -> None:
    """Write each sheet's tables and heatmaps (and Gene IDs) to its folder."""
    for output_dir, (rules_tsv, compiled) in sheets.items():
        output_dir.mkdir(parents=True, exist_ok=True)
        write_tables(results[output_dir]["dfs"], output_dir)
        if save_needed_features:
            (output_dir / "needed_features.txt").write_text(
                "\n".join(compiled.needed_features_raw)
            )
        write_heatmaps(Dashboard(**results[output_dir], output_dir=output_dir))
        logger.info(f"Visualized {rules_tsv} in {output_dir}")


def serve(results: dict[Path, dict], port: int, save: bool) -> None:
    """Serve the live dashboard of the (one) sheet; each browser session gets
    its own Dashboard, saving downloads to the sheet's folder if `save`."""
    [(output_dir, args)] = results.items()
    output_dir = output_dir if save else None
    pn.serve(lambda: Dashboard(**args, output_dir=output_dir), port=port)


def write_tables(dfs: dict[str, dict[str, pl.DataFrame]], output_dir: Path) -> None:
    """The vertical and wide tables, per Evaluation Level."""
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
        # Abundances are part of the row key, so float noise in their last bits
        # (from summing in a different order) would split a rule's row.
        wide = df.with_columns(
            pl.col(c).round_sig_figs(ABUNDANCE_SIG_FIGS)
            for c in ("mean_sample_abundance", "sample_abundance")
            if c in df.columns
        )
        wide.pivot(
            on=[y_col],
            values=val_cols,
            # drop taxonomy because taxonomy and genome have a relationship, so either we would
            # need to concatenate columns awkwardly, or drop it, or make a taxonomy column which would
            # blow up the dataframe
            index=[
                col
                for col in wide.columns
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


def title(rules_tsv: Path) -> str:
    """The heatmap title, from the sheet's file name: rules_gut.tsv is GUT."""
    rule_type = (
        Path(rules_tsv)
        .stem.lower()
        .removeprefix("rule_")
        .removeprefix("rules_")
        .removeprefix("distill_")
        .removesuffix("_rules")
        .removesuffix("_rule")
    ).upper()
    return f"DRAM {rule_type} VISUALIZATION"


def write_heatmaps(dashboard: Dashboard) -> None:
    """Save every heatmap the dashboard can show, as HTML."""
    column_options = dashboard.column_options.options
    current_y = dashboard.y_col
    current_color_col = dashboard.color_col
    for y_col, c_cols in column_options.items():
        for c_col in c_cols:
            if y_col == current_y and c_col == current_color_col:
                continue
            dashboard.y_col = y_col
            dashboard.color_col = c_col
            dashboard.update_plot()
            dashboard.download_heatmap()
