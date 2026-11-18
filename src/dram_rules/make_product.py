"""The visualize refinement: heatmaps and tables of each rule's Presence, Coverage
and abundance per genome (and per sample, with a mapping file)."""

from __future__ import annotations

import logging
import time
from dataclasses import replace
from pathlib import Path

import panel as pn
import polars as pl

from dram_rules.apps.heatmap import Dashboard
from dram_rules.definitions import (
    DEFAULT_LABEL_COLUMN,
)
from dram_rules.hits import Hits
from dram_rules.processing.process_annotations import (
    build_tax_edge_df_pl,
    build_taxonomy_df_pl,
    build_tree_pl,
)
from dram_rules.rule_parser.rules import (
    CompiledRules,
    get_val_cols,
)

logger = logging.getLogger("dram.viz")

# Significant figures of abundances in the wide tables' row key.
ABUNDANCE_SIG_FIGS = 10

pn.config.reuse_sessions = True
pn.config.global_loading_spinner = True


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
    compiled = CompiledRules.from_rules(
        rules_path=path,
        label_col=label_col,
        alias_col=alias_col,
        group_col=group_col,
        common_rules_path=common_rules_tsv,
        function_files=function_files,
    )
    label = compiled.label_col
    # Downstream (heatmaps, tables) reads rule names from "name".
    df = compiled.df
    if label != DEFAULT_LABEL_COLUMN:
        df = df.drop(DEFAULT_LABEL_COLUMN, strict=False).rename(
            {label: DEFAULT_LABEL_COLUMN}
        )
    if "long_name" not in df.columns:
        df = df.with_columns(long_name=pl.col(DEFAULT_LABEL_COLUMN))
    df = df.with_columns(pl.col("long_name").fill_null(pl.col(DEFAULT_LABEL_COLUMN)))
    return replace(compiled, df=df, label_col=DEFAULT_LABEL_COLUMN)


def visualize(
    sheets: dict[Path, tuple[Path, CompiledRules]], hits: Hits
) -> dict[Path, dict]:
    """Evaluate each sheet, `{output folder: (sheet path, compiled sheet)}`,
    for its heatmaps and tables, at the sample level too when `hits` has a
    mapping file. The taxonomy is loaded once and shared. Returns each
    folder's `Dashboard` arguments; nothing is written."""
    s = time.time()
    genes = hits.genes
    anno_columns = genes.collect_schema().names()
    extra_cols = [
        col for col in ["Completeness", "Contamination"] if col in anno_columns
    ]
    meta = None
    if extra_cols:
        meta = (
            genes.unpivot(
                index="genome", on=extra_cols, variable_name=DEFAULT_LABEL_COLUMN
            )
            .unique()
            .collect()
        )

    tax_df = tax_tree_data = None
    if "taxonomy" in anno_columns:
        tax_df = build_taxonomy_df_pl(genes, "genome")
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
        if hits.samples is not None:
            dfs["sample"] = hits.evaluate(compiled, "sample")
        dfs["genome"] = hits.evaluate(compiled, "genome", additional_cols=["long_name"])
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
    """The heatmap title, from the sheet's file name: gut.tsv or rules_gut.tsv is GUT."""
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
