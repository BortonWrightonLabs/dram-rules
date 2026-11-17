"""The summarize refinement: Gene Counts per genome for Summarize Sheets, written as
Summary Collections (one workbook with a worksheet per Group, plus a TSV).

Ported from DRAM2's `bin/distill.py` (ADR 0002: dram-rules is standalone)."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Literal

import polars as pl
from xlsxwriter import Workbook

from dram_rules.hits import Hits
from dram_rules.rule_parser.rules import CompiledRules, check_gene_level

logger = logging.getLogger("dram.rules")

# Sheet columns a Summary Collection keeps, after the label and before the Gene Counts.
SUMMARY_COLUMNS = ["gene_description", "module", "header", "subheader"]
# A Collection's Group column (the worksheet names); distill's name for it.
GROUP_COLUMN = "topic_ecosystem"
SORT_COLUMNS = ["header", "subheader", "module"]  # then the label, as distill sorts

# How Gene Abundance aggregates over a rule's genes in each sample.
Abundance = Literal["mean", "sum"]
ABUNDANCES = ("mean", "sum")

SUMMARIZE_SHEETS_DIR = Path(__file__).parent / "data/summarize_sheets"
# Built-in Summarize Sheets by their --topics/--ecos names, in distill's order.
TOPICS = {
    "assim": "assimilation_and_cofactor_metabolism.tsv",
    "cell": "cellular_machinery.tsv",
    "energy": "energy_acquisition_bioenergetics.tsv",
    "env": "environmental_interaction_and_adaptation.tsv",
}
ECOSYSTEMS = {
    "ag": "distill_ag.tsv",
    "eng_sys": "distill_eng_sys.tsv",
    "bgc": "distill_bgc.tsv",
    "gut": "distill_gut.tsv",
    "aquatic": "distill_aquatic.tsv",
}


def summary_collections(
    topics: list[str], ecos: list[str], custom: list[Path]
) -> dict[str, list[Path]]:
    """Each Summary Collection's name and sheets: the chosen built-in topic
    and ecosystem sheets, and the custom sheets."""
    return {
        "topic_summary": [
            SUMMARIZE_SHEETS_DIR / f for t, f in TOPICS.items() if t in topics
        ],
        "ecosystem_summary": [
            SUMMARIZE_SHEETS_DIR / f for e, f in ECOSYSTEMS.items() if e in ecos
        ],
        "custom_summary": custom,
    }


def load_summarize_sheet(
    path: Path,
    common_rules_tsv: Path,
    function_files: tuple[str, ...],
    label_col: str | None = None,
    alias_col: str | None = None,
    group_col: str | None = None,
) -> CompiledRules:
    """Compile a Summarize Sheet; coverage settings in it are an error. The
    label, alias and Group columns are looked up unless given."""
    compiled = CompiledRules.from_rules(
        rules_path=path,
        label_col=label_col,
        alias_col=alias_col,
        group_col=group_col,
        common_rules_path=common_rules_tsv,
        function_files=function_files,
    )
    check_gene_level(compiled)
    return compiled


def summarize_sheet(
    compiled: CompiledRules, hits: Hits, abundance: Abundance = "mean"
) -> pl.DataFrame | None:
    """One row per rule, in sheet order: the label, the Group (as GROUP_COLUMN),
    the SUMMARY_COLUMNS the sheet has, then one Gene Count column per genome.
    With a mapping file, then one `<abundance>_abundance_<sample>` Gene
    Abundance column per sample. Rules left out by `include_when` are absent;
    None if no rule is left."""
    label, group = compiled.label_col, compiled.group_col
    dfs = hits.evaluate(compiled, "gene")
    if not dfs:
        return None
    rules = pl.concat(dfs.values())
    counts = rules.select(label, pl.col("gene_count").struct.unnest())
    if "abundance" in rules.columns:
        columns = {f"{abundance}_abundance_{s}": s for s in hits.samples}
        sums = [
            pl.col("abundance").struct.field(s).alias(c) for c, s in columns.items()
        ]
        if abundance == "mean":
            genes = pl.lit(rules["gene_count"].struct.unnest().sum_horizontal())
            # A rule that holds for no gene gets 0, not NaN.
            sums = [pl.when(genes > 0).then(s / genes).otherwise(0.0) for s in sums]
        counts = counts.with_columns(rules.select(*sums))
    return (
        compiled.df.drop_nulls(label)
        .select(
            label,
            pl.col(group).alias(GROUP_COLUMN),
            *(c for c in SUMMARY_COLUMNS if c in compiled.df.columns),
        )
        .join(counts, on=label, how="inner", maintain_order="left")
    )


def summarize(
    collections: dict[str, list[CompiledRules]],
    hits: Hits,
    abundance: Abundance = "mean",
) -> dict[str, pl.DataFrame | None]:
    """Each Summary Collection's table: the `summarize_sheet` rows of its
    sheets, one sheet after another. The label column keeps the sheets' name,
    or is `name` when they disagree. None when no rule is left. `abundance`
    picks how Gene Abundance aggregates, when Hits was loaded with a mapping file."""
    tables = {}
    for name, sheets in collections.items():
        frames = [
            df
            for c in sheets
            if (df := summarize_sheet(c, hits, abundance)) is not None
        ]
        if len({df.columns[0] for df in frames}) > 1:
            frames = [df.rename({df.columns[0]: "name"}) for df in frames]
        tables[name] = pl.concat(frames, how="diagonal_relaxed") if frames else None
    return tables


def write_collection(df: pl.DataFrame | None, name: str, out_dir: Path) -> None:
    """Write a `summarize` table as `<name>.xlsx`, a worksheet per Group, and
    `<name>.tsv`. Nothing is written for None."""
    if df is None:
        return
    label = df.columns[0]
    out_dir.mkdir(parents=True, exist_ok=True)
    df.write_csv(out_dir / f"{name}.tsv", separator="\t")
    sort = [c for c in SORT_COLUMNS if c in df.columns] + [label]
    with Workbook(out_dir / f"{name}.xlsx") as wb:
        for (group,), frame in df.group_by(GROUP_COLUMN, maintain_order=True):
            frame.sort(sort).drop(GROUP_COLUMN).write_excel(
                workbook=wb, worksheet=group, autofit=True
            )
    logger.info(f"Saved {out_dir / name}.xlsx and .tsv")
