"""The summarize step: Gene Counts per genome for Summarize Sheets, written as
Summary Collections (one workbook with a worksheet per Group, plus a TSV).

Ported from DRAM2's `bin/distill.py` (ADR 0002: dram-rules is standalone)."""

from __future__ import annotations

import logging
from pathlib import Path

import polars as pl
from xlsxwriter import Workbook

from dram_rules.rule_parser.rules import CompiledRules, evaluate

logger = logging.getLogger("dram.rules")

# Sheet columns a Summary Collection keeps, after the label and before the Gene Counts.
SUMMARY_COLUMNS = ["gene_description", "module", "header", "subheader"]
# A Collection's Group column (the worksheet names); distill's name for it.
GROUP_COLUMN = "topic_ecosystem"
SORT_COLUMNS = ["header", "subheader", "module"]  # then the label, as distill sorts


def add_trna(annotations: pl.LazyFrame, trna_path: Path) -> pl.LazyFrame:
    """Add tRNAscan rows (`input_fasta`, `query_id`, `gene_id`, as DRAM collects
    them) to the annotations, with the tRNA's Gene ID in a `trna_id` column.

    Each tRNA is its own gene. Its `query_id` is its contig, which tRNAs on one
    contig share, so the row number is appended to keep them apart."""
    trna = pl.read_csv(trna_path, separator="\t", infer_schema_length=None)
    if trna.columns == ["NULL"] or trna.is_empty():  # DRAM's file for "no tRNAs"
        return annotations
    trna = trna.lazy().select(
        genome=pl.col("input_fasta").cast(pl.String),
        query_id=pl.format(
            "{}_trna_{}",
            pl.col("query_id").cast(pl.String).str.strip_chars(),
            pl.int_range(pl.len()),
        ),
        trna_id=pl.col("gene_id"),
    )
    return pl.concat([annotations, trna], how="diagonal_relaxed")


def summarize_sheet(
    compiled: CompiledRules,
    hits: pl.DataFrame,
    annotations: pl.LazyFrame,
    genomes: list[str],
) -> pl.DataFrame | None:
    """One row per rule, in sheet order: the label, the Group (as GROUP_COLUMN),
    the SUMMARY_COLUMNS the sheet has, then one Gene Count column per genome.
    Rules left out by `include_when` are absent; None if no rule is left."""
    label, group = compiled.label_col, compiled.group_col
    dfs = evaluate(compiled, hits, annotations, "gene", genomes)
    if not dfs:
        return None
    counts = pl.concat(dfs.values()).pivot(
        on="genome", index=label, values="gene_count"
    )
    return (
        compiled.df.drop_nulls(label)
        .select(
            label,
            pl.col(group).alias(GROUP_COLUMN),
            *(c for c in SUMMARY_COLUMNS if c in compiled.df.columns),
        )
        .join(counts, on=label, how="inner", maintain_order="left")
    )


def write_collection(
    sheets: list[pl.DataFrame | None], name: str, out_dir: Path
) -> None:
    """Write `summarize_sheet` frames as one Summary Collection: `<name>.xlsx`,
    a worksheet per Group, and `<name>.tsv`. The label column keeps the sheets'
    name, or is `name` when they disagree. Nothing is written without sheets."""
    sheets = [df for df in sheets if df is not None]
    if not sheets:
        return
    labels = {df.columns[0] for df in sheets}
    if len(labels) > 1:
        sheets = [df.rename({df.columns[0]: "name"}) for df in sheets]
    label = sheets[0].columns[0]
    df = pl.concat(sheets, how="diagonal_relaxed")
    out_dir.mkdir(parents=True, exist_ok=True)
    df.write_csv(out_dir / f"{name}.tsv", separator="\t")
    sort = [c for c in SORT_COLUMNS if c in df.columns] + [label]
    with Workbook(out_dir / f"{name}.xlsx") as wb:
        for (group,), frame in df.group_by(GROUP_COLUMN, maintain_order=True):
            frame.sort(sort).drop(GROUP_COLUMN).write_excel(
                workbook=wb, worksheet=group, autofit=True
            )
    logger.info(f"Saved {out_dir / name}.xlsx and .tsv")
