"""Hits: every Hit in one run's annotations (tRNAs included), loaded once and
evaluated by every Refinement at any Evaluation Level."""

from __future__ import annotations

import logging
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path

import polars as pl

from dram_rules.definitions import BACKUP_FASTA_COLUMNS, DEFAULT_FASTA_COLUMN
from dram_rules.processing.process_annotations import load_mapping, scan_tsv
from dram_rules.rule_parser.rules import (
    ID_PARSERS,
    CompiledRules,
    EvaluationLevel,
    evaluate,
    prepare_present_map_df,
)

logger = logging.getLogger("dram.rules")


@dataclass(frozen=True)
class Hits:
    """Build with `Hits.load`; evaluate sheets with `evaluate`."""

    hits: pl.DataFrame  # query_id, genome, hit: one row per Hit
    rows: pl.LazyFrame  # annotation and tRNA rows, for rules that read rows
    genes: pl.LazyFrame  # annotation rows only, for Meta and taxonomy
    genomes: list[str]  # annotation and tRNA genomes, sorted
    samples: list[str] | None = None  # the mapping file's samples, if given
    genome_hits: pl.DataFrame | None = None  # hits + mean_sample_abundance
    sample_hits: pl.DataFrame | None = None  # hits + sample, sample_abundance
    sample_rows: pl.LazyFrame | None = None  # rows, once per sample holding them

    @classmethod
    def load(
        cls,
        annotations: Path,
        sheets: Iterable[CompiledRules],
        *,
        fasta_column: str = DEFAULT_FASTA_COLUMN,
        trna: Path | None = None,
        mapping: Path | None = None,
    ) -> Hits:
        """The Hits of every Gene ID the `sheets` need. The genome column is
        `fasta_column` or the first of its backups; `trna` is a tRNAscan file
        and `mapping` a mapping file of abundances per gene or per genome.
        Raises ValueError for annotations without a genome column, or a
        mapping file that matches neither genes nor genomes."""
        sheets = list(sheets)
        genes = scan_tsv(annotations, infer_schema_length=10_000)
        names = genes.collect_schema().names()
        candidates = [fasta_column, *BACKUP_FASTA_COLUMNS]
        column = next((c for c in candidates if c in names), None)
        if column is None:
            raise ValueError(
                f"Annotations have none of the genome columns {candidates}"
            )
        genes = genes.rename({column: "genome"}).with_columns(
            pl.col("genome").cast(pl.String)
        )

        needed = set().union(*(s.needed_features for s in sheets))
        hits = prepare_present_map_df(
            genes,
            count_col="query_id",
            group_col="genome",
            besthit_cols=list(ID_PARSERS),
            needed_features=needed,
        )
        trna_rows = read_trna(trna) if trna else None
        if trna_rows is not None:
            trna_hits = prepare_present_map_df(
                trna_rows,
                count_col="query_id",
                group_col="genome",
                besthit_cols=[],
                needed_features=needed,
            )
            hits = pl.concat([hits, trna_hits], how="vertical_relaxed")
        logger.info(f"Loaded the annotations' Hits once, for {len(hits)} hits")

        # Rules that read rows would otherwise read the file again for every sheet.
        # After the Hits: building them from the in-memory copy peaked higher.
        if any(s.reads_rows for s in sheets):
            genes = genes.collect().lazy()
        rows = genes
        if trna_rows is not None:
            rows = pl.concat([genes, trna_rows], how="diagonal_relaxed")
        genomes = rows.select("genome").unique().sort("genome").collect()
        genomes = genomes["genome"].to_list()
        if mapping is None:
            return cls(hits, rows, genes, genomes)

        pairs = genes.select("query_id", "genome").collect()
        means, abundances, samples = load_mapping(mapping, pairs)
        key = means.columns[0]  # query_id or genome
        return cls(
            hits,
            rows,
            genes,
            genomes,
            samples=samples,
            genome_hits=hits.join(means, on=key),
            sample_hits=hits.join(abundances, on=key),
            # A gene belongs to every sample where its abundance isn't 0.
            sample_rows=rows.join(abundances.lazy().select(key, "sample"), on=key),
        )

    def evaluate(
        self,
        compiled: CompiledRules,
        level: EvaluationLevel,
        additional_cols: list[str] | None = None,
    ) -> dict[str, pl.DataFrame]:
        """`rules.evaluate` at `level`, with abundance at the genome and sample
        levels when a mapping file was given. The sample level needs one."""
        if level == "sample":
            if self.samples is None:
                raise ValueError("The sample level needs a mapping file")
            return evaluate(
                compiled,
                self.sample_hits,
                self.sample_rows,
                "sample",
                self.samples,
                additional_cols=additional_cols,
                abundance_col="sample_abundance",
            )
        if level == "genome" and self.samples is not None:
            return evaluate(
                compiled,
                self.genome_hits,
                self.rows,
                "genome",
                self.genomes,
                additional_cols=additional_cols,
                abundance_col="mean_sample_abundance",
            )
        return evaluate(
            compiled,
            self.hits,
            self.rows,
            level,
            self.genomes,
            additional_cols=additional_cols,
        )


def read_trna(trna_path: Path) -> pl.LazyFrame | None:
    """tRNAscan rows (`input_fasta`, `query_id`, `gene_id`, as DRAM collects
    them) as annotation rows: `genome`, `query_id`, and the tRNA's Gene ID in a
    `trna_id` column. None for DRAM's file for "no tRNAs".

    Each tRNA is its own gene. Its `query_id` is its contig, which tRNAs on one
    contig share, so the row number is appended to keep them apart."""
    trna = pl.read_csv(trna_path, separator="\t", infer_schema_length=None)
    if trna.columns == ["NULL"] or trna.is_empty():
        return None
    return trna.lazy().select(
        genome=pl.col("input_fasta").cast(pl.String),
        query_id=pl.format(
            "{}_trna_{}",
            pl.col("query_id").cast(pl.String).str.strip_chars(),
            pl.int_range(pl.len()),
        ),
        trna_id=pl.col("gene_id"),
    )
