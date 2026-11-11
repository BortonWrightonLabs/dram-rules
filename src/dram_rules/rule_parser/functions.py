"""Built-in rule functions.

These use the same API as user functions (see dram_rules.rules_api), so this
file doubles as a set of examples. `not`, `&` and `|` are part of the rule
language itself (see rules.py), since they work on row filters as well as
per-sample values.
"""

from __future__ import annotations

from typing import Literal

import numpy as np
import polars as pl

from dram_rules.rules_api import (
    Column,
    Group,
    Op,
    PerSample,
    RowFilter,
    Rows,
    rule_function,
)


@rule_function
def percent(n: int, group: Group) -> PerSample:
    """At least `n` percent of the group's parts are present."""
    return group.present.mean(axis=1) >= n / 100.0


@rule_function
def at_least(
    n: int, group: Group, mode: Literal["PRESENCE", "COUNTS"] = "PRESENCE"
) -> PerSample:
    """PRESENCE: at least `n` of the group's parts are present.
    COUNTS: the parts' hit counts add up to at least `n`."""
    hits = group.counts if mode == "COUNTS" else group.present
    return hits.sum(axis=1) >= n


@rule_function
def tax(rows: Rows, label: str, regex: bool = False) -> PerSample:
    """Any of the sample's rows has a taxonomy containing `label` (plain text
    unless `regex=true`). Without a taxonomy column, no sample matches."""
    if "taxonomy" not in rows.df.columns:
        return np.zeros(len(rows.samples), dtype=bool)
    matches = pl.col("taxonomy").cast(pl.String).str.contains(label, literal=not regex)
    return rows.per_sample(matches.any(), fill=False)


@rule_function
def filter_contains(col: Column, text: str, regex: bool = False) -> RowFilter:
    """Rows whose `col` contains `text` (plain text unless `regex=true`)."""
    return pl.col(col).cast(pl.String).str.contains(text, literal=not regex)


@rule_function
def filter_compare(col: Column, op: Op, value: float) -> RowFilter:
    """Rows whose numeric `col` compares to `value` with `op`, e.g. `>=`."""
    return op(pl.col(col).cast(pl.Float64), value)


@rule_function
def count_rows(rows: Rows, op: Op, n: float) -> PerSample:
    """The sample's number of rows compares to `n` with `op`."""
    return op(rows.per_sample(pl.len(), fill=0), n)


@rule_function
def sum_column(rows: Rows, col: Column, op: Op, value: float) -> PerSample:
    """The total of `col` over the sample's rows compares to `value` with `op`."""
    return op(rows.per_sample(pl.col(col).cast(pl.Float64).sum(), fill=0.0), value)
