from __future__ import annotations

import logging

import pandas as pd
import polars as pl

from dram_rules.definitions import DEFAULT_FASTA_COLUMN

logger = logging.getLogger("dram.viz")


def build_taxonomy_df_pl(
    annotations_df: pd.DataFrame, groupby_column=DEFAULT_FASTA_COLUMN
) -> pl.DataFrame:
    cols = [groupby_column, "taxonomy"]
    tax_df = annotations_df[cols].unique()
    tax_df = tax_df.rename({groupby_column: "genome"})
    return tax_df


def build_tax_edge_df_pl(
    tax_df,
) -> tuple[pl.DataFrame, pl.DataFrame]:
    """Edges of the taxonomy tree, one per (parent lineage, child lineage).

    Nodes are whole lineage prefixes ("d__A;p__B"), not rank names, so empty
    ranks shared by many lineages ("f__") don't merge their children. A
    lineage's first rank has a null source; lineages of any length are kept.
    """
    parts = pl.col("taxonomy").str.split(";")
    depth = tax_df.select(parts.list.len().max()).item() or 0
    if not depth:  # no genome has a taxonomy
        return pl.DataFrame(schema={"source": pl.String, "target": pl.String}), tax_df

    def edges_at(i: int) -> pl.DataFrame:
        source = parts.list.head(i).list.join(";") if i else pl.lit(None, pl.String)
        return tax_df.filter(parts.list.len() > i).select(
            source=source, target=parts.list.head(i + 1).list.join(";")
        )

    tax_edge_df = pl.concat([edges_at(i) for i in range(depth)]).unique(
        maintain_order=True
    )
    return tax_edge_df, tax_df


def build_tree_pl(
    edge_df,
    source_col: str = "source",
    target_col: str = "target",
    id_cb=None,
):
    """
    Builds a tree structure from an edge DataFrame (see `build_tax_edge_df_pl`).

    Parameters:
    - edge_df (DataFrame): Edges between ";"-joined lineages; roots have a null source.
    - source_col (str): The name of the column in edge_df that represents the source nodes. Default is "source".
    - target_col (str): The name of the column in edge_df that represents the target nodes. Default is "target".
    - id_cb (callable): id_cb(source, label, parent_id) gives a child's id; a root's id is its label.

    Returns:
    - tree_data (list): A list of dictionaries representing the tree structure.

    """

    if id_cb is None:
        # sensible default: deterministic id
        def id_cb(source, target, parent_id):
            return f"{parent_id}->{target}"

    adj_df = edge_df.group_by(source_col, maintain_order=True).agg(
        pl.col(target_col).alias("targets")
    )
    # Keyed by the source lineage; roots are under the null source.
    adjacency = dict(zip(adj_df[source_col].to_list(), adj_df["targets"].to_list()))

    def recurse_tree(source, parent_id=None):
        children = []
        for target in adjacency.get(source, []):
            label = target.rsplit(";", 1)[-1]
            node_id = label if parent_id is None else id_cb(source, label, parent_id)
            children.append(
                {
                    "label": label,
                    "items": recurse_tree(target, parent_id=node_id),
                    "id": node_id,
                }
            )
        return children

    return recurse_tree(None)
