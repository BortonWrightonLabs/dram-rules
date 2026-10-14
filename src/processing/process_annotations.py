from __future__ import annotations

import logging

import pandas as pd
import polars as pl

from src.definitions import DEFAULT_FASTA_COLUMN

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
    # regex = "".join([regex for regex in TAXONOMY_RANKS_REGEX.values()])  # regex that might be useful later
    # tree = tax_df["taxonomy"].str.extractall(regex)

    ranks = ["domain", "phylum", "class", "order", "family", "genus", "species"]

    tree_df = tax_df.with_columns(
        pl.col("taxonomy")
        .str.split(";")
        .list.to_struct(fields=ranks)
        .alias("taxonomy_struct")
    ).unnest("taxonomy_struct")

    # generate successive rank edges dynamically
    tax_edge_df = pl.concat(
        [
            tree_df.select(
                pl.col(ranks[i]).alias("source"), pl.col(ranks[i + 1]).alias("target")
            )
            for i in range(len(ranks) - 1)
        ],
        how="vertical",
    ).unique()

    return tax_edge_df, tax_df


def build_tree_pl(
    edge_df,
    source_col: str = "source",
    target_col: str = "target",
    id_cb=None,
):
    """
    Builds a tree structure from an edge DataFrame.

    Parameters:
    - edge_df (DataFrame): The edge DataFrame containing the source and target nodes.
    - source_col (str): The name of the column in edge_df that represents the source nodes. Default is "source".
    - target_col (str): The name of the column in edge_df that represents the target nodes. Default is "target".
    - id_cb (callable): A callback function that generates unique IDs for the tree nodes. Default is None.

    Returns:
    - tree_data (list): A list of dictionaries representing the tree structure.

    """

    if id_cb is None:
        # sensible default: deterministic id
        def id_cb(source, target, parent_id):
            return f"{parent_id}->{target}"

    # --- build adjacency list: source -> list of targets ---
    # This is the key optimization vs filtering inside recursion.
    adj_df = (
        edge_df.select([pl.col(source_col), pl.col(target_col)])
        .group_by(source_col, maintain_order=True)
        .agg(pl.col(target_col).alias("targets"))
    )

    adjacency = dict(zip(adj_df[source_col].to_list(), adj_df["targets"].to_list()))

    # --- roots = sources that never appear as targets ---
    roots = (
        edge_df.select(pl.col(source_col))
        .filter(
            ~pl.col(source_col).is_in(
                edge_df.select(pl.col(target_col)).to_series().implode()
            )
        )
        .unique()
        .to_series()
        .to_list()
    )

    def recurse_tree(source, parent_id=None):
        children = []
        for target in adjacency.get(source, []):
            node_id = id_cb(source, target, parent_id)
            children.append(
                {
                    "label": target,
                    "items": recurse_tree(target, parent_id=node_id),
                    "id": node_id,
                }
            )
        return children

    tree_data = [
        {
            "id": root,
            "label": root,
            "items": recurse_tree(root, parent_id=root),
        }
        for root in roots
    ]

    return tree_data


def build_tax_tree_selected_recurse(tax_tree, flat_tree=None):
    if flat_tree is None:
        flat_tree = []
    for node in tax_tree:
        flat_tree.append(node["id"])
        if "items" in node:
            build_tax_tree_selected_recurse(node["items"], flat_tree)
    return flat_tree
