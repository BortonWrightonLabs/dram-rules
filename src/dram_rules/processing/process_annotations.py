from __future__ import annotations

import logging

import polars as pl

from dram_rules.definitions import DEFAULT_FASTA_COLUMN

logger = logging.getLogger("dram.viz")

# Non-sample columns that feature-count style mapping files can carry.
MAPPING_META_COLUMNS = ["Chr", "Start", "End", "Strand", "Length", "KO", "Description"]


def scan_tsv(path, *, drop_empty_rows: bool = False, **kwargs) -> pl.LazyFrame:
    """Lazily read a TSV. `kwargs` go to `pl.scan_csv`."""
    lf = pl.scan_csv(path, separator="\t", **kwargs)
    return lf.filter(~pl.all_horizontal(pl.all().is_null())) if drop_empty_rows else lf


def load_mapping(
    path, genes: pl.DataFrame
) -> tuple[pl.DataFrame, pl.DataFrame, list[str]]:
    """Read a mapping file of abundances per gene or per genome.

    Its first column holds gene names (`query_id`) or genome names; that
    becomes the key column, named `query_id` or `genome`, of both outputs:
    - means: key, mean_sample_abundance (over all samples, missing as 0)
    - abundances: key, sample, sample_abundance; non-zero abundances only
    Sample columns are every numeric column except MAPPING_META_COLUMNS.
    `genes` is the annotations' (query_id, genome) pairs.
    """
    # Types come from the first rows only, which picks the sample columns
    # fast; every sample column is then read as Float64, so a decimal after
    # those rows still parses and a stray string fails loudly.
    schema = scan_tsv(path, infer_schema_length=10_000).collect_schema()
    first, *rest = schema.names()
    samples = [
        c for c in rest if c not in MAPPING_META_COLUMNS and schema[c].is_numeric()
    ]
    wide = scan_tsv(
        path,
        infer_schema_length=10_000,
        schema_overrides={first: pl.String, **{c: pl.Float64 for c in samples}},
    ).select(first, *samples)
    means, abundances = pl.collect_all(
        [
            wide.select(
                first,
                mean_sample_abundance=pl.mean_horizontal(pl.col(samples).fill_null(0)),
            ),
            wide.unpivot(
                index=first,
                variable_name="sample",
                value_name="sample_abundance",
            ).filter(pl.col("sample_abundance").fill_null(0) != 0),
        ],
        engine="streaming",
    )

    names = means[first]
    if names.is_duplicated().any():
        raise ValueError(
            f"Mapping file has repeated names in its first column '{first}'"
        )
    if names.is_in(genes["query_id"].implode()).all():
        key = "query_id"
        if genes.unique()["query_id"].is_duplicated().any():
            raise ValueError("A query_id in the annotations belongs to several genomes")
    elif names.is_in(genes["genome"].implode()).all():
        key = "genome"
    else:
        raise ValueError(
            "First Column in Mapping file can't be mapped to raw annotation file either to query_id column (gene level) or input_fasta (fasta file name, stand in for genome)"
        )
    return means.rename({first: key}), abundances.rename({first: key}), samples


def build_taxonomy_df_pl(
    annotations_df: pl.LazyFrame | pl.DataFrame, groupby_column=DEFAULT_FASTA_COLUMN
) -> pl.DataFrame:
    cols = [groupby_column, "taxonomy"]
    tax_df = annotations_df.lazy().select(cols).unique().collect()
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
