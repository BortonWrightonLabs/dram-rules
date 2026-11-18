from __future__ import annotations

import logging

import polars as pl

from dram_rules.definitions import DEFAULT_FASTA_COLUMN, SAMPLE_COLUMN_PREFIX

logger = logging.getLogger("dram.viz")

# Non-sample columns that feature-count style mapping files can carry.
MAPPING_META_COLUMNS = ["Chr", "Start", "End", "Strand", "Length", "KO", "Description"]


def scan_tsv(path, *, drop_empty_rows: bool = False, **kwargs) -> pl.LazyFrame:
    """Lazily read a TSV. `kwargs` go to `pl.scan_csv`."""
    lf = pl.scan_csv(path, separator="\t", **kwargs)
    return lf.filter(~pl.all_horizontal(pl.all().is_null())) if drop_empty_rows else lf


def read_mapping(path) -> tuple[pl.LazyFrame, str, list[str]]:
    """A mapping file's first column (names, as text) and its Sample columns.

    Returns the wide table, the first column's name and the Samples: every
    numeric column except MAPPING_META_COLUMNS, as Float64, in file order.
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
    return wide, first, samples


def mapping_key(names: pl.Series, genes: pl.DataFrame) -> str:
    """ "query_id" when every mapping name is a gene, else "genome" when every
    one is a genome. `genes` is the annotations' (query_id, genome) pairs.
    Raises ValueError for repeated names or names matching neither."""
    if names.is_duplicated().any():
        raise ValueError(
            f"Mapping file has repeated names in its first column '{names.name}'"
        )
    if names.is_in(genes["query_id"].implode()).all():
        if genes.unique()["query_id"].is_duplicated().any():
            raise ValueError("A query_id in the annotations belongs to several genomes")
        return "query_id"
    if names.is_in(genes["genome"].implode()).all():
        return "genome"
    raise ValueError(
        "First Column in Mapping file can't be mapped to raw annotation file either to query_id column (gene level) or input_fasta (fasta file name, stand in for genome)"
    )


def join_mapping(
    annotations: pl.DataFrame,
    path,
    gene_column: str = "query_id",
    genome_column: str = DEFAULT_FASTA_COLUMN,
) -> pl.DataFrame:
    """The annotations unchanged, with a mapping file's Samples joined on as
    Sample Columns (`sample:<name>`, in file order) after every other column.

    The mapping joins on `gene_column` or `genome_column`, as `load_mapping`
    decides; rows the mapping lacks get nulls. Raises ValueError as
    `mapping_key` does, or when a Sample Column is already in the annotations.
    """
    wide, first, samples = read_mapping(path)
    wide = wide.collect()
    genes = annotations.select(
        query_id=pl.col(gene_column).cast(pl.String),
        genome=pl.col(genome_column).cast(pl.String),
    )
    key = {"query_id": gene_column, "genome": genome_column}[
        mapping_key(wide[first], genes)
    ]
    columns = {s: SAMPLE_COLUMN_PREFIX + s for s in samples}
    if clash := [c for c in columns.values() if c in annotations.columns]:
        raise ValueError(f"The annotations already have Sample Columns {clash}")
    wide = wide.rename({first: key, **columns}).with_columns(
        pl.col(key).cast(annotations.schema[key])
    )
    return annotations.join(wide, on=key, how="left", maintain_order="left")


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
    wide, first, samples = read_mapping(path)
    means, abundances = means_and_abundances(wide, first, samples)
    key = mapping_key(means[first], genes)
    return means.rename({first: key}), abundances.rename({first: key}), samples


def means_and_abundances(
    wide: pl.LazyFrame, key: str, samples: list[str]
) -> tuple[pl.DataFrame, pl.DataFrame]:
    """`load_mapping`'s means and abundances from a wide table of `key` and
    one column per Sample."""
    return tuple(
        pl.collect_all(
            [
                wide.select(
                    key,
                    mean_sample_abundance=pl.mean_horizontal(
                        pl.col(samples).fill_null(0)
                    ),
                ),
                wide.unpivot(
                    index=key,
                    on=samples,
                    variable_name="sample",
                    value_name="sample_abundance",
                ).filter(pl.col("sample_abundance").fill_null(0) != 0),
            ],
            engine="streaming",
        )
    )


def sample_columns(names: list[str]) -> list[str]:
    """The Sample Columns among annotation column `names`."""
    return [c for c in names if c.startswith(SAMPLE_COLUMN_PREFIX)]


def load_sample_columns(
    genes: pl.LazyFrame, columns: list[str]
) -> tuple[pl.DataFrame, pl.DataFrame, list[str]]:
    """`load_mapping`'s outputs from annotations carrying Sample `columns`
    (genome column already named `genome`).

    Genome-level when every gene of each genome has the same abundances in
    every Sample, so a genome's tRNA rows still belong to its Samples;
    otherwise gene-level. Rows with no abundance in any Sample are left out,
    as genes or genomes missing from a mapping file are.
    """
    samples = [c.removeprefix(SAMPLE_COLUMN_PREFIX) for c in columns]
    wide = genes.select(
        "query_id",
        "genome",
        *(pl.col(c).cast(pl.Float64).alias(s) for c, s in zip(columns, samples)),
    ).filter(pl.any_horizontal(pl.col(samples).is_not_null()))
    per_genome = (
        wide.group_by("genome")
        .agg(pl.col(samples).n_unique().max())
        .select(pl.max_horizontal(samples).max())
        .collect()
        .item()
    )
    if per_genome in (None, 1):
        wide, key = wide.select("genome", *samples).unique("genome"), "genome"
    else:
        wide, key = wide.select("query_id", *samples), "query_id"
    means, abundances = means_and_abundances(wide, key, samples)
    return means, abundances, samples


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
