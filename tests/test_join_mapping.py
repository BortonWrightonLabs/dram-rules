"""`join_mapping`: a mapping file attached to annotations as Sample Columns."""

import polars as pl
import pytest

from dram_rules.processing.process_annotations import join_mapping

ANNOTATIONS = pl.DataFrame(
    {
        "query_id": ["bin1_1", "bin1_2", "bin2_1"],
        "input_fasta": ["bin1", "bin1", "bin2"],
        "kegg_id": ["K00399", None, "K00399"],
        "rank": ["B", "E", "B"],
    }
)


def write(tmp_path, text):
    path = tmp_path / "mapping.tsv"
    path.write_text(text)
    return path


def test_gene_mapping_adds_sample_columns_on_the_right(tmp_path):
    # featureCounts layout: meta columns and a text column are not Samples;
    # bin1_2 has no row, so its Sample Columns are empty.
    mapping = write(
        tmp_path,
        "Geneid\tChr\tLength\tNotes\tS1\tS2\n"
        "bin1_1\tc1\t900\tx\t3\t0\n"
        "bin2_1\tc2\t300\ty\t1.5\t7\n",
    )
    joined = join_mapping(ANNOTATIONS, mapping)
    expected = ANNOTATIONS.with_columns(
        pl.Series("sample:S1", [3.0, None, 1.5]),
        pl.Series("sample:S2", [0.0, None, 7.0]),
    )
    assert joined.equals(expected)


def test_genome_mapping_gives_each_gene_its_genome_abundance(tmp_path):
    # bin2 has no row: its gene's Sample Column is empty.
    mapping = write(tmp_path, "bin\tS1\nbin1\t10\n")
    joined = join_mapping(ANNOTATIONS, mapping)
    expected = ANNOTATIONS.with_columns(pl.Series("sample:S1", [10.0, 10.0, None]))
    assert joined.equals(expected)


def test_names_matching_neither_genes_nor_genomes_fail(tmp_path):
    mapping = write(tmp_path, "Geneid\tS1\nbin1_1\t1\nbin9_9\t2\n")
    with pytest.raises(ValueError, match="can't be mapped"):
        join_mapping(ANNOTATIONS, mapping)


def test_a_sample_column_already_in_the_annotations_fails(tmp_path):
    annotations = ANNOTATIONS.with_columns(pl.lit(1.0).alias("sample:S1"))
    mapping = write(tmp_path, "Geneid\tS1\nbin1_1\t1\n")
    with pytest.raises(
        ValueError, match=r"already have Sample Columns \['sample:S1'\]"
    ):
        join_mapping(annotations, mapping)
