"""Hits: loading the annotations once and evaluating sheets at any level."""

from pathlib import Path

import polars as pl
import pytest

from dram_rules.hits import Hits
from dram_rules.rule_parser.rules import CompiledRules

ANNOTATIONS = """\
query_id\tinput_fasta\tkegg_id\tCompleteness
bin1_1\tbin1\tK00399\t90
bin1_2\tbin1\t\t90
bin2_1\tbin2\tK00399\t80
"""
MAPPING = "Geneid\tS1\tS2\nbin1_1\t1\t0\nbin1_2\t3\t4\nbin2_1\t5\t6\n"
# bin3 has tRNAs but no annotations.
TRNA = """\
input_fasta\tquery_id\tgene_id
bin1\tbin1_contig1\tAla (AGC)
bin3\tbin3_contig1\tVal (GAC)
"""


def sheet(*rules: str) -> CompiledRules:
    return CompiledRules.from_rules(
        rules=pl.DataFrame(
            {"name": [f"r{i}" for i in range(len(rules))], "rule": rules, "group": "G"}
        ).lazy()
    )


def write(tmp_path: Path, name: str, text: str) -> Path:
    (tmp_path / name).write_text(text)
    return tmp_path / name


def load(tmp_path: Path, *sheets: CompiledRules, annotations=ANNOTATIONS, **files):
    paths = {k: write(tmp_path, f"{k}.tsv", v) for k, v in files.items()}
    return Hits.load(write(tmp_path, "anno.tsv", annotations), sheets, **paths)


def present(dfs: dict[str, pl.DataFrame], unit: str, rule: str = "r0") -> dict:
    """{genome or sample: Presence} of `rule`."""
    df = pl.concat(dfs.values()).filter(pl.col("name") == rule)
    return dict(zip(df[unit], df["present"]))


def in_memory(frame: pl.LazyFrame) -> bool:
    return "SCAN" not in frame.explain().upper()


def test_rows_are_read_once_into_memory_when_a_rule_reads_them(tmp_path):
    hits = load(tmp_path, sheet("count_rows(>=, 1)"), sheet("K00399"), mapping=MAPPING)
    assert in_memory(hits.rows)
    assert in_memory(hits.sample_rows)


def test_rows_stay_on_disk_when_no_rule_reads_them(tmp_path):
    hits = load(tmp_path, sheet("K00399"), mapping=MAPPING)
    assert not in_memory(hits.rows)
    assert not in_memory(hits.sample_rows)


def test_the_sample_level_needs_a_mapping_file(tmp_path):
    compiled = sheet("K00399")
    with pytest.raises(ValueError, match="mapping file"):
        load(tmp_path, compiled).evaluate(compiled, "sample")


def test_sample_level_with_a_mapping_file(tmp_path):
    compiled = sheet("K00399")
    hits = load(tmp_path, compiled, mapping=MAPPING)
    assert hits.samples == ["S1", "S2"]
    # bin1_1 (S1 only) and bin2_1 hold K00399.
    assert present(hits.evaluate(compiled, "sample"), "sample") == {
        "S1": True,
        "S2": True,
    }


def test_trna_rows_reach_every_level_but_not_the_gene_rows(tmp_path):
    compiled = sheet("`Ala (AGC)`", "`Val (GAC)`")
    hits = load(tmp_path, compiled, trna=TRNA)
    assert hits.genomes == ["bin1", "bin2", "bin3"]
    genome = hits.evaluate(compiled, "genome")
    assert present(genome, "genome") == {"bin1": True, "bin2": False, "bin3": False}
    assert present(genome, "genome", "r1") == {
        "bin1": False,
        "bin2": False,
        "bin3": True,
    }
    gene_counts = pl.concat(hits.evaluate(compiled, "gene").values())
    assert gene_counts["gene_count"].sum() == 2
    # Meta and taxonomy read the gene rows, which have no tRNA rows.
    genes = hits.genes.collect()
    assert genes["genome"].unique().sort().to_list() == ["bin1", "bin2"]
    assert genes["Completeness"].null_count() == 0


def test_genome_column_falls_back_to_fasta(tmp_path):
    compiled = sheet("K00399")
    annotations = ANNOTATIONS.replace("input_fasta", "fasta", 1)
    hits = load(tmp_path, compiled, annotations=annotations)
    assert hits.genomes == ["bin1", "bin2"]


def test_annotations_without_a_genome_column_are_an_error(tmp_path):
    annotations = ANNOTATIONS.replace("input_fasta", "bin", 1)
    with pytest.raises(ValueError, match="genome columns"):
        load(tmp_path, sheet("K00399"), annotations=annotations)
