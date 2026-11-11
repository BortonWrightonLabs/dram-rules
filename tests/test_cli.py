"""End-to-end runs of the `dram_viz` command (make_product.main)."""

from pathlib import Path

import polars as pl
from click.testing import CliRunner

from dram_rules.make_product import main

ANNOTATIONS = """\
query_id\tinput_fasta\tkegg_id\tfegenie_id\theme_regulatory_motif_count
bin1_1\tbin1\tK00399\t\t5
bin1_2\tbin1\t\tCyc2_repCluster1\t0
bin2_1\tbin2\tK00399\t\t1
"""

RULES = """\
name\tgroup\trule
Mcr\tG1\tK00399
Cyc2\tG1\tCyc2_repCluster1
"""

MAPPING = """\
Geneid\tS1\tS2
bin1_1\t1.0\t2.0
bin1_2\t3.0\t4.0
bin2_1\t5.0\t6.0
"""


def run(
    tmp_path: Path,
    *args: str,
    annotations: str = ANNOTATIONS,
    rules: str = RULES,
    mapping: str | None = MAPPING,
):
    """Run the CLI on small inputs written to `tmp_path`; return the output dir."""
    (tmp_path / "anno.tsv").write_text(annotations)
    (tmp_path / "rules.tsv").write_text(rules)
    (tmp_path / "common.tsv").write_text("name\talias\trule\n")
    out = tmp_path / "out"
    cli_args = [
        "-a", str(tmp_path / "anno.tsv"),
        "-r", str(tmp_path / "rules.tsv"),
        "-c", str(tmp_path / "common.tsv"),
        "-o", str(out),
        "--save_dataframes",
    ]  # fmt: skip
    if mapping is not None:
        (tmp_path / "mapping.tsv").write_text(mapping)
        cli_args += ["-m", str(tmp_path / "mapping.tsv")]
    result = CliRunner().invoke(main, [*cli_args, *args])
    assert result.exit_code == 0, result.output + repr(result.exception)
    return out


def read(out: Path, kind: str) -> pl.DataFrame:
    return pl.read_csv(out / f"df_vertical_{kind}.tsv", separator="\t")


def value(df: pl.DataFrame, col: str, **where) -> float:
    """The single value of `col` in the row matching `where`."""
    row = df.filter(*(pl.col(k) == v for k, v in where.items()))
    return row.get_column(col).item()


def test_mixed_case_gene_id_gets_its_abundance(tmp_path):
    out = run(tmp_path)
    # Cyc2_repCluster1 is only on bin1_2: S1 = 3.0, S2 = 4.0.
    sample = read(out, "sample")
    assert value(sample, "sample_abundance", name="Cyc2", sample="S1") == 3.0
    assert value(sample, "sample_abundance", name="Cyc2", sample="S2") == 4.0
    genome = read(out, "genome")
    assert value(genome, "mean_sample_abundance", name="Cyc2", genome="bin1") == 3.5


def test_mapping_decimal_after_many_integer_rows_is_kept(tmp_path):
    # Polars guesses column types from the first rows; 150 integer rows
    # before the decimal one must not turn 3.5 into 0.
    padding = [f"pad_{i}" for i in range(150)]
    annotations = ANNOTATIONS + "".join(f"{g}\tbin3\t\t\t0\n" for g in padding)
    mapping = (
        "Geneid\tS1\tS2\n"
        + "".join(f"{g}\t0\t0\n" for g in padding)
        + "bin1_1\t1\t2\nbin1_2\t3.5\t4\nbin2_1\t5\t6\n"
    )
    out = run(tmp_path, annotations=annotations, mapping=mapping)
    sample = read(out, "sample")
    assert value(sample, "sample_abundance", name="Cyc2", sample="S1") == 3.5
