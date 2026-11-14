"""End-to-end runs of `dram-rules --sum` (the summarize step)."""

from pathlib import Path

import polars as pl
from click.testing import CliRunner

from dram_rules.cli import main

# bin1 has two mcrA genes and a cyc2 gene; bin2 has one mcrA gene.
ANNOTATIONS = """\
query_id\tinput_fasta\tkegg_id\tfegenie_id
bin1_1\tbin1\tK00399\t
bin1_2\tbin1\tK00399\tCyc2_repCluster1
bin1_3\tbin1\t\tCyc2_repCluster1
bin2_1\tbin2\tK00399\t
"""

SHEET = """\
gene_id\tgene_description\tmodule\ttopic_ecosystem\tnotes
K00399\tmcrA\tMethanogenesis\tEnergy\tdropped
Cyc2_repCluster1\tcyc2\tIron\tEnergy\tdropped
K99999\tnothing\tNone\tOther\tdropped
"""


def invoke(tmp_path: Path, *args: str, annotations: str = ANNOTATIONS):
    (tmp_path / "anno.tsv").write_text(annotations)
    (tmp_path / "common.tsv").write_text("name\talias\trule\n")
    cli_args = [
        "-a", str(tmp_path / "anno.tsv"),
        "-c", str(tmp_path / "common.tsv"),
        "-o", str(tmp_path / "out"),
    ]  # fmt: skip
    return CliRunner().invoke(main, [*cli_args, *args])


def run(tmp_path: Path, *args: str, **kwargs) -> Path:
    result = invoke(tmp_path, *args, **kwargs)
    assert result.exit_code == 0, result.output + repr(result.exception)
    return tmp_path / "out" / "summarize"


def write(tmp_path: Path, name: str, text: str) -> str:
    (tmp_path / name).write_text(text)
    return str(tmp_path / name)


def counts(df: pl.DataFrame, label_col: str = "gene_id") -> dict:
    """{label: {genome: Gene Count}}."""
    return {
        row[label_col]: {g: row[g] for g in ("bin1", "bin2")}
        for row in df.iter_rows(named=True)
    }


def test_custom_sheet_gene_counts(tmp_path):
    sheet = write(tmp_path, "a.tsv", SHEET)
    out = run(tmp_path, "--sum", "--sum_custom", sheet)
    df = pl.read_csv(out / "custom_summary.tsv", separator="\t")
    assert counts(df) == {
        "K00399": {"bin1": 2, "bin2": 1},
        "Cyc2_repCluster1": {"bin1": 2, "bin2": 0},
        "K99999": {"bin1": 0, "bin2": 0},
    }


def test_workbook_has_a_sorted_worksheet_per_group(tmp_path):
    sheet = write(tmp_path, "a.tsv", SHEET)
    out = run(tmp_path, "--sum", "--sum_custom", sheet)
    book = pl.read_excel(out / "custom_summary.xlsx", sheet_id=0)
    assert list(book) == ["Energy", "Other"]
    energy = book["Energy"]
    # Only the summary columns: no `notes`, and the Group is the worksheet.
    assert energy.columns == ["gene_id", "gene_description", "module", "bin1", "bin2"]
    # Sorted by module (distill's order: header, subheader, module, label).
    assert energy["gene_id"].to_list() == ["Cyc2_repCluster1", "K00399"]
    assert energy.row(1) == ("K00399", "mcrA", "Methanogenesis", 2, 1)


def test_sheets_with_different_label_columns_use_name(tmp_path):
    a = write(tmp_path, "a.tsv", SHEET)
    b = write(tmp_path, "b.tsv", "name\ttopic_ecosystem\trule\nMcr\tEnergy\tK00399\n")
    out = run(tmp_path, "--sum", "--sum_custom", f"{a},{b}")
    df = pl.read_csv(out / "custom_summary.tsv", separator="\t")
    assert counts(df, "name")["Mcr"] == {"bin1": 2, "bin2": 1}
    assert counts(df, "name")["K00399"] == {"bin1": 2, "bin2": 1}
    book = pl.read_excel(out / "custom_summary.xlsx", sheet_id=0)
    assert book["Energy"].columns[0] == "name"


# tRNAscan rows as DRAM collects them: query_id is the contig, so two tRNAs
# on one contig share it.
TRNA = """\
input_fasta\tquery_id\ttRNA #\tbegin\tend\ttype\tcodon\tscore\tnote\tgene_id
bin1\tbin1_contig1 \t1\t10\t80\tAla\tAGC\t70.0\t\tAla (AGC)
bin1\tbin1_contig1 \t2\t90\t160\tAla\tAGC\t70.0\t\tAla (AGC)
bin2\tbin2_contig1\t1\t10\t80\tVal\tGAC\t70.0\t\tVal (GAC)
"""

TRNA_SHEET = """\
gene_id\tgene_description\ttopic_ecosystem
`Ala (AGC)`\tAlanine (AGC)\ttRNAs
`Val (GAC)`\tValine (GAC)\ttRNAs
K00399\tmcrA\tEnergy
"""


def test_trna_rows_get_gene_counts(tmp_path):
    sheet = write(tmp_path, "a.tsv", TRNA_SHEET)
    trna = write(tmp_path, "trna.tsv", TRNA)
    out = run(tmp_path, "--sum", "--sum_custom", sheet, "--trna_path", trna)
    df = pl.read_csv(out / "custom_summary.tsv", separator="\t")
    assert counts(df) == {
        "`Ala (AGC)`": {"bin1": 2, "bin2": 0},
        "`Val (GAC)`": {"bin1": 0, "bin2": 1},
        "K00399": {"bin1": 2, "bin2": 1},
    }


def test_coverage_options_are_an_error_naming_the_rule(tmp_path):
    sheet = write(
        tmp_path,
        "a.tsv",
        "name\ttopic_ecosystem\trule\trule_options\n"
        "Mcr\tEnergy\tK00399\t\n"
        'Methanogenesis\tCH4\tK00399, Cyc2_repCluster1\t{"coverage": "steps"}\n',
    )
    result = invoke(tmp_path, "--sum", "--sum_custom", sheet)
    assert result.exit_code != 0
    assert "Methanogenesis" in str(result.exception)
    assert "coverage" in str(result.exception)
    assert not (tmp_path / "out" / "summarize").exists()


def test_null_trna_file_adds_nothing(tmp_path):
    # DRAM writes a file holding only "NULL" when annotate found no tRNAs.
    sheet = write(tmp_path, "a.tsv", SHEET)
    trna = write(tmp_path, "trna.tsv", "NULL\n")
    out = run(tmp_path, "--sum", "--sum_custom", sheet, "--trna_path", trna)
    df = pl.read_csv(out / "custom_summary.tsv", separator="\t")
    assert counts(df)["K00399"] == {"bin1": 2, "bin2": 1}


NEEDS_TAXONOMY = '{"include_when": {"all_columns_present": ["taxonomy"]}}'


def test_rules_left_out_by_include_when_are_absent(tmp_path):
    # The annotations have no taxonomy, so Group "Taxa" loses its only rule
    # and sheet b loses every rule.
    a = write(
        tmp_path,
        "a.tsv",
        "gene_id\ttopic_ecosystem\trule_options\n"
        "K00399\tEnergy\t\n"
        f"Cyc2_repCluster1\tTaxa\t{NEEDS_TAXONOMY}\n",
    )
    b = write(
        tmp_path,
        "b.tsv",
        f"name\ttopic_ecosystem\trule_options\nK00399\tEnergy\t{NEEDS_TAXONOMY}\n",
    )
    out = run(tmp_path, "--sum", "--sum_custom", f"{a},{b}")
    df = pl.read_csv(out / "custom_summary.tsv", separator="\t")
    # Sheet b added nothing, so the label column stays a's `gene_id`.
    assert counts(df) == {"K00399": {"bin1": 2, "bin2": 1}}
    assert list(pl.read_excel(out / "custom_summary.xlsx", sheet_id=0)) == ["Energy"]


def test_no_collection_when_every_rule_is_left_out(tmp_path):
    b = write(
        tmp_path,
        "b.tsv",
        f"gene_id\ttopic_ecosystem\trule_options\nK00399\tEnergy\t{NEEDS_TAXONOMY}\n",
    )
    out = run(tmp_path, "--sum", "--sum_custom", b)
    assert not (out / "custom_summary.tsv").exists()
    assert not (out / "custom_summary.xlsx").exists()
