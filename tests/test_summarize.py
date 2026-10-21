import polars as pl
import pytest
from click.testing import CliRunner

from dram_rules.make_summarize import (
    convert_rules,
    load_definitions,
    main,
    summarize_annotations,
)


@pytest.fixture
def definitions_dir(tmp_path):
    directory = tmp_path / "definitions"
    directory.mkdir()
    (directory / "summarize_genes.tsv").write_text(
        "gene_id\tgene_description\nA\tAlpha [EC 1.2.3.4]\nSELF\tself\n"
    )
    (directory / "fam-substrate-mapping.tsv").write_text(
        "# disclosure\nFamily\tEC_Number\tName\n"
        "GH9\t 3.2.1.4 \tCellulase\nCBM1\t\tBinding module\n"
    )
    (directory / "dbcan.fam-activities.tsv").write_text(
        "# disclosure\nGH9\tHydrolase\nCBM1\tBinding\textra ignored field\n"
    )
    return directory


@pytest.fixture
def rule_files(tmp_path):
    rules = tmp_path / "rules.tsv"
    rules.write_text(
        "name\trule\tmodule\tgroup\ttopic_ecosystem\tnotes\n"
        "First\t@Shared | GH9\tM1\tG\tag\tone\n"
        "Second\tA | B\tM2\tG\tbgc\ttwo\n"
        "First\tC\tM3\tH\tgut\tthree\n"
    )
    common = tmp_path / "common.tsv"
    common.write_text("name\talias\trule\n\tShared\tA | B\nUnrelated\t\tUNRELATED\n")
    return rules, common


def test_definition_sources_and_override(definitions_dir, tmp_path, caplog):
    custom = tmp_path / "custom.tsv"
    custom.write_text("gene_id\tgene_description\na\tCustom Alpha\nB\tBeta\n")
    definitions = load_definitions(definitions_dir, custom)
    assert definitions == {
        "A": "Custom Alpha",
        "B": "Beta",
        "SELF": "",
        "GH9": "Hydrolase",
        "GH9;EC:3.2.1.4": "Cellulase",
        "CBM1": "Binding",
    }
    assert "Skipped 1 rows" in caplog.text
    assert load_definitions(definitions_dir)["A"] == "Alpha [EC:1.2.3.4]"


def test_builtin_duplicates_fail_before_override(definitions_dir, tmp_path):
    with (definitions_dir / "summarize_genes.tsv").open("a") as handle:
        handle.write("a\tAlpha\nGH9\tHydrolase\n")
    custom = tmp_path / "custom.tsv"
    custom.write_text("gene_id\tgene_description\na\tReplacement\n")
    with pytest.raises(ValueError) as error:
        load_definitions(definitions_dir, custom)
    assert "\nA\t" in str(error.value)
    assert "\nGH9\t" in str(error.value)
    assert "dbcan.fam-activities.tsv" in str(error.value)


def test_custom_duplicates_fail(definitions_dir, tmp_path):
    custom = tmp_path / "custom.tsv"
    custom.write_text("gene_id\tgene_description\na\tSame\nA\tSame\n")
    with pytest.raises(ValueError, match="Duplicate definition"):
        load_definitions(definitions_dir, custom)


def test_rule_conversion_metadata_and_aliases(rule_files):
    rules, common = rule_files
    sheet = convert_rules(rules, {"A": "Alpha"}, common, all_metadata=True)
    assert sheet.columns == [
        "gene_id",
        "gene_description",
        "module",
        "group",
        "topic_ecosystem",
        "notes",
        "name",
    ]
    assert sheet["gene_id"].to_list() == ["A", "B", "C", "GH9"]
    alpha = sheet.filter(pl.col("gene_id") == "A").row(0, named=True)
    assert alpha == {
        "gene_id": "A",
        "gene_description": "Alpha",
        "module": "M1,M2",
        "group": "G",
        "topic_ecosystem": "ag,bgc",
        "notes": "one,two",
        "name": "First,Second",
    }
    assert sheet.filter(pl.col("gene_id") == "C")["name"].to_list() == ["First"]
    assert convert_rules(rules, {}, common).columns == [
        "gene_id",
        "gene_description",
        "topic_ecosystem",
        "name",
    ]


def test_strict_conversion_reports_all_missing_including_cazy(rule_files):
    rules, common = rule_files
    with pytest.raises(ValueError) as error:
        convert_rules(rules, {"A": "Alpha", "B": ""}, common, strict=True)
    assert str(error.value).splitlines() == [
        "Missing gene descriptions:",
        "B",
        "C",
        "GH9",
    ]


@pytest.mark.parametrize("genome_column", ["genome", "input_fasta"])
def test_counts_distinct_genes_across_databases_and_zero_genomes(
    rule_files, genome_column
):
    rules, common = rule_files
    sheet = convert_rules(rules, {}, common)
    raw = pl.DataFrame(
        {
            "query_id": ["q1", "q1", "q2", "q3", "q1", "q4"],
            genome_column: ["g1", "g1", "g1", "g1", "g2", "g3"],
            "kegg_id": ["a", "a", "A/B", "-", "B", None],
            "kofam_id": ["A", "A", "A", "C", "B", None],
            "dbcan_id": [None, None, None, "GH9", None, None],
        }
    )
    result = summarize_annotations(sheet, raw)
    assert result["gene_id"].to_list() == ["A", "B", "C", "GH9"]
    assert result["g1"].to_list() == [2, 1, 1, 1]
    assert result["g2"].to_list() == [0, 1, 0, 0]
    assert result["g3"].to_list() == [0, 0, 0, 0]
    assert result.schema["g1"] == pl.Int64


def test_cazy_combination_count():
    sheet = pl.DataFrame(
        {"gene_id": ["GH9;EC:3.2.1.4"], "gene_description": ["Cellulase"]}
    )
    raw = pl.DataFrame(
        {
            "query_id": ["q1", "q2"],
            "genome": ["g1", "g1"],
            "dbcan_id": ["GH9", "GH9"],
            "dbcan_sub_ec": ["3.2.1.4:10", "3.2.1.4:20"],
        }
    )
    assert summarize_annotations(sheet, raw)["g1"].to_list() == [2]


def test_cli_output_and_strict_failure(definitions_dir, tmp_path):
    rules = tmp_path / "example.tsv"
    rules.write_text("name\trule\ttopic_ecosystem\nFirst\tA | B\tag\n")
    common = tmp_path / "common.tsv"
    common.write_text("name\talias\trule\n")
    annotations = tmp_path / "annotations.tsv"
    annotations.write_text("query_id\tinput_fasta\tkegg_id\nq1\tg1\tA\n")
    output = tmp_path / "output"
    args = [
        "-a",
        str(annotations),
        "-r",
        str(rules),
        "-c",
        str(common),
        "--definitions-dir",
        str(definitions_dir),
        "-o",
        str(output),
    ]
    runner = CliRunner()
    strict = runner.invoke(main, [*args, "--strict"])
    assert strict.exit_code != 0
    assert "Missing gene descriptions:\nB" in strict.output
    assert not output.exists()
    success = runner.invoke(main, [*args, "--save-converted-sheet", "--save-lookup"])
    assert success.exit_code == 0, success.output
    assert sorted(path.name for path in output.iterdir()) == [
        "example_summarize_rules.tsv",
        "example_summary.tsv",
        "master_definitions.tsv",
    ]
    result = pl.read_csv(output / "example_summary.tsv", separator="\t")
    assert result["g1"].to_list() == [1, 0]
