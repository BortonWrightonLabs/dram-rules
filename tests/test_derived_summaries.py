"""Derived Summarize Sheets (ADR 0004): ecosystem summaries built from rule
sheets and the Description Tables."""

from pathlib import Path

import polars as pl
import pytest
from click.testing import CliRunner

from dram_rules import cli, make_summarize
from dram_rules.rule_parser.rules import RuleError

# bin1: two mcrA genes (one also cyc2), a GH18 chitinase, a cyc2-only gene.
ANNOTATIONS = """\
query_id\tinput_fasta\tkegg_id\tfegenie_id\tdbcan_id\tdbcan_sub_ec
bin1_1\tbin1\tK00399\t\t\t
bin1_2\tbin1\tK00399\tCyc2_repCluster1\t\t
bin1_3\tbin1\t\tCyc2_repCluster1\t\t
bin1_4\tbin1\t\t\tGH18\t3.2.1.14:2
bin2_1\tbin2\tK00399\t\t\t
"""

COMMON = "name\talias\trule\n\tiron\tCyc2_repCluster1\n"

# Rules in sheet order; K00399 is in two rules, modules and Groups.
SHEET = """\
name\trule\tmodule\tgroup
Methanogenesis\tK00399 & not(K11111)\tEnergy\tMethane
Iron, oxidation\t@iron | k00399\tMetals\tIron
Chitin\tGH18;EC:3.2.1.14\tCarbon\tPolymer
Chitin, again\tGH18;EC:3.2.1.14 | K22222\tEnergy\tPolymer
"""

DESCRIPTIONS = (
    "gene_id\tdescription\n"
    "K00399\tmcrA\n"
    "k11111\tnegated\n"
    "Cyc2_repCluster1\tcyc2\n"
    "K22222\tsomething\n"
)


@pytest.fixture
def tables(tmp_path, monkeypatch):
    """The built-in Description Tables, replaced with small fixture tables."""
    dbcan = tmp_path / "dbcan.tsv"
    dbcan.write_text(
        "# Copyright disclosure\nGH18\t  chitinase\t(EC 3.2.1.14) ;\nCE20\t  \n"
    )
    genes = tmp_path / "gene_descriptions.tsv"
    genes.write_text(DESCRIPTIONS)
    monkeypatch.setattr(make_summarize, "DBCAN_FAMILIES", dbcan)
    monkeypatch.setattr(make_summarize, "GENE_DESCRIPTIONS", genes)
    topics = tmp_path / "topics"
    topics.mkdir()
    (topics / "t1.tsv").write_text(TOPIC_1)
    (topics / "t2.tsv").write_text(TOPIC_2)
    monkeypatch.setattr(make_summarize, "SUMMARIZE_SHEETS_DIR", topics)
    monkeypatch.setattr(make_summarize, "TOPICS", {"t1": "t1.tsv", "t2": "t2.tsv"})


# Topic Summarize Sheets, also read as Description Tables.
TOPIC_1 = """\
\ufefftopic_ecosystem\tgene_id\tgene_description
Energy\tK1|K2\tshared
Energy\tK3\tK3
Energy\tK00399\tloses to gene_descriptions.tsv
Energy\t`Ala (AGC)`\tAlanine tRNA
"""
TOPIC_2 = "topic_ecosystem\tgene_id\tgene_description\nEnergy\tK1\tlater sheet\n"


def write(tmp_path: Path, name: str, text: str) -> Path:
    (tmp_path / name).write_text(text)
    return tmp_path / name


def derive(tmp_path, sheet: str = SHEET, custom=()):
    compiled = make_summarize.load_rule_sheet(
        write(tmp_path, "s.tsv", sheet), write(tmp_path, "common.tsv", COMMON), ()
    )
    return make_summarize.derive_sheet(
        compiled, make_summarize.load_descriptions(list(custom))
    )


def test_one_row_per_gene_id_listing_its_rules_modules_and_groups(tmp_path, tables):
    df, missing = derive(tmp_path)
    assert missing == []
    assert df.columns == make_summarize.DERIVED_COLUMNS
    rows = {r[0]: r[1:] for r in df.iter_rows()}
    # Matched regardless of case; listed as first written, with every use.
    assert rows["K00399"] == (
        "mcrA",
        "Methanogenesis | Iron, oxidation",
        "Energy | Metals",
        "Methane | Iron",
    )
    # Negated Gene IDs and Gene IDs from aliases are rows too.
    assert rows["K11111"][:2] == ("negated", "Methanogenesis")
    assert rows["Cyc2_repCluster1"][:2] == ("cyc2", "Iron, oxidation")


def test_family_ec_gene_ids_are_described_from_the_sheet(tmp_path, tables):
    df, _ = derive(tmp_path)
    [description] = df.filter(gene_id="GH18;EC:3.2.1.14")["gene_description"]
    assert description == "Carbon; Polymer; Chitin | Energy; Polymer; Chitin, again"


def test_rows_are_ordered_by_module_then_group_then_first_appearance(tmp_path, tables):
    df, _ = derive(tmp_path)
    # Energy (Methane, then Polymer), then Metals, then Carbon; a Gene ID in
    # several sits under its first module and Group.
    assert df["gene_id"].to_list() == [
        "K00399",
        "K11111",
        "K22222",
        "Cyc2_repCluster1",
        "GH18;EC:3.2.1.14",
    ]


def test_missing_descriptions_and_families_are_reported(tmp_path, tables):
    sheet = "name\trule\tmodule\tgroup\nA\tK99999 | CE20;EC:3.1.1.1 | GH99;EC:1\tM\tG\n"
    df, missing = derive(tmp_path, sheet)
    # CE20 is in the dbCAN table with no description, so it counts as missing.
    assert missing == ["K99999", "family CE20", "family GH99"]
    assert df.filter(gene_id="K99999")["gene_description"].to_list() == [None]


def test_topic_summarize_sheets_describe_gene_ids(tmp_path, tables):
    descriptions = make_summarize.load_descriptions([])
    # Each of a row's `|`-joined Gene IDs; the first sheet wins a conflict.
    assert descriptions["K1"] == descriptions["K2"] == "shared"
    # A description that only repeats the Gene ID is none.
    assert "K3" not in descriptions
    assert descriptions["ALA (AGC)"] == "Alanine tRNA"
    assert descriptions["K00399"] == "mcrA"


def test_custom_tables_win_over_built_in_ones(tmp_path, tables):
    custom = write(
        tmp_path, "mine.tsv", "﻿gene_id\tdescription\tnotes\nK00399\tmine\tx\n"
    )
    descriptions = make_summarize.load_descriptions([custom])
    assert descriptions["K00399"] == "mine"
    # dbCAN descriptions lose their padding and inner tabs.
    assert descriptions["GH18"] == "chitinase (EC 3.2.1.14) ;"
    assert "CE20" not in descriptions


def test_a_custom_table_without_the_columns_is_an_error(tmp_path, tables):
    bad = write(tmp_path, "bad.tsv", "id\tdesc\nK00399\tx\n")
    with pytest.raises(RuleError, match="bad.tsv"):
        make_summarize.load_descriptions([bad])


def test_strict_lists_every_missing_id_by_sheet(tmp_path, tables):
    common = write(tmp_path, "common.tsv", COMMON)
    sheets = {
        name: make_summarize.load_rule_sheet(
            write(tmp_path, f"{name}.tsv", f"name\trule\tgroup\nR\t{ids}\tG\n"),
            common,
            (),
        )
        for name, ids in [("a", "K1 | K00399"), ("b", "K2"), ("c", "K00399")]
    }
    with pytest.raises(RuleError) as error:
        make_summarize.derive_collection(sheets, {"K00399": "mcrA"}, strict=True)
    assert "  a: K1\n  b: K2" in str(error.value)
    assert "c:" not in str(error.value)


def test_not_strict_warns_once(tmp_path, tables, caplog):
    compiled = make_summarize.load_rule_sheet(
        write(tmp_path, "a.tsv", "name\trule\tgroup\nR\tK1 | K2\tG\n"),
        write(tmp_path, "common.tsv", COMMON),
        (),
    )
    with caplog.at_level("WARNING", logger="dram.rules"):
        [df] = make_summarize.derive_collection({"a": compiled}, {}, strict=False)
    assert len(caplog.records) == 1
    assert "2 Gene IDs" in caplog.text
    assert df["sheet"].to_list() == ["a", "a"]


@pytest.fixture
def ecosystems(tmp_path, monkeypatch, tables):
    """Every built-in ecosystem sheet replaced with SHEET."""
    sheet = write(tmp_path, "eco.tsv", SHEET)
    monkeypatch.setattr(cli, "ECOSYSTEMS", dict.fromkeys(cli.ECOSYSTEMS, sheet))


def run(tmp_path, *args, ok=True):
    anno = write(tmp_path, "anno.tsv", ANNOTATIONS)
    common = write(tmp_path, "common.tsv", COMMON)
    out = tmp_path / "out"
    result = CliRunner().invoke(
        cli.main, ["-a", str(anno), "-o", str(out), "-c", str(common), *args]
    )
    if ok:
        assert result.exit_code == 0, result.output + repr(result.exception)
    return result, out


def test_ecosystem_summary_has_a_worksheet_per_sheet(tmp_path, ecosystems):
    _, out = run(tmp_path, "--sum", "--topics", "none", "--ecos", "gut,bgc")
    df = pl.read_csv(out / "summarize/ecosystem_summary.tsv", separator="\t")
    assert df.columns == [
        *make_summarize.DERIVED_COLUMNS,
        "sheet",
        "bin1",
        "bin2",
    ]
    # The derived columns follow every Collection's order.
    derived = make_summarize.DERIVED_COLUMNS[1:]
    assert derived == [c for c in make_summarize.COLUMN_ORDER if c in derived]
    gut = df.filter(sheet="gut")
    counts = {r["gene_id"]: (r["bin1"], r["bin2"]) for r in gut.iter_rows(named=True)}
    assert counts == {
        "K00399": (2, 1),
        "K11111": (0, 0),
        "K22222": (0, 0),
        "Cyc2_repCluster1": (2, 0),
        "GH18;EC:3.2.1.14": (1, 0),
    }
    book = pl.read_excel(out / "summarize/ecosystem_summary.xlsx", sheet_id=0)
    assert list(book) == ["gut", "bgc"]
    # Rows keep the derived order, not distill's sort.
    assert book["gut"]["gene_id"].to_list() == list(counts)


def test_built_in_sheets_are_strict(tmp_path, ecosystems):
    write(tmp_path, "gene_descriptions.tsv", DESCRIPTIONS.replace("K22222", "K0"))
    result, out = run(tmp_path, "--sum", "--topics", "none", "--ecos", "gut", ok=False)
    assert result.exit_code != 0
    assert "gut: K22222" in str(result.exception)
    assert not out.exists()


def test_viz_custom_sheets_get_a_custom_derived_summary(tmp_path, ecosystems):
    viz = write(tmp_path, "mine.tsv", "name\trule\tgroup\nR\tK00399 | K77777\tG\n")
    _, out = run(
        tmp_path,
        "--sum",
        "--topics",
        "none",
        "--ecos",
        "none",
        "--viz_custom",
        str(viz),
    )
    assert sorted(p.name for p in (out / "summarize").iterdir()) == [
        "custom_derived_summary.tsv",
        "custom_derived_summary.xlsx",
    ]
    df = pl.read_csv(out / "summarize/custom_derived_summary.tsv", separator="\t")
    # Not strict: K77777 has no description, and the run still succeeds.
    assert df.select("gene_id", "gene_description", "sheet").rows() == [
        ("K00399", "mcrA", "mine"),
        ("K77777", None, "mine"),
    ]


def test_descriptions_option_adds_custom_tables(tmp_path, ecosystems):
    viz = write(tmp_path, "mine.tsv", "name\trule\tgroup\nR\tK77777\tG\n")
    mine = write(tmp_path, "d.tsv", "gene_id\tdescription\nK77777\tmine\n")
    _, out = run(
        tmp_path, "--sum", "--topics", "none", "--ecos", "none",
        "--viz_custom", str(viz), "--descriptions", str(mine),
    )  # fmt: skip
    df = pl.read_csv(out / "summarize/custom_derived_summary.tsv", separator="\t")
    assert df["gene_description"].to_list() == ["mine"]


def test_viz_custom_without_summarize_writes_no_summary(tmp_path, ecosystems):
    viz = write(tmp_path, "mine.tsv", "name\trule\tgroup\nR\tK00399\tG\n")
    _, out = run(tmp_path, "--viz", "--ecos", "none", "--viz_custom", str(viz))
    assert not (out / "summarize").exists()
    assert (out / "visualize" / "mine").is_dir()


def test_built_in_ecosystem_sheets_and_tables_are_packaged():
    assert set(cli.ECOSYSTEMS) == {"ag", "aquatic", "bgc", "eng_sys", "gut"}
    assert make_summarize.DBCAN_FAMILIES.is_file()
    assert make_summarize.GENE_DESCRIPTIONS.is_file()


def test_every_built_in_gene_id_has_a_description():
    sheets = {
        name: make_summarize.load_rule_sheet(path, cli.COMMON_RULES, ())
        for name, path in cli.ECOSYSTEMS.items()
    }
    make_summarize.derive_collection(
        sheets, make_summarize.load_descriptions([]), strict=True
    )
