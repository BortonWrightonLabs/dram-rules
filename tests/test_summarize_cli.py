"""The summarize refinement: `summarize()` directly, and `dram-rules --sum` end to end."""

import shutil
from pathlib import Path

import polars as pl
import pytest
from click.testing import CliRunner

from dram_rules import make_summarize
from dram_rules.cli import COMMON_RULES, main
from dram_rules.hits import Hits

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


# The built-in topic Summarize Sheets, by their --topics names.
BUILT_IN = {
    "assim": "assimilation_and_cofactor_metabolism.tsv",
    "cell": "cellular_machinery.tsv",
    "energy": "energy_acquisition_bioenergetics.tsv",
    "env": "environmental_interaction_and_adaptation.tsv",
}


@pytest.fixture
def ecosystems(tmp_path, monkeypatch):
    """Every built-in ecosystem sheet replaced with one describable rule, so
    strict ecosystem summaries pass (tests/test_derived_summaries.py covers
    deriving them)."""
    from dram_rules import cli

    sheet = tmp_path / "eco.tsv"
    sheet.write_text("name\tgroup\trule\nMcr\tG\tK00399\n")
    monkeypatch.setattr(cli, "ECOSYSTEMS", dict.fromkeys(cli.ECOSYSTEMS, sheet))
    monkeypatch.setattr(
        make_summarize, "load_descriptions", lambda custom: {"K00399": "mcrA"}
    )


# Summarize only the custom sheets, skipping the slow-to-load built-in ones.
CUSTOM_ONLY = ("--sum", "--topics", "none", "--ecos", "none")


def invoke(tmp_path: Path, *args: str, annotations: str = ANNOTATIONS):
    (tmp_path / "anno.tsv").write_text(annotations)
    cli_args = [
        "-a", str(tmp_path / "anno.tsv"),
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


def summarize(tmp_path: Path, *sheets: str, trna: str | None = None) -> pl.DataFrame:
    """The `summarize` table of one Summary Collection of `sheets` (TSV text)."""
    common = write(tmp_path, "common.tsv", "name\talias\trule\n")
    compiled = [
        make_summarize.load_summarize_sheet(
            Path(write(tmp_path, f"s{i}.tsv", text)), Path(common), ()
        )
        for i, text in enumerate(sheets)
    ]
    hits = Hits.load(
        Path(write(tmp_path, "anno.tsv", ANNOTATIONS)),
        compiled,
        trna=Path(write(tmp_path, "trna.tsv", trna)) if trna else None,
    )
    return make_summarize.summarize({"c": compiled}, hits)["c"]


def test_custom_sheet_gene_counts(tmp_path):
    assert counts(summarize(tmp_path, SHEET)) == {
        "K00399": {"bin1": 2, "bin2": 1},
        "Cyc2_repCluster1": {"bin1": 2, "bin2": 0},
        "K99999": {"bin1": 0, "bin2": 0},
    }


def test_workbook_has_a_sorted_worksheet_per_group(tmp_path):
    sheet = write(tmp_path, "a.tsv", SHEET)
    out = run(tmp_path, *CUSTOM_ONLY, "--sum_custom", sheet)
    book = pl.read_excel(out / "custom_summary.xlsx", sheet_id=0)
    assert list(book) == ["Energy", "Other"]
    energy = book["Energy"]
    # Only the summary columns: no `notes`, and the Group is the worksheet.
    assert energy.columns == ["gene_id", "gene_description", "module", "bin1", "bin2"]
    # Sorted by module (distill's order: header, subheader, module, label).
    assert energy["gene_id"].to_list() == ["Cyc2_repCluster1", "K00399"]
    assert energy.row(1) == ("K00399", "mcrA", "Methanogenesis", 2, 1)


def test_sheets_with_different_label_columns_use_name(tmp_path):
    df = summarize(
        tmp_path, SHEET, "name\ttopic_ecosystem\trule\nMcr\tEnergy\tK00399\n"
    )
    assert df.columns[0] == "name"
    assert counts(df, "name")["Mcr"] == {"bin1": 2, "bin2": 1}
    assert counts(df, "name")["K00399"] == {"bin1": 2, "bin2": 1}


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
    assert counts(summarize(tmp_path, TRNA_SHEET, trna=TRNA)) == {
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
    df = summarize(tmp_path, SHEET, trna="NULL\n")
    assert counts(df)["K00399"] == {"bin1": 2, "bin2": 1}


NEEDS_TAXONOMY = '{"include_when": {"all_columns_present": ["taxonomy"]}}'


def test_rules_left_out_by_include_when_are_absent(tmp_path):
    # The annotations have no taxonomy, so Group "Taxa" loses its only rule
    # and sheet b loses every rule.
    a = (
        "gene_id\ttopic_ecosystem\trule_options\n"
        "K00399\tEnergy\t\n"
        f"Cyc2_repCluster1\tTaxa\t{NEEDS_TAXONOMY}\n"
    )
    b = f"name\ttopic_ecosystem\trule_options\nK00399\tEnergy\t{NEEDS_TAXONOMY}\n"
    df = summarize(tmp_path, a, b)
    # Sheet b added nothing, so the label column stays a's `gene_id`.
    assert counts(df) == {"K00399": {"bin1": 2, "bin2": 1}}
    assert df["topic_ecosystem"].to_list() == ["Energy"]


def test_no_collection_when_every_rule_is_left_out(tmp_path):
    sheet = (
        f"gene_id\ttopic_ecosystem\trule_options\nK00399\tEnergy\t{NEEDS_TAXONOMY}\n"
    )
    assert summarize(tmp_path, sheet) is None
    # The command then writes no workbook.
    out = run(tmp_path, *CUSTOM_ONLY, "--sum_custom", write(tmp_path, "b.tsv", sheet))
    assert not (out / "custom_summary.tsv").exists()
    assert not (out / "custom_summary.xlsx").exists()


def worksheets(path: Path) -> set[str]:
    return set(pl.read_excel(path, sheet_id=0))


def groups(*keys: str) -> set[str]:
    """The Groups (topic_ecosystem values) of the built-in sheets `keys`."""
    sheets = [
        make_summarize.load_summarize_sheet(
            make_summarize.SUMMARIZE_SHEETS_DIR / BUILT_IN[key], COMMON_RULES, ()
        )
        for key in keys
    ]
    return {g for s in sheets for g in s.df[s.group_col].drop_nulls()}


def test_topics_and_ecos_write_their_own_collections(tmp_path, ecosystems):
    out = run(tmp_path, "--sum", "--topics", "energy,env", "--ecos", "gut")
    assert worksheets(out / "topic_summary.xlsx") == groups("energy", "env")
    assert worksheets(out / "ecosystem_summary.xlsx") == {"gut"}
    df = pl.read_csv(out / "ecosystem_summary.tsv", separator="\t")
    assert counts(df)["K00399"] == {"bin1": 2, "bin2": 1}
    assert (out / "topic_summary.tsv").exists()
    assert not (out / "custom_summary.tsv").exists()
    assert not (out / "custom_summary.xlsx").exists()


def test_defaults_are_all_topics_and_ecosystems(tmp_path, ecosystems):
    from dram_rules.cli import ECOSYSTEMS

    out = run(tmp_path, "--sum")
    topics = groups("assim", "cell", "energy", "env")
    assert worksheets(out / "topic_summary.xlsx") == topics
    assert worksheets(out / "ecosystem_summary.xlsx") == set(ECOSYSTEMS)


def test_all_is_every_topic_and_ecosystem_sheet(tmp_path, ecosystems):
    from dram_rules.cli import ECOSYSTEMS

    out = run(tmp_path, "--sum", "--topics", "all,energy", "--ecos", "all,gut")
    assert worksheets(out / "topic_summary.xlsx") == groups(*BUILT_IN)
    assert worksheets(out / "ecosystem_summary.xlsx") == set(ECOSYSTEMS)


def test_no_built_in_sheet_has_a_reserved_name():
    from dram_rules.cli import ECOSYSTEMS, TOPICS, VISUALIZE_ONLY_ECOS

    # all and none are split_choices' words; default and legacy name other sheets.
    assert not {"all", "none", *VISUALIZE_ONLY_ECOS} & {*ECOSYSTEMS, *TOPICS}


def test_summarize_skips_visualize_only_ecos_with_a_warning(
    tmp_path, ecosystems, caplog
):
    with caplog.at_level("WARNING", logger="dram.rules"):
        out = run(tmp_path, "--viz", "--sum", "--ecos", "gut,legacy")
    assert sorted(p.name for p in (out.parent / "visualize").iterdir()) == [
        "gut",
        "legacy",
    ]
    assert worksheets(out / "ecosystem_summary.xlsx") == {"gut"}
    [warning] = [r for r in caplog.records if r.name == "dram.rules"]
    assert "legacy" in warning.getMessage()


def test_summarize_skips_default_rather_than_reading_it_as_gut(tmp_path, caplog):
    with caplog.at_level("WARNING", logger="dram.rules"):
        out = run(tmp_path, "--sum", "--topics", "none", "--ecos", "default")
    assert not out.exists()
    [warning] = [r for r in caplog.records if r.name == "dram.rules"]
    assert "default" in warning.getMessage()


def test_viz_custom_sheets_sharing_a_name_are_an_error(tmp_path):
    for folder in ("a", "b"):
        (tmp_path / folder).mkdir()
        write(tmp_path, f"{folder}/x.tsv", "name\tgroup\trule\nMcr\tG\tK00399\n")
    sheets = f"{tmp_path / 'a/x.tsv'},{tmp_path / 'b/x.tsv'}"
    result = invoke(tmp_path, *CUSTOM_ONLY, "--viz_custom", sheets)
    assert result.exit_code != 0
    assert "share a worksheet in custom_derived_summary: ['x']" in result.output
    assert not (tmp_path / "out").exists()


def test_unknown_topic_is_an_error(tmp_path):
    result = invoke(tmp_path, "--sum", "--topics", "energy,metals")
    assert result.exit_code != 0
    assert "metals" in result.output


@pytest.mark.parametrize("none", ["", "none"])
def test_no_topics_and_ecos_leave_only_custom_sheets(tmp_path, none):
    sheet = write(tmp_path, "a.tsv", SHEET)
    out = run(
        tmp_path, "--sum", "--topics", none, "--ecos", none, "--sum_custom", sheet
    )
    assert sorted(p.name for p in out.iterdir()) == [
        "custom_summary.tsv",
        "custom_summary.xlsx",
    ]


def test_built_in_sheets_are_packaged():
    data = Path(make_summarize.__file__).parent / "data/summarize_sheets"
    assert all((data / name).is_file() for name in BUILT_IN.values())


def test_column_overrides_read_custom_sheets_not_built_in(tmp_path, ecosystems):
    sheet = write(
        tmp_path, "a.tsv",
        "my_label\tmy_alias\tmy_group\trule\n"
        "\tmcr\t\tK00399\n"
        "Mcr\t\tEnergy\t@mcr\n",
    )  # fmt: skip
    out = run(
        tmp_path, "--sum", "--topics", "none", "--ecos", "gut", "--sum_custom", sheet,
        "-lc", "my_label", "-acc", "my_alias", "-gc", "my_group",
    )  # fmt: skip
    df = pl.read_csv(out / "custom_summary.tsv", separator="\t")
    assert counts(df, "my_label") == {"Mcr": {"bin1": 2, "bin2": 1}}
    assert list(pl.read_excel(out / "custom_summary.xlsx", sheet_id=0)) == ["Energy"]
    # The built-in sheet still uses the lookup (name, group).
    eco = pl.read_csv(out / "ecosystem_summary.tsv", separator="\t")
    assert counts(eco)["K00399"] == {"bin1": 2, "bin2": 1}


@pytest.mark.parametrize("alias_column", ["aliases", "parent"])
def test_backup_alias_columns_hold_aliases(tmp_path, alias_column):
    sheet = write(
        tmp_path, "a.tsv",
        f"gene_id\t{alias_column}\ttopic_ecosystem\trule\n"
        "\tmcr\t\tK00399\n"
        "Mcr\t\tEnergy\t@mcr\n",
    )  # fmt: skip
    out = run(
        tmp_path, "--sum", "--topics", "none", "--ecos", "none", "--sum_custom", sheet
    )
    df = pl.read_csv(out / "custom_summary.tsv", separator="\t")
    assert counts(df) == {"Mcr": {"bin1": 2, "bin2": 1}}


def test_a_sheet_without_a_group_column_names_the_sheet(tmp_path):
    sheet = write(tmp_path, "nogroup.tsv", "gene_id\tgene_description\nK00399\tmcrA\n")
    result = invoke(
        tmp_path, "--sum", "--topics", "none", "--ecos", "none", "--sum_custom", sheet
    )
    assert result.exit_code != 0
    assert "nogroup.tsv" in str(result.exception)
    assert "No Group column" in str(result.exception)


# A visualize/traits sheet, for runs combining refinements.
VIZ_SHEET = "name\tgroup\trule\nMcr\tG1\tK00399\nCyc2\tG1\tCyc2_repCluster1\n"
MAPPING = "Geneid\tS1\nbin1_1\t1\nbin1_2\t2\nbin1_3\t3\nbin2_1\t4\n"


def test_a_bad_sheet_in_any_refinement_stops_the_run_before_any_output(tmp_path):
    viz = write(tmp_path, "viz.tsv", VIZ_SHEET)
    bad = write(tmp_path, "bad.tsv", SHEET + "K00399\tagain\tM\tEnergy\t\n")
    result = invoke(
        tmp_path, "--viz", "--ecos", "none", "--viz_custom", viz,
        "--sum", "--topics", "none", "--sum_custom", bad,
    )  # fmt: skip
    assert result.exit_code != 0
    assert "bad.tsv" in str(result.exception)
    assert not (tmp_path / "out").exists()


def test_a_failure_while_evaluating_leaves_no_outputs(tmp_path):
    # Compiles, but reading a missing column fails at evaluate (ADR 0001).
    sheet = write(
        tmp_path, "a.tsv",
        "name\ttopic_ecosystem\trule\n"
        'Bad\tEnergy\tfilter_contains(no_such_column, "x") -> count_rows(>=, 1)\n',
    )  # fmt: skip
    viz = write(tmp_path, "viz.tsv", VIZ_SHEET)
    result = invoke(
        tmp_path, "--viz", "--ecos", "none", "--viz_custom", viz,
        "--sum", "--sum_custom", sheet,
    )  # fmt: skip
    assert result.exit_code != 0
    assert "no_such_column" in str(result.exception)
    assert not (tmp_path / "out").exists()


def test_refinements_share_one_load_of_hits(tmp_path, monkeypatch):
    from dram_rules import hits

    calls = []
    prepare = hits.prepare_present_map_df
    monkeypatch.setattr(
        hits,
        "prepare_present_map_df",
        lambda *a, **kw: calls.append(kw["needed_features"]) or prepare(*a, **kw),
    )
    viz = write(tmp_path, "viz.tsv", VIZ_SHEET)
    traits = write(tmp_path, "traits.tsv", "name\tgroup\trule\nOnly\tG\tK11111\n")
    sheet = write(tmp_path, "a.tsv", SHEET)
    run(
        tmp_path, "--viz", "--ecos", "none", "--viz_custom", viz,
        "--traits", "--traits_custom", traits,
        *CUSTOM_ONLY, "--sum_custom", sheet,
    )  # fmt: skip
    # One build of the Hits, for the Gene IDs of every refinement's sheets.
    [needed] = calls
    assert {"K00399", "CYC2_REPCLUSTER1", "K11111", "K99999"} <= needed


def test_summarize_ignores_the_mapping_file(tmp_path, caplog):
    sheet = write(tmp_path, "a.tsv", SHEET)
    mapping = write(tmp_path, "m.tsv", MAPPING)
    with caplog.at_level("INFO", logger="dram.rules"):
        out = run(tmp_path, *CUSTOM_ONLY, "--sum_custom", sheet, "-m", mapping)
    assert "ignores --mapping" in caplog.text
    df = pl.read_csv(out / "custom_summary.tsv", separator="\t")
    assert counts(df)["K00399"] == {"bin1": 2, "bin2": 1}
    assert counts(df)["Cyc2_repCluster1"] == {"bin1": 2, "bin2": 0}


def files(folder: Path) -> dict[str, pl.DataFrame]:
    """Every TSV under `folder` by relative path, rows sorted."""
    return {
        str(p.relative_to(folder)): (df := pl.read_csv(p, separator="\t")).sort(
            df.columns
        )
        for p in sorted(folder.rglob("*.tsv"))
    }


def test_each_refinement_in_a_combined_run_matches_running_it_alone(tmp_path):
    viz = write(tmp_path, "viz.tsv", VIZ_SHEET)
    trna = write(tmp_path, "trna.tsv", TRNA)
    sheet = write(tmp_path, "a.tsv", TRNA_SHEET)
    mapping = write(tmp_path, "m.tsv", MAPPING)
    refinements = {
        "visualize": ["--viz"],
        "traits": ["--traits", "--traits_custom", viz],
        "summarize": ["--sum", "--sum_custom", sheet, "--trna_path", trna],
    }
    # --viz_custom also gives summarize a custom_derived_summary, so every run
    # passes it; only --viz makes its heatmaps.
    shared = ["--ecos", "none", "--topics", "none", "-m", mapping, "--viz_custom", viz]
    out = run(
        tmp_path, *shared, *(a for args in refinements.values() for a in args)
    ).parent
    combined = files(out)
    for refinement, args in refinements.items():
        shutil.rmtree(out)
        run(tmp_path, *shared, *args)
        alone = files(out)
        assert alone.keys() == {
            k for k in combined if k.startswith(refinement + "/")
        }, refinement
        for key, df in alone.items():
            assert df.equals(combined[key]), key
