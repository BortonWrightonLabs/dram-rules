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


# The built-in Summarize Sheets, by their --topics/--ecos names.
BUILT_IN = {
    "assim": "assimilation_and_cofactor_metabolism.tsv",
    "cell": "cellular_machinery.tsv",
    "energy": "energy_acquisition_bioenergetics.tsv",
    "env": "environmental_interaction_and_adaptation.tsv",
    "ag": "distill_ag.tsv",
    "eng_sys": "distill_eng_sys.tsv",
    "bgc": "distill_bgc.tsv",
    "gut": "distill_gut.tsv",
    "aquatic": "distill_aquatic.tsv",
}


@pytest.fixture(autouse=True)
def fake_built_in_sheets(tmp_path, monkeypatch):
    """Stand-ins for the built-in sheets (which fail the duplicate check until
    they are deduplicated). Each counts K00399, in a Group named by its key."""
    sheets = tmp_path / "built_in"
    sheets.mkdir()
    for key, name in BUILT_IN.items():
        (sheets / name).write_text(f"gene_id\ttopic_ecosystem\nK00399\t{key}\n")
    monkeypatch.setattr(make_summarize, "SUMMARIZE_SHEETS_DIR", sheets)


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


def summarize(
    tmp_path: Path,
    *sheets: str,
    trna: str | None = None,
    annotations: str = ANNOTATIONS,
    mapping: str | None = None,
    **kwargs,
) -> pl.DataFrame:
    """The `summarize` table of one Summary Collection of `sheets` (TSV text).
    `kwargs` go to `make_summarize.summarize`."""
    common = write(tmp_path, "common.tsv", "name\talias\trule\n")
    compiled = [
        make_summarize.load_summarize_sheet(
            Path(write(tmp_path, f"s{i}.tsv", text)), Path(common), ()
        )
        for i, text in enumerate(sheets)
    ]
    hits = Hits.load(
        Path(write(tmp_path, "anno.tsv", annotations)),
        compiled,
        trna=Path(write(tmp_path, "trna.tsv", trna)) if trna else None,
        mapping=Path(write(tmp_path, "map.tsv", mapping)) if mapping else None,
    )
    return make_summarize.summarize({"c": compiled}, hits, **kwargs)["c"]


# The spec's worked example: bin2_2 carries no KO, so K12345's genes are
# bin1_1, bin1_2 and bin2_1.
EXAMPLE_ANNOTATIONS = """\
query_id\tinput_fasta\tkegg_id
bin1_1\tbin1\tK12345
bin1_2\tbin1\tK12345
bin2_1\tbin2\tK12345
bin2_2\tbin2\t
"""
EXAMPLE_MAPPING = """\
Geneid\tsample1\tsample2
bin1_1\t0\t1
bin1_2\t1\t1
bin2_1\t0\t2
bin2_2\t2\t2
"""
EXAMPLE_SHEET = "gene_id\ttopic_ecosystem\nK12345\tEnergy\n"


def example(tmp_path: Path, **kwargs) -> dict:
    """The worked example's one summary row, without its Group."""
    df = summarize(
        tmp_path,
        EXAMPLE_SHEET,
        annotations=EXAMPLE_ANNOTATIONS,
        mapping=EXAMPLE_MAPPING,
        **kwargs,
    )
    return df.drop("topic_ecosystem").row(0, named=True)


def test_gene_abundance_sums_each_sample_over_the_rules_genes(tmp_path):
    assert example(tmp_path, abundance="sum") == {
        "gene_id": "K12345",
        "bin1": 2,
        "bin2": 1,
        "sum_abundance_sample1": 1.0,
        "sum_abundance_sample2": 4.0,
    }


def test_gene_abundance_averages_over_the_rules_genes_by_default(tmp_path):
    row = example(tmp_path)
    assert row["mean_abundance_sample1"] == pytest.approx(1 / 3)
    assert row["mean_abundance_sample2"] == pytest.approx(4 / 3)
    assert "sum_abundance_sample1" not in row


# g1 carries K1 and K2, g2 carries K1; g3 holds no rule. All in genome G.
ONCE_ANNOTATIONS = """\
query_id\tinput_fasta\tkegg_id
g1\tG\tK1
g1\tG\tK2
g2\tG\tK1
g3\tG\tK9
"""
ONCE_SHEET = "name\ttopic_ecosystem\trule\nEither\tE\tK1 | K2\n"


@pytest.mark.parametrize(
    "mapping, expected_sum, expected_mean",
    [
        # Once per gene: 10 + 4, where Sample Abundance (once per Hit) gives 24.
        ("Geneid\tS\ng1\t10\ng2\t4\ng3\t100\n", 14.0, 7.0),
        # Per genome: both of G's genes take G's 10.
        ("genome\tS\nG\t10\n", 20.0, 10.0),
    ],
)
def test_gene_abundance_counts_each_gene_once(
    tmp_path, mapping, expected_sum, expected_mean
):
    def gene_abundance(abundance):
        df = summarize(
            tmp_path,
            ONCE_SHEET,
            annotations=ONCE_ANNOTATIONS,
            mapping=mapping,
            abundance=abundance,
        )
        return df[f"{abundance}_abundance_S"].item()

    assert gene_abundance("sum") == expected_sum
    assert gene_abundance("mean") == expected_mean


def test_unmapped_genes_count_as_zero_in_the_mean(tmp_path):
    # The mapping file lists only g1, so g2 counts as 0 and stays a gene of the rule.
    df = summarize(
        tmp_path,
        ONCE_SHEET + "Nothing\tE\tK7\n",
        annotations=ONCE_ANNOTATIONS,
        mapping="Geneid\tS\ng1\t10\n",
    )
    assert df["mean_abundance_S"].to_list() == [5.0, 0.0]


def test_sample_columns_follow_the_mapping_file_after_the_genomes(tmp_path):
    # A sample named like a genome doesn't collide with its Gene Count column.
    df = summarize(
        tmp_path,
        SHEET,
        mapping="Geneid\tzeta\tbin1\nbin1_1\t1\t2\n",
        abundance="sum",
    )
    assert df.columns[-4:] == [
        "bin1",
        "bin2",
        "sum_abundance_zeta",
        "sum_abundance_bin1",
    ]
    assert counts(df)["K00399"] == {"bin1": 2, "bin2": 1}
    assert df["sum_abundance_bin1"].to_list() == [2.0, 0.0, 0.0]


def test_custom_sheet_gene_counts(tmp_path):
    assert counts(summarize(tmp_path, SHEET)) == {
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
    out = run(tmp_path, "--sum", "--sum_custom", write(tmp_path, "b.tsv", sheet))
    assert not (out / "custom_summary.tsv").exists()
    assert not (out / "custom_summary.xlsx").exists()


def worksheets(path: Path) -> list[str]:
    return list(pl.read_excel(path, sheet_id=0))


def test_topics_and_ecos_write_their_own_collections(tmp_path):
    out = run(tmp_path, "--sum", "--topics", "energy,env", "--ecos", "gut")
    assert worksheets(out / "topic_summary.xlsx") == ["energy", "env"]
    assert worksheets(out / "ecosystem_summary.xlsx") == ["gut"]
    df = pl.read_csv(out / "ecosystem_summary.tsv", separator="\t")
    assert counts(df) == {"K00399": {"bin1": 2, "bin2": 1}}
    assert (out / "topic_summary.tsv").exists()
    assert not (out / "custom_summary.tsv").exists()
    assert not (out / "custom_summary.xlsx").exists()


def test_defaults_are_all_topics_and_ecosystems(tmp_path):
    out = run(tmp_path, "--sum")
    topics = ["assim", "cell", "energy", "env"]
    assert worksheets(out / "topic_summary.xlsx") == topics
    ecos = ["ag", "eng_sys", "bgc", "gut", "aquatic"]
    assert worksheets(out / "ecosystem_summary.xlsx") == ecos


@pytest.mark.parametrize("eco", ["legacy", "default"])
def test_visualize_only_ecos_are_an_error(tmp_path, eco):
    result = invoke(tmp_path, "--sum", "--ecos", f"gut,{eco}")
    assert result.exit_code != 0
    assert eco in result.output
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


def test_column_overrides_read_custom_sheets_not_built_in(tmp_path):
    sheet = write(
        tmp_path, "a.tsv",
        "my_label\tmy_alias\tmy_group\trule\n"
        "\tmcr\t\tK00399\n"
        "Mcr\t\tEnergy\t@mcr\n",
    )  # fmt: skip
    out = run(
        tmp_path, "--sum", "--topics", "energy", "--ecos", "none", "--sum_custom", sheet,
        "-lc", "my_label", "-acc", "my_alias", "-gc", "my_group",
    )  # fmt: skip
    df = pl.read_csv(out / "custom_summary.tsv", separator="\t")
    assert counts(df, "my_label") == {"Mcr": {"bin1": 2, "bin2": 1}}
    assert list(pl.read_excel(out / "custom_summary.xlsx", sheet_id=0)) == ["Energy"]
    # The built-in sheet still uses the lookup (gene_id, topic_ecosystem).
    topic = pl.read_csv(out / "topic_summary.tsv", separator="\t")
    assert counts(topic) == {"K00399": {"bin1": 2, "bin2": 1}}


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
        "--traits", "--traits_custom", traits, "-c", str(COMMON_RULES),
        "--sum", "--sum_custom", sheet,
    )  # fmt: skip
    # One build of the Hits, for the Gene IDs of every refinement's sheets.
    [needed] = calls
    assert {"K00399", "CYC2_REPCLUSTER1", "K11111", "K99999"} <= needed


@pytest.mark.parametrize(
    "args, column, expected",
    [
        ([], "mean_abundance_S1", [7 / 3, 2.5, 0.0]),
        (["--abundance", "sum"], "sum_abundance_S1", [7.0, 5.0, 0.0]),
    ],
)
def test_mapping_adds_gene_abundance_to_summaries(tmp_path, args, column, expected):
    sheet = write(tmp_path, "a.tsv", SHEET)
    mapping = write(tmp_path, "m.tsv", MAPPING)
    out = run(tmp_path, "--sum", "--sum_custom", sheet, "-m", mapping, *args)
    df = pl.read_csv(out / "custom_summary.tsv", separator="\t")
    assert counts(df)["K00399"] == {"bin1": 2, "bin2": 1}
    assert df.columns[-1] == column
    assert df[column].to_list() == pytest.approx(expected)
    book = pl.read_excel(out / "custom_summary.xlsx", sheet_id=0)
    assert book["Energy"].columns[-1] == column


def test_abundance_without_mapping_changes_nothing(tmp_path):
    sheet = write(tmp_path, "a.tsv", SHEET)
    out = run(tmp_path, "--sum", "--sum_custom", sheet)
    plain = (out / "custom_summary.tsv").read_bytes()
    out = run(tmp_path, "--sum", "--sum_custom", sheet, "--abundance", "sum")
    assert (out / "custom_summary.tsv").read_bytes() == plain


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
        "visualize": ["--viz", "--viz_custom", viz],
        "traits": ["--traits", "--traits_custom", viz],
        "summarize": ["--sum", "--sum_custom", sheet, "--trna_path", trna],
    }
    # The built-in traits sheet needs the real Common Alias Sheet.
    shared = [
        "--ecos",
        "none",
        "--topics",
        "none",
        "-m",
        mapping,
        "-c",
        str(COMMON_RULES),
    ]
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
