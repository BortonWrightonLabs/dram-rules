"""End-to-end runs of `dram-rules --viz` (the visualize step)."""

from pathlib import Path

import polars as pl
import pytest
from click.testing import CliRunner

from dram_rules.cli import TRAITS_SHEET, main

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
    """Run `--viz` on one custom sheet written to `tmp_path`; return its folder."""
    (tmp_path / "rules.tsv").write_text(rules)
    invoke_ok(
        tmp_path,
        "--ecos", "none", "--viz_custom", str(tmp_path / "rules.tsv"),
        *args,
        annotations=annotations,
        mapping=mapping,
    )  # fmt: skip
    return tmp_path / "out" / "visualize" / "rules"


def invoke(
    tmp_path: Path,
    *args: str,
    annotations: str = ANNOTATIONS,
    mapping: str | None = None,
    out: str = "out",
    empty_common: bool = True,
    step: str = "--viz",
):
    """`dram-rules --viz` on small inputs written to `tmp_path`. Built-in
    sheets need the default Common Alias Sheet (`empty_common=False`)."""
    (tmp_path / "anno.tsv").write_text(annotations)
    cli_args = ["-a", str(tmp_path / "anno.tsv"), "-o", str(tmp_path / out), step]
    if empty_common:
        (tmp_path / "common.tsv").write_text("name\talias\trule\n")
        cli_args += ["-c", str(tmp_path / "common.tsv")]
    if mapping is not None:
        (tmp_path / "mapping.tsv").write_text(mapping)
        cli_args += ["-m", str(tmp_path / "mapping.tsv")]
    return CliRunner().invoke(main, [*cli_args, *args])


def invoke_ok(tmp_path: Path, *args: str, out: str = "out", **kwargs) -> Path:
    result = invoke(tmp_path, *args, out=out, **kwargs)
    assert result.exit_code == 0, result.output + repr(result.exception)
    return tmp_path / out / "visualize"


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


def test_genome_mapping_drops_text_columns_and_counts_zeros_in_mean(tmp_path):
    # Keyed by genome: each gene takes its genome's abundance. `Notes` is text,
    # so it isn't a sample.
    mapping = "genome\tNotes\tS1\tS2\nbin1\tx\t2\t0\nbin2\ty\t5\t6\n"
    out = run(tmp_path, mapping=mapping)
    sample = read(out, "sample")
    assert set(sample["sample"]) == {"S1", "S2"}
    assert value(sample, "sample_abundance", name="Cyc2", sample="S1") == 2.0
    assert value(sample, "present", name="Cyc2", sample="S2") is False
    genome = read(out, "genome")
    assert value(genome, "mean_sample_abundance", name="Cyc2", genome="bin1") == 1.0


def test_row_reading_rule_with_mapping(tmp_path):
    rules = (
        "name\tgroup\trule\n"
        "Heme\tG1\tfilter_compare(heme_regulatory_motif_count, >=, 4)"
        " -> count_rows(>=, 1)\n"
    )
    # bin1_1 is the only gene with >= 4 heme motifs, and only S1 has it.
    mapping = "Geneid\tS1\tS2\nbin1_1\t1\t0\nbin1_2\t3\t4\nbin2_1\t5\t6\n"
    out = run(tmp_path, rules=rules, mapping=mapping)
    sample = read(out, "sample")
    assert value(sample, "present", name="Heme", sample="S1") is True
    assert value(sample, "present", name="Heme", sample="S2") is False


def test_rules_named_in_label_column(tmp_path):
    out = run(tmp_path, rules=RULES.replace("name\t", "label\t", 1), mapping=None)
    assert set(read(out, "genome")["name"]) == {"Mcr", "Cyc2"}
    assert (out / "heatmap_genome_presence_coverage.html").exists()


def test_rules_named_by_label_column_option(tmp_path):
    rules = RULES.replace("name\t", "trait\t", 1)
    out = run(tmp_path, "-lc", "trait", rules=rules, mapping=None)
    assert set(read(out, "genome")["name"]) == {"Mcr", "Cyc2"}
    assert (out / "heatmap_genome_presence_coverage.html").exists()


# Labelled, aliased and grouped by columns the lookup doesn't know.
CUSTOM_COLUMNS = """\
my_label\tmy_alias\tmy_group\trule
\tmcr\t\tK00399
Mcr\t\tG1\t@mcr
Cyc2\t\tG2\tCyc2_repCluster1
"""
OVERRIDES = ("-lc", "my_label", "-acc", "my_alias", "-gc", "my_group")


@pytest.mark.parametrize("step", ["--viz", "--traits"])
def test_column_overrides_read_custom_sheets_not_built_in(tmp_path, step):
    (tmp_path / "mine.tsv").write_text(CUSTOM_COLUMNS)
    custom = "--viz_custom" if step == "--viz" else "--traits_custom"
    # The built-in gut and traits sheets have no my_label column: they keep the lookup.
    invoke_ok(
        tmp_path, "--ecos", "gut", custom, str(tmp_path / "mine.tsv"), *OVERRIDES,
        step=step, empty_common=False,
    )  # fmt: skip
    folder = tmp_path / "out" / step.strip("-").replace("viz", "visualize")
    genome = read(folder / "mine", "genome")
    assert value(genome, "present", name="Mcr", genome="bin2") is True
    assert value(genome, "my_group", name="Cyc2", genome="bin1") == "G2"
    built_in = folder / ("gut" if step == "--viz" else "default")
    assert (built_in / "df_vertical_genome.tsv").exists()


def test_override_column_missing_from_a_custom_sheet_names_it(tmp_path):
    (tmp_path / "rules.tsv").write_text(RULES)
    result = invoke(
        tmp_path, "--ecos", "none", "--viz_custom", str(tmp_path / "rules.tsv"),
        "-lc", "my_label",
    )  # fmt: skip
    assert result.exit_code != 0
    assert "my_label" in str(result.exception)
    assert "rules.tsv" in str(result.exception)


@pytest.mark.parametrize("old", ["-l", "-ac", "-g"])
def test_all_sheet_column_overrides_are_gone(tmp_path, old):
    result = invoke(tmp_path, "--ecos", "gut", old, "name")
    assert result.exit_code == 2
    # -ac reads as -a c, then "name" is an extra argument.
    assert ("No such option" if old != "-ac" else "extra argument") in result.output


SUMMARIZE_SHEET = """\
gene_id\tgene_description\ttopic_ecosystem
K00399\tmcrA\tMethanogenesis
Cyc2_repCluster1\tcyc2\tIron
"""


def test_summarize_sheet_draws_a_heatmap(tmp_path):
    out = run(tmp_path, rules=SUMMARIZE_SHEET, mapping=None)
    df = read(out, "genome")
    assert value(df, "present", genome="bin1", name="K00399") is True
    assert value(df, "present", genome="bin2", name="Cyc2_repCluster1") is False
    assert list(out.glob("*.html"))


def test_repeated_label_names_the_sheet(tmp_path):
    (tmp_path / "rules.tsv").write_text(RULES + "Mcr\tG1\tK00400\n")
    result = invoke(
        tmp_path, "--ecos", "none", "--viz_custom", str(tmp_path / "rules.tsv")
    )
    assert "rules.tsv" in str(result.exception)
    assert "Mcr" in str(result.exception)


# Mcr is in both sheets, with different expressions.
SHEET_A = "name\tgroup\trule\nMcr\tG1\tK00399\nCyc2\tG1\tCyc2_repCluster1\n"
SHEET_B = "name\tgroup\trule\nMcr\tG2\tCyc2_repCluster1\n"


def tables(folder: Path) -> dict[str, pl.DataFrame]:
    """Each table in a sheet's folder, rows sorted."""
    return {
        p.name: (df := pl.read_csv(p, separator="\t")).sort(df.columns)
        for p in sorted(folder.glob("*.tsv"))
    }


def same_tables(a: Path, b: Path) -> bool:
    ta, tb = tables(a), tables(b)
    return ta.keys() == tb.keys() and all(ta[k].equals(tb[k]) for k in ta)


def test_each_sheet_matches_a_single_sheet_run(tmp_path):
    (tmp_path / "a.tsv").write_text(SHEET_A)
    (tmp_path / "b.tsv").write_text(SHEET_B)
    both = invoke_ok(
        tmp_path, "--ecos", "none", "--viz_custom", f"{tmp_path / 'a.tsv'},{tmp_path / 'b.tsv'}",
        mapping=MAPPING,
    )  # fmt: skip
    for sheet in ("a", "b"):
        alone = invoke_ok(
            tmp_path, "--ecos", "none", "--viz_custom", str(tmp_path / f"{sheet}.tsv"),
            mapping=MAPPING, out=f"out_{sheet}",
        )  # fmt: skip
        assert sorted(p.name for p in (both / sheet).iterdir()) == sorted(
            p.name for p in (alone / sheet).iterdir()
        )
        assert tables(both / sheet).keys() == {
            "df_genome.tsv",
            "df_sample.tsv",
            "df_vertical_genome.tsv",
            "df_vertical_sample.tsv",
        }
        assert same_tables(both / sheet, alone / sheet), sheet
    # Each sheet used its own Mcr: only bin1 has Cyc2_repCluster1.
    genome = pl.read_csv(both / "b" / "df_vertical_genome.tsv", separator="\t")
    assert value(genome, "present", name="Mcr", genome="bin2") is False
    genome = pl.read_csv(both / "a" / "df_vertical_genome.tsv", separator="\t")
    assert value(genome, "present", name="Mcr", genome="bin2") is True


def test_ecos_choose_built_in_sheets(tmp_path):
    out = invoke_ok(tmp_path, "--ecos", "gut,ag,default", empty_common=False)
    assert sorted(p.name for p in out.iterdir()) == ["ag", "default", "gut"]
    assert same_tables(out / "default", out / "gut")
    assert not same_tables(out / "ag", out / "gut")


def test_no_ecos_is_default_which_is_gut(tmp_path):
    out = invoke_ok(tmp_path, empty_common=False)
    assert [p.name for p in out.iterdir()] == ["default"]
    gut = invoke_ok(tmp_path, "--ecos", "gut", empty_common=False, out="gut")
    assert same_tables(out / "default", gut / "gut")


def test_sheets_sharing_a_folder_are_an_error(tmp_path):
    (tmp_path / "gut.tsv").write_text(RULES)
    result = invoke(
        tmp_path, "--ecos", "gut", "--viz_custom", str(tmp_path / "gut.tsv")
    )
    assert result.exit_code != 0
    assert "gut" in result.output
    assert not (tmp_path / "out").exists()


def test_needed_features_go_in_each_sheets_folder(tmp_path):
    (tmp_path / "a.tsv").write_text(SHEET_A)
    (tmp_path / "b.tsv").write_text(SHEET_B)
    out = invoke_ok(
        tmp_path,
        "-sf",
        "--ecos",
        "none",
        "--viz_custom",
        f"{tmp_path / 'a.tsv'},{tmp_path / 'b.tsv'}",
    )
    assert set((out / "a" / "needed_features.txt").read_text().split()) == {
        "K00399",
        "Cyc2_repCluster1",
    }
    assert (out / "b" / "needed_features.txt").read_text().split() == [
        "Cyc2_repCluster1"
    ]


def test_genome_tables_carry_taxonomy(tmp_path):
    # Completeness and Contamination go to the dashboard, not the tables.
    annotations = (
        "query_id\tinput_fasta\tkegg_id\ttaxonomy\tCompleteness\tContamination\n"
        "bin1_1\tbin1\tK00399\td__Archaea;p__Halo\t90.5\t1.2\n"
        "bin2_1\tbin2\t\td__Bacteria;p__P\t80\t3\n"
    )
    genome = read(run(tmp_path, annotations=annotations, mapping=None), "genome")
    assert value(genome, "taxonomy", name="Mcr", genome="bin1") == "d__Archaea;p__Halo"
    assert value(genome, "present", name="Mcr", genome="bin2") is False


def test_custom_sheets_add_to_the_default(tmp_path):
    (tmp_path / "my_rules.tsv").write_text(SHEET_B)
    out = invoke_ok(
        tmp_path, "--viz_custom", str(tmp_path / "my_rules.tsv"), empty_common=False
    )
    assert sorted(p.name for p in out.iterdir()) == ["default", "my_rules"]


def test_sheet_without_group_column_errors_before_any_output(tmp_path):
    (tmp_path / "a.tsv").write_text(SHEET_A)
    (tmp_path / "b.tsv").write_text("name\trule\nMcr\tK00399\n")
    result = invoke(
        tmp_path, "--ecos", "none", "--viz_custom",
        f"{tmp_path / 'a.tsv'},{tmp_path / 'b.tsv'}",
    )  # fmt: skip
    assert "b.tsv" in str(result.exception)
    assert "Group" in str(result.exception)
    assert not (tmp_path / "out").exists()


@pytest.fixture
def served(monkeypatch):
    """Record `pn.serve` calls instead of starting a server."""
    from dram_rules import make_product

    calls = []
    monkeypatch.setattr(
        make_product.pn, "serve", lambda app, port: calls.append((app(), port))
    )
    return calls


def test_dashboard_serves_the_one_sheet_and_still_writes(tmp_path, served):
    out = invoke_ok(
        tmp_path, "--ecos", "gut", "--dashboard", "--port", "5123", empty_common=False
    )
    [(dashboard, port)] = served
    assert port == 5123
    assert dashboard.title == "DRAM GUT VISUALIZATION"
    assert list((out / "gut").glob("*.html"))
    assert (out / "gut" / "df_vertical_genome.tsv").exists()


def test_dashboard_with_two_sheets_errors_before_evaluating(tmp_path, served):
    result = invoke(tmp_path, "--ecos", "gut,ag", "--dashboard", empty_common=False)
    assert result.exit_code != 0
    assert "--dashboard" in result.output
    assert not served
    assert not (tmp_path / "out").exists()


def test_no_save_serves_and_writes_nothing(tmp_path, served, monkeypatch):
    (tmp_path / "anno.tsv").write_text(ANNOTATIONS)
    monkeypatch.chdir(tmp_path)
    result = CliRunner().invoke(
        main,
        ["-a", "anno.tsv", "--viz", "--ecos", "gut", "--dashboard", "--no_save"],
    )
    assert result.exit_code == 0, result.output + repr(result.exception)
    assert len(served) == 1
    assert sorted(p.name for p in tmp_path.iterdir()) == ["anno.tsv"]


def test_no_save_without_dashboard_errors(tmp_path):
    result = invoke(tmp_path, "--ecos", "gut", "--no_save", empty_common=False)
    assert result.exit_code != 0
    assert "--no_save" in result.output


def test_output_dir_is_required_without_no_save(tmp_path):
    (tmp_path / "anno.tsv").write_text(ANNOTATIONS)
    result = CliRunner().invoke(main, ["-a", str(tmp_path / "anno.tsv"), "--viz"])
    assert result.exit_code != 0
    assert "--output_dir" in result.output


def test_traits_match_visualizing_the_traits_sheet(tmp_path):
    out = invoke_ok(tmp_path, step="--traits", empty_common=False)
    traits = tmp_path / "out" / "traits"
    assert not out.exists()
    assert (traits / "df_vertical_genome.tsv").exists()
    assert list(traits.glob("*.html"))
    viz = invoke_ok(
        tmp_path, "--ecos", "none", "--viz_custom", str(TRAITS_SHEET),
        empty_common=False, out="viz",
    )  # fmt: skip
    assert same_tables(traits, viz / "traits_rules")


def test_traits_custom_adds_sheets_beside_default(tmp_path):
    (tmp_path / "extra.tsv").write_text(SHEET_B)
    invoke_ok(
        tmp_path, "--ecos", "none", "--traits", "--traits_custom", str(tmp_path / "extra.tsv"),
        empty_common=False,
    )  # fmt: skip
    traits = tmp_path / "out" / "traits"
    assert sorted(p.name for p in traits.iterdir()) == ["default", "extra"]
    assert (traits / "default" / "df_vertical_genome.tsv").exists()
    assert (traits / "extra" / "df_vertical_genome.tsv").exists()


def test_viz_and_traits_in_one_run(tmp_path):
    out = invoke_ok(tmp_path, "--ecos", "gut", "--traits", empty_common=False)
    assert [p.name for p in out.iterdir()] == ["gut"]
    assert (tmp_path / "out" / "traits" / "df_vertical_genome.tsv").exists()


def test_dashboard_with_viz_and_traits_errors(tmp_path, served):
    result = invoke(
        tmp_path, "--ecos", "gut", "--traits", "--dashboard", empty_common=False
    )
    assert result.exit_code != 0
    assert not served


def test_port_without_dashboard_errors(tmp_path):
    result = invoke(tmp_path, "--ecos", "gut", "--port", "5123", empty_common=False)
    assert result.exit_code != 0
    assert "--port" in result.output
    assert not (tmp_path / "out").exists()


@pytest.fixture
def annotation_sources(monkeypatch):
    """Where each genome- or sample-level evaluation got its annotations:
    `(level, "file")` for a scan, `(level, "memory")` otherwise."""
    from dram_rules import make_product

    sources = []
    evaluate = make_product.evaluate

    def spy(compiled, hits, annotations, level, *args, **kwargs):
        plan = annotations.lazy().explain().upper()
        sources.append((level, "file" if "SCAN" in plan else "memory"))
        return evaluate(compiled, hits, annotations, level, *args, **kwargs)

    monkeypatch.setattr(make_product, "evaluate", spy)
    return sources


ROWS_SHEET = "name\tgroup\trule\nAny\tG1\tcount_rows(>=, 1)\nMcr\tG1\tK00399\n"


def test_sheets_reading_rows_share_one_copy_in_memory(tmp_path, annotation_sources):
    (tmp_path / "a.tsv").write_text(ROWS_SHEET)
    (tmp_path / "b.tsv").write_text(SHEET_B)
    invoke_ok(
        tmp_path, "--ecos", "none", "--viz_custom", f"{tmp_path / 'a.tsv'},{tmp_path / 'b.tsv'}",
    )  # fmt: skip
    assert annotation_sources == [("genome", "memory")] * 2


def test_sample_level_reads_rows_from_memory_too(tmp_path, annotation_sources):
    (tmp_path / "a.tsv").write_text(ROWS_SHEET)
    (tmp_path / "b.tsv").write_text(SHEET_B)
    invoke_ok(
        tmp_path, "--ecos", "none", "--viz_custom", f"{tmp_path / 'a.tsv'},{tmp_path / 'b.tsv'}",
        mapping=MAPPING,
    )  # fmt: skip
    assert (
        sorted(annotation_sources)
        == [("genome", "memory")] * 2 + [("sample", "memory")] * 2
    )


def test_annotations_stay_on_disk_when_no_rule_reads_rows(tmp_path, annotation_sources):
    run(tmp_path)  # with the mapping file
    assert sorted(annotation_sources) == [("genome", "file"), ("sample", "file")]
