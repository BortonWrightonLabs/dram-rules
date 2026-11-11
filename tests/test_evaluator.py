import gc
import textwrap
import weakref
from pathlib import Path

import numpy as np
import polars as pl
import pytest

from dram_rules.rule_parser.rules import (
    CompiledRules,
    Evaluator,
    RuleError,
    evaluate_cycles,
    evaluate_rules,
    evaluate_rules_on_anno,
    get_val_cols,
)
from dram_rules.rules_api import FunctionRegistry, RuleFunctionError, default_functions

SAMPLES = ["g1", "g2"]


def compile_rule(rule: str, **kwargs) -> CompiledRules:
    return CompiledRules.from_rules(
        rules=pl.DataFrame({"name": ["r"], "rule": [rule]}).lazy(), **kwargs
    )


def evaluate(rule, present_map=None, annotations=None, samples=SAMPLES, **kwargs):
    if annotations is None:
        annotations = pl.DataFrame({"genome": samples})
    out = evaluate_rules(
        compile_rule(rule, **kwargs),
        samples,
        present_map or {},
        annotations=annotations,
        sample_col="genome",
    )
    return out["r"].to_list()


ANNOTATIONS = pl.DataFrame(
    {
        "genome": ["g1", "g1", "g1", "g2", "g2"],
        "desc": ["sulfur a", "sulfur other", "sulfur c", "sulfur b", "sulfur other"],
        "score": [1.0, 1.0, 1.0, 1.0, 1.0],
        "taxonomy": ["d__Bacteria", "d__Bacteria", "d__Bacteria", "d__Archaea", "d__Archaea"],
    }
)


# --- pipes and row filters --------------------------------------------------


def test_not_of_a_filter_after_another_filter():
    # sulfur AND NOT other: g1 keeps 2 rows, g2 keeps 1.
    rule = (
        "filter_contains(desc, sulfur) -> not(filter_contains(desc, other))"
        " -> count_rows(>=, 2)"
    )
    assert evaluate(rule, annotations=ANNOTATIONS) == [True, False]


def test_not_of_a_filter_at_start_of_pipe():
    # The pattern used by hemeCount in rules_common.tsv.
    rule = "not(filter_contains(desc, other)) -> count_rows(>=, 2)"
    assert evaluate(rule, annotations=ANNOTATIONS) == [True, False]


def test_filters_combined_with_and_before_pipe():
    rule = 'filter_contains(desc, "sulfur") & not(filter_contains(desc, other)) -> count_rows(>=, 2)'
    assert evaluate(rule, annotations=ANNOTATIONS) == [True, False]


def test_filters_combined_with_or_before_pipe():
    rule = 'filter_contains(desc, " a") | filter_contains(desc, " b") -> count_rows(==, 1)'
    assert evaluate(rule, annotations=ANNOTATIONS) == [True, True]


def test_tax_respects_pipe_filters():
    # Only g2 has a "sulfur b" row, and g2 is Archaea.
    rule = 'filter_contains(desc, "sulfur b") -> tax(Archaea)'
    assert evaluate(rule, annotations=ANNOTATIONS) == [False, True]


def test_count_rows_sees_samples_with_no_rows_as_zero():
    # No row of g1 or g2 contains "zzz": both have 0 rows, and 0 < 1.
    rule = 'filter_contains(desc, "zzz") -> count_rows(<, 1)'
    assert evaluate(rule, annotations=ANNOTATIONS) == [True, True]


def test_sum_column():
    annotations = ANNOTATIONS.with_columns(score=pl.Series([1.0, 2.0, 3.0, 10.0, 0.5]))
    assert evaluate("sum_column(score, >=, 7)", annotations=annotations) == [False, True]


def test_filter_compare():
    annotations = ANNOTATIONS.with_columns(score=pl.Series([1.0, 2.0, 3.0, 10.0, 0.5]))
    rule = "filter_compare(score, >, 1.5) -> count_rows(>=, 2)"
    assert evaluate(rule, annotations=annotations) == [True, False]


def test_matching_is_plain_text_unless_regex():
    annotations = pl.DataFrame({"genome": ["g1", "g2"], "desc": ["EC:1.2.3.4", "EC:1x2x3x4"]})
    plain = 'filter_contains(desc, "1.2.3.4") -> count_rows(>=, 1)'
    regex = 'filter_contains(desc, "1.2.3.4", regex=true) -> count_rows(>=, 1)'
    assert evaluate(plain, annotations=annotations) == [True, False]
    assert evaluate(regex, annotations=annotations) == [True, True]


# --- gene groups -------------------------------------------------------------


@pytest.mark.parametrize("mode", ["COUNTS", "counts"])
def test_at_least_counts_adds_up_hits(mode):
    # K1 hit 3 times in g1; K2 absent everywhere.
    present_map = {"K1": np.array([3.0, 1.0]), "K2": np.array([0.0, 0.0])}
    assert evaluate(f"at_least(2, [K1 | K2], mode={mode})", present_map) == [True, False]


@pytest.mark.parametrize("rule", ["at_least(2, [K1 | K2])", "at_least(2, [K1, K2], mode=presence)"])
def test_at_least_presence_counts_distinct_parts(rule):
    present_map = {"K1": np.array([3.0, 1.0]), "K2": np.array([0.0, 1.0])}
    assert evaluate(rule, present_map) == [False, True]


def test_percent():
    present_map = {"K1": np.array([1.0, 1.0]), "K2": np.array([0.0, 1.0])}
    assert evaluate("percent(100, [K1, K2])", present_map) == [False, True]
    assert evaluate("percent(50, [K1 | K2])", present_map) == [True, True]


# --- coverage (rule_options) ---------------------------------------------------


def coverage_output(rule: str, mode: str, present_map) -> pl.DataFrame:
    compiled = CompiledRules.from_rules(
        rules=pl.DataFrame(
            {
                "name": ["r"],
                "rule": [rule],
                "rule_options": [f'{{"coverage":"{mode}"}}'],
                "group": ["g"],
            }
        ).lazy()
    )
    dfs = evaluate_cycles(
        compiled,
        SAMPLES,
        present_map,
        annotations=pl.DataFrame({"genome": SAMPLES}),
        sample_col="genome",
    )
    return dfs["g"]


def test_coverage_steps_counts_steps_not_hits():
    # Two steps; K1 has 3 hits in g1 but is still one step.
    out = coverage_output(
        "K1, [K2 | K3]", "steps", {"K1": np.array([3.0, 0.0]), "K2": np.array([0.0, 0.0])}
    )
    assert out["steps"].to_list() == [2, 2]
    assert out["steps_present"].to_list() == [1, 0]
    assert out["coverage_percentage"].to_list() == [0.5, 0.0]


def test_coverage_subunits_counts_top_level_and_parts():
    out = coverage_output(
        "K1 & [K2 | K3] & K4",
        "subunits",
        {"K1": np.array([1.0, 0.0]), "K3": np.array([1.0, 1.0])},
    )
    assert out["subunits"].to_list() == [3, 3]
    assert out["subunits_present"].to_list() == [2, 1]


def test_coverage_rule_still_works_as_present_absent():
    present_map = {"K1": np.array([1.0, 1.0]), "K2": np.array([1.0, 0.0])}
    assert evaluate("K1, K2", present_map) == [True, False]


def test_unknown_coverage_mode_is_an_error():
    with pytest.raises(RuleError, match="coverage must be one of"):
        CompiledRules.from_rules(
            rules=pl.DataFrame(
                {"name": ["r"], "rule": ["K1"], "rule_options": ['{"coverage":"pathway"}']}
            ).lazy()
        )


def test_group_mixing_coverage_settings_is_an_error():
    with pytest.raises(RuleError, match="same rule_options coverage"):
        CompiledRules.from_rules(
            rules=pl.DataFrame(
                {
                    "name": ["a", "b"],
                    "rule": ["K1, K2", "K3"],
                    "rule_options": ['{"coverage":"steps"}', None],
                    "group": ["g", "g"],
                }
            ).lazy()
        )


# --- checks when rules load ------------------------------------------------------


@pytest.mark.parametrize(
    "rule, message",
    [
        ("count_rows(>=, 1) -> count_rows(>=, 1)", "per-sample and can't be used before '->'"),
        ("K1 -> count_rows(>=, 1)", "per-sample and can't be used before '->'"),
        ("filter_contains(desc, x)", "row filter; row filters can only be used before '->'"),
        ("K1 & filter_contains(desc, x)", "mixes row filters and per-sample values"),
        ("filter_contains(desc, x) -> percent(50, [K1])", "doesn't read annotation rows"),
        ("filter_contains(desc, x) -> filter_contains(desc, y)", "can't end a pipe"),
        ("percent(fifty, [K1])", "expects an integer"),
        ("at_least(2, [K1], mode=SOME)", "expects one of PRESENCE, COUNTS"),
        ("count_rows(ge, 1)", "use >="),
        ("at_least(2, [K1], nope=1)", "unknown keyword"),
        ("at_least(2, [K1], COUNTS, 4)", "too many arguments"),
        ("at_least(2)", "missing argument 'group'"),
        ("at_lest(2, [K1])", "Did you mean at_least"),
        ("at_least(n=2, [K1])", "positional argument"),
    ],
)
def test_rule_errors_name_the_problem(rule, message):
    with pytest.raises(RuleError, match="Rule 'r'|in rule:") as error:
        compile_rule(rule)
    assert message in str(error.value)


def compile_sheet(names, aliases, rules):
    return CompiledRules.from_rules(
        rules=pl.DataFrame({"name": names, "alias": aliases, "rule": rules}).lazy()
    )


def test_errors_quote_aliases_as_written():
    # @F is a row filter, so it can't be a whole rule. The message should say
    # "@F" (what the sheet says), not the expanded filter_contains(...).
    with pytest.raises(RuleError) as error:
        compile_sheet(["r", None], [None, "F"], ["@F", "filter_contains(desc, x)"])
    assert str(error.value).startswith("Rule 'r': @F is a row filter")


def test_errors_inside_aliases_name_the_alias_chain():
    with pytest.raises(RuleError) as error:
        compile_sheet(
            ["r", None, None],
            [None, "Outer", "Inner"],
            ["@Outer", "not(@Inner)", "percent(fifty, [K1])"],
        )
    assert str(error.value).startswith(
        "Alias '@Outer': in @Inner: percent() argument 'n' expects an integer, got fifty"
    )


def test_gene_ids_are_collected_from_sub_rules_only():
    compiled = compile_rule("at_least(2, [K1 | K2], mode=COUNTS) & percent(50, [K3, K4])")
    assert compiled.features_by_rules["r"] == {"K1", "K2", "K3", "K4"}


# --- user functions ------------------------------------------------------------


def write_functions(tmp_path: Path, code: str) -> Path:
    path = tmp_path / "my_functions.py"
    path.write_text(textwrap.dedent(code))
    return path


def test_user_functions_file(tmp_path):
    path = write_functions(
        tmp_path,
        """
        import polars as pl
        from dram_rules.rules_api import Column, Group, Op, PerSample, RowFilter, Rows, rule_function

        @rule_function
        def majority(group: Group) -> PerSample:
            return group.present.mean(axis=1) > 0.5

        @rule_function
        def filter_high(col: Column, value: float = 1.5) -> RowFilter:
            return pl.col(col) > value

        @rule_function
        def mean_column(rows: Rows, col: Column, op: Op, value: float) -> PerSample:
            return op(rows.per_sample(pl.col(col).mean(), fill=0.0), value)
        """,
    )
    present_map = {"K1": np.array([1.0, 1.0]), "K2": np.array([1.0, 0.0]), "K3": np.array([0.0, 0.0])}
    assert evaluate("majority([K1, K2, K3])", present_map, function_files=[path]) == [True, False]

    annotations = ANNOTATIONS.with_columns(score=pl.Series([1.0, 2.0, 3.0, 10.0, 0.5]))
    rule = "filter_high(score) -> mean_column(score, >=, 3)"
    assert evaluate(rule, annotations=annotations, function_files=[path]) == [False, True]

    # Loading a file registers its functions for that run only.
    assert "majority" not in default_functions
    with pytest.raises(RuleError, match="Unknown function 'majority'"):
        compile_rule("majority([K1, K2, K3])")


def test_user_function_clashing_with_builtin(tmp_path):
    path = write_functions(
        tmp_path,
        """
        from dram_rules.rules_api import Group, PerSample, rule_function

        @rule_function
        def percent(group: Group) -> PerSample:
            return group.present.any(axis=1)
        """,
    )
    with pytest.raises(RuleFunctionError, match="already defined"):
        compile_rule("K1", function_files=[path])


@pytest.mark.parametrize(
    "code, message",
    [
        (
            """
            import polars as pl
            from dram_rules.rules_api import RowFilter, rule_function

            @rule_function
            def high(col: str) -> RowFilter:
                return pl.col(col) > 1
            """,
            "rename it to 'filter_high'",
        ),
        (
            """
            from dram_rules.rules_api import PerSample, rule_function

            @rule_function
            def filter_count(n: int) -> PerSample:
                return n
            """,
            "must be named 'filter_...'",
        ),
        (
            """
            from dram_rules.rules_api import PerSample, rule_function

            @rule_function
            def untyped(n) -> PerSample:
                return n
            """,
            "needs a type hint",
        ),
        ("x = 1\n", "No @rule_function functions"),
    ],
)
def test_bad_user_functions_fail_when_loaded(tmp_path, code, message):
    path = write_functions(tmp_path, code)
    with pytest.raises(RuleFunctionError, match=message):
        compile_rule("K1", function_files=[path])


def test_user_function_cannot_change_shared_values(tmp_path):
    path = write_functions(
        tmp_path,
        """
        from dram_rules.rules_api import PerSample, rule_function

        @rule_function
        def clip_in_place(x: PerSample) -> PerSample:
            x[x > 1] = 1
            return x > 0
        """,
    )
    present_map = {"K1": np.array([3.0, 0.0])}
    with pytest.raises(RuleError, match="read-only"):
        evaluate("clip_in_place(K1)", present_map, function_files=[path])
    assert present_map["K1"].tolist() == [3.0, 0.0]


@pytest.mark.parametrize(
    "code, message",
    [
        (
            """
            from dram_rules.rules_api import PerSample, RowFilter, rule_function

            @rule_function
            def either(x: PerSample) -> PerSample | RowFilter:
                return x
            """,
            "return type hint must be PerSample or RowFilter",
        ),
        (
            """
            from dram_rules.rules_api import PerSample, rule_function

            @rule_function(name="not")
            def my_not(x: PerSample) -> PerSample:
                return ~x
            """,
            "'not' is part of the rule language",
        ),
    ],
)
def test_rejected_function_definitions(tmp_path, code, message):
    path = write_functions(tmp_path, code)
    with pytest.raises(RuleFunctionError, match=message):
        compile_rule("K1", function_files=[path])


@pytest.mark.parametrize(
    "body, rule, message",
    [
        (
            "def high(col: Column) -> PerSample:\n    return pl.col(col) > 1",
            "high(score)",
            (
                "returned a Polars expression, which is a row filter. If it is a filter,"
                " give it the return hint RowFilter and name it 'filter_high'"
            ),
        ),
        (
            "def filter_high(col: Column) -> RowFilter:\n    return np.array([True, False])",
            "filter_high(score) -> count_rows(>=, 1)",
            "returned ndarray of shape (2,). A row filter is a Polars expression",
        ),
        (
            "def anyof(group: Group) -> PerSample:\n    return group.present.any()",
            "anyof([K1, K2])",
            "returned a single value (np.True_); expected one value per sample (an array of length 2)",
        ),
        (
            "def words(group: Group) -> PerSample:\n    return np.array(['yes', 'no'])",
            "words([K1, K2])",
            "returned values of type <U3; expected numbers or True/False",
        ),
    ],
)
def test_wrong_results_are_explained(tmp_path, body, rule, message):
    path = write_functions(
        tmp_path,
        "import numpy as np\nimport polars as pl\n"
        "from dram_rules.rules_api import Column, Group, PerSample, RowFilter, rule_function\n"
        "@rule_function\n" + body + "\n",
    )
    annotations = pl.DataFrame({"genome": SAMPLES, "score": [2.0, 0.0]})
    with pytest.raises(RuleError) as error:
        evaluate(rule, {"K1": np.array([1.0, 0.0])}, annotations=annotations, function_files=[path])
    assert str(error.value).startswith("Rule 'r': ")
    assert "my_functions.py:4" in str(error.value)
    assert message in str(error.value)


def test_not_takes_the_kind_of_its_argument():
    present_map = {"K1": np.array([1.0, 0.0])}
    assert evaluate("not(K1)", present_map) == [False, True]
    rule = "not(filter_contains(desc, other)) -> count_rows(>=, 2)"
    assert evaluate(rule, annotations=ANNOTATIONS) == [True, False]
    with pytest.raises(RuleError, match="exactly one argument"):
        compile_rule("not(K1, K2)")


def test_given_registry_is_used_even_if_empty():
    with pytest.raises(RuleError, match="Unknown function 'percent'"):
        compile_rule("percent(50, [K1])", functions=FunctionRegistry())


def test_user_function_returning_wrong_shape(tmp_path):
    path = write_functions(
        tmp_path,
        """
        import numpy as np
        from dram_rules.rules_api import Group, PerSample, rule_function

        @rule_function
        def broken(group: Group) -> PerSample:
            return np.array([True])
        """,
    )
    with pytest.raises(RuleError, match="expected one value per sample"):
        evaluate("broken([K1])", function_files=[path])


# --- missing columns and rule_options ------------------------------------------


def test_sample_without_annotation_rows():
    annotations = pl.DataFrame({"genome": ["g1"], "score": [5.0], "taxonomy": ["d__Bacteria"]})
    assert evaluate("sum_column(score, >=, 1)", annotations=annotations) == [True, False]
    assert evaluate("count_rows(>=, 1)", annotations=annotations) == [True, False]
    assert evaluate("tax(Bacteria)", annotations=annotations) == [True, False]


def test_tax_without_taxonomy_column_is_false():
    assert evaluate("tax(Bacteria)", annotations=ANNOTATIONS.drop("taxonomy")) == [False, False]


@pytest.mark.parametrize(
    "rule",
    [
        "filter_contains(nope, x) -> count_rows(>=, 1)",
        "filter_compare(nope, >=, 1) -> count_rows(>=, 1)",
        "sum_column(nope, >=, 1)",
    ],
)
def test_missing_column_raises(rule):
    with pytest.raises(RuleError, match="include_when"):
        evaluate(rule, annotations=ANNOTATIONS)


def test_rule_options_skip_rule_in_evaluate_rules():
    compiled = CompiledRules.from_rules(
        rules=pl.DataFrame(
            {
                "name": ["needs_desc", "always"],
                "rule": ["filter_contains(nope, x) -> count_rows(>=, 1)", "K1"],
                "rule_options": ['{"include_when":{"all_columns_present":["nope"]}}', None],
            }
        ).lazy()
    )
    out = evaluate_rules(
        compiled,
        SAMPLES,
        {"K1": np.array([1.0, 0.0])},
        annotations=ANNOTATIONS,
        sample_col="genome",
    )
    assert out.columns == ["genome", "always"]


def test_invalid_rule_options_json_names_the_rule():
    with pytest.raises(RuleError, match="bad_rule"):
        CompiledRules.from_rules(
            rules=pl.DataFrame(
                {"name": ["bad_rule"], "rule": ["K1"], "rule_options": ["{not json"]}
            ).lazy()
        )


@pytest.mark.parametrize("has_kegg_description", [True, False])
def test_traits_run_with_and_without_kegg_description(has_kegg_description):
    data = Path(__file__).parents[1] / "src/dram_rules/data"
    annotations = pl.read_csv(Path(__file__).parent / "data/basic_anno.tsv", separator="\t")
    if not has_kegg_description:
        annotations = annotations.drop("kegg_description")

    out = evaluate_rules_on_anno(
        count_col="query_id",
        group_col="input_fasta",
        annotations=annotations,
        rules_path=data / "traits_rules.tsv",
        common_rules_path=data / "rules_common.tsv",
    )
    has_columns = has_kegg_description and "heme_regulatory_motif_count" in annotations.columns
    assert ("probable-iron-reduction" in out.columns) == has_columns


# --- misc ---------------------------------------------------------------------------


def test_get_val_cols_raises_when_no_value_column():
    with pytest.raises(RuleError):
        get_val_cols(pl.DataFrame({"genome": ["g1"]}))


def test_evaluator_is_garbage_collected():
    ev = Evaluator(
        samples=SAMPLES,
        present_map={"K1": np.array([1.0, 0.0])},
        sample_col="genome",
        annotations=pl.DataFrame({"genome": SAMPLES}),
    )
    ev.evaluate(compile_rule("K1").rules["r"])
    ref = weakref.ref(ev)
    del ev
    gc.collect()
    assert ref() is None
