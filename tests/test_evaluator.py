import gc
import weakref
from pathlib import Path

import numpy as np
import polars as pl
import pytest

from dram_rules.rule_parser.rules import (
    CompiledRules,
    Evaluator,
    RuleError,
    evaluate_rules,
    evaluate_rules_on_anno,
    get_val_cols,
)

SAMPLES = ["g1", "g2"]


def compile_rule(rule: str, **kwargs) -> CompiledRules:
    return CompiledRules.from_rules(
        rules=pl.DataFrame({"name": ["r"], "rule": [rule]}).lazy(), **kwargs
    )


def evaluate(rule, present_map=None, annotations=None, samples=SAMPLES):
    if annotations is None:
        annotations = pl.DataFrame({"genome": samples})
    out = evaluate_rules(
        compile_rule(rule),
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
    }
)


def test_not_of_a_filter_after_another_filter():
    # sulfur AND NOT other: g1 keeps 2 rows, g2 keeps 1.
    rule = (
        "filter_contains(desc, sulfur) -> not(filter_contains(desc, other))"
        " -> column_count_values(score, ge, 0, ge, 2)"
    )
    assert evaluate(rule, annotations=ANNOTATIONS) == [True, False]


def test_not_of_a_filter_at_start_of_pipe():
    # The pattern used by hemeCount in rules_common.tsv.
    rule = (
        "not(filter_contains(desc, other)) -> column_count_values(score, ge, 0, ge, 2)"
    )
    assert evaluate(rule, annotations=ANNOTATIONS) == [True, False]


@pytest.mark.parametrize("mode", ["COUNTS", "counts"])
def test_at_least_counts_adds_up_hits(mode):
    # K1 hit 3 times in g1; K2 absent everywhere.
    present_map = {"K1": np.array([3.0, 1.0]), "K2": np.array([0.0, 0.0])}
    assert evaluate(f"at_least(2, {mode}, [K1 | K2])", present_map) == [True, False]


@pytest.mark.parametrize("mode", ["PRESENCE", "presence"])
def test_at_least_presence_counts_distinct_parts(mode):
    present_map = {"K1": np.array([3.0, 1.0]), "K2": np.array([0.0, 1.0])}
    assert evaluate(f"at_least(2, {mode}, [K1 | K2])", present_map) == [False, True]


def test_path_steps_counts_steps_not_hits():
    # Two steps; K1 has 3 hits in g1 but is still one step.
    compiled = compile_rule("path_steps(K1, K2)", allow_visualize_functions=True)
    ev = Evaluator(
        samples=SAMPLES,
        present_map={"K1": np.array([3.0, 0.0]), "K2": np.array([0.0, 0.0])},
        sample_col="genome",
        annotations=pl.DataFrame({"genome": SAMPLES}),
    )
    out = ev.eval_bool(compiled.rules["r"])
    assert out["steps_present"].to_list() == [1, 0]
    assert out["coverage_percentage"].to_list() == [0.5, 0.0]


@pytest.mark.parametrize(
    "rule",
    [
        "column_sum_values(score, ge, 1)",
        "column_count_values(score, ge, 0, ge, 1)",
        "tax(Bacteria)",
    ],
)
def test_sample_without_annotation_rows_is_false(rule):
    annotations = pl.DataFrame(
        {"genome": ["g1"], "score": [5.0], "taxonomy": ["d__Bacteria"]}
    )
    assert evaluate(rule, annotations=annotations, samples=["g1", "g2"]) == [
        True,
        False,
    ]


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
    ev.eval_bool(compile_rule("K1").rules["r"])
    ref = weakref.ref(ev)
    del ev
    gc.collect()
    assert ref() is None


@pytest.mark.parametrize(
    "rule",
    [
        "filter_contains(nope, x) -> column_count_values(score, ge, 0, ge, 1)",
        "filter_compare(nope, ge, 1) -> column_count_values(score, ge, 0, ge, 1)",
    ],
)
def test_filter_on_missing_column_raises(rule):
    with pytest.raises(RuleError, match="include_when"):
        evaluate(rule, annotations=ANNOTATIONS)


def test_rule_options_skip_rule_in_evaluate_rules():
    compiled = CompiledRules.from_rules(
        rules=pl.DataFrame(
            {
                "name": ["needs_desc", "always"],
                "rule": [
                    "filter_contains(nope, x) -> column_count_values(score, ge, 0, ge, 1)",
                    "K1",
                ],
                "rule_options": [
                    '{"include_when":{"all_columns_present":["nope"]}}',
                    None,
                ],
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
    annotations = pl.read_csv(
        Path(__file__).parent / "data/basic_anno.tsv", separator="\t"
    )
    if not has_kegg_description:
        annotations = annotations.drop("kegg_description")

    out = evaluate_rules_on_anno(
        count_col="query_id",
        group_col="input_fasta",
        annotations=annotations,
        rules_path=data / "traits_rules.tsv",
        common_rules_path=data / "rules_common.tsv",
    )
    assert ("probable-iron-reduction" in out.columns) == has_kegg_description
