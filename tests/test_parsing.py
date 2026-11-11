from pathlib import Path

import polars as pl
import pytest

from dram_rules.rule_parser.rules import (
    CompiledRules,
    ID_PARSERS,
    Name,
    RuleError,
    evaluate_rules_on_anno,
    load_rules,
    parse_ids,
    prepare_present_map_df,
)

DATA_DIR = Path(__file__).parent / "data"

BASIC_NEEDS = {
    "K00399",
    "K00400",
    "K00401",
    "K00402",
    "K03421",
    "K03422",
    "K02591",
    "K02586",
    "K02588",
}

COMPOUND_NEEDS = {
    "K00400",
    "K00399",
    "K03422",
    "K02591",
    "K13380",
    "K02588",
    "K00401",
    "K00402",
    "K00333",
    "K03421",
    "K02586",
    "K14126",
}

# What each rule in the test sheets means with its aliases written out in full.
NOT_B4 = "not([K02591 | K02586 | K02588])"
BASIC_EXPANDED = {
    "A": f"[K00399 | K00400 | K00401 | K00402 | K03421 | K03422] & {NOT_B4}",
}

COMPOUND_EXPANDED = {
    "A": BASIC_EXPANDED["A"],
    "B": "percent(50, [K00399, [K00401 & K00402 & [K00333 | K03421] & [K03422 | K13380]],"
    " K00400, K00401, K00402, K03421, K03422]) | [K02591 | K02586 | K02588]",
    "C": 'not(filter_contains(kegg_description, "nitrate reductase"))'
    " -> filter_compare(heme_regulatory_motif_count, >=, 2) -> count_rows(>=, 2)",
    "D": "percent(50, [K00399, K14126, K02588])",
}


def compile_text(rule: str):
    """Compile one rule written without aliases."""
    return CompiledRules.from_rules(
        rules=pl.DataFrame({"name": ["r"], "rule": [rule]}).lazy()
    ).rules["r"]


def test_dbcan_sub_ec_parsing():
    annotations = pl.DataFrame(
        {
            "dbcan_sub_ec": [
                "2.4.1.25:4;3.2.1.68:75;3.2.1.-:1",
                "1.2.3.4:2",
            ]
        }
    )

    parsed = annotations.select(parse_ids("dbcan_sub_ec")).to_series().to_list()

    # Parsers only find the ids; the "EC:" prefix is normalized at match time.
    assert parsed == [
        ["2.4.1.25", "3.2.1.68", "3.2.1.-"],
        ["1.2.3.4"],
    ]


@pytest.mark.parametrize(
    "column, annotation, expected",
    [
        # Versioned Pfam ids lose their version; several hits per cell.
        (
            "pfam_hits",
            "Pkinase [PF00069.28]; PK_Tyr [PF07714.20]",
            ["PF00069", "PF07714"],
        ),
        ("pfam_id", "PF00001.23", ["PF00001"]),
        # CAZy ECs written either "(EC x)" or "EC:x".
        (
            "cazy_hits",
            "GH5 (EC 3.2.1.4); GH9 (EC 3.2.1.-)",
            ["EC:3.2.1.4", "EC:3.2.1.-"],
        ),
        ("cazy_hits", "GH5 [EC:3.2.1.4]", ["EC:3.2.1.4"]),
        ("methyl_id", "M1 some desc, M2 other", ["M1", "M2"]),
        ("unknown_db_id", "abc123", ["abc123"]),
        ("unknown_db_acc", "XYZ0001.4", ["XYZ0001"]),
    ],
)
def test_parse_ids(column, annotation, expected):
    annotations = pl.DataFrame({column: [annotation]})
    assert annotations.select(parse_ids(column)).item().to_list() == expected


@pytest.mark.parametrize(
    "column, annotation, feature",
    [
        ("pfam_hits", "Pkinase [PF00069.28]", "PF00069"),
        ("cazy_hits", "GH5 (EC 3.2.1.4)", "EC:3.2.1.4"),
        ("cazy_hits", "GH5 (EC 3.2.1.4)", "3.2.1.4"),
    ],
)
def test_known_formats_match_rule_ids(column, annotation, feature):
    annotations = pl.DataFrame({"query_id": ["gene1"], column: [annotation]})

    parsed = prepare_present_map_df(
        annotations,
        count_col="query_id",
        besthit_cols=list(ID_PARSERS),
        needed_features={feature},
    )

    assert parsed.rows() == [("gene1", feature)]


@pytest.mark.parametrize(
    "column, annotation, feature",
    [
        ("custom_ec", "3.2.1.212", "EC:3.2.1.212"),
        ("custom_ec", "EC:3.2.1.212", "3.2.1.212"),
        ("custom_ec", "ec:3.2.1.-", "EC:3.2.1.-"),
        ("custom_ec", "GH107;3.2.1.212", "GH107;EC:3.2.1.212"),
        ("dbcan_sub_ec", "3.2.1.212:4", "EC:3.2.1.212"),
        ("dbcan_sub_ec", "3.2.1.212:4", "3.2.1.212"),
    ],
)
def test_ec_prefix_fallback(column, annotation, feature):
    annotations = pl.DataFrame({"query_id": ["gene1"], column: [annotation]})

    parsed = prepare_present_map_df(
        annotations,
        count_col="query_id",
        besthit_cols=list(ID_PARSERS),
        needed_features={feature},
    )

    assert parsed.rows() == [("gene1", feature)]
    assert annotations[column][0] == annotation


def test_dbcan_compound_ec_fallback():
    annotations = pl.DataFrame(
        {
            "query_id": ["gene1", "gene2", "gene3", "gene4"],
            "dbcan_id": ["GH107", "GH168", "GH107", "GH107"],
            "dbcan_sub_ec": ["3.2.1.212:4", "3.2.1.-:1", "-", None],
            "custom_ec": ["EC:3.2.1.212", None, None, None],
        }
    )

    parsed = prepare_present_map_df(
        annotations,
        count_col="query_id",
        besthit_cols=list(ID_PARSERS),
        needed_features={
            "GH107;EC:3.2.1.212",
            "GH107;3.2.1.212",
            "GH168;EC:3.2.1.-",
            "EC:3.2.1.212",
            "3.2.1.212",
        },
    )

    assert set(parsed.iter_rows()) == {
        ("gene1", "GH107;EC:3.2.1.212"),
        ("gene1", "GH107;3.2.1.212"),
        ("gene2", "GH168;EC:3.2.1.-"),
        ("gene1", "EC:3.2.1.212"),
        ("gene1", "3.2.1.212"),
    }
    assert parsed.height == 5


def test_ec_prefix_fallback_preserves_other_identifiers():
    annotations = pl.DataFrame(
        {
            "query_id": ["gene1", "gene2", "gene3"],
            "custom_id": ["SPEC:123", "K00001", None],
        }
    )

    parsed = prepare_present_map_df(
        annotations,
        count_col="query_id",
        besthit_cols=list(ID_PARSERS),
        needed_features={"SPEC:123", "SP123", "K00001"},
    )

    assert set(parsed.iter_rows()) == {("gene1", "SPEC:123"), ("gene2", "K00001")}


def test_custom_accession_columns():
    annotations = pl.DataFrame(
        {
            "query_id": ["gene1", "gene2", "gene3", "gene4", "gene5"],
            "custom_hmm_acc": [
                "PF00001.23",
                "PF00002",
                "custom.domain.12",
                "custom.domain",
                None,
            ],
            "custom_id": ["ID.2", None, None, None, None],
        }
    )

    parsed = prepare_present_map_df(
        annotations,
        count_col="query_id",
        besthit_cols=list(ID_PARSERS),
        needed_features={"PF00001", "PF00002", "CUSTOM.DOMAIN", "ID.2"},
    )

    assert set(parsed.iter_rows()) == {
        ("gene1", "PF00001"),
        ("gene1", "ID.2"),
        ("gene2", "PF00002"),
        ("gene3", "CUSTOM.DOMAIN"),
        ("gene4", "CUSTOM.DOMAIN"),
    }
    assert annotations["custom_hmm_acc"][0] == "PF00001.23"


@pytest.mark.parametrize(
    "square, round_",
    [
        ("K12345 & [K54321 | K98765]", "K12345 & (K54321 | K98765)"),
        ("[K12345, K12346, K12347]", "(K12345, K12346, K12347)"),
        ("[K12345]", "(K12345)"),
        ("percent(50, [K12345, K12346])", "percent(50, (K12345, K12346))"),
        ("not([K12345 | K12346])", "not((K12345 | K12346))"),
        ("[K12345, [K12346 | K12347]]", "(K12345, [K12346 | K12347])"),
        ("[K12345, [K12346 | K12347]]", "[K12345, (K12346 | K12347)]"),
        ("percent(50, [K12345, K12346])", "percent(50, @CalvinCycle)"),
    ],
)
def test_group_delimiters(square, round_):
    compiled = CompiledRules.from_rules(
        rules=pl.DataFrame(
            {
                "name": ["square", "round", ""],
                "alias": ["", "", "CalvinCycle"],
                "rule": [square, round_, "(K12345, K12346)"],
            }
        ).lazy()
    )

    assert compiled.rules["square"] == compiled.rules["round"]
    assert compiled.features_by_rules["square"] == compiled.features_by_rules["round"]


@pytest.mark.parametrize(
    "expression",
    ["(K12345 | K12346]", "[K12345, K12346)", "(K12345", "K12345)"],
)
def test_mismatched_group_delimiters(expression):
    with pytest.raises(RuleError):
        CompiledRules.from_rules(
            rules=pl.DataFrame({"name": ["invalid"], "rule": [expression]}).lazy()
        )


def test_alias_only_rules_preserve_nulls():
    definitions, rules, frame = load_rules(
        rules=pl.DataFrame(
            {
                "name": [None],
                "alias": ["Example"],
                "rule": ["K12345"],
                "description": [None],
            }
        ).lazy()
    )

    assert definitions == {"Example": Name("K12345")}
    assert rules == {}
    assert frame.row(0, named=True)["description"] is None


def test_rules_with_all_null_aliases():
    definitions, rules, _ = load_rules(
        rules=pl.DataFrame(
            {"name": ["Example"], "alias": [None], "rule": ["K12345"]}
        ).lazy()
    )

    assert definitions == {}
    assert rules == {"Example": Name("K12345")}


def test_backup_group_col_fills_nulls_without_mutating_defaults():
    from dram_rules.definitions import BACKUP_GROUPBY_COLUMNS

    before = list(BACKUP_GROUPBY_COLUMNS)
    rules = pl.DataFrame(
        {
            "name": ["A", "B"],
            "rule": ["K1", "K2"],
            "topic_ecosystem": ["eco", None],
            "module": ["mod_a", "mod_b"],
        }
    ).lazy()

    for _ in range(2):
        compiled = CompiledRules.from_rules(rules=rules)
        # topic_ecosystem is the group col; its null is filled from module.
        assert compiled.df["topic_ecosystem"].to_list() == ["eco", "mod_b"]
    assert BACKUP_GROUPBY_COLUMNS == before


def test_basic_compiled_rules():
    compiled = CompiledRules.from_rules(rules_path=DATA_DIR / "basic_rules.tsv")
    assert compiled.needed_features == BASIC_NEEDS

    assert compiled.rules == {k: compile_text(v) for k, v in BASIC_EXPANDED.items()}


def test_basic_rules_pass_annotation():

    df = evaluate_rules_on_anno(
        rules_path=DATA_DIR / "basic_rules.tsv",
        annotations_path=DATA_DIR / "basic_anno.tsv",
        count_col="query_id",
        group_col="input_fasta",
    )

    true_df = pl.DataFrame(
        {
            "input_fasta": ["bin1", "bin2", "bin3", "bin4", "bin5", "bin6"],
            "A": [True, False, False, False, False, False],
        }
    )
    assert df.equals(true_df)


def test_compoun_compiled_rules():
    compiled = CompiledRules.from_rules(rules_path=DATA_DIR / "compound_rules.tsv")
    assert compiled.needed_features == COMPOUND_NEEDS

    assert compiled.rules == {k: compile_text(v) for k, v in COMPOUND_EXPANDED.items()}


def test_compound_rules_pass_annotation():

    df = evaluate_rules_on_anno(
        rules_path=DATA_DIR / "compound_rules.tsv",
        annotations_path=DATA_DIR / "basic_anno.tsv",
        count_col="query_id",
        group_col="input_fasta",
    )

    true_df = pl.DataFrame(
        {
            "input_fasta": ["bin1", "bin2", "bin3", "bin4", "bin5", "bin6"],
            "A": [True, False, False, False, False, False],
            "B": [False, False, False, False, False, True],
            "C": [True, False, False, False, False, False],
            "D": [True, False, False, False, False, False],
        }
    )
    assert df.equals(true_df)
