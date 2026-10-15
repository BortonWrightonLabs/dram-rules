from pathlib import Path

import polars as pl

from dram_rules.rule_parser.rules import (
    And,
    Call,
    CompiledRules,
    ID_EXPR_DICT,
    Name,
    Or,
    PipeChain,
    Steps,
    String,
    evaluate_rules_on_anno,
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

BASIC_PARSE = {
    "A": And(
        parts=(
            Or(
                parts=(
                    Name(value="K00399", qualifier=None),
                    Name(value="K00400", qualifier=None),
                    Name(value="K00401", qualifier=None),
                    Name(value="K00402", qualifier=None),
                    Name(value="K03421", qualifier=None),
                    Name(value="K03422", qualifier=None),
                )
            ),
            Call(
                value="not",
                args=(
                    Or(
                        parts=(
                            Name(value="K02591", qualifier=None),
                            Name(value="K02586", qualifier=None),
                            Name(value="K02588", qualifier=None),
                        )
                    ),
                ),
            ),
        )
    )
}

COMPOUND_PARSE = {
    "A": And(
        parts=(
            Or(
                parts=(
                    Name(value="K00399", qualifier=None),
                    Name(value="K00400", qualifier=None),
                    Name(value="K00401", qualifier=None),
                    Name(value="K00402", qualifier=None),
                    Name(value="K03421", qualifier=None),
                    Name(value="K03422", qualifier=None),
                )
            ),
            Call(
                value="not",
                args=(
                    Or(
                        parts=(
                            Name(value="K02591", qualifier=None),
                            Name(value="K02586", qualifier=None),
                            Name(value="K02588", qualifier=None),
                        )
                    ),
                ),
            ),
        )
    ),
    "B": Or(
        parts=(
            Call(
                value="percent",
                args=(
                    Name(value="50", qualifier=None),
                    Steps(
                        parts=(
                            Name(value="K00399", qualifier=None),
                            And(
                                parts=(
                                    Name(value="K00401", qualifier=None),
                                    Name(value="K00402", qualifier=None),
                                    Or(
                                        parts=(
                                            Name(value="K00333", qualifier=None),
                                            Name(value="K03421", qualifier=None),
                                        )
                                    ),
                                    Or(
                                        parts=(
                                            Name(value="K03422", qualifier=None),
                                            Name(value="K13380", qualifier=None),
                                        )
                                    ),
                                )
                            ),
                            Name(value="K00400", qualifier=None),
                            Name(value="K00401", qualifier=None),
                            Name(value="K00402", qualifier=None),
                            Name(value="K03421", qualifier=None),
                            Name(value="K03422", qualifier=None),
                        )
                    ),
                ),
            ),
            Or(
                parts=(
                    Name(value="K02591", qualifier=None),
                    Name(value="K02586", qualifier=None),
                    Name(value="K02588", qualifier=None),
                )
            ),
        )
    ),
    "C": PipeChain(
        calls=(
            Call(
                value="not",
                args=(
                    Call(
                        value="filter_contains",
                        args=(
                            Name(value="kegg_description", qualifier=None),
                            String(value="nitrate reductase"),
                        ),
                    ),
                ),
            ),
            Call(
                value="column_count_values",
                args=(
                    Name(value="heme_regulatory_motif_count", qualifier=None),
                    Name(value="ge", qualifier=None),
                    Name(value="2", qualifier=None),
                    Name(value="ge", qualifier=None),
                    Name(value="2", qualifier=None),
                ),
            ),
        )
    ),
    "D": Call(
        value="percent",
        args=(
            Name(value="50", qualifier=None),
            Steps(
                parts=(
                    Name(value="K00399", qualifier=None),
                    Name(value="K14126", qualifier=None),
                    Name(value="K02588", qualifier=None),
                )
            ),
        ),
    ),
}


def test_dbcan_sub_ec_parsing():
    annotations = pl.DataFrame(
        {
            "dbcan_sub_ec": [
                "2.4.1.25:4;3.2.1.68:75;3.2.1.-:1",
                "1.2.3.4:2",
            ]
        }
    )

    parsed = annotations.select(ID_EXPR_DICT["dbcan_sub_ec"]).to_series().to_list()

    assert parsed == [
        ["EC:2.4.1.25", "EC:3.2.1.68", "EC:3.2.1.-"],
        ["EC:1.2.3.4"],
    ]


def test_basic_compiled_rules():
    compiled = CompiledRules.from_rules(rules_path=DATA_DIR / "basic_rules.tsv")
    assert compiled.needed_features == BASIC_NEEDS

    assert compiled.rules == BASIC_PARSE


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

    assert compiled.rules == COMPOUND_PARSE


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
