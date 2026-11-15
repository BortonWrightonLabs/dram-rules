import polars as pl
import pytest

from dram_rules.rule_parser import rules
from dram_rules.rule_parser.rules import (
    ID_PARSERS,
    CompiledRules,
    prepare_present_map_df,
)


@pytest.fixture
def evaluate_annotations():
    """Evaluate a rule sheet file (put in one group, "grp") per genome on
    annotations keyed by input_fasta."""

    def evaluate_annotations(rules_path, annotations: pl.DataFrame, **kwargs):
        annotations = annotations.rename({"input_fasta": "genome"})
        compiled = CompiledRules.from_rules(
            rules=pl.scan_csv(rules_path, separator="\t").with_columns(
                group=pl.lit("grp")
            ),
            **kwargs,
        )
        hits = prepare_present_map_df(
            annotations,
            count_col="query_id",
            group_col="genome",
            besthit_cols=list(ID_PARSERS),
            needed_features=compiled.needed_features,
        )
        genomes = annotations["genome"].unique().sort().to_list()
        return rules.evaluate(compiled, hits, annotations, "genome", genomes)

    return evaluate_annotations
