"""The `dram-rules` command: load annotations once, then run each enabled step."""

from __future__ import annotations

from pathlib import Path

import click
import polars as pl

from dram_rules.definitions import BACKUP_FASTA_COLUMNS, DEFAULT_FASTA_COLUMN
from dram_rules.make_summarize import (
    ECOSYSTEMS,
    TOPICS,
    add_trna,
    load_summarize_sheet,
    summarize_sheet,
    summary_collections,
    write_collection,
)
from dram_rules.processing.process_annotations import scan_tsv
from dram_rules.rule_parser.rules import (
    ID_PARSERS,
    CompiledRules,
    RuleError,
    prepare_present_map_df,
)

COMMON_RULES = Path(__file__).parent / "data/rules_common.tsv"
VISUALIZE_ONLY_ECOS = ("legacy", "default")


def split_commas(ctx, param, value: str | None) -> list[Path]:
    return [Path(p) for p in value.split(",") if p] if value else []


def split_choices(*choices: str):
    """Click callback: a comma-separated list of `choices` ("" or "none" is an
    empty list), or None when the option is not given."""

    def callback(ctx, param, value: str | None) -> list[str] | None:
        if value is None:
            return None
        choice = click.Choice([*choices, "none"])
        values = [choice.convert(v, param, ctx) for v in value.split(",") if v]
        return [v for v in values if v != "none"]

    return callback


def compile_sheets(paths, load, *args) -> dict[Path, CompiledRules]:
    """`{path: load(path, *args)}`, naming the sheet in any RuleError."""
    compiled = {}
    for path in paths:
        try:
            compiled[path] = load(path, *args)
        except RuleError as e:
            raise RuleError(f"Rule sheet {path}: {e}") from e
    return compiled


@click.command()
@click.option("--annotations", "-a", type=Path, required=True, help="Annotations TSV")
@click.option(
    "--fasta_column",
    "-f",
    default=DEFAULT_FASTA_COLUMN,
    help=f"Annotations column of genome names; falls back to {BACKUP_FASTA_COLUMNS}",
)
@click.option(
    "--output_dir", "-o", type=Path, default=Path.cwd(), help="Output directory"
)
@click.option(
    "--common_rules_tsv",
    "-c",
    type=click.Path(exists=True),
    default=COMMON_RULES,
    help="Common Alias Sheet: aliases every rule sheet can use",
)
@click.option(
    "--functions",
    "-fn",
    "function_files",
    type=click.Path(exists=True, dir_okay=False),
    multiple=True,
    help="Python file of extra rule functions. Can be given more than once.",
)
@click.option(
    "--mapping",
    "-m",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="Gene or genome abundances per sample: a Geneid (matching the annotations' "
    "query_id) or genome column, then one column per sample (visualize)",
)
@click.option(
    "--ecos",
    callback=split_choices(*ECOSYSTEMS, *VISUALIZE_ONLY_ECOS),
    help='Built-in ecosystem sheets, comma separated; none or "" for none. '
    "Visualize uses default (gut) and summarize all five when omitted; "
    "legacy and default are visualize-only.",
)
@click.option("--visualize", "--viz", "-v", is_flag=True, help="Run the visualize step")
@click.option(
    "--viz_custom",
    callback=split_commas,
    help="Custom visualize sheets, comma separated",
)
@click.option(
    "--save_needed_features",
    "-sf",
    is_flag=True,
    help="Write each visualize sheet's Gene IDs to needed_features.txt in its folder",
)
@click.option("--summarize", "--sum", "-s", is_flag=True, help="Run the summarize step")
@click.option(
    "--topics",
    default="default",
    callback=split_choices(*TOPICS, "default"),
    help='Built-in topic Summarize Sheets, comma separated; none or "" for none. '
    "default is all four.",
)
@click.option(
    "--sum_custom",
    callback=split_commas,
    help="Custom Summarize Sheets, comma separated",
)
@click.option(
    "--trna_path",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="tRNAscan output from DRAM annotate, so tRNA rules get Gene Counts (summarize)",
)
def main(
    annotations,
    fasta_column,
    output_dir,
    common_rules_tsv,
    function_files,
    mapping,
    ecos,
    visualize,
    viz_custom,
    save_needed_features,
    summarize,
    topics,
    sum_custom,
    trna_path,
):
    """Visualize, score traits of, and summarize DRAM annotations."""
    if not (visualize or summarize):
        raise click.UsageError("Choose at least one step: --visualize, --summarize")

    # Every sheet is compiled before the annotations load, so a bad sheet
    # stops the run early.
    viz_sheets = {}  # output folder: (sheet path, compiled sheet)
    if visualize:
        from dram_rules import make_product  # imports panel, which is slow

        viz_ecos = ["default"] if ecos is None else ecos
        named = [
            *((e, make_product.RULES_SYSTEMS[e]) for e in viz_ecos),
            *((p.stem, p) for p in viz_custom),
        ]
        names = [name for name, _ in named]
        if shared := sorted({n for n in names if names.count(n) > 1}):
            raise click.UsageError(f"Visualize sheets would share a folder: {shared}")
        compiled = compile_sheets(
            dict.fromkeys(p for _, p in named),
            make_product.load_sheet,
            common_rules_tsv,
            function_files,
        )
        viz_sheets = {
            output_dir / "visualize" / name: (p, compiled[p]) for name, p in named
        }

    if summarize:
        topics = list(TOPICS) if "default" in topics else topics
        sum_ecos = list(ECOSYSTEMS) if ecos is None else ecos
        if visualize_only := [e for e in sum_ecos if e in VISUALIZE_ONLY_ECOS]:
            raise click.UsageError(
                f"--ecos {visualize_only}: visualize only, not summarize"
            )
        collections = summary_collections(topics, sum_ecos, sum_custom)
        sheets = compile_sheets(
            dict.fromkeys(p for paths in collections.values() for p in paths),
            load_summarize_sheet,
            common_rules_tsv,
            function_files,
        )

    raw_anno = scan_tsv(annotations, infer_schema_length=10_000)
    fasta_column = next(
        (
            c
            for c in [fasta_column, *BACKUP_FASTA_COLUMNS]
            if c in raw_anno.collect_schema().names()
        ),
        None,
    )
    if fasta_column is None:
        raise click.UsageError(
            f"Annotations have none of the genome columns {[DEFAULT_FASTA_COLUMN, *BACKUP_FASTA_COLUMNS]}"
        )
    raw_anno = raw_anno.rename({fasta_column: "genome"}).with_columns(
        pl.col("genome").cast(pl.String)
    )

    if visualize:
        make_product.visualize(viz_sheets, raw_anno, mapping, save_needed_features)

    if summarize:
        # shortcut: summarize builds its own Hits, as tRNA rows are summarize-only; ticket 09 shares one load.
        sum_anno = add_trna(raw_anno, trna_path) if trna_path else raw_anno
        hits = prepare_present_map_df(
            sum_anno,
            count_col="query_id",
            group_col="genome",
            besthit_cols=list(ID_PARSERS),
            needed_features=set().union(*(c.needed_features for c in sheets.values())),
        )
        genomes = sum_anno.select("genome").unique().sort("genome").collect()
        # Every sheet is evaluated before any is written.
        summaries = {
            name: [
                summarize_sheet(sheets[p], hits, sum_anno, genomes["genome"].to_list())
                for p in paths
            ]
            for name, paths in collections.items()
        }
        for name, frames in summaries.items():
            write_collection(frames, name, output_dir / "summarize")


if __name__ == "__main__":
    main()
