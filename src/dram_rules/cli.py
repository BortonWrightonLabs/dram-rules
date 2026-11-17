"""The `dram-rules` command: load annotations once, then run each enabled refinement."""

from __future__ import annotations

import logging
from pathlib import Path

import click

from dram_rules import make_summarize
from dram_rules.definitions import BACKUP_FASTA_COLUMNS, DEFAULT_FASTA_COLUMN
from dram_rules.hits import Hits
from dram_rules.make_summarize import (
    ABUNDANCES,
    ECOSYSTEMS,
    TOPICS,
    load_summarize_sheet,
    summary_collections,
    write_collection,
)
from dram_rules.rule_parser.rules import (
    CompiledRules,
    RuleError,
)

COMMON_RULES = Path(__file__).parent / "data/rules_common.tsv"
TRAITS_SHEET = Path(__file__).parent / "data/traits_rules.tsv"
VISUALIZE_ONLY_ECOS = ("legacy", "default")

logger = logging.getLogger("dram.rules")


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


def sheet_folders(step_dir: Path, named) -> list[tuple[Path, Path]]:
    """`[(step_dir / name, sheet path)]` for each `(name, sheet path)`; two
    sheets sharing a folder are an error."""
    names = [name for name, _ in named]
    if shared := sorted({n for n in names if names.count(n) > 1}):
        raise click.UsageError(f"Sheets would share a folder in {step_dir}: {shared}")
    return [(step_dir / name, path) for name, path in named]


def compile_sheets(
    paths, load, *args, custom=(), columns=None
) -> dict[Path, CompiledRules]:
    """`{path: load(path, *args)}`, naming the sheet in any RuleError. Sheets
    in `custom` are read with the `columns` overrides (label_col, alias_col,
    group_col), which they must have."""
    compiled = {}
    for path in paths:
        overrides = (columns or {}) if path in custom else {}
        try:
            compiled[path] = load(path, *args, **overrides)
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
    "--output_dir",
    "-o",
    type=Path,
    help="Output directory; optional only with --no_save",
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
    "--custom_label_column",
    "-lc",
    help="Label column of every custom sheet; built-in sheets look it up",
)
@click.option(
    "--custom_alias_column",
    "-acc",
    help="Alias column of every custom sheet; built-in sheets look it up",
)
@click.option(
    "--custom_group_column",
    "-gc",
    help="Group column of every custom sheet; built-in sheets look it up",
)
@click.option(
    "--mapping",
    "-m",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="Gene or genome abundances per sample: a Geneid (matching the annotations' "
    "query_id) or genome column, then one column per sample",
)
@click.option(
    "--abundance",
    type=click.Choice(ABUNDANCES),
    default="mean",
    show_default=True,
    help="With --mapping, how summarize aggregates each rule's abundance in a "
    "sample over the genes where it holds",
)
@click.option(
    "--ecos",
    callback=split_choices(*ECOSYSTEMS, *VISUALIZE_ONLY_ECOS),
    help='Built-in ecosystem sheets, comma separated; none or "" for none. '
    "Visualize uses default (gut) and summarize all five when omitted; "
    "legacy and default are visualize-only.",
)
@click.option(
    "--visualize", "--viz", "-v", is_flag=True, help="Run the visualize refinement"
)
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
@click.option("--traits", "-t", is_flag=True, help="Run the traits refinement")
@click.option(
    "--traits_custom",
    callback=split_commas,
    help="Custom traits sheets, comma separated, added to the built-in one",
)
@click.option(
    "--dashboard",
    "-d",
    is_flag=True,
    help="Serve the live dashboard; the run must have exactly one visualize or "
    "traits sheet",
)
@click.option(
    "--port",
    type=click.INT,
    help="Port to serve the dashboard on (default 5006); only with --dashboard",
)
@click.option(
    "--no_save",
    is_flag=True,
    help="Write nothing (no HTML, no tables); only with --dashboard",
)
@click.option(
    "--summarize", "--sum", "-s", is_flag=True, help="Run the summarize refinement"
)
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
    help="tRNAscan output from DRAM annotate; its tRNAs count as genes in every refinement",
)
def main(
    annotations,
    fasta_column,
    output_dir,
    common_rules_tsv,
    function_files,
    custom_label_column,
    custom_alias_column,
    custom_group_column,
    mapping,
    abundance,
    ecos,
    visualize,
    viz_custom,
    save_needed_features,
    traits,
    traits_custom,
    dashboard,
    port,
    no_save,
    summarize,
    topics,
    sum_custom,
    trna_path,
):
    """Visualize, score traits of, and summarize DRAM annotations."""
    if not (visualize or traits or summarize):
        raise click.UsageError(
            "Choose at least one refinement: --visualize, --traits, --summarize"
        )
    if port is not None and not dashboard:
        raise click.UsageError("--port is only valid with --dashboard")
    if no_save and not dashboard:
        raise click.UsageError("--no_save is only valid with --dashboard")
    if no_save and summarize:
        raise click.UsageError(
            "--no_save would discard the summarize refinement's outputs"
        )
    if output_dir is None and not no_save:
        raise click.UsageError("Missing option '--output_dir' / '-o'")
    # Nothing is written with --no_save; the folders only name the sheets.
    out = output_dir or Path.cwd()
    custom = {*viz_custom, *traits_custom, *sum_custom}
    columns = {
        k: v
        for k, v in [
            ("label_col", custom_label_column),
            ("alias_col", custom_alias_column),
            ("group_col", custom_group_column),
        ]
        if v
    }

    # Every sheet is compiled before the annotations load, so a bad sheet
    # stops the run early. Visualize and traits sheets are the same kind.
    folders = []  # (output folder, sheet path)
    if visualize or traits:
        from dram_rules import make_product  # imports panel, which is slow

    if visualize:
        viz_ecos = ["default"] if ecos is None else ecos
        folders += sheet_folders(
            out / "visualize",
            [
                *((e, make_product.RULES_SYSTEMS[e]) for e in viz_ecos),
                *((p.stem, p) for p in viz_custom),
            ],
        )
    if traits:
        if traits_custom:
            folders += sheet_folders(
                out / "traits",
                [("default", TRAITS_SHEET), *((p.stem, p) for p in traits_custom)],
            )
        else:  # the common case stays flat
            folders.append((out / "traits", TRAITS_SHEET))
    if dashboard and len(folders) != 1:
        raise click.UsageError(
            f"--dashboard needs exactly one visualize or traits sheet, not {len(folders)}"
        )
    heatmap_sheets = {}  # output folder: (sheet path, compiled sheet)
    if folders:
        compiled = compile_sheets(
            dict.fromkeys(p for _, p in folders),
            make_product.load_sheet,
            common_rules_tsv,
            function_files,
            custom=custom,
            columns=columns,
        )
        heatmap_sheets = {folder: (p, compiled[p]) for folder, p in folders}

    sheets = {}
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
            custom=custom,
            columns=columns,
        )

    every_sheet = [*(c for _, c in heatmap_sheets.values()), *sheets.values()]
    try:
        # One load of the annotations' Hits, for every Gene ID of every refinement.
        hits = Hits.load(
            annotations,
            every_sheet,
            fasta_column=fasta_column,
            trna=trna_path,
            mapping=mapping,
        )
    except ValueError as e:
        raise click.UsageError(str(e)) from e

    # Every refinement evaluates before any writes, so a failure leaves no outputs.
    visualized = make_product.visualize(heatmap_sheets, hits) if heatmap_sheets else {}
    if summarize:
        summaries = make_summarize.summarize(
            {name: [sheets[p] for p in paths] for name, paths in collections.items()},
            hits,
            abundance,
        )

    if visualized and not no_save:
        make_product.save_visualized(heatmap_sheets, visualized, save_needed_features)
    if summarize:
        for name, df in summaries.items():
            write_collection(df, name, output_dir / "summarize")
    if dashboard:
        make_product.serve(visualized, 5006 if port is None else port, not no_save)


if __name__ == "__main__":
    main()
