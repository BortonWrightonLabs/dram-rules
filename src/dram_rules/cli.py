"""The `dram-rules` command: load annotations once, then run each enabled refinement."""

from __future__ import annotations

import logging
from collections.abc import Iterable
from io import BytesIO
from pathlib import Path

import click
from xlsxwriter import Workbook
from xlsxwriter.exceptions import DuplicateWorksheetName, InvalidWorksheetName

from dram_rules import make_summarize
from dram_rules.definitions import BACKUP_FASTA_COLUMNS, DEFAULT_FASTA_COLUMN
from dram_rules.hits import Hits
from dram_rules.make_summarize import (
    SUMMARIZE_SHEETS_DIR,
    TOPICS,
    load_rule_sheet,
    load_summarize_sheet,
    write_collection,
)
from dram_rules.rule_parser.rules import (
    CompiledRules,
    RuleError,
)

COMMON_RULES = Path(__file__).parent / "data/rules_common.tsv"
TRAITS_SHEET = Path(__file__).parent / "data/traits_rules.tsv"
# The built-in ecosystem sheets by --ecos name: every sheet in the folder.
ECOSYSTEMS = {
    p.stem: p for p in sorted((Path(__file__).parent / "data/ecosystems").glob("*.tsv"))
}
VISUALIZE_ONLY_ECOS = {
    "legacy": Path(__file__).parent / "data/rules.tsv",
    "default": ECOSYSTEMS["gut"],
}

logger = logging.getLogger("dram.rules")


def split_commas(ctx, param, value: str | None) -> list[Path]:
    return [Path(p) for p in value.split(",") if p] if value else []


def split_choices(*choices: str, every: Iterable[str]):
    """Click callback: a comma-separated list of `choices`, where "all" is
    each of `every` and "" or "none" is an empty list, repeats removed; or None
    when the option is not given."""

    def callback(ctx, param, value: str | None) -> list[str] | None:
        if value is None:
            return None
        choice = click.Choice([*choices, "all", "none"])
        values = [choice.convert(v, param, ctx) for v in value.split(",") if v]
        expanded = (n for v in values for n in (every if v == "all" else [v]))
        return list(dict.fromkeys(n for n in expanded if n != "none"))

    return callback


def check_unique(names: list[str], place: str) -> None:
    """Two sheets with one of `names` would share `place`: an error."""
    if shared := sorted({n for n in names if names.count(n) > 1}):
        raise click.UsageError(f"Sheets would share {place}: {shared}")


def check_worksheets(names: list[str], collection: str) -> None:
    """`names` that can't all be worksheets of one workbook (too long, a
    character Excel forbids, or shared regardless of case): an error naming
    `collection`. xlsxwriter checks them, so its rules aren't copied here."""
    try:
        with Workbook(BytesIO()) as wb:
            for name in names:
                wb.add_worksheet(name)
    except (DuplicateWorksheetName, InvalidWorksheetName) as e:
        raise click.UsageError(f"{collection}: {e}") from e


def sheet_folders(step_dir: Path, named) -> list[tuple[Path, Path]]:
    """`[(step_dir / name, sheet path)]` for each `(name, sheet path)`; two
    sheets sharing a folder are an error."""
    check_unique([name for name, _ in named], f"a folder in {step_dir}")
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
    "query_id) or genome column, then one column per sample (visualize)",
)
@click.option(
    "--ecos",
    callback=split_choices(*ECOSYSTEMS, *VISUALIZE_ONLY_ECOS, every=ECOSYSTEMS),
    help=f"Ecosystem Sheets ({', '.join(ECOSYSTEMS)}), comma separated; all for "
    'every one, none or "" for none. Visualize uses default (gut) and '
    "summarize all when omitted; summarize skips legacy and default "
    "(visualize-only) with a warning. Summarize derives each sheet's ecosystem_summary "
    "worksheet from it.",
)
@click.option(
    "--visualize", "--viz", "-v", is_flag=True, help="Run the visualize refinement"
)
@click.option(
    "--viz_custom",
    callback=split_commas,
    help="Custom visualize sheets, comma separated. With --summarize, each "
    "also gets a derived worksheet in custom_derived_summary.",
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
    callback=split_choices(*TOPICS, "default", every=TOPICS),
    help=f"Built-in topic Summarize Sheets ({', '.join(TOPICS)}), comma separated; "
    'all for every one, none or "" for none. default (the default) is also '
    "every one.",
)
@click.option(
    "--sum_custom",
    callback=split_commas,
    help="Custom Summarize Sheets, comma separated",
)
@click.option(
    "--descriptions",
    callback=split_commas,
    help="Custom Description Tables (gene_id and description columns), comma "
    "separated; they win over the built-in ones in derived summaries",
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
    descriptions,
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
                *((e, {**ECOSYSTEMS, **VISUALIZE_ONLY_ECOS}[e]) for e in viz_ecos),
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

    sheets, derived_sheets = {}, {}
    if summarize:
        topics = list(TOPICS) if "default" in topics else topics
        sum_ecos = list(ECOSYSTEMS) if ecos is None else ecos
        if visualize_only := [e for e in sum_ecos if e in VISUALIZE_ONLY_ECOS]:
            logger.warning("Summarize skips visualize-only --ecos %s", visualize_only)
            sum_ecos = [e for e in sum_ecos if e not in VISUALIZE_ONLY_ECOS]
        check_worksheets([p.stem for p in viz_custom], "custom_derived_summary")
        check_worksheets([p.stem for p in sum_custom], "custom_summary")
        # Each Collection's sheets by worksheet name.
        collections = {
            "topic_summary": {
                t: SUMMARIZE_SHEETS_DIR / f for t, f in TOPICS.items() if t in topics
            },
            "custom_summary": {p.stem: p for p in sum_custom},
        }
        sheets = compile_sheets(
            dict.fromkeys(p for named in collections.values() for p in named.values()),
            load_summarize_sheet,
            common_rules_tsv,
            function_files,
            custom=custom,
            columns=columns,
        )
        # Derived Summarize Sheets: built-in ones strict, custom ones not.
        derived = {
            "ecosystem_summary": ({e: ECOSYSTEMS[e] for e in sum_ecos}, True),
            "custom_derived_summary": ({p.stem: p for p in viz_custom}, False),
        }
        derived_sheets = compile_sheets(
            dict.fromkeys(p for named, _ in derived.values() for p in named.values()),
            load_rule_sheet,
            common_rules_tsv,
            function_files,
            custom=custom,
            columns=columns,
        )
        lookup = make_summarize.load_descriptions(descriptions)
        derived_frames = {
            name: make_summarize.derive_collection(
                {n: derived_sheets[p] for n, p in named.items()}, lookup, strict
            )
            for name, (named, strict) in derived.items()
        }

    every_sheet = [
        *(c for _, c in heatmap_sheets.values()),
        *sheets.values(),
        *derived_sheets.values(),
    ]
    try:
        # One load of the annotations' Hits, for every Gene ID of every refinement.
        hits = Hits.load(
            annotations,
            every_sheet,
            fasta_column=fasta_column,
            trna=trna_path,
            # Summarize ignores the mapping file, so only heatmaps load it.
            mapping=mapping if heatmap_sheets else None,
        )
    except ValueError as e:
        raise click.UsageError(str(e)) from e

    # Every refinement evaluates before any writes, so a failure leaves no outputs.
    visualized = make_product.visualize(heatmap_sheets, hits) if heatmap_sheets else {}
    if summarize:
        if mapping:
            logger.info("Summarize ignores --mapping for now")
        summaries = make_summarize.summarize(
            {
                name: {n: sheets[p] for n, p in named.items()}
                for name, named in collections.items()
            },
            hits,
        )
        derived_summaries = {
            name: make_summarize.count_derived(frames, hits)
            for name, frames in derived_frames.items()
        }

    if visualized and not no_save:
        make_product.save_visualized(heatmap_sheets, visualized, save_needed_features)
    if summarize:
        for name, df in summaries.items():
            write_collection(df, name, output_dir / "summarize")
        for name, df in derived_summaries.items():
            write_collection(df, name, output_dir / "summarize", sort=False)
    if dashboard:
        make_product.serve(visualized, 5006 if port is None else port, not no_save)


if __name__ == "__main__":
    main()
