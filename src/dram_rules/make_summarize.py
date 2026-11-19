"""The summarize refinement: Gene Counts per genome for Summarize Sheets, written as
Summary Collections (one workbook with a worksheet per sheet, plus a TSV). Derived
Summarize Sheets are built from rule sheets and the Description Tables (ADR 0004).

Ported from DRAM2's `bin/distill.py` (ADR 0002: dram-rules is standalone)."""

from __future__ import annotations

import logging
from pathlib import Path

import polars as pl
from xlsxwriter import Workbook

from dram_rules.hits import Hits
from dram_rules.rule_parser.rules import (
    CompiledRules,
    Name,
    RuleError,
    check_gene_level,
    walk,
)

logger = logging.getLogger("dram.rules")

# Sheet columns a hand-written Summary Collection keeps, besides the label and Group.
SUMMARY_COLUMNS = ["gene_description", "module", "header", "subheader"]
# A hand-written Collection's Group column; distill's name for it.
GROUP_COLUMN = "topic_ecosystem"
# Every Collection's columns after the label, those present in this order, coarse
# to fine: a topic sheet's Group holds its modules, a rule sheet's module its Groups.
COLUMN_ORDER = [
    "gene_description",
    "rule",
    GROUP_COLUMN,
    "module",
    "group",
    "header",
    "subheader",
]
# A Collection table's worksheet names: the sheet each row comes from.
SHEET_COLUMN = "sheet"
# Hand-written worksheets' sort; then the label, as distill sorts.
SORT_COLUMNS = [GROUP_COLUMN, "header", "subheader", "module"]

# A Derived Summarize Sheet's columns before the Gene Counts, in COLUMN_ORDER;
# `module` and `group` come from the rule sheet's module and Group columns.
DERIVED_COLUMNS = ["gene_id", "gene_description", "rule", "module", "group"]
# Joins a derived row's several rules, modules, Groups or descriptions; rule
# names contain commas.
LIST_SEPARATOR = " | "

SUMMARIZE_SHEETS_DIR = Path(__file__).parent / "data/summarize_sheets"
# Built-in topic Summarize Sheets by their --topics names, in distill's order.
TOPICS = {
    "assim": "assimilation_and_cofactor_metabolism.tsv",
    "cell": "cellular_machinery.tsv",
    "energy": "energy_acquisition_bioenergetics.tsv",
    "env": "environmental_interaction_and_adaptation.tsv",
}
DEFINITIONS_DIR = Path(__file__).parent / "data/definitions"
# An unedited dbCAN release: headerless, `#` comments, family then description.
DBCAN_FAMILIES = DEFINITIONS_DIR / "dbcan.fam-activities.tsv"
# Hand curated, gene_id and description; wins over the dbCAN table.
GENE_DESCRIPTIONS = DEFINITIONS_DIR / "gene_descriptions.tsv"


def load_rule_sheet(
    path: Path,
    common_rules_tsv: Path,
    function_files: tuple[str, ...],
    label_col: str | None = None,
    alias_col: str | None = None,
    group_col: str | None = None,
) -> CompiledRules:
    """Compile a rule sheet. The label, alias and Group columns are looked up
    unless given."""
    return CompiledRules.from_rules(
        rules_path=path,
        label_col=label_col,
        alias_col=alias_col,
        group_col=group_col,
        common_rules_path=common_rules_tsv,
        function_files=function_files,
    )


def load_summarize_sheet(path: Path, *args, **columns) -> CompiledRules:
    """`load_rule_sheet` for a Summarize Sheet; coverage settings in it are
    an error."""
    compiled = load_rule_sheet(path, *args, **columns)
    check_gene_level(compiled)
    return compiled


def _described_rows(path: Path, description_col: str) -> list[tuple[str, str]]:
    """A table's (gene_id, description) rows with both non-empty, stripped."""
    df = pl.read_csv(path, separator="\t", infer_schema=False)
    df = df.rename({c: c.lstrip("\ufeff").strip() for c in df.columns})
    if missing := {"gene_id", description_col} - set(df.columns):
        raise RuleError(f"Description Table {path}: no {sorted(missing)} column")
    rows = df.select("gene_id", description_col).iter_rows()
    return [
        (g.strip(), d.strip()) for g, d in rows if g and d and g.strip() and d.strip()
    ]


def load_descriptions(custom: list[Path]) -> dict[str, str]:
    """The Description Tables as one lookup, keyed by uppercase Gene ID: the
    dbCAN families, then the built-in topic Summarize Sheets, then
    `gene_descriptions.tsv`, then each `custom` table (gene_id, description),
    each winning over the ones before."""
    descriptions = {}
    for line in DBCAN_FAMILIES.read_text().splitlines():
        family, _, description = line.partition("\t")
        # Some dbCAN lines hold extra tabs inside the description.
        description = " ".join(description.split())
        if not line.startswith("#") and description:
            descriptions[family.strip().upper()] = description
    # A topic sheet's gene_id can be several Gene IDs joined with `|`, each
    # taking the row's description. Where the sheets disagree, the first in
    # TOPICS order wins.
    topics = {}
    for name in TOPICS.values():
        for gene_ids, description in _described_rows(
            SUMMARIZE_SHEETS_DIR / name, "gene_description"
        ):
            for gene_id in gene_ids.split("|"):
                key = gene_id.strip(" `").upper()
                # A description that only repeats the Gene ID describes nothing.
                if key and description.upper() != key:
                    topics.setdefault(key, description)
    descriptions |= topics
    for path in [GENE_DESCRIPTIONS, *custom]:
        for gene_id, description in _described_rows(path, "description"):
            descriptions[gene_id.upper()] = description
    return descriptions


def _join(values) -> str | None:
    """`values` without empties or repeats, in order, joined with LIST_SEPARATOR."""
    return LIST_SEPARATOR.join(dict.fromkeys(v for v in values if v)) or None


def derive_sheet(
    compiled: CompiledRules, descriptions: dict[str, str]
) -> tuple[pl.DataFrame, list[str]]:
    """The Derived Summarize Sheet of a rule sheet, and the Gene IDs and
    families it lacks descriptions for.

    One row (DERIVED_COLUMNS) per Gene ID any rule looks up, aliases expanded
    and negated ones included, matched regardless of case. A Family;EC Gene ID
    (one with a `;`) is described as `<module>; <group>; <name>` of each rule
    using it, and its family must have a description; any other Gene ID is
    looked up whole. Rows are ordered by module, then Group, each in sheet
    order, then by where the Gene ID first appears; a Gene ID with several
    sits under its first."""
    label, group = compiled.label_col, compiled.group_col
    sheet = compiled.df.drop_nulls(label)
    places = {
        row[label]: (row.get("module"), row.get(group))
        for row in sheet.iter_rows(named=True)
    }
    uses = {}  # uppercase Gene ID: [(as written, rule, module, group)]
    for rule, expr in compiled.rules.items():
        names = (n.value for n in walk(expr) if isinstance(n, Name))
        for gene_id in dict.fromkeys(names):
            uses.setdefault(gene_id.upper(), []).append((gene_id, rule, *places[rule]))
    modules = {m: i for i, m in enumerate(dict.fromkeys(m for m, _ in places.values()))}
    groups = {g: i for i, g in enumerate(dict.fromkeys(g for _, g in places.values()))}
    ordered = sorted(
        uses.values(), key=lambda u: (modules[u[0][2]], groups[u[0][3]])
    )  # stable, so ties keep first appearance

    rows, missing = [], []
    for found in ordered:
        gene_id = found[0][0]
        if ";" in gene_id:
            family = gene_id.split(";")[0]
            if family.upper() not in descriptions:
                missing.append(f"family {family}")
            description = _join(
                "; ".join(str(v) for v in (m, g, r) if v is not None)
                for _, r, m, g in found
            )
        else:
            description = descriptions.get(gene_id.upper())
            if description is None:
                missing.append(gene_id)
        rows.append(
            (
                gene_id,
                description,
                *(
                    _join(str(u[i]) for u in found if u[i] is not None)
                    for i in (1, 2, 3)
                ),
            )
        )
    df = pl.DataFrame(
        rows, schema=dict.fromkeys(DERIVED_COLUMNS, pl.String), orient="row"
    )
    return df, list(dict.fromkeys(missing))


def derive_collection(
    sheets: dict[str, CompiledRules], descriptions: dict[str, str], strict: bool
) -> list[pl.DataFrame]:
    """The Derived Summarize Sheet of each rule sheet in `sheets` (by
    worksheet name), its name in SHEET_COLUMN. A Gene ID or family without
    a description raises a RuleError listing them all by sheet when `strict`,
    and is otherwise left undescribed with one warning."""
    frames, missing = [], {}
    for name, compiled in sheets.items():
        df, lacking = derive_sheet(compiled, descriptions)
        frames.append(df.with_columns(pl.lit(name).alias(SHEET_COLUMN)))
        if lacking:
            missing[name] = lacking
    if missing and strict:
        listed = "\n".join(f"  {n}: {', '.join(ids)}" for n, ids in missing.items())
        raise RuleError(
            "Gene IDs and families without a description; add them to a "
            f"Description Table:\n{listed}"
        )
    if missing:
        ids = [i for lacking in missing.values() for i in lacking]
        logger.warning(
            f"{len(ids)} Gene IDs and families in {sorted(missing)} have no "
            f"description, e.g. {', '.join(ids[:5])}"
        )
    return frames


def count_derived(frames: list[pl.DataFrame], hits: Hits) -> pl.DataFrame | None:
    """`frames` (Derived Summarize Sheets) one after another, each Gene ID
    with one Gene Count column per genome: the genome's genes annotated with
    it. None when there are no rows."""
    if not frames or not (df := pl.concat(frames)).height:
        return None
    key = pl.col("gene_id").str.to_uppercase()
    counts = (
        hits.hits.filter(pl.col("hit").is_in(df.select(key).to_series().implode()))
        .group_by("hit", "genome")
        .agg(gene_count=pl.col("query_id").n_unique())
        .pivot(on="genome", index="hit", values="gene_count")
    )
    df = df.with_columns(hit=key).join(
        counts, on="hit", how="left", maintain_order="left"
    )
    return df.select(
        *DERIVED_COLUMNS,
        SHEET_COLUMN,
        *(
            pl.col(g).fill_null(0) if g in counts.columns else pl.lit(0).alias(g)
            for g in hits.genomes
        ),
    )


def summarize_sheet(
    compiled: CompiledRules, hits: Hits, name: str
) -> pl.DataFrame | None:
    """One row per rule, in sheet order: the label, the Group (as GROUP_COLUMN)
    and the SUMMARY_COLUMNS the sheet has in COLUMN_ORDER, `name` (as
    SHEET_COLUMN), then one Gene Count column per genome.
    Rules left out by `include_when` are absent; None if no rule is left."""
    label, group = compiled.label_col, compiled.group_col
    dfs = hits.evaluate(compiled, "gene")
    if not dfs:
        return None
    counts = pl.concat(dfs.values()).pivot(
        on="genome", index=label, values="gene_count"
    )
    kept = {
        GROUP_COLUMN: pl.col(group),
        **{c: pl.col(c) for c in SUMMARY_COLUMNS if c in compiled.df.columns},
    }
    return (
        compiled.df.drop_nulls(label)
        .select(
            label,
            *(kept[c].alias(c) for c in COLUMN_ORDER if c in kept),
            pl.lit(name).alias(SHEET_COLUMN),
        )
        .join(counts, on=label, how="inner", maintain_order="left")
    )


def summarize(
    collections: dict[str, dict[str, CompiledRules]], hits: Hits
) -> dict[str, pl.DataFrame | None]:
    """Each Summary Collection's table: the `summarize_sheet` rows of its
    sheets (by worksheet name), one sheet after another. The label column keeps the sheets' name,
    or is `name` when they disagree. None when no rule is left."""
    tables = {}
    for name, sheets in collections.items():
        frames = [
            df
            for n, c in sheets.items()
            if (df := summarize_sheet(c, hits, n)) is not None
        ]
        if len({df.columns[0] for df in frames}) > 1:
            frames = [df.rename({df.columns[0]: "name"}) for df in frames]
        tables[name] = pl.concat(frames, how="diagonal_relaxed") if frames else None
    return tables


def write_collection(
    df: pl.DataFrame | None, name: str, out_dir: Path, sort: bool = True
) -> None:
    """Write a Summary Collection's table as `<name>.xlsx`, a worksheet per
    SHEET_COLUMN value (sorted by SORT_COLUMNS and the label unless not `sort`), and
    `<name>.tsv`. Nothing is written for None."""
    if df is None:
        return
    label = df.columns[0]
    out_dir.mkdir(parents=True, exist_ok=True)
    df.write_csv(out_dir / f"{name}.tsv", separator="\t")
    order = [c for c in SORT_COLUMNS if c in df.columns] + [label] if sort else []
    with Workbook(out_dir / f"{name}.xlsx") as wb:
        for (sheet,), frame in df.group_by(SHEET_COLUMN, maintain_order=True):
            (frame.sort(order) if order else frame).drop(SHEET_COLUMN).write_excel(
                workbook=wb, worksheet=sheet, autofit=True
            )
    logger.info(f"Saved {out_dir / name}.xlsx and .tsv")
