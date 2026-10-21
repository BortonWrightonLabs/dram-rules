from __future__ import annotations

import csv
import logging
import re
from collections import defaultdict
from pathlib import Path

import click
import polars as pl

from dram_rules.cli_utils import get_column_name
from dram_rules.definitions import BACKUP_FASTA_COLUMNS, DEFAULT_FASTA_COLUMN
from dram_rules.rule_parser.rules import (
    ID_EXPR_DICT,
    CompiledRules,
    RuleError,
    build_present_map,
    prepare_present_map_df,
)

logger = logging.getLogger(__name__)
DATA_DIR = Path(__file__).parent / "data"
COMMON_RULES = DATA_DIR / "rules_common.tsv"
RULES_SYSTEMS = {
    name: DATA_DIR / filename
    for name, filename in {
        "default": "rules.tsv",
        "ag": "rules_ag.tsv",
        "bgc": "rules_bgc.tsv",
        "eng_sys": "rules_eng_sys.tsv",
        "gut": "rules_gut.tsv",
        "traits": "traits_rules.tsv",
    }.items()
}


def read_tsv(path: Path, required: tuple[str, ...]) -> pl.DataFrame:
    frame = pl.read_csv(path, separator="\t", infer_schema=False, comment_prefix="#")
    missing = set(required) - set(frame.columns)
    if missing:
        raise ValueError(f"{path}: missing columns: {', '.join(sorted(missing))}")
    return frame


def clean_description(gene_id: str, description: str | None) -> str:
    description = (description or "").strip()
    if re.fullmatch(r"[;:\s]*", description) or description.upper() == gene_id.upper():
        return ""
    return re.sub(
        r"\bEC\s*:?\s*([0-9-]+\.[0-9-]+\.[0-9-]+\.[0-9-]+)",
        r"EC:\1",
        description,
    )


def validate_definitions(
    records: list[tuple[str, str | None, str]],
) -> dict[str, str]:
    sources: dict[str, list[str]] = defaultdict(list)
    descriptions = {}
    for gene_id, description, source in records:
        gene_id = (gene_id or "").strip()
        if not gene_id:
            raise ValueError(f"Empty gene_id in {source}")
        key = gene_id.upper()
        sources[key].append(source)
        descriptions[key] = clean_description(gene_id, description)
    duplicates = sorted(key for key, locations in sources.items() if len(locations) > 1)
    if duplicates:
        details = "\n".join(
            f"{key}\t{', '.join(sorted(set(sources[key])))}" for key in duplicates
        )
        raise ValueError(
            f"Duplicate definition gene_ids (case-insensitive):\n{details}"
        )
    return descriptions


def load_definitions(
    definitions_dir: Path, lookup_tsv: Path | None = None
) -> dict[str, str]:
    records = []
    path = definitions_dir / "summarize_genes.tsv"
    for gene_id, description in (
        read_tsv(path, ("gene_id", "gene_description"))
        .select("gene_id", "gene_description")
        .iter_rows()
    ):
        records.append((gene_id, description, str(path)))

    path = definitions_dir / "fam-substrate-mapping.tsv"
    mapping = read_tsv(path, ("Family", "EC_Number", "Name"))
    skipped = 0
    for family, ec, description in mapping.select(
        "Family", "EC_Number", "Name"
    ).iter_rows():
        family, ec = (family or "").strip(), (ec or "").strip()
        if not family or not ec:
            skipped += 1
            continue
        records.append((f"{family};EC:{ec}", description, str(path)))
    if skipped:
        logger.warning(
            "Skipped %s rows without Family or EC_Number in %s", skipped, path
        )

    path = definitions_dir / "dbcan.fam-activities.tsv"
    with path.open(newline="") as handle:
        for row in csv.reader(handle, delimiter="\t"):
            if not row or not row[0].strip() or row[0].lstrip().startswith("#"):
                continue
            if len(row) < 2:
                raise ValueError(f"{path}: expected at least two columns for {row[0]}")
            records.append((row[0], row[1], str(path)))

    definitions = validate_definitions(records)
    if lookup_tsv is not None:
        overrides = read_tsv(lookup_tsv, ("gene_id", "gene_description"))
        definitions.update(
            validate_definitions(
                [
                    (gene, description, str(lookup_tsv))
                    for gene, description in overrides.select(
                        "gene_id", "gene_description"
                    ).iter_rows()
                ]
            )
        )
    return definitions


def convert_rules(
    rules_tsv: Path,
    definitions: dict[str, str],
    common_rules_tsv: Path = COMMON_RULES,
    all_metadata: bool = False,
    strict: bool = False,
) -> pl.DataFrame:
    rules = pl.read_csv(
        rules_tsv, separator="\t", infer_schema=False, truncate_ragged_lines=True
    ).filter(~pl.all_horizontal(pl.all().is_null()))
    if not {"name", "rule"}.issubset(rules.columns):
        raise ValueError(f"{rules_tsv}: rules require name and rule columns")

    labels = [
        f"__summarize_{index}" if (name or "").strip() else ""
        for index, name in enumerate(rules.get_column("name"))
    ]
    compiled = CompiledRules.from_rules(
        rules=rules.with_columns(pl.Series("name", labels)).lazy(),
        common_rules_path=common_rules_tsv,
        allow_visualize_functions=True,
    )
    excluded = {"name", "rule", "alias", "gene_id", "gene_description"}
    metadata = [
        column
        for column in rules.columns
        if column not in excluded and (all_metadata or column == "topic_ecosystem")
    ]
    values: dict[str, dict[str, set[str]]] = {}
    spellings: dict[str, set[str]] = defaultdict(set)
    for label, row in zip(labels, rules.iter_rows(named=True)):
        if not label:
            continue
        for feature in compiled.features_by_rules[label]:
            key = feature.upper()
            spellings[key].add(feature)
            associated = values.setdefault(
                key, {column: set() for column in [*metadata, "name"]}
            )
            for column in associated:
                value = (row[column] or "").strip()
                if value:
                    associated[column].add(value)

    missing = sorted(min(spellings[key]) for key in values if not definitions.get(key))
    if strict and missing:
        raise ValueError("Missing gene descriptions:\n" + "\n".join(missing))
    columns = ["gene_id", "gene_description", *metadata, "name"]
    rows = [
        {
            "gene_id": min(spellings[key]),
            "gene_description": definitions.get(key, ""),
            **{column: ",".join(sorted(items)) for column, items in associated.items()},
        }
        for key, associated in values.items()
    ]
    return pl.DataFrame(rows, schema={column: pl.String for column in columns}).sort(
        "gene_id"
    )


def summarize_annotations(
    sheet: pl.DataFrame,
    raw_annotations: pl.DataFrame,
    fasta_column: str = DEFAULT_FASTA_COLUMN,
) -> pl.DataFrame:
    fasta_column = get_column_name(
        fasta_column,
        ["genome", *BACKUP_FASTA_COLUMNS],
        name_of_data="annotations",
        df=raw_annotations,
    )
    if "query_id" not in raw_annotations.columns:
        raise ValueError("Annotations require a query_id column")
    if fasta_column == "query_id":
        raise ValueError("The genome column must differ from query_id")
    raw_annotations = raw_annotations.with_columns(
        pl.col(fasta_column).cast(pl.String).alias("genome")
    )
    if fasta_column != "genome":
        raw_annotations = raw_annotations.drop(fasta_column)
    for column in ("genome", "query_id"):
        if raw_annotations.filter(
            pl.col(column).is_null()
            | (pl.col(column).cast(pl.String).str.strip_chars() == "")
        ).height:
            raise ValueError(f"Annotations contain empty {column} values")
    genomes = sorted(raw_annotations.get_column("genome").unique().to_list())
    collisions = set(genomes) & set(sheet.columns)
    if collisions:
        raise ValueError(
            f"Genome names conflict with metadata columns: {', '.join(sorted(collisions))}"
        )
    needed = {gene.upper() for gene in sheet.get_column("gene_id")}
    identifier_columns = [
        column
        for column in raw_annotations.columns
        if column not in {"query_id", "genome"}
        and (column in ID_EXPR_DICT or column.endswith("_id"))
    ]
    if identifier_columns and needed:
        annotations = prepare_present_map_df(
            df=raw_annotations,
            count_col="query_id",
            group_col="genome",
            besthit_cols=identifier_columns,
            needed_features=needed,
        )
        counts = build_present_map(
            df=annotations,
            sample_col="genome",
            needed_features=needed,
            sample_names=genomes,
        )
        return sheet.with_columns(
            pl.Series(
                genome,
                [int(counts[gene.upper()][index]) for gene in sheet["gene_id"]],
                dtype=pl.Int64,
            )
            for index, genome in enumerate(genomes)
        )
    return sheet.with_columns(
        pl.lit(0, dtype=pl.Int64).alias(genome) for genome in genomes
    )


@click.command()
@click.option(
    "--annotations",
    "-a",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--fasta-column",
    "--fasta_column",
    "-f",
    default=DEFAULT_FASTA_COLUMN,
    show_default=True,
)
@click.option(
    "--output-dir",
    "--output_dir",
    "-o",
    type=Path,
    default=Path("."),
    show_default=True,
)
@click.option(
    "--rules-system", "--rules_system", "-rs", type=click.Choice(list(RULES_SYSTEMS))
)
@click.option(
    "--rules-tsv",
    "--rules_tsv",
    "-r",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--common-rules-tsv",
    "--common_rules_tsv",
    "-c",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=COMMON_RULES,
)
@click.option(
    "--definitions-dir",
    type=click.Path(exists=True, file_okay=False, path_type=Path),
    default=DATA_DIR / "definitions",
    show_default=True,
)
@click.option(
    "--lookup-tsv",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="Custom gene_id/gene_description TSV overriding validated built-in definitions.",
)
@click.option(
    "--strict",
    is_flag=True,
    help="Fail if any required feature lacks a description, including CAZy features.",
)
@click.option(
    "--all-metadata",
    is_flag=True,
    help="Carry all rule metadata columns; default: topic_ecosystem and name.",
)
@click.option(
    "--save-converted-sheet",
    is_flag=True,
    help="Also save the feature sheet without genome counts.",
)
@click.option(
    "--save-lookup",
    is_flag=True,
    help="Also save the combined definition lookup after custom overrides.",
)
def main(
    annotations,
    fasta_column,
    output_dir,
    rules_system,
    rules_tsv,
    common_rules_tsv,
    definitions_dir,
    lookup_tsv,
    strict,
    all_metadata,
    save_converted_sheet,
    save_lookup,
):
    """Convert visualization rules into per-feature gene counts for each genome."""
    if rules_system and rules_tsv:
        raise click.UsageError("Supply either --rules-system or --rules-tsv")
    rules_tsv = rules_tsv or RULES_SYSTEMS[rules_system or "default"]
    try:
        definitions = load_definitions(definitions_dir, lookup_tsv)
        sheet = convert_rules(
            rules_tsv, definitions, common_rules_tsv, all_metadata, strict
        )
        raw_annotations = pl.read_csv(
            annotations,
            separator="\t",
            infer_schema_length=None,
            schema_overrides={"query_id": pl.String, fasta_column: pl.String},
        )
        summary = summarize_annotations(sheet, raw_annotations, fasta_column)
    except (ValueError, OSError, RuleError, pl.exceptions.PolarsError) as error:
        raise click.ClickException(str(error)) from error

    output_dir.mkdir(parents=True, exist_ok=True)
    output = output_dir / f"{rules_tsv.stem}_summary.tsv"
    summary.write_csv(output, separator="\t")
    if save_converted_sheet:
        sheet.write_csv(
            output_dir / f"{rules_tsv.stem}_summarize_rules.tsv", separator="\t"
        )
    if save_lookup:
        pl.DataFrame(
            sorted(definitions.items()),
            schema={"gene_id": pl.String, "gene_description": pl.String},
            orient="row",
        ).write_csv(output_dir / "master_definitions.tsv", separator="\t")
    click.echo(f"Wrote {summary.height} gene IDs to {output}")


if __name__ == "__main__":
    main()
