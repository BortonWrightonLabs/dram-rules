from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Dict, Tuple, List, Optional, Set, Iterable
import operator
import os
from pathlib import Path
import functools
import json

import networkx as nx
import numpy as np
import polars as pl
from lark import Lark, Transformer, LarkError
from dram_rules.definitions import (
    BACKUP_GROUPBY_COLUMNS,
    DEFAULT_GROUPBY_COLUMN,
)

logger = logging.getLogger("drams.rules")

OP_TO_EXPR = {
    "gt": operator.gt,
    "ge": operator.ge,
    "lt": operator.lt,
    "le": operator.le,
    "eq": operator.eq,
    "ne": operator.ne,
}
ALLOWED_CMPOPS = set(OP_TO_EXPR.keys())
AT_LEAST_MODES = set(["PRESENCE", "COUNTS"])


class RuleError(Exception):
    pass


# ---------------------------------------------------------------------------
# Parsing annotation columns into identifiers
#
# A parser takes a column name and returns an expression giving the list of
# identifiers in each cell. Parsers only *find* the ids; case and "EC:" prefixes
# are normalized once, later, by `_match_key`, so "3.2.1.4", "EC:3.2.1.4" and
# "ec:3.2.1.4" all match each other.
# ---------------------------------------------------------------------------


def _text(col: str) -> pl.Expr:
    return pl.col(col).cast(pl.String)


def split_on(sep: str):
    """Parser for cells holding ids separated by `sep`, e.g. "K00001,K00002"."""
    return lambda col: _text(col).str.split(sep)


def find_all(pattern: str):
    """Parser whose ids are every match of the regex `pattern` in the cell.
    Note: Polars returns the *whole* match and ignores capture groups."""
    return lambda col: _text(col).str.extract_all(pattern)


def whole_cell(col: str) -> pl.Expr:
    """Parser for cells that are a single id."""
    return pl.concat_list(_text(col))


def without_version(col: str) -> pl.Expr:
    """Parser for a single accession, dropping its version: PF00001.23 -> PF00001."""
    return pl.concat_list(_text(col).str.replace(r"\.\d+$", ""))


def first_word_of_each(sep: str):
    """Parser for "<id> <description>" entries separated by `sep`,
    e.g. "M1 some desc, M2 other" -> ["M1", "M2"]."""
    return lambda col: (
        _text(col)
        .str.split(sep)
        .list.eval(pl.element().str.strip_chars().str.split(" ").list.first())
    )


def dbcan_ecs(col: str) -> pl.Expr:
    """dbCAN sub EC numbers with hit counts, "3.2.1.4:4;3.2.1.-:1" -> ["3.2.1.4", "3.2.1.-"].
    A lone "-" means no EC."""
    return _text(col).replace("-", None).str.replace_all(r":\d+", "").str.split(";")


def cazy_ecs(col: str) -> pl.Expr:
    """EC numbers inside CAZy hit descriptions, written "(EC 3.2.1.4)" or "EC:3.2.1.4"."""
    return _text(col).str.replace_all(r"EC\s+", "EC:").str.extract_all(r"EC:[\d.\-]+")


EC_NUMBER = r"EC:\d+\.\d+\.\d+\.\d+"
PFAM = r"PF\d{5}"

# Known annotation columns and how to read them. Any other column is still
# searched: see `parse_ids` for the fallbacks.
ID_PARSERS = {
    "camper_id": whole_cell,
    "fegenie_id": whole_cell,
    "sulfur_id": whole_cell,
    "kegg_genes_id": whole_cell,
    "kofam_genes_id": whole_cell,
    "ko_id": split_on(","),
    "kegg_id": split_on("/"),
    "kofam_id": split_on("/"),
    "peptidase_family": split_on(";"),
    "merops_family": split_on(";"),
    "dbcan_id": split_on(";"),
    "dbcan_sub_id": split_on(";"),
    "dbcan_sub_substrate": split_on(";"),
    "camper_EC": find_all(EC_NUMBER),
    "kegg_EC": find_all(EC_NUMBER),
    "kofam_EC": find_all(EC_NUMBER),
    "cazy_subfam_ec": find_all(r"[\d.\-]+"),
    "pfam_hits": find_all(PFAM),
    "pfam_id": find_all(PFAM),
    "cazy_hits": cazy_ecs,
    "dbcan_sub_ec": dbcan_ecs,
    "methyl_id": first_word_of_each(","),
}

# Fallbacks for columns not listed above, by name suffix, then for anything else.
SUFFIX_PARSERS = {"_acc": without_version}
DEFAULT_PARSER = whole_cell


def parse_ids(col: str) -> pl.Expr:
    """Expression turning annotation column `col` into a list of identifiers."""
    if col in ID_PARSERS:
        return ID_PARSERS[col](col)
    for suffix, parser in SUFFIX_PARSERS.items():
        if col.endswith(suffix):
            return parser(col)
    return DEFAULT_PARSER(col)


# Extra identifiers built by joining ids from several columns with ";".
# Every id from the first column is paired with every id from the next,
# e.g. dbcan_id "GH5;GH9" + dbcan_sub_ec "3.2.1.4" -> "GH5;3.2.1.4", "GH9;3.2.1.4".
COMPUTED_COLUMNS: Dict[str, Tuple[str, ...]] = {
    "dbcan_id_ec": ("dbcan_id", "dbcan_sub_ec"),
}

# Annotation columns that never hold rule identifiers, so they are not searched.
NON_ID_ANNOTATION_COLUMNS = {
    "genome",
    "scaffold",
    "gene_number",
    "start_position",
    "stop_position",
    "strandedness",
    "rank",
    "heme_regulatory_motif_count",
}


def _match_key(hit: pl.Expr) -> pl.Expr:
    """Key for matching annotation ids to rule ids: uppercase, then drop "EC:"
    at the start or after a ";". ${1} keeps the separator, so GH107;EC:3.2.1.212
    becomes GH107;3.2.1.212 rather than GH1073.2.1.212."""
    return hit.str.to_uppercase().str.replace_all(r"(^|;)EC:", "${1}")


CALL_FUNCTIONS = {
    "not",
    "percent",
    "at_least",
    "tax",
    "column_count_values",
    "column_sum_values",
    "filter_contains",
    "filter_compare",
}

VISUALIZE_FUNCTIONS = {
    "path_steps",
    "path_subunits",
}

# AST node definitions and Grammar definition


@dataclass(frozen=True)
class Expr:
    def validate(self):
        pass

    def pprint(self):
        """Pretty print the expression for debugging."""
        from pprint import pprint

        pprint(self)


@dataclass(frozen=True)
class Name(Expr):
    value: str
    qualifier: str | None = None


@dataclass(frozen=True)
class Alias(Expr):
    value: str


@dataclass(frozen=True)
class Number(Expr):
    value: float


@dataclass(frozen=True)
class String(Expr):
    value: str


@dataclass(frozen=True)
class And(Expr):
    parts: Tuple[Expr, ...]


@dataclass(frozen=True)
class Or(Expr):
    parts: Tuple[Expr, ...]


@dataclass(frozen=True)
class Steps(Expr):
    parts: Tuple[Expr, ...]  # list-of-parts (coverage/count units)


@dataclass(frozen=True)
class Call(Expr):
    value: str
    args: Tuple[Expr, ...]
    allow_visualize_functions: bool = False

    def validate(self):
        # Assume visualization function args are valid since they can
        # take a wide variety of args
        if self.value in VISUALIZE_FUNCTIONS and self.allow_visualize_functions:
            return
        if self.value in VISUALIZE_FUNCTIONS and not self.allow_visualize_functions:
            raise RuleError(
                f"Visualize function {self.value} cannot be used in rules when allow_visualize_functions is False."
                " Check you aren't using a visualization rule function like path_steps() or path_subunits() in your rules,"
                " when doing non-visualization rule evaluation such as traits."
            )

        if self.value not in CALL_FUNCTIONS:
            raise RuleError(f"Unknown function: {self.value}")
        n_args = len(self.args)
        match self.value:
            case "not" | "tax":
                if n_args != 1:
                    raise RuleError(f"{self.value} expects 1 arg. args: {self.args}")
            case "percent":
                if n_args != 2:
                    raise RuleError(
                        f"{self.value}(n, Group) expects 2 args. args: {self.args}"
                    )
            case "at_least":
                if n_args != 3:
                    raise RuleError(
                        f"{self.value}(n, mode, Group) expects 3 args. args: {self.args}"
                    )
            case "column_count_values":
                if n_args != 5:
                    raise RuleError(
                        "column_count_values(column, val_op, val_threshold, count_op, count_thrreshold)"
                        f" expects 5 args. args: {self.args}"
                    )
                for op in (self.args[1], self.args[3]):
                    op = _as_str(op)
                    if op not in ALLOWED_CMPOPS:
                        raise RuleError(
                            f"column_count_values ops must be one of {sorted(ALLOWED_CMPOPS)}; got {op}"
                        )
            case "column_sum_values":
                if n_args != 3:
                    raise RuleError(
                        f"column_sum_values(column, op, threshold) expects 3 args. args{self.args}"
                    )
                op = _as_str(self.args[1])
                if op not in ALLOWED_CMPOPS:
                    raise RuleError(
                        f"column_sum_values op must be one of {sorted(ALLOWED_CMPOPS)}; got {op}"
                    )
            case "filter_contains":
                if n_args != 2:
                    raise RuleError(
                        f"filter_contains(column, value) expects 2. args: {self.args}"
                    )
            case "filter_compare":
                if n_args != 3:
                    raise RuleError(
                        f"filter_compare(column, op, threshold) expects 3 args. args: {self.args}"
                    )
            case _:
                raise RuleError(
                    f"Unable to parse function in rules. Function: {self.value}"
                )


@dataclass(frozen=True)
class Pipe(Expr):
    left: Expr
    right: Expr


@dataclass(frozen=True)
class PipeChain(Expr):
    calls: tuple[Call, ...]  # or Expr, depending on what you permit

    def validate(self):
        for call in self.calls:
            if not isinstance(call, Call):
                raise RuleError(f"Expected Call in PipeChain, got {call}")


def _as_str(e: Expr) -> str:
    try:
        return e.value
    except AttributeError:
        raise RuleError(f"Could not convert node to string, got {e}")


def _as_int(e: Expr) -> int:
    if isinstance(e, Number) and float(e.value).is_integer():
        return int(e.value)
    if isinstance(e, Name) and e.value.isdigit():
        return int(e.value)
    raise RuleError(f"Expected integer literal, got {e}")


def _as_bool(e: Expr) -> bool:
    if isinstance(e, (Name, String)):
        name = e.value.lower()
        if name in {"true", "1", "yes", "t"}:
            return True
        if name in {"false", "0", "no", "f"}:
            return False
    raise RuleError(f"Expected boolean literal, got {e}")


def _as_float(e: Expr) -> float:
    if isinstance(e, Number):
        return float(e.value)
    try:
        return float(e.value)
    except (ValueError, AttributeError):
        raise RuleError(
            f"Could not convert node to gloat. Expected numeric literal, got {e}"
        )


class ASTTransformer(Transformer):
    def __init__(self, allow_visualize_functions=False):
        self.allow_visualize_functions = allow_visualize_functions
        super().__init__()

    def simple_name(self, items):
        return Name(value=str(items[0]), qualifier=None)

    def qualified_name(self, items):
        qualifier, value = items
        return Name(value=str(value), qualifier=str(qualifier))

    def quoted_name(self, items):
        val_index = 0 if (len(items) == 1) else 1
        items[val_index] = (
            items[val_index][1:-1].replace(r"\`", "`").replace(r"\\", "\\")
        )  # remove backticks
        if len(items) == 2:
            return self.qualified_name(items)
        return self.simple_name(items)

    def alias_name(self, items):
        return Alias(value=str(items[0]))

    def quoted_alias_name(self, items):
        value = str(items[0])
        value = value[1:-1].replace(r"\`", "`").replace(r"\\", "\\")
        return Alias(value=value)

    def number(self, items):
        return Number(float(str(items[0])))

    def string(self, items):
        s = str(items[0])
        if len(s) >= 2 and s[0] == s[-1] and s[0] in ("'", '"'):
            s = s[1:-1]
        return String(s)

    def and_(self, items):
        return And(parts=tuple(items))

    def or_(self, items):
        return Or(parts=tuple(items))

    def group(self, items):
        # square brackets are just grouping, not a node
        assert len(items) == 1, "Grouping isn't supposed to change number of items"
        return items[0]

    def step_(self, items):
        return Steps(tuple(items))

    def call(self, items):
        name = str(items[0])
        args = tuple(items[1:])
        return Call(
            name, args, allow_visualize_functions=self.allow_visualize_functions
        )

    def pipe_(self, items):
        left, right = items

        def to_list(x):
            if isinstance(x, PipeChain):
                return list(x.calls)
            return [x]

        calls = tuple(to_list(left) + to_list(right))
        return PipeChain(calls)


@functools.cache
def _get_parser(allow_visualize_functions: bool) -> Lark:
    """Build the LALR parser once per flag value; building the tables isn't free."""
    grammar = (Path(__file__).parent / "rules.lark").read_text()
    return Lark(
        grammar,
        parser="lalr",
        transformer=ASTTransformer(allow_visualize_functions=allow_visualize_functions),
    )


def _parse_rule_expr(parser: Lark, expr_str: str) -> Expr:
    """Parse one rule string, turning Lark errors into friendlier RuleErrors."""
    try:
        return parser.parse(expr_str)
    except LarkError as e:
        n_left_b = expr_str.count("[")
        n_right_b = expr_str.count("]")
        if n_left_b != n_right_b:
            raise RuleError(
                f"Possible mismatched brackets in rule expression: {expr_str}"
            ) from e
        if expr_str.count("&") > 0 and expr_str.count("|") > 0:
            # both used; likely missing brackets
            raise RuleError(
                f"Possible ambiguous use of '&' or '|' without surrounding brackets `[ ]`. Check that you don't have a situation like `A | B & C` or `A & B | C`. These are generally not allowed because they can be ambiguous. These should be written as `[A | B] & C` or `A | [B & C]`: {expr_str}."
            ) from e
        if hasattr(e, "get_context"):
            logger.error(e.get_context(expr_str))
            raise RuleError(
                f"Error parsing rule expression for unknown reason. Possible hints could be found from the lark parsing context: {e.get_context(expr_str)}"
            ) from e
        else:
            raise RuleError(
                f"Error parsing rule expression for unknown reason. Possible hints could be found in the above exceptions from the parsing library Lark: {expr_str}"
            ) from e


def _scan_rules_tsv(path: str | Path) -> pl.LazyFrame:
    """Lazily read a rules TSV, dropping rows that are entirely empty."""
    return pl.scan_csv(
        path,
        separator="\t",
        infer_schema_length=None,
        truncate_ragged_lines=True,
    ).filter(~pl.all_horizontal(pl.all().is_null()))


@dataclass(frozen=True)
class CompiledRules:
    rules: Dict[str, Expr]  # top-level rules with macros expanded
    needed_features: Set[str]
    needed_features_raw: Set[str]
    features_by_rules: Dict[str, Set[str]]
    trees_by_rules: Dict[str, nx.DiGraph]
    df: pl.DataFrame
    # Parsed `rule_options` JSON per rule name; rules without options are absent.
    rule_options: Dict[str, dict] = field(default_factory=dict)

    @classmethod
    def from_rules(
        cls,
        *args,
        label_col: str = "name",
        group_col=None,
        group_col_backups=BACKUP_GROUPBY_COLUMNS,
        **kwargs,
    ) -> CompiledRules:
        definitions, rules, df = load_rules(*args, label_col=label_col, **kwargs)

        defs_expanded = {
            k: expand_macros(v, definitions) for k, v in definitions.items()
        }

        # Expand rules using expanded defs
        # and we still need to add needed features from rules
        features_by_rules = {k: set() for k in rules.keys()}
        trees_by_rules = {k: nx.DiGraph() for k in rules.keys()}
        rules_expanded = {
            k: expand_macros(
                v,
                defs_expanded,
                needed_features=features_by_rules[k],
                graph=trees_by_rules[k],
            )
            for k, v in rules.items()
        }
        needed_features = set().union(*[s for s in features_by_rules.values()])

        # Pick the group column: the default if present, else the first backup present.
        # (Don't mutate group_col_backups: its default is the shared module-level list.)
        if group_col is None:
            group_col = next(
                (c for c in [DEFAULT_GROUPBY_COLUMN, *group_col_backups] if c in df),
                None,
            )
        if group_col:
            # Fill nulls in the group column from the remaining backups, in order.
            backups = [c for c in group_col_backups if c in df and c != group_col]
            if backups:
                df = df.with_columns(pl.coalesce(group_col, *backups).alias(group_col))

        return cls(
            rules=rules_expanded,
            needed_features={item.upper() for item in needed_features},
            needed_features_raw=needed_features,
            features_by_rules=features_by_rules,
            trees_by_rules=trees_by_rules,
            df=df,
            rule_options=_parse_rule_options(df, label_col),
        )


def _parse_rule_options(df: pl.DataFrame, label_col: str) -> Dict[str, dict]:
    """Read the optional `rule_options` column (JSON per rule) into a dict by rule name."""
    if "rule_options" not in df.columns:
        return {}
    options = {}
    rows = df.select(label_col, "rule_options").drop_nulls()
    for name, raw in zip(rows[label_col], rows["rule_options"]):
        try:
            options[name] = json.loads(raw)
        except json.JSONDecodeError as e:
            raise RuleError(
                f"Invalid rule_options JSON for rule '{name}': {raw}"
            ) from e
    return options


def load_rules(
    rules_path: str | Path = None,
    rules: pl.LazyFrame = None,
    label_col: str = "name",
    alias_col: str = "alias",
    rules_col: str = "rule",
    allow_visualize_functions: bool = False,
    common_rules: pl.LazyFrame = None,
    common_rules_path: str | Path = None,
    check_for_dupes: bool = True,
) -> Tuple[Dict[str, Expr], Dict[str, Expr], pl.DataFrame]:
    """
    Assumes TSV has columns at least: name, rule
    Convention used here:
      - if `name` is non-empty: this is an OUTPUT RULE whose expression is in `rule`
      - if `name` is empty and `alias` is non-empty: this is a DEFINITION macro: alias := rule
    If your file uses a slightly different convention, adjust this function only.

    Returns (definitions, rules, df). `df` is the cleaned rules table; its rule
    column holds the raw rule strings, the parsed expressions are in the dicts.
    """
    assert (rules_path is not None) != (rules is not None), (
        "Either rules_path or rules DataFrame must be provided, but not both."
    )
    lf = _scan_rules_tsv(rules_path) if rules_path else rules

    if common_rules_path:
        clf = _scan_rules_tsv(common_rules_path)
    else:
        clf = common_rules

    cols = set(lf.collect_schema().names())
    required = {label_col, rules_col}
    if not required.issubset(cols):
        raise ValueError(
            f"rules TSV missing required columns {required}. Found: {sorted(cols)}"
        )

    if not (alias_col and alias_col in cols):
        alias_col = alias_col or "alias"
        lf = lf.with_columns(pl.lit(None, dtype=pl.String).alias(alias_col))

    if clf is not None:
        # Common rules always use name/alias/rule; map them onto this sheet's columns.
        clf = clf.with_columns(
            pl.col("rule").alias(rules_col),
            pl.col("name").alias(label_col),
            pl.col("alias").alias(alias_col),
        )
        lf = pl.concat([lf, clf], how="diagonal_relaxed")

    # Clean up text lazily, then materialize exactly once. Every collect() re-runs the
    # whole plan (CSV scan, concat, ...), so we build everything below from this one frame.
    df = (
        lf.with_columns(pl.col(label_col, alias_col, rules_col).cast(pl.String))
        # Strip stray whitespace users may have left around labels, aliases and rules,
        .with_columns(pl.col(label_col, alias_col, rules_col).str.strip_chars())
        # then treat empty cells (including whitespace-only ones) as missing.
        .with_columns(pl.col(pl.String).replace("", None))
        .collect()
    )

    if check_for_dupes:
        labels = df.get_column(label_col).drop_nulls()
        duplicate_labels = labels.filter(labels.is_duplicated()).unique().sort()
        if len(duplicate_labels):
            raise RuleError(
                f"Rule column '{label_col}' needs to be unique."
                f" Duplicate labels: {duplicate_labels.to_list()}"
            )

    # Parse each distinct rule string once (aliases are often reused across sheets).
    parser = _get_parser(allow_visualize_functions)
    parsed = {
        s: _parse_rule_expr(parser, s)
        for s in df.get_column(rules_col).drop_nulls().unique()
    }

    def rules_keyed_by(key_col: str) -> Dict[str, Expr]:
        keyed = df.select(key_col, rules_col).drop_nulls(key_col)
        return {
            key: parsed.get(rule) for key, rule in zip(keyed[key_col], keyed[rules_col])
        }

    return rules_keyed_by(alias_col), rules_keyed_by(label_col), df


def expand_macros(
    expr: Expr,
    definitions: Dict[str, Expr],
    needed_features: Set[str] = None,
    graph: Optional[nx.DiGraph] = None,
) -> Expr:
    """Expand recursively macros in expr using definitions"""
    # TODO: make the expanded aliases get parsed as a group
    memo: Dict[
        tuple[Expr, bool],
        tuple[Expr, set[Expr], set[Expr]],
    ] = {}
    stack: list[str] = []
    skip_needed = False
    if needed_features is None:
        skip_needed = True

    def recurse(expr: Expr, add_name_to_needed: bool = True) -> Expr:
        memo_key = (expr, add_name_to_needed)
        if memo_key in memo:
            return memo[memo_key]

        entries, exits = None, None

        # We first check if it is an Alias to extract out
        if isinstance(expr, Alias):
            name = expr.value
            if name not in definitions:
                raise RuleError(f"Undefined alias: @{name}")
            if name in stack:
                raise RuleError(f"Cycle detected: {' -> '.join(stack + [name])}")
            stack.append(name)
            out, entries, exits = recurse(
                definitions[name],
                add_name_to_needed=add_name_to_needed,
            )
            stack.pop()
        elif isinstance(expr, Name):
            name = expr.value
            if name in definitions:
                raise RuleError(
                    "Identity found in rule that matches alias name without alias syntax."
                    f"Name in question: `{name}`"
                    "\nGene ids in rules and aliases need to have different and unique names."
                    " Check if you meant to reference an alias with an `@` symbol."
                    " Or change your alias name to be unique from your gene id."
                )
            out = expr
            # This will slightly
            if add_name_to_needed and not skip_needed:
                needed_features.add(name)
            entries, exits = {expr}, {expr}
        elif isinstance(expr, (Number, String)):
            out = expr
        elif isinstance(expr, And):
            is_atomic = all(not isinstance(p, Steps) for p in expr.parts)
            and_entries, and_exits = set(), set()
            if is_atomic:
                and_entries, and_exits = {expr}, {expr}
            outs = []
            for part in expr.parts:
                o, entries, exits = recurse(part)
                if not is_atomic:
                    and_entries |= set(entries)
                    and_exits |= set(exits)
                outs.append(o)
            entries, exits = and_entries, and_exits

            out = And(parts=tuple(outs))
        elif isinstance(expr, Or):
            is_atomic = all(not isinstance(p, Steps) for p in expr.parts)
            or_entries, or_exits = set(), set()
            if is_atomic:
                or_entries, or_exits = {expr}, {expr}
            outs = []
            for part in expr.parts:
                o, entries, exits = recurse(part)
                if not is_atomic:
                    or_entries |= set(entries)
                    or_exits |= set(exits)
                outs.append(o)
            entries, exits = or_entries, or_exits

            out = Or(parts=tuple(outs))
        elif isinstance(expr, Steps):
            prev_exits = None
            all_entries = None
            all_exits = None
            parts = []
            for part in expr.parts:
                out, entries, exits = recurse(part)
                parts.append(out)
                if all_entries is None:
                    all_entries = set(entries)

                if graph is not None and prev_exits is not None:
                    # connect every previous exit to every current entry
                    for u in prev_exits:
                        for v in entries:
                            graph.add_edge(u, v, kind="seq")

                prev_exits = set(exits)
                all_exits = set(exits)
            entries, exits = all_entries, all_exits
            out = Steps(tuple(parts))
        elif isinstance(expr, Call):
            fn, args = expr.value, expr.args
            call_rec = []
            # Determine which args should have their Names added to needed_features
            # since some call args are just numbers, ops, and columns.
            # needed_features are the genes we need to look for in the annotations.
            for i, arg in enumerate(args):
                if fn in {"at_least"}:
                    add_name = i == 2
                elif fn == "percent":
                    add_name = 1 == i
                elif fn == "not":
                    add_name = 0 == i
                elif fn in {
                    "column_count_values",
                    "column_sum_values",
                    "filter_contains",
                    "filter_compare",
                }:
                    add_name = False
                # unknown function; be conservative. Worse case is we count extra needed features
                # such as numbers or ops, which is not ideal, but quite unlikely to be
                # a gene id. So almost always harmless.
                else:
                    add_name = True
                o, _, _ = recurse(arg, add_name_to_needed=add_name)
                call_rec.append(o)
            out = Call(
                expr.value,
                tuple(call_rec),
                allow_visualize_functions=expr.allow_visualize_functions,
            )
        elif isinstance(expr, PipeChain):
            calls = list(expr.calls)
            for i, call in enumerate(expr.calls):
                o, _, _ = recurse(call)
                calls[i] = o
            out = PipeChain(calls=tuple(calls))
        else:
            raise TypeError(expr)
        # Some nodes need validation and have to be done after rules are expanded
        try:
            out.validate()
        except Exception:
            print(
                f"Error validating expression after macro expansion. Expression: {out}"
            )
            print(expr)
            raise
        ret = out, entries or set(), exits or set()
        memo[memo_key] = ret
        return ret

    return recurse(expr)[0]


def prepare_present_map_df(
    df: pl.DataFrame,
    count_col: str,
    besthit_cols: List[str],
    needed_features: Set[str],
    additional_cols: List[str] = None,
    group_col: str = None,
):
    """Long table of (index cols, hit, additional cols): one row per gene and
    rule identifier it was annotated with, keeping only identifiers in
    `needed_features`.

    Each annotation column becomes its own lazy branch (split into ids, then
    exploded), so columns are processed in parallel and never multiply each
    other's rows. Everything is collected once at the end.
    """
    additional_cols = additional_cols or []
    index_cols = [count_col, group_col] if group_col else [count_col]
    keep_cols = [*index_cols, *additional_cols]
    hit_col = "hit"

    hit_cols = [col for col in besthit_cols if col in df.columns]
    hit_cols += [
        col
        for col in df.columns
        if col
        and col not in hit_cols
        and col not in keep_cols
        and col not in NON_ID_ANNOTATION_COLUMNS
    ]

    lf = df.lazy()

    def hits_from(*cols: str) -> pl.LazyFrame:
        """Explode each column's ids; with several columns, pair every id of one with every id of the next."""
        out = lf.select(*keep_cols, *(parse_ids(col).alias(col) for col in cols))
        for col in cols:
            out = out.explode(col, empty_as_null=True)
        if len(cols) == 1:
            return out.rename({cols[0]: hit_col})
        return out.select(*keep_cols, pl.concat_str(cols, separator=";").alias(hit_col))

    branches = [hits_from(col) for col in hit_cols]
    branches += [
        hits_from(*sources)
        for sources in COMPUTED_COLUMNS.values()
        if all(col in hit_cols for col in sources)
    ]

    # Rule identifiers, keyed the same way as annotation identifiers below.
    needed = pl.LazyFrame(
        {hit_col: sorted(needed_features)}, schema={hit_col: pl.String}
    ).with_columns(_match_id=_match_key(pl.col(hit_col)))

    return (
        pl.concat(branches, how="vertical_relaxed")
        .drop_nulls()
        .with_columns(_match_id=_match_key(pl.col(hit_col)))
        .drop(hit_col)
        # Keep only identifiers the rules ask for. The join restores `hit` using
        # the rule's spelling, including any EC: prefix.
        .join(needed, on="_match_id", how="inner")
        .select(*index_cols, hit_col, *additional_cols)
        # Multiple columns can match the same rule identifier; count each gene once.
        .unique()
        .collect()
    )


def build_present_map(
    df: pl.DataFrame,
    sample_col: str,
    needed_features: Set[str],
    sample_names: list[str],
    value_col: str = None,
) -> Tuple[List[str], Dict[str, np.ndarray]]:
    """Build present_map of needed gene_ids from annotations DataFrame"""
    hit_col = "hit"

    df: pl.DataFrame = df.drop_nulls(subset=sample_col)

    sample_index = {s: i for i, s in enumerate(sample_names)}
    n = len(sample_names)
    present_map: Dict[str, np.ndarray] = {f: np.zeros(n) for f in needed_features}

    # Group by hit_id and set booleans
    for hit_id, sub in df.group_by(hit_col):
        hit = hit_id[0]
        arr = present_map.get(hit)
        if arr is None:
            continue

        if value_col:
            # for name, value in sub.select([sample_col, value_col]).iter_rows():
            for name, value in (
                sub.group_by([sample_col]).agg(pl.col(value_col).sum()).iter_rows()
            ):
                arr[sample_index[name]] = value
        else:
            for name, count in sub[sample_col].value_counts().iter_rows():
                arr[sample_index[name]] = count
    return present_map


class Evaluator:
    def __init__(
        self,
        samples: List[str],
        present_map: Dict[str, np.ndarray],
        sample_col: str,
        annotations: Optional[pl.DataFrame],
    ):
        self.samples = samples
        self.present_map = present_map
        self.annotations = annotations
        self.sample_col = sample_col
        self._order_df = pl.DataFrame(
            {
                self.sample_col: self.samples,
                "_order": range(len(self.samples)),
            }
        )

        self._all_false = np.zeros(len(samples), dtype=bool)

        # Per-instance caches. `functools.cache` on the methods themselves would be
        # one cache shared by the whole class, keeping every Evaluator (and its
        # annotations DataFrame) alive for the life of the process.
        self.eval_bool = functools.cache(self._eval_bool)
        self.eval_cycle = functools.cache(self._eval_cycle)

    def _eval_bool(self, expr: Expr, reduce_outer_and=True) -> np.ndarray:
        out = None
        if isinstance(expr, Alias):
            raise RuleError(
                "Alias was not resolved properly in parsing. Please file a github bug report."
            )

        if isinstance(expr, Name):
            out = self.present_map.get(expr.value.upper(), self._all_false)

        if isinstance(expr, And):
            out = np.stack([self.eval_bool(part) for part in expr.parts], axis=1)
            if reduce_outer_and:
                out = np.all(out, axis=1)

        if isinstance(expr, Or):
            try:
                out = np.any(
                    np.stack([self.eval_bool(part) for part in expr.parts], axis=1),
                    axis=1,
                )
            except ValueError as e:
                raise RuleError(f"Error evaluating OR expression: {expr}") from e

        if isinstance(expr, PipeChain):
            kwargs = {"df": self.annotations, "masks": []}
            for call in expr.calls:
                out = self.eval_call(call, **kwargs)
                if isinstance(out, (pl.DataFrame, pl.LazyFrame)):
                    kwargs["df"] = out
                if isinstance(out, pl.Expr):
                    kwargs["masks"].append(out)

        if isinstance(expr, Call):
            out = self.eval_call(expr)

        if isinstance(expr, Number):
            # numeric alone isn't boolean; treat as error
            raise RuleError(f"Literal used where boolean expected: {expr}")

        if isinstance(expr, Steps):
            out = self.eval_cycle(expr)

        if out is not None:
            return out

        raise RuleError(
            "Something failed to parse properly."
            f" Tried to evaluate truth value of {expr}"
        )

    def _eval_cycle(
        self, expr: Steps | list[Expr], simplify=True, **kwargs
    ) -> np.ndarray:
        """Evaluate each part of `expr` into one column per part.

        Gene names keep their hit counts (not just True/False) so that
        `at_least(n, COUNTS, ...)` can count; callers that need presence
        convert with `.astype(bool)`. With `simplify`, reduce to "all parts present".
        """
        if isinstance(expr, Iterable):
            parts = [self.eval_bool(p, **kwargs) for p in expr]
        elif isinstance(expr, Steps):
            parts = [self.eval_bool(p, **kwargs) for p in expr.parts]
        else:
            raise RuleError(
                f"Expected Steps or Iterable of Expr for cycle evaluation, got {expr}"
            )
        mat = (
            np.stack(parts, axis=1)
            if parts
            else np.zeros((len(self.samples), 0), dtype=bool)
        )
        if simplify:
            mat = mat.all(axis=1)
        return mat

    def len_cycle(self, expr: Steps) -> int:
        return len(expr.parts)

    def eval_call(
        self, call: Call, df: pl.DataFrame = None, masks: list[pl.Expr] = None
    ) -> np.ndarray:
        fn = call.value
        args = call.args
        kwargs = {}
        if df is not None:
            kwargs["df"] = df
        if masks:
            kwargs["masks"] = masks

        match fn:
            case "not":
                # Negates its own argument. Inside a pipe, a negated filter is
                # combined (AND) with any filters before it.
                return self.not_(self.eval_bool(args[0]))
            case "percent":
                expr = args[1]
                if isinstance(expr, Or):
                    expr = Steps(parts=expr.parts)
                return self.percent(
                    _as_int(args[0]), self.eval_cycle(expr, simplify=False)
                )
            case "at_least":
                expr = args[2]
                if isinstance(expr, Or):
                    expr = Steps(parts=expr.parts)
                return self.at_least(
                    _as_int(args[0]),
                    _as_str(args[1]),
                    self.eval_cycle(expr, simplify=False),
                )
            case "tax":
                return self.tax(_as_str(args[0]), **kwargs)
            case "column_count_values":
                return self.column_count_values(
                    col=_as_str(args[0]),
                    val_op=_as_str(args[1]),
                    val_thr=_as_float(args[2]),
                    count_op=_as_str(args[3]),
                    count_thr=_as_float(args[4]),
                    **kwargs,
                )
            case "column_sum_values":
                return self.column_sum_values(
                    col=_as_str(args[0]),
                    op=_as_str(args[1]),
                    thr=_as_float(args[2]),
                    **kwargs,
                )
            case "filter_contains":
                return self.filter_contains(
                    col=_as_str(args[0]), val=_as_str(args[1]), **kwargs
                )
            case "filter_compare":
                return self.filter_compare(
                    col=_as_str(args[0]),
                    op=_as_str(args[1]),
                    thr=_as_float(args[2]),
                    **kwargs,
                )
            case "path_steps":
                return self.path_steps(args)

            case "path_subunits":
                return self.path_subunits(args)
            case _:
                raise RuleError(f"Unable to parse function in rules. Function: {fn}")

    def _per_sample(self, df: pl.DataFrame, col: str) -> np.ndarray:
        """Boolean result column `col` of a per-sample `df`, as an array in
        `self.samples` order. Samples missing from `df` are False."""
        return (
            self._order_df.join(df, on=self.sample_col, how="left")
            .sort("_order")
            .get_column(col)
            .fill_null(False)
            .to_numpy()
        )

    def eval_filter_dec(func):
        """Decorator to evaluate filter functions with optional df and masks
        Applies masks to the specified column before calling the function."""

        @functools.wraps(func)
        def wrapper(
            self,
            col,
            *args,
            df: pl.DataFrame = None,
            masks: list[pl.Expr] = None,
            **kwargs,
        ):
            df = self.annotations if df is None else df
            masks = masks or []
            if masks:
                df = df.with_columns(
                    pl.when(pl.lit(True).and_(*masks))
                    .then(pl.col(col))
                    .otherwise(pl.lit(None))
                )
            return func(self, col, *args, df=df, **kwargs)

        return wrapper

    # Call functions
    @staticmethod
    def not_(x: np.ndarray | pl.Expr) -> np.ndarray | pl.Expr:
        """Negate a per-sample result, or a row filter (from filter_* functions)."""
        if isinstance(x, pl.Expr):
            return ~x
        return np.logical_not(x)

    @staticmethod
    def percent(n: int, x: np.ndarray) -> np.ndarray:
        x = x.astype(bool)
        thr = n / 100.0
        if x.shape[1] == 0:
            return np.zeros(x.shape[0], dtype=bool)
        cov = x.mean(axis=1)
        return cov >= thr

    @staticmethod
    def at_least(k: int, mode: str, x: np.ndarray) -> np.ndarray:
        """PRESENCE: at least k of the parts are present.
        COUNTS: the parts' hit counts add up to at least k."""
        mode = mode.upper()
        if mode not in AT_LEAST_MODES:
            raise ValueError(
                f"Unsupported at_least mode: {mode}. Use one of {AT_LEAST_MODES}"
            )
        if mode == "PRESENCE":
            x = x.astype(bool)
        return x.sum(axis=1) >= k

    def tax(self, tax_label, *, df: pl.DataFrame = None):
        df = self.annotations if df is None else df
        if "taxonomy" not in df.columns:
            return np.zeros(len(self.samples), dtype=bool)
        return self._per_sample(
            df.group_by(self.sample_col).agg(
                pl.col("taxonomy").str.contains(tax_label).any()
            ),
            "taxonomy",
        )

    @eval_filter_dec
    def column_count_values(
        self,
        col: str,
        val_op: str,
        val_thr: float,
        count_op: str,
        count_thr: float,
        *,
        df: pl.DataFrame,
    ) -> np.ndarray:
        df = self.annotations if df is None else df

        if col not in df.columns:
            return np.zeros(len(self.samples), dtype=bool)
        if val_op not in OP_TO_EXPR:
            raise ValueError(
                f"Unsupported value operation={val_op!r} to compare values for column {col}. Use one of {sorted(OP_TO_EXPR)}"
            )
        val_cmp_fn = OP_TO_EXPR[val_op]
        if count_op not in OP_TO_EXPR:
            raise ValueError(
                f"Unsupported count operation={count_op!r} to compare value counts for column {col}. Use one of {sorted(OP_TO_EXPR)}"
            )
        count_cmp_fn = OP_TO_EXPR[count_op]

        df = (
            df.group_by(self.sample_col)
            # first do the val cmp fn on the column values to get a bool of which rows
            # pass the value threshold (e.g., per row: col_val >= val_thr)
            # Then sum those booleans to get a count of how many rows per sample pass
            .agg(
                count_cmp_fn(
                    val_cmp_fn(pl.col(col).cast(float), val_thr).sum(), count_thr
                )
            )
        )
        return self._per_sample(df, col)

    @eval_filter_dec
    def column_sum_values(
        self, col: str, op: str, thr: float, *, df: pl.DataFrame = None
    ) -> np.ndarray:
        df = self.annotations if df is None else df

        if col not in df.columns:
            return np.zeros(len(self.samples), dtype=bool)
        try:
            cmp_fn = OP_TO_EXPR[op]
        except KeyError:
            raise ValueError(f"Unsupported op={op!r}. Use one of {sorted(OP_TO_EXPR)}")

        df = df.group_by(self.sample_col).agg(
            cmp_fn(pl.col(col).cast(float).sum(), thr)
        )
        return self._per_sample(df, col)

    @staticmethod
    def _require_column(df: pl.DataFrame, col: str, fn: str) -> None:
        if col not in df.columns:
            raise RuleError(
                f"Missing column '{col}' for {fn}(): the annotations don't have it."
                " If the rule should be skipped when the column is missing, give it"
                f' rule_options {{"include_when":{{"all_columns_present":["{col}"]}}}}'
            )

    def filter_contains(
        self, col: str, val: str, df: pl.DataFrame = None, **kwargs
    ) -> pl.Expr:
        df = self.annotations if df is None else df
        self._require_column(df, col, "filter_contains")
        return pl.col(col).str.contains(val)

    def filter_compare(
        self, col: str, op: str, thr: float, df: pl.DataFrame = None, **kwargs
    ) -> pl.Expr:
        df = self.annotations if df is None else df
        self._require_column(df, col, "filter_compare")
        try:
            cmp_fn = OP_TO_EXPR[op]
        except KeyError:
            raise ValueError(f"Unsupported op={op!r}. Use one of {sorted(OP_TO_EXPR)}")
        return cmp_fn(pl.col(col).cast(float), thr)

    def path_steps(self, expr: Steps | list[Expr], **kwargs) -> pl.DataFrame:
        # Count steps present, not gene hits within them.
        cycle = self.eval_cycle(expr, simplify=False, **kwargs).astype(bool)

        df = pl.DataFrame(
            {
                "steps": cycle.shape[-1],
                # sum across all axis except the first one. We don't know the number of dim
                # since path_subunits doesn't reduce the and
                "steps_present": cycle.sum(
                    axis=tuple(i for i in range(len(cycle.shape)) if i > 0)
                ),
                self.sample_col: self.samples,
            }
        )
        df = df.with_columns(
            coverage_percentage=pl.col("steps_present") / pl.col("steps")
        )
        return df

    def path_subunits(self, expr: Steps | list[Expr]) -> pl.DataFrame:
        df = self.path_steps(expr, reduce_outer_and=False)
        df = df.rename({"steps": "subunits", "steps_present": "subunits_present"})
        return df


def should_include(options: dict, available_columns: set[str]) -> bool:
    condition = options.get("include_when", {})

    present = set(condition.get("all_columns_present", []))
    absent = set(condition.get("all_columns_absent", []))

    return present <= available_columns and absent.isdisjoint(available_columns)


def get_val_cols(df, allow_multiple=False):
    cols = []
    if "coverage_percentage" in df.columns:
        cols.append("coverage_percentage")
    if "present" in df.columns:
        cols.append("present")
    if "value" in df.columns:
        cols.append("value")
    if not cols:
        raise RuleError(
            f"Rules parsing problem, can't determine check_col for dataframe rules: {df}. File github issue"
        )
    # it doesn't matter how many in cols if we allow_multiple and len > 0, so we return
    if allow_multiple:
        return cols
    # if we don't allow multiple, and we have multiple, we raise RuleError
    if len(cols) > 1:
        raise RuleError(
            "In rule sheet, all rules inside a single group must use the same evaluation style."
            "\nThat is, if any in a group uses `path_subunits`, all must use it"
            "\n same for `path_steps`. Or none can use them."
            "You can't mix and match inside a single `group`."
        )
    # if we expect a single return, we extract it and return
    return cols[0]


def evaluate_rules(
    compiled: CompiledRules,
    samples: List[str],
    present_map: Dict[str, np.ndarray],
    annotations: Optional[pl.DataFrame] = None,
    sample_col: Optional[str] = None,
    bool_output: bool = True,
) -> pl.DataFrame:
    ev = Evaluator(
        samples=samples,
        present_map=present_map,
        sample_col=sample_col,
        annotations=annotations,
    )

    # Skip rules whose rule_options say they don't apply to these annotations.
    available_columns = set(annotations.columns) if annotations is not None else set()
    rule_names = sorted(
        rn
        for rn in compiled.rules
        if should_include(compiled.rule_options.get(rn, {}), available_columns)
    )
    data = {}
    for rn in rule_names:
        out = ev.eval_bool(compiled.rules[rn])
        if bool_output:
            out = out.astype(bool)
        data[rn] = out

    df = pl.DataFrame(data)
    df = df.with_columns(pl.Series(sample_col, samples)).select(
        [sample_col] + rule_names
    )

    return df


def evaluate_cycles(
    compiled: CompiledRules,
    samples: List[str],
    present_map: Dict[str, np.ndarray],
    annotations: Optional[pl.DataFrame],
    group_col: Optional[str] = "group",
    label_col: Optional[str] = "name",
    sample_col: Optional[str] = None,
    additional_cols: Optional[List[str]] = None,
    anno_df: pl.DataFrame = None,
    value_col: str = None,
) -> pl.DataFrame:
    assert (anno_df is None and value_col is None) or (
        anno_df is not None and value_col is not None
    ), "anno_df and value_col must both be None or neither be None"
    ev = Evaluator(
        samples=samples,
        present_map=present_map,
        sample_col=sample_col,
        annotations=annotations,
    )
    # if group_col in set(compiled.lf.collect_schema().names()):
    dfs = {}
    additional_cols = [pl.col(c) for c in additional_cols] if additional_cols else []
    if compiled.df.drop_nulls(subset=pl.col(label_col))[group_col].null_count() > 0:
        raise RuleError(
            "Group column, and backup group columns contain no entries! Make sure each row for the group column contains values."
            f"Rule names with null group columns: {', '.join(val for val in compiled.df.drop_nulls(subset=pl.col(label_col)).filter(pl.col(group_col).is_null())[label_col])}"
        )
    for group, frame in compiled.df.drop_nulls(subset=[label_col]).group_by(
        group_col, maintain_order=True
    ):
        group = group[0]
        dfs[group] = []
        for rn in frame.select(pl.col(label_col)).to_series():
            if rn is None:
                continue
            # for rn, expr in compiled.rules.items():

            if "rule_options" in frame.columns:
                rule_options = (
                    frame.filter(pl.col(label_col) == rn).item(0, "rule_options")
                    or "{}"
                )
                if not should_include(
                    json.loads(rule_options), set(annotations.columns)
                ):
                    continue

            expr = compiled.rules[rn]
            out = ev.eval_bool(expr)
            if isinstance(out, np.ndarray):
                out = pl.DataFrame(
                    {"present": out.astype(bool), sample_col: ev.samples}
                )
            out = out.with_columns(pl.lit(rn).alias(label_col))
            features = sorted({f.upper() for f in compiled.features_by_rules[rn]})
            hits = [[] for _ in samples]

            for feature in features:
                counts = present_map.get(feature)
                if counts is not None:
                    for i in np.flatnonzero(counts):
                        hits[i].append(feature)

            # ensure they are sorted the same way by sample col in case eval_bool returns sample_col in differnt ordering
            hit_df = pl.DataFrame(
                {
                    sample_col: samples,
                    "gene_ids_present": [", ".join(ids) for ids in hits],
                }
            )
            out = out.join(
                hit_df, on=sample_col, how="left", validate="1:1", maintain_order="left"
            )
            if anno_df is not None:
                check_col = get_val_cols(out)

                out = (
                    out.join(
                        # We grab the the anno_df that only match the hits
                        # then we group by the sample col
                        # and se sum them over the value column
                        (
                            anno_df.filter(
                                pl.col("hit").is_in(compiled.features_by_rules[rn])
                            )
                            .group_by([sample_col])
                            .agg(pl.col(value_col).sum())
                        ),
                        on=sample_col,
                        how="left",
                        validate="1:1",
                    )
                    # we want to scale the summed value columns by the check column (presence/absence or coverage %)
                    .with_columns(
                        pl.col(value_col).fill_null(strategy="zero") * pl.col(check_col)
                    )
                )
            dfs[group].append(out)

        assert all(isinstance(df, pl.DataFrame) for df in dfs[group]), (
            f"All rules in group {group} should evaluate to the same type, but got different types: {[type(df) for df in dfs[group]]}"
        )
        if len(dfs[group]):
            df = pl.concat(dfs[group])
            df = df.join(
                compiled.df.select(
                    pl.col(label_col), pl.col(group_col), *additional_cols
                ),
                on=label_col,
            )
            dfs[group] = df
        else:
            del dfs[group]

    return dfs


def evaluate_rules_on_anno(
    count_col: str,
    annotations_path: os.PathLike = None,
    annotations: pl.DataFrame = None,
    group_col: str | None = None,
    *args,
    **kwargs,
):
    compiled = CompiledRules.from_rules(*args, check_for_dupes=False, **kwargs)
    print(f"Need {len(compiled.needed_features)} features")
    assert (annotations_path is not None) != (annotations is not None), (
        "Exactly one of annotations_path or annotations DataFrame must be provided."
    )
    if annotations_path:
        try:
            annotations = pl.read_csv(
                annotations_path, separator="\t", infer_schema_length=10_000
            )
        except Exception:
            annotations = pl.read_csv(
                annotations_path, separator="\t", infer_schema_length=None
            )

    sample_col = group_col if group_col else count_col
    sample_names = (
        annotations.select(sample_col).unique().sort(sample_col).to_series().to_list()
    )
    df = prepare_present_map_df(
        df=annotations,
        count_col=count_col,
        besthit_cols=list(ID_PARSERS),
        needed_features=compiled.needed_features,
        group_col=group_col,
    )

    present_map = build_present_map(
        df=df,
        sample_col=sample_col,
        needed_features=compiled.needed_features,
        sample_names=sample_names,
    )

    df = evaluate_rules(
        compiled,
        sample_names,
        present_map,
        annotations=annotations,
        sample_col=sample_col,
    )
    return df
