from __future__ import annotations

import logging
from dataclasses import dataclass
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


ID_EXPR_DICT = {
    "camper_id": pl.col("camper_id").cast(pl.Utf8).cast(pl.List(pl.Utf8)),
    "fegenie_id": pl.col("fegenie_id").cast(pl.Utf8).cast(pl.List(pl.Utf8)),
    "sulfur_id": pl.col("sulfur_id").cast(pl.Utf8).cast(pl.List(pl.Utf8)),
    "kegg_genes_id": pl.col("kegg_genes_id").cast(pl.Utf8).cast(pl.List(pl.Utf8)),
    "kofam_genes_id": pl.col("kofam_genes_id").cast(pl.Utf8).cast(pl.List(pl.Utf8)),
    "ko_id": pl.col("ko_id").str.split(","),
    "kegg_id": pl.col("kegg_id").str.split("/"),
    "kofam_id": pl.col("kofam_id").str.split("/"),
    "peptidase_family": pl.col("peptidase_family").str.split(";"),
    "merops_family": pl.col("merops_family").str.split(";"),
    "camper_EC": pl.col("camper_EC").str.extract_all(r"\[EC:(\d+\.\d+\.\d+\.\d+)\]"),
    "kegg_EC": pl.col("kegg_EC").str.extract_all(r"\[EC:(\d+\.\d+\.\d+\.\d+)\]"),
    "kofam_EC": pl.col("kofam_EC").str.extract_all(r"\[EC:(\d+\.\d+\.\d+\.\d+)\]"),
    "pfam_hits": pl.col("pfam_hits").str.extract_all(r"\[(PF\d{5})\.\d+\]"),
    "pfam_id": pl.col("pfam_id").str.extract_all(r"\[(PF\d{5})\.\d+\]"),
    "cazy_hits": (
        pl.col("cazy_hits")
        .cast(pl.Utf8)
        # capture just the EC number part inside "(EC ...)"
        .str.extract_all(r"\(EC\s+([\d\.\-]+)\)")
        # prefix each captured number with "EC:"
        .list.eval(pl.concat_str([pl.lit("EC:"), pl.element()]))
    ),
    "cazy_subfam_ec": (
        pl.col("cazy_subfam_ec")
        .cast(pl.Utf8)
        .str.extract_all(r"([\d\.\-]+)")
        .list.eval(pl.concat_str([pl.lit("EC:"), pl.element()]))
    ),
    "dbcan_sub_ec": (
        pl.col("dbcan_sub_ec")
        .cast(pl.Utf8)
        .replace({"-": None})
        .str.split(";")
        .list.eval(
            pl.concat_str([pl.lit("EC:"), pl.element().str.split(":").list.first()])
        )
    ),
    "dbcan_id": pl.col("dbcan_id").cast(pl.Utf8).str.split(";"),
    "dbcan_sub_id": pl.col("dbcan_sub_id").cast(pl.Utf8).str.split(";"),
    "dbcan_sub_substrate": pl.col("dbcan_sub_substrate").cast(pl.Utf8).str.split(";"),
    "methyl_id": (
        pl.col("methyl_id")
        .cast(pl.Utf8)
        .str.split(",")
        .list.eval(pl.element().str.strip_chars().str.split(" ").list.first())
    ),
    "DEFAULT": lambda col: pl.col(col).cast(pl.Utf8).cast(pl.List(pl.Utf8)),
}

COMPUTED_COLUMNS = {
    "dbcan_id_ec": pl.concat_str(
        [pl.col("dbcan_id"), pl.col("dbcan_sub_ec")], separator=";"
    )
}

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


@dataclass(frozen=True)
class CompiledRules:
    rules: Dict[str, Expr]  # top-level rules with macros expanded
    needed_features: Set[str]
    needed_features_raw: Set[str]
    features_by_rules: Dict[str, Set[str]]
    trees_by_rules: Dict[str, nx.DiGraph]
    df: pl.DataFrame

    @classmethod
    def from_rules(cls, *args, **kwargs) -> CompiledRules:
        definitions, rules, lf = load_rules(*args, **kwargs)

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

        return cls(
            rules=rules_expanded,
            needed_features={item.upper() for item in needed_features},
            needed_features_raw=needed_features,
            features_by_rules=features_by_rules,
            trees_by_rules=trees_by_rules,
            df=lf.collect(),
        )


def load_rules(
    rules_path: str | Path = None,
    rules: pl.LazyFrame = None,
    label_col: str = "name",
    alias_col: str = "alias",
    rules_col: str = "rule",
    allow_visualize_functions: bool = False,
    common_rules: pl.LazyFrame = None,
    common_rules_path: str | Path = None,
) -> Tuple[Dict[str, Expr], Dict[str, Expr]]:
    """
    Assumes TSV has columns at least: name, rule
    Convention used here:
      - if `name` is non-empty: this is an OUTPUT RULE whose expression is in `rule`
      - if `name` is empty and `alias` is non-empty: this is a DEFINITION macro: alias := rule
    If your file uses a slightly different convention, adjust this function only.
    """
    assert (rules_path is not None) != (rules is not None), (
        "Either rules_path or rules DataFrame must be provided, but not both."
    )
    if rules_path:
        lf = (
            pl.scan_csv(
                rules_path,
                separator="\t",
                infer_schema_length=None,
                truncate_ragged_lines=True,
            )
            .filter(~pl.all_horizontal(pl.all().is_null()))
            .fill_null("")
        )
    else:
        lf = rules.fill_null("")

    if common_rules_path:
        clf = (
            pl.scan_csv(
                common_rules_path,
                separator="\t",
                infer_schema_length=None,
                truncate_ragged_lines=True,
            )
            .filter(~pl.all_horizontal(pl.all().is_null()))
            .fill_null("")
        )
    elif common_rules is not None:
        clf = common_rules.fill_null("")
    else:
        clf = None

    # Normalize columns
    cols = set(lf.collect_schema().names())
    required = {label_col, rules_col}
    if not required.issubset(cols):
        raise ValueError(
            f"rules TSV missing required columns {required}. Found: {lf.columns}"
        )

    has_alias_col = alias_col and (alias_col in cols)

    with open(Path(__file__).parent.absolute() / "rules.lark") as f:
        parser = Lark(
            f,
            parser="lalr",
            transformer=ASTTransformer(
                allow_visualize_functions=allow_visualize_functions
            ),
        )

    def parse_rule_expr(expr_str: str) -> Expr:
        """We use a closure to capture the parser instance since
        polars map_elements doesn't let us pass extra args."""
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

    lf = lf.with_columns(
        [
            pl.col(rules_col)
            .str.strip_chars()
            .map_elements(parse_rule_expr, return_dtype=pl.Object)
        ]
    )

    lf = lf.with_columns(pl.col(pl.String).replace("", None))

    if not has_alias_col:
        if not alias_col:
            alias_col = "alias"
        lf = lf.with_columns(pl.lit(None).cast(pl.String).alias(alias_col))

    if clf is not None:
        clf = clf.with_columns(
            [
                (
                    pl.col("rule")
                    .str.strip_chars()
                    .map_elements(parse_rule_expr, return_dtype=pl.Object)
                    .alias(rules_col)
                ),
                pl.col("name").alias(label_col),
                pl.col("alias").alias(alias_col),
            ],
        )

        lf = pl.concat([lf, clf], how="diagonal_relaxed")
    # We strip white space around labels incase users left spaces at the beg or end of cell
    # which we don't want to match to
    lf = lf.with_columns(
        pl.col(label_col).str.strip_chars(), pl.col(alias_col).str.strip_chars()
    )

    rules = {
        a: b
        for a, b in lf.filter(~pl.col(label_col).is_null())
        .select([pl.col(label_col), pl.col(rules_col)])
        .collect()
        .iter_rows()
    }

    definitions = {
        a: b
        for a, b in lf.filter(~pl.col(alias_col).is_null())
        .select([pl.col(alias_col), pl.col(rules_col)])
        .collect()
        .iter_rows()
    }

    return definitions, rules, lf


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
    additional_cols = additional_cols or []
    index_cols = [count_col, group_col] if group_col else [count_col]
    besthit_cols = [col for col in besthit_cols if col in df.columns]
    besthit_cols.extend(
        [
            col
            for col in df.columns
            if (col.endswith("_id") and col not in besthit_cols and col != "query_id")
        ]
    )
    df = df.select(index_cols + besthit_cols + additional_cols)
    for col in besthit_cols:
        if col not in ID_EXPR_DICT:
            explode_col = ID_EXPR_DICT["DEFAULT"](col).alias(col)
        else:
            explode_col = ID_EXPR_DICT[col].alias(col)
        df = (
            df.with_columns(explode_col)
            .explode(col)
            .with_columns(pl.col(col).str.to_uppercase())
        )

    for col in COMPUTED_COLUMNS:
        try:
            df = df.with_columns(COMPUTED_COLUMNS[col].alias(col))
            besthit_cols.append(col)
        except pl.exceptions.ColumnNotFoundError:
            pass

    # unpivot to long (sample, hit)
    hit_col = "hit"
    df = (
        df.unpivot(
            index=index_cols + additional_cols,
            on=besthit_cols,
            variable_name="db",
            value_name=hit_col,
        )
        .drop("db")
        .drop_nulls()
        .filter(pl.col(hit_col).is_in(list(needed_features)))
        .select([*index_cols, hit_col, *additional_cols])
        .unique(subset=[*index_cols, hit_col, *additional_cols])
    )
    return df


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

    @functools.cache
    def eval_bool(self, expr: Expr, reduce_outer_and=True) -> np.ndarray:
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

    @functools.cache
    def eval_cycle(
        self, expr: Steps | list[Expr], simplify=True, **kwargs
    ) -> np.ndarray:
        if isinstance(expr, Iterable):
            parts = [self.eval_bool(p, **kwargs) for p in expr]
        elif isinstance(expr, Steps):
            parts = [self.eval_bool(p, **kwargs) for p in expr.parts]
        else:
            raise RuleError(
                f"Expected Steps or Iterable of Expr for cycle evaluation, got {expr}"
            )
        mat = (
            np.stack(
                parts, axis=1, dtype=bool, casting="unsafe"
            )  # unsafe here because we want float -> bool
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
                if "masks" in kwargs and kwargs["masks"]:
                    return self.not_(masks=kwargs["masks"])
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

    def _sort_df_to_ordered_df(self, df):
        return (
            df.join(self._order_df, on=self.sample_col, how="left")
            .sort("_order")
            .drop("_order")
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
    def not_(x: np.ndarray = None) -> np.ndarray | pl.Expr:
        if isinstance(x, list) and len(x) > 0 or isinstance(x, pl.Expr):
            if isinstance(x, pl.Expr):
                x = [x]
            masks = x
            if len(masks) == 1:
                mask = masks[0]
            else:
                mask = masks[0].and_(*[mask for mask in masks[1:]])
            return ~mask
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
        if mode.upper() not in AT_LEAST_MODES:
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
        return (
            self._sort_df_to_ordered_df(
                df.group_by(self.sample_col).agg(
                    pl.col("taxonomy").str.contains(tax_label).any()
                )
            )
            .select("taxonomy")
            .to_series()
            .to_numpy()
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

        df = self._sort_df_to_ordered_df(df)
        return df.select(pl.col(col)).to_series().to_numpy()

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
        df = self._sort_df_to_ordered_df(df)
        return df.select(pl.col(col)).to_series().to_numpy()

    def filter_contains(
        self, col: str, val: str, df: pl.DataFrame = None, **kwargs
    ) -> pl.DataFrame:
        df = self.annotations if df is None else df
        if col not in df:
            return np.zeros(len(self.samples), dtype=bool)
        return pl.col(col).str.contains(val)

    def filter_compare(
        self, col: str, op: str, thr: float, df: pl.DataFrame = None, **kwargs
    ) -> pl.DataFrame:
        df = self.annotations if df is None else df
        if col not in df:
            raise RuleError(f"Missing column '{col}' for filter_compare()")
        try:
            cmp_fn = OP_TO_EXPR[op]
        except KeyError:
            raise ValueError(f"Unsupported op={op!r}. Use one of {sorted(OP_TO_EXPR)}")
        return cmp_fn(pl.col(col).cast(float), thr)

    def path_steps(self, expr: Steps | list[Expr], **kwargs) -> pl.DataFrame:
        cycle = self.eval_cycle(expr, simplify=False, **kwargs)

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


def get_val_col(df):
    if "coverage_percentage" in df.columns:
        return "coverage_percentage"
    elif "present" in df.columns:
        return "present"
    elif "value" in df.columns:
        return "value"
    raise AssertionError(
        f"Rules parsing problem, can't determine check_col for dataframe rules: {df}. File github issue"
    )


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

    rule_names = sorted(compiled.rules.keys())
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
    # We drop nones here to make sure any rules without groups just aren't included
    for group, frame in compiled.df.drop_nulls(subset=[group_col]).group_by(
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
            if anno_df is not None:
                check_col = get_val_col(out)

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
    compiled = CompiledRules.from_rules(*args, **kwargs)
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
        besthit_cols=list(ID_EXPR_DICT.keys()),
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
