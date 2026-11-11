from __future__ import annotations

import contextlib
import functools
import json
import logging
import operator
import os
from collections.abc import Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, ClassVar

import numpy as np
import polars as pl
from lark import Lark, LarkError, Transformer

from dram_rules.definitions import (
    BACKUP_GROUPBY_COLUMNS,
    DEFAULT_GROUPBY_COLUMN,
)
from dram_rules.rules_api import (
    ArgKind,
    FunctionRegistry,
    Group,
    Op,
    Param,
    Rows,
    RuleError,
    RuleFunction,
    ValueKind,
    default_functions,
)

logger = logging.getLogger("drams.rules")

COVERAGE_MODES = ("steps", "subunits")


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
COMPUTED_COLUMNS: dict[str, tuple[str, ...]] = {
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


# ---------------------------------------------------------------------------
# Rule syntax tree
#
# Every node prints as rule text (`str(node)`), for messages, and lists its
# child nodes (`children()`), so tree walks are written once (`walk`).
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class Expr:
    """A node of a parsed rule."""

    def children(self) -> tuple[Expr, ...]:
        return ()

    def with_children(self, children: tuple[Expr, ...]) -> Expr:
        """A copy with new child nodes (same order as `children()`)."""
        return self

    def pprint(self):
        """Pretty print the expression for debugging."""
        from pprint import pprint

        pprint(self)


def walk(expr: Expr):
    """`expr` and every node below it."""
    yield expr
    for child in expr.children():
        yield from walk(child)


def _bracketed(part: Expr) -> str:
    return f"[{part}]" if isinstance(part, (Combination, PipeChain)) else str(part)


@dataclass(frozen=True)
class Name(Expr):
    """A gene ID. As a function argument it can also be a bare word, such as a
    column name, a number or a choice like COUNTS, depending on the parameter."""

    value: str
    qualifier: str | None = None

    def __str__(self):
        return f"{self.qualifier}@{self.value}" if self.qualifier else self.value


@dataclass(frozen=True)
class Alias(Expr):
    """`@name`, before expansion."""

    value: str

    def __str__(self):
        return f"@{self.value}"


@dataclass(frozen=True)
class AliasRef(Expr):
    """An expanded `@name`: its definition, plus where it came from, so
    messages can name the alias the user wrote. Removed when compiling."""

    name: str
    body: Expr

    def children(self):
        return (self.body,)

    def with_children(self, children):
        return AliasRef(self.name, children[0])

    def __str__(self):
        return f"@{self.name}"


@dataclass(frozen=True)
class Number(Expr):
    value: float

    def __str__(self):
        return f"{self.value:g}"


@dataclass(frozen=True)
class String(Expr):
    value: str

    def __str__(self):
        return json.dumps(self.value)


@dataclass(frozen=True)
class OpLit(Expr):
    """A comparison operator written as a function argument, e.g. `>=`."""

    symbol: str

    def __str__(self):
        return self.symbol


@dataclass(frozen=True)
class Not(Expr):
    """`not(x)`: negates a per-sample value, or a row filter. Part of the
    language (like `&` and `|`), so it takes the kind of its argument."""

    part: Expr

    def children(self):
        return (self.part,)

    def with_children(self, children):
        return Not(children[0])

    def __str__(self):
        return f"not({self.part})"


@dataclass(frozen=True)
class Combination(Expr):
    """Parts joined by one operator: `a & b`, `a | b`, or steps `a, b`."""

    parts: tuple[Expr, ...]

    separator: ClassVar[str]
    needs_all: ClassVar[bool]  # True: all parts must be present; False: any part

    def children(self):
        return self.parts

    def with_children(self, children):
        return type(self)(tuple(children))

    def __str__(self):
        return self.separator.join(_bracketed(p) for p in self.parts)


class And(Combination):
    separator = " & "
    needs_all = True


class Or(Combination):
    separator = " | "
    needs_all = False


class Steps(Combination):
    """`a, b, c`: present when every step is present."""

    separator = ", "
    needs_all = True


@dataclass(frozen=True)
class Call(Expr):
    """A function call as written, before it is matched to its definition."""

    value: str
    args: tuple[Expr, ...] = ()
    kwargs: tuple[tuple[str, Expr], ...] = ()

    def children(self):
        return (*self.args, *(v for _, v in self.kwargs))

    def with_children(self, children):
        n = len(self.args)
        kwargs = tuple((k, v) for (k, _), v in zip(self.kwargs, children[n:]))
        return Call(self.value, tuple(children[:n]), kwargs)

    def __str__(self):
        items = [str(a) for a in self.args] + [f"{k}={v}" for k, v in self.kwargs]
        return f"{self.value}({', '.join(items)})"


@dataclass(frozen=True)
class PipeChain(Expr):
    """`source -> call -> ... -> call`: row filters, then a per-sample summary."""

    parts: tuple[Expr, ...]

    def children(self):
        return self.parts

    def with_children(self, children):
        return PipeChain(tuple(children))

    def __str__(self):
        return " -> ".join(_bracketed(p) for p in self.parts)


@dataclass(frozen=True)
class BoundCall(Expr):
    """A call matched to its function, with its arguments converted: one value
    per parameter of `function` (a tuple for `*args` parameters). Literal
    arguments are Python values; sub-rule arguments are compiled `Expr`s."""

    function: RuleFunction
    args: tuple[Any, ...]
    kind: ValueKind  # what the call returns: ValueKind.PER_SAMPLE or ValueKind.ROW_FILTER

    def children(self):
        found = []
        for param, value in zip(self.function.params, self.args):
            if param.kind.is_expression:
                found += value if param.variadic else (value,)
        return tuple(found)

    def __str__(self):
        def text(value):
            if isinstance(value, bool):
                return "true" if value else "false"
            return str(value)

        items = []
        for param, value in zip(self.function.params, self.args):
            if param.variadic:
                items += [text(v) for v in value]
            elif param.required and not param.keyword_only:
                items.append(text(value))
            else:
                items.append(f"{param.name}={text(value)}")
        return f"{self.function.name}({', '.join(items)})"


@dataclass(frozen=True)
class _KeywordArg:
    name: str
    value: Expr


class ASTTransformer(Transformer):
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
        return String(json.loads(str(items[0])))

    def op(self, items):
        return OpLit(str(items[0]))

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

    def kwarg(self, items):
        name, value = items
        return _KeywordArg(str(name), value)

    def call(self, items):
        name = str(items[0])
        args, kwargs = [], []
        for item in items[1:]:
            if item is None:  # empty argument list
                continue
            if isinstance(item, _KeywordArg):
                kwargs.append((item.name, item.value))
            elif kwargs:
                raise RuleError(
                    f"{name}(): positional argument {item} after keyword arguments"
                )
            else:
                args.append(item)
        if name == "not":
            if len(args) != 1 or kwargs:
                raise RuleError("not() takes exactly one argument")
            return Not(args[0])
        return Call(name, tuple(args), tuple(kwargs))

    def pipe_(self, items):
        return PipeChain(tuple(items))


@functools.cache
def _get_parser() -> Lark:
    """Build the LALR parser once; building the tables isn't free."""
    grammar = (Path(__file__).parent / "rules.lark").read_text()
    return Lark(grammar, parser="lalr", transformer=ASTTransformer())


def _parse_rule_expr(parser: Lark, expr_str: str) -> Expr:
    """Parse one rule string, turning Lark errors into friendlier RuleErrors."""
    try:
        return parser.parse(expr_str)
    except RuleError as e:  # raised by ASTTransformer while parsing
        raise RuleError(f"{e} in rule: {expr_str}") from None
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


# ---------------------------------------------------------------------------
# Compiling: match calls to functions, convert arguments, check kinds
# ---------------------------------------------------------------------------


def builtin_functions() -> FunctionRegistry:
    """The built-in rule functions, registered by importing their module."""
    import dram_rules.rule_parser.functions  # noqa: F401

    return default_functions


_TRUE_WORDS = {"true", "yes", "t", "1"}
_FALSE_WORDS = {"false", "no", "f", "0"}
_OLD_OPERATORS = {"ge": ">=", "gt": ">", "le": "<=", "lt": "<", "eq": "==", "ne": "!="}


def _word(expr: Expr) -> str | None:
    """A bare word or quoted string argument, as text."""
    if isinstance(expr, Name) and expr.qualifier is None:
        return expr.value
    if isinstance(expr, String):
        return expr.value
    return None


def _number(expr: Expr) -> float | None:
    if isinstance(expr, Number):
        return expr.value
    if isinstance(expr, Name) and expr.qualifier is None:
        try:
            return float(expr.value)
        except ValueError:
            return None
    return None


class _ArgumentError(RuleError):
    """A problem with a call's own arguments (not with a sub-rule inside them)."""


class RuleCompiler:
    """Turns parsed rules (aliases already expanded) into compiled rules.

    Each call is matched to its function and its arguments converted by
    parameter type. Every piece gets a kind, per-sample value or row filter,
    and the kinds are checked: row filters only before `->`, per-sample values
    everywhere else.
    """

    def __init__(self, functions: FunctionRegistry):
        self.functions = functions

    def compile_rule(self, expr: Expr) -> Expr:
        """Compile a whole rule, which must give a per-sample value."""
        return self.expect(expr, ValueKind.PER_SAMPLE)

    def compile(self, expr: Expr) -> tuple[Expr, ValueKind]:
        match expr:
            case Name():
                return expr, ValueKind.PER_SAMPLE
            case Alias(value=name):
                raise RuleError(f"Undefined alias: @{name}")
            case Number() | String() | OpLit():
                raise RuleError(
                    f"{expr} is a literal value; it can only be a function argument"
                )
            case AliasRef(name=name, body=body):
                # Say which alias an error is in: the text in the message is
                # from that alias's definition, not from the rule itself.
                try:
                    return self.compile(body)
                except RuleError as e:
                    raise RuleError(f"in @{name}: {e}") from None
            case Not(part=part):  # same kind as its argument
                compiled, kind = self.compile(part)
                return Not(compiled), kind
            case Steps(parts=parts):  # steps are always per-sample
                parts = tuple(self.expect(p, ValueKind.PER_SAMPLE) for p in parts)
                return Steps(parts), ValueKind.PER_SAMPLE
            case Combination(parts=parts):  # & and |: all parts of one kind
                compiled = [self.compile(p) for p in parts]
                kinds = {kind for _, kind in compiled}
                if len(kinds) > 1:
                    raise RuleError(
                        f"'{expr}' mixes row filters and per-sample values"
                        f" with '{expr.separator.strip()}'"
                    )
                return expr.with_children(tuple(e for e, _ in compiled)), kinds.pop()
            case Call():
                return self.bind(expr)
            case PipeChain():
                return self.compile_pipe(expr), ValueKind.PER_SAMPLE
        raise TypeError(f"Unexpected rule node: {expr!r}")

    def expect(self, expr: Expr, kind: ValueKind) -> Expr:
        compiled, got = self.compile(expr)
        if got is not kind:
            if got is ValueKind.ROW_FILTER:
                raise RuleError(
                    f"{expr} is a row filter; row filters can only be used before '->'"
                )
            raise RuleError(
                f"{expr} is per-sample and can't be used before '->';"
                " only row filters (filter_...) can"
            )
        return compiled

    def compile_pipe(self, pipe: PipeChain) -> PipeChain:
        *filters, last = pipe.parts
        compiled = [self.expect(f, ValueKind.ROW_FILTER) for f in filters]
        summary, kind = self.bind(last)
        if kind is not ValueKind.PER_SAMPLE:
            raise RuleError(
                f"{last} is a row filter, so it can't end a pipe; the last call"
                " must summarize rows per sample, e.g. count_rows(>=, 1)"
            )
        if summary.function.rows_param is None:
            raise RuleError(
                f"{summary.function.name}() doesn't read annotation rows, so it can't be used after '->'"
            )
        return PipeChain((*compiled, summary))

    def bind(self, call: Call) -> tuple[BoundCall, ValueKind]:
        function = self.functions.get(call.value)
        try:
            values = self._bind_arguments(function, call)
        except _ArgumentError as e:
            raise RuleError(f"{e} (expected {function.signature()})") from None
        return BoundCall(function, tuple(values), function.returns), function.returns

    def _bind_arguments(self, function: RuleFunction, call: Call):
        positional = list(call.args)
        keywords: dict[str, Expr] = {}
        for name, value in call.kwargs:
            if name in keywords:
                raise _ArgumentError(f"{function.name}() got '{name}' twice")
            keywords[name] = value

        values = []
        for param in function.params:
            if param.variadic:
                raws, positional = positional, []
                values.append(tuple(self._convert(function, param, raw) for raw in raws))
                continue
            if positional and not param.keyword_only:
                raw = positional.pop(0)
                if param.name in keywords:
                    raise _ArgumentError(f"{function.name}() got two values for '{param.name}'")
            elif param.name in keywords:
                raw = keywords.pop(param.name)
            elif param.required:
                raise _ArgumentError(f"{function.name}() is missing argument '{param.name}'")
            else:
                values.append(param.default)
                continue
            values.append(self._convert(function, param, raw))

        if positional:
            raise _ArgumentError(
                f"{function.name}() got too many arguments: {', '.join(str(a) for a in positional)}"
            )
        if keywords:
            valid = ", ".join(p.name for p in function.params) or "none"
            raise _ArgumentError(
                f"{function.name}() got unknown keyword(s) {', '.join(keywords)}; valid: {valid}"
            )
        return values

    def _convert(self, function: RuleFunction, param: Param, raw: Expr) -> Any:
        """Convert one argument for `param`: a Python value, or a compiled sub-rule."""

        def wrong(expected: str, suggestion: str = ""):
            return _ArgumentError(
                f"{function.name}() argument '{param.name}' expects {expected}, got {raw}"
                + (f"; {suggestion}" if suggestion else "")
            )

        match param.kind:
            case ArgKind.INT:
                number = _number(raw)
                if number is None or not float(number).is_integer():
                    raise wrong("an integer")
                return int(number)
            case ArgKind.FLOAT:
                number = _number(raw)
                if number is None:
                    raise wrong("a number")
                return float(number)
            case ArgKind.BOOL:
                word = (_word(raw) or "").lower()
                if word in _TRUE_WORDS:
                    return True
                if word in _FALSE_WORDS:
                    return False
                raise wrong("true or false")
            case ArgKind.STR | ArgKind.COLUMN:
                word = _word(raw)
                if word is None:
                    raise wrong(param.kind.value)
                return word
            case ArgKind.CHOICE:
                word = (_word(raw) or "").upper()
                for choice in param.choices:
                    if choice.upper() == word:
                        return choice
                raise wrong("one of " + ", ".join(param.choices))
            case ArgKind.OP:
                if isinstance(raw, OpLit):
                    return Op(raw.symbol)
                old = _OLD_OPERATORS.get((_word(raw) or "").lower())
                raise wrong("a comparison operator (>=, >, <=, <, ==, !=)", f"use {old}" if old else "")
            case ArgKind.PER_SAMPLE | ArgKind.GROUP:
                return self.expect(raw, ValueKind.PER_SAMPLE)
            case ArgKind.ROW_FILTER:
                return self.expect(raw, ValueKind.ROW_FILTER)
        raise TypeError(param.kind)


def gene_ids(expr: Expr) -> set[str]:
    """The gene IDs a compiled rule looks up. Compiling turns literal
    arguments (column names, numbers, ...) into values, so every `Name` left
    is a gene ID."""
    return {node.value for node in walk(expr) if isinstance(node, Name)}


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
    rules: dict[str, Expr]  # compiled top-level rules, aliases expanded
    needed_features: set[str]  # gene IDs any rule looks up, uppercase
    needed_features_raw: set[str]  # the same, as written in the rules
    features_by_rules: dict[str, set[str]]
    df: pl.DataFrame
    # Parsed `rule_options` JSON per rule name; rules without options are absent.
    rule_options: dict[str, dict] = field(default_factory=dict)
    functions: FunctionRegistry = field(default_factory=builtin_functions)

    @classmethod
    def from_rules(
        cls,
        *args,
        label_col: str = "name",
        group_col=None,
        group_col_backups=BACKUP_GROUPBY_COLUMNS,
        functions: FunctionRegistry | None = None,
        function_files: Sequence[str | Path] = (),
        **kwargs,
    ) -> CompiledRules:
        """Load, expand and compile rule sheets (arguments as for `load_rules`).

        `functions` are the functions rules may call (default: the built-ins);
        `function_files` are Python files of extra `@rule_function`s to add.
        """
        if functions is None:
            functions = builtin_functions()
        functions = functions.with_files(function_files)
        definitions, rules, df = load_rules(*args, label_col=label_col, **kwargs)
        compiler = RuleCompiler(functions)

        # Check each alias on its own, so errors in unused aliases are found too.
        for name, expr in definitions.items():
            try:
                compiler.compile(expand_macros(expr, definitions))
            except RuleError as e:
                raise RuleError(f"Alias '@{name}': {e}") from None

        rules_compiled, features_by_rules = {}, {}
        for name, expr in rules.items():
            try:
                rules_compiled[name] = compiler.compile_rule(expand_macros(expr, definitions))
                features_by_rules[name] = gene_ids(rules_compiled[name])
                _check_not_alias_names(features_by_rules[name], definitions)
            except RuleError as e:
                raise RuleError(f"Rule '{name}': {e}") from None
        needed_features = set().union(*features_by_rules.values())

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

        rule_options = _parse_rule_options(df, label_col)
        _check_coverage_options(rule_options, df, label_col, group_col)

        return cls(
            rules=rules_compiled,
            needed_features={item.upper() for item in needed_features},
            needed_features_raw=needed_features,
            features_by_rules=features_by_rules,
            df=df,
            rule_options=rule_options,
            functions=functions,
        )


def _check_not_alias_names(genes: set[str], definitions: dict[str, Expr]) -> None:
    clashes = sorted(genes & definitions.keys())
    if clashes:
        raise RuleError(
            f"Gene ID(s) {', '.join(clashes)} match alias names but are written without '@'."
            " Gene IDs and alias names must differ; did you mean"
            f" {', '.join('@' + c for c in clashes)}?"
        )


def coverage_mode(rule_options: dict) -> str | None:
    """The rule's `coverage` option: "steps", "subunits", or None (present/absent)."""
    return rule_options.get("coverage")


def _check_coverage_options(
    rule_options: dict[str, dict], df: pl.DataFrame, label_col: str, group_col: str | None
) -> None:
    for name, options in rule_options.items():
        mode = coverage_mode(options)
        if mode is not None and mode not in COVERAGE_MODES:
            raise RuleError(
                f"Rule '{name}': rule_options coverage must be one of"
                f" {', '.join(COVERAGE_MODES)}, got {mode!r}"
            )
    if not group_col or group_col not in df.columns:
        return
    # A heatmap group shows one kind of value, so its rules must agree.
    named = df.select(label_col, group_col).drop_nulls(label_col)
    for (group,), frame in named.group_by(group_col, maintain_order=True):
        by_mode: dict[str, list[str]] = {}
        for label in frame[label_col]:
            mode = coverage_mode(rule_options.get(label, {})) or "present/absent"
            by_mode.setdefault(mode, []).append(label)
        if len(by_mode) > 1:
            details = "; ".join(f"{m}: {', '.join(ls[:5])}" for m, ls in by_mode.items())
            raise RuleError(
                f"Rules in group '{group}' must all use the same rule_options coverage."
                f" Found {details}"
            )


def _parse_rule_options(df: pl.DataFrame, label_col: str) -> dict[str, dict]:
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
    common_rules: pl.LazyFrame = None,
    common_rules_path: str | Path = None,
    check_for_dupes: bool = True,
) -> tuple[dict[str, Expr], dict[str, Expr], pl.DataFrame]:
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
    parser = _get_parser()
    parsed = {
        s: _parse_rule_expr(parser, s)
        for s in df.get_column(rules_col).drop_nulls().unique()
    }

    def rules_keyed_by(key_col: str) -> dict[str, Expr]:
        keyed = df.select(key_col, rules_col).drop_nulls(key_col)
        return {
            key: parsed.get(rule) for key, rule in zip(keyed[key_col], keyed[rules_col])
        }

    return rules_keyed_by(alias_col), rules_keyed_by(label_col), df


def expand_macros(expr: Expr, definitions: dict[str, Expr]) -> Expr:
    """Replace every `@alias` in `expr` by its definition, recursively. Each
    expansion is kept as an `AliasRef`, so messages can name the alias."""
    stack: list[str] = []  # aliases being expanded, to detect cycles

    def expand(expr: Expr) -> Expr:
        if not isinstance(expr, Alias):
            return expr.with_children(tuple(expand(c) for c in expr.children()))
        name = expr.value
        if name not in definitions:
            raise RuleError(f"Undefined alias: @{name}")
        if name in stack:
            raise RuleError(f"Cycle detected: {' -> '.join(stack + [name])}")
        stack.append(name)
        body = expand(definitions[name])
        stack.pop()
        return AliasRef(name, body)

    return expand(expr)


def prepare_present_map_df(
    df: pl.DataFrame,
    count_col: str,
    besthit_cols: list[str],
    needed_features: set[str],
    additional_cols: list[str] = None,
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
    needed_features: set[str],
    sample_names: list[str],
    value_col: str = None,
) -> tuple[list[str], dict[str, np.ndarray]]:
    """Build present_map of needed gene_ids from annotations DataFrame"""
    hit_col = "hit"

    df: pl.DataFrame = df.drop_nulls(subset=sample_col)

    sample_index = {s: i for i, s in enumerate(sample_names)}
    n = len(sample_names)
    present_map: dict[str, np.ndarray] = {f: np.zeros(n) for f in needed_features}

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
    """Evaluates compiled rules for a set of samples.

    Args:
        samples: sample names, in output order.
        present_map: hit counts per (uppercase) gene ID, one value per sample.
        sample_col: name of the sample column in `annotations`.
        annotations: annotation rows, for functions that read rows (`Rows`).
    """

    def __init__(
        self,
        samples: list[str],
        present_map: dict[str, np.ndarray],
        sample_col: str,
        annotations: pl.DataFrame | None,
    ):
        self.samples = samples
        self.present_map = present_map
        self.annotations = annotations
        self.sample_col = sample_col
        self._all_false = np.zeros(len(samples), dtype=bool)
        self._all_rows = (
            Rows(annotations, sample_col, samples) if annotations is not None else None
        )
        # Per-instance cache. `functools.cache` on the method itself would be one
        # cache shared by the whole class, keeping every Evaluator (and its
        # annotations DataFrame) alive for the life of the process.
        self.evaluate = functools.cache(self._evaluate)

    def _evaluate(self, expr: Expr) -> np.ndarray:
        """Per-sample value of a compiled rule: hit counts or True/False."""
        match expr:
            case Name(value=value):
                return self.present_map.get(value.upper(), self._all_false)
            case Not(part=part):
                return np.logical_not(self.evaluate(part))
            case Combination(parts=parts):
                stacked = np.stack([self.evaluate(p) for p in parts], axis=1)
                return np.all(stacked, axis=1) if expr.needs_all else np.any(stacked, axis=1)
            case BoundCall():
                rows = self._rows(expr) if expr.function.rows_param else None
                return self._call(expr, rows)
            case PipeChain(parts=parts):
                rows = self._rows(expr)
                for row_filter in parts[:-1]:
                    rows = rows.filter(self.row_filter(row_filter))
                return self._call(parts[-1], rows)
        raise RuleError(f"Can't evaluate {expr} per sample")

    def _rows(self, expr: Expr) -> Rows:
        """All annotation rows, for `expr` that reads them."""
        if self._all_rows is None:
            raise RuleError(f"{expr} reads annotation rows, but no annotations were given")
        return self._all_rows

    def row_filter(self, expr: Expr) -> pl.Expr:
        """The Polars expression of a compiled row filter."""
        match expr:
            case BoundCall():
                return self._call(expr, None)
            case Not(part=part):
                return ~self.row_filter(part)
            case And(parts=parts) | Or(parts=parts):
                combine = operator.and_ if expr.needs_all else operator.or_
                return functools.reduce(combine, (self.row_filter(p) for p in parts))
        raise RuleError(f"{expr} is not a row filter")

    def group(self, expr: Expr) -> Group:
        """The parts of a gene group: `[a, b]`, `[a | b]`, or a single part."""
        parts = expr.parts if isinstance(expr, (Steps, Or)) else (expr,)
        counts = np.stack([np.asarray(self.evaluate(p), dtype=float) for p in parts], axis=1)
        return Group(_read_only(counts), tuple(str(p) for p in parts))

    def coverage(self, expr: Expr, mode: str) -> pl.DataFrame:
        """How much of a rule is present per sample (rule_options `coverage`).

        "steps" counts the rule's top-level steps (`a, b, c`); "subunits" its
        top-level `&` parts. A rule of another shape counts as one part. (The
        mode was checked when the rules were loaded.)
        """
        top = Steps if mode == "steps" else And
        parts = expr.parts if isinstance(expr, top) else (expr,)
        present = np.stack(
            [np.asarray(self.evaluate(p)).astype(bool) for p in parts], axis=1
        )
        df = pl.DataFrame(
            {
                mode: present.shape[1],
                f"{mode}_present": present.sum(axis=1),
                self.sample_col: self.samples,
            }
        )
        return df.with_columns(
            coverage_percentage=pl.col(f"{mode}_present") / pl.col(mode)
        )

    def _call(self, call: BoundCall, rows: Rows | None):
        function = call.function
        values = {
            param.name: (
                tuple(self._argument(call, param, v) for v in value)
                if param.variadic
                else self._argument(call, param, value)
            )
            for param, value in zip(function.params, call.args)
        }
        if function.rows_param:
            values[function.rows_param] = rows
        try:
            result = function(values)
        except RuleError:
            raise
        except Exception as e:
            raise RuleError(
                f"{call}: {function.name}() failed ({function.source}):"
                f" {type(e).__name__}: {e}"
            ) from e
        return self._checked(call, result)

    def _argument(self, call: BoundCall, param: Param, value: Any) -> Any:
        # Per-sample values are shared (cached, and the caller's present_map),
        # so functions get read-only views: an in-place change can't leak.
        match param.kind:
            case ArgKind.PER_SAMPLE:
                return _read_only(self.evaluate(value))
            case ArgKind.GROUP:
                return self.group(value)
            case ArgKind.ROW_FILTER:
                return self.row_filter(value)
            case ArgKind.COLUMN:
                if self.annotations is None or value not in self.annotations.columns:
                    raise RuleError(
                        f"Missing column '{value}' for {call.function.name}():"
                        " the annotations don't have it. If the rule should be skipped"
                        " when the column is missing, give it rule_options"
                        f' {{"include_when":{{"all_columns_present":["{value}"]}}}}'
                    )
                return value
        return value

    def _checked(self, call: BoundCall, result: Any) -> Any:
        """Check a function's result against its return hint, and say what's
        wrong (and how to fix the function) if they don't match."""
        function = call.function
        n = len(self.samples)
        hint = "RowFilter" if function.returns is ValueKind.ROW_FILTER else "PerSample"
        declared = f"{function.name}() has the return hint {hint} ({function.source}) but returned"
        if function.returns is ValueKind.ROW_FILTER:
            if isinstance(result, pl.Expr):
                return result
            raise RuleError(
                f"{declared} {_describe(result)}. A row filter is a Polars expression"
                " giving True/False per annotation row, e.g. pl.col(col) > 1."
            )
        if isinstance(result, pl.Expr):
            raise RuleError(
                f"{declared} a Polars expression, which is a row filter. If it is a filter,"
                f" give it the return hint RowFilter and name it 'filter_{function.name}';"
                f" otherwise return one value per sample (an array of length {n})."
            )
        values = np.asarray(result)
        if values.shape != (n,):
            raise RuleError(
                f"{declared} {_describe(result)}; expected one value per sample"
                f" (an array of length {n})."
            )
        if values.dtype.kind not in "biuf":
            raise RuleError(
                f"{declared} values of type {values.dtype}; expected numbers or True/False."
            )
        return values


@contextlib.contextmanager
def _in_rule(name: str):
    """Prefix errors raised while evaluating a rule with the rule's name. The
    original cause (e.g. an error inside a user function) stays attached."""
    try:
        yield
    except RuleError as e:
        raise RuleError(f"Rule '{name}': {e}") from (e.__cause__ or e)


def _describe(value: Any) -> str:
    """A short description of a function's result, for messages."""
    if isinstance(value, pl.Expr):
        return "a Polars expression"
    if np.isscalar(value) or value is None:
        return f"a single value ({value!r})"
    shape = getattr(value, "shape", None)
    if shape is not None:
        return f"{type(value).__name__} of shape {shape}"
    return type(value).__name__


def _read_only(values: np.ndarray) -> np.ndarray:
    view = np.asarray(values).view()
    view.flags.writeable = False
    return view


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
            "In rule sheet, all rules inside a single group must use the same"
            " rule_options coverage (steps, subunits, or none)."
        )
    # if we expect a single return, we extract it and return
    return cols[0]


def evaluate_rules(
    compiled: CompiledRules,
    samples: list[str],
    present_map: dict[str, np.ndarray],
    annotations: pl.DataFrame | None = None,
    sample_col: str | None = None,
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
        with _in_rule(rn):
            out = ev.evaluate(compiled.rules[rn])
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
    samples: list[str],
    present_map: dict[str, np.ndarray],
    annotations: pl.DataFrame | None,
    group_col: str | None = "group",
    label_col: str | None = "name",
    sample_col: str | None = None,
    additional_cols: list[str] | None = None,
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
    available_columns = set(annotations.columns) if annotations is not None else set()
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
        for rn in frame[label_col]:
            options = compiled.rule_options.get(rn, {})
            if not should_include(options, available_columns):
                continue

            expr = compiled.rules[rn]
            mode = coverage_mode(options)
            with _in_rule(rn):
                if mode:
                    out = ev.coverage(expr, mode)
                else:
                    out = pl.DataFrame(
                        {"present": ev.evaluate(expr).astype(bool), sample_col: ev.samples}
                    )
            out = out.with_columns(pl.lit(rn).alias(label_col))
            features = sorted({f.upper() for f in compiled.features_by_rules[rn]})
            hits = [[] for _ in samples]

            for feature in features:
                counts = present_map.get(feature)
                if counts is not None:
                    for i in np.flatnonzero(counts):
                        hits[i].append(feature)

            # join on the sample column, so the order of `out` doesn't matter
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
