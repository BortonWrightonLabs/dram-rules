from __future__ import annotations

import contextlib
import functools
import json
import logging
import operator
from collections.abc import Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, ClassVar, Literal

import numpy as np
import polars as pl
from lark import Lark, LarkError, Transformer

from dram_rules.definitions import (
    BACKUP_ALIAS_COLUMNS,
    BACKUP_GROUPBY_COLUMNS,
    BACKUP_LABEL_COLUMNS,
    DEFAULT_ALIAS_COLUMN,
    DEFAULT_GROUPBY_COLUMN,
    DEFAULT_LABEL_COLUMN,
)
from dram_rules.processing.process_annotations import scan_tsv
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

EC_NUMBER = r"EC:\d+\.\d+\.\d+\.\d+"
PFAM = r"PF\d{5}"

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

_TRUE_WORDS = {"true", "yes", "t", "1"}
_FALSE_WORDS = {"false", "no", "f", "0"}
_OLD_OPERATORS = {"ge": ">=", "gt": ">", "le": "<=", "lt": "<", "eq": "==", "ne": "!="}

# Functions removed from the rule language, and what to write instead.
_REMOVED_FUNCTIONS = {
    "path_steps": "write the steps as the rule itself (A, B, C) and give it"
    ' rule_options {"coverage":"steps"}',
    "path_subunits": "write the subunits as the rule itself (A & B & C) and give it"
    ' rule_options {"coverage":"subunits"}',
    "column_count_values": "use filter_compare(col, >=, value) -> count_rows(>=, count),"
    " e.g. column_count_values(col,ge,4,ge,3) becomes"
    " filter_compare(col, >=, 4) -> count_rows(>=, 3)",
    "column_sum_values": "use sum_column(col, op, value), e.g. sum_column(col, >=, 10)",
}

EvaluationLevel = Literal["gene", "genome", "sample"]
EVALUATION_LEVELS = ("gene", "genome", "sample")


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
    kind: (
        ValueKind  # what the call returns: ValueKind.PER_SAMPLE or ValueKind.ROW_FILTER
    )

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
        if call.value in _REMOVED_FUNCTIONS and call.value not in self.functions:
            raise RuleError(
                f"{call.value}() is no longer part of the rule language;"
                f" {_REMOVED_FUNCTIONS[call.value]}"
            )
        function = self.functions.get(call.value)
        try:
            values = self._bind_arguments(function, call)
        except _ArgumentError as e:
            raise RuleError(f"{e} (expected {function.signature()})") from None
        return BoundCall(function, tuple(values), function.returns), function.returns

    def _bind_arguments(self, function: RuleFunction, call: Call):
        self._check_choices_given_by_keyword(function, call)
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
                values.append(
                    tuple(self._convert(function, param, raw) for raw in raws)
                )
                continue
            if positional and not param.keyword_only:
                raw = positional.pop(0)
                if param.name in keywords:
                    raise _ArgumentError(
                        f"{function.name}() got two values for '{param.name}'"
                    )
            elif param.name in keywords:
                raw = keywords.pop(param.name)
            elif param.required:
                raise _ArgumentError(
                    f"{function.name}() is missing argument '{param.name}'"
                )
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

    def _check_choices_given_by_keyword(
        self, function: RuleFunction, call: Call
    ) -> None:
        """Catch a choice word passed positionally to an optional choice
        parameter that comes later, e.g. the earlier form
        at_least(2, PRESENCE, [...]). Otherwise the word would be read as the
        argument in that position (here a gene named PRESENCE)."""
        for position, raw in enumerate(call.args):
            word = _word(raw)
            if word is None:
                continue
            for param in function.params[position + 1 :]:
                if (
                    param.kind is ArgKind.CHOICE
                    and not param.required
                    and any(word.upper() == choice.upper() for choice in param.choices)
                ):
                    raise _ArgumentError(
                        f"{function.name}(): {word} is a value for '{param.name}', which is"
                        f" optional and comes later; write it as {param.name}={word}"
                        " after the other arguments"
                    )

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
                raise wrong(
                    "a comparison operator (>=, >, <=, <, ==, !=)",
                    f"use {old}" if old else "",
                )
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
    return scan_tsv(
        path, infer_schema_length=None, truncate_ragged_lines=True, drop_empty_rows=True
    )


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
    label_col: str = DEFAULT_LABEL_COLUMN  # the sheet's label column, as found
    group_col: str | None = None  # the sheet's Group column, as found

    @property
    def reads_rows(self) -> bool:
        """Whether any rule calls a function that reads annotation rows."""
        return any(
            isinstance(node, BoundCall) and node.function.rows_param is not None
            for expr in self.rules.values()
            for node in walk(expr)
        )

    @classmethod
    def from_rules(
        cls,
        rules_path: str | Path | None = None,
        rules: pl.LazyFrame | None = None,
        label_col: str | None = None,
        alias_col: str | None = None,
        group_col=None,
        group_col_backups=BACKUP_GROUPBY_COLUMNS,
        functions: FunctionRegistry | None = None,
        function_files: Sequence[str | Path] = (),
        **kwargs,
    ) -> CompiledRules:
        """Load, expand and compile a rule sheet (other arguments as for `load_rules`).

        The label, alias and Group columns are the first of their default and
        backups the sheet has (`definitions.py`), unless given; a given column
        the sheet lacks, or no Group column at all, is a RuleError.
        `functions` are the functions rules may call (default: the built-ins);
        `function_files` are Python files of extra `@rule_function`s to add.
        """
        if functions is None:
            functions = builtin_functions()
        functions = functions.with_files(function_files)
        sheet = _scan_rules_tsv(rules_path) if rules_path else rules
        columns = sheet.collect_schema().names()
        given = [label_col, alias_col, group_col]
        if missing := [c for c in given if c and c not in columns]:
            raise RuleError(f"No column {missing}, given as a custom column")
        label_col = label_col or find_label_column(columns)
        alias_col = alias_col or next(
            (c for c in [DEFAULT_ALIAS_COLUMN, *BACKUP_ALIAS_COLUMNS] if c in columns),
            DEFAULT_ALIAS_COLUMN,
        )
        definitions, rules, df = load_rules(
            rules=sheet, label_col=label_col, alias_col=alias_col, **kwargs
        )
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
                rules_compiled[name] = compiler.compile_rule(
                    expand_macros(expr, definitions)
                )
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
        if group_col is None:
            raise RuleError(
                "No Group column: expected one of"
                f" {[DEFAULT_GROUPBY_COLUMN, *group_col_backups]}"
            )
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
            label_col=label_col,
            group_col=group_col,
        )


def find_label_column(columns: Sequence[str]) -> str:
    """The rule sheet's label column: the first of `name`, `label`, `gene_id` present."""
    for col in [DEFAULT_LABEL_COLUMN, *BACKUP_LABEL_COLUMNS]:
        if col in columns:
            return col
    raise RuleError(
        f"Rule sheet has no label column: expected one of"
        f" {[DEFAULT_LABEL_COLUMN, *BACKUP_LABEL_COLUMNS]}, found {list(columns)}"
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
    rule_options: dict[str, dict],
    df: pl.DataFrame,
    label_col: str,
    group_col: str | None,
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
            details = "; ".join(
                f"{m}: {', '.join(ls[:5])}" for m, ls in by_mode.items()
            )
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
    rules_path: str | Path | None = None,
    rules: pl.LazyFrame = None,
    label_col: str = "name",
    alias_col: str = "alias",
    rules_col: str = "rule",
    common_rules: pl.LazyFrame = None,
    common_rules_path: str | Path | None = None,
) -> tuple[dict[str, Expr], dict[str, Expr], pl.DataFrame]:
    """
    Assumes TSV has at least a label column (`label_col`, e.g. name).
    Convention used here:
      - if `name` is non-empty: this is an OUTPUT RULE whose expression is in `rule`;
        with no `rule` column, or an empty `rule`, the label itself is the rule
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
    if label_col not in cols:
        raise ValueError(
            f"Rule sheet is missing its label column {label_col!r}. Found: {sorted(cols)}"
        )
    if rules_col not in cols:
        lf = lf.with_columns(pl.lit(None, dtype=pl.String).alias(rules_col))

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
        # A rule with no expression is its own label (e.g. a Gene ID).
        .with_columns(pl.col(rules_col).fill_null(pl.col(label_col)))
        .collect()
    )

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
    df: pl.DataFrame | pl.LazyFrame,
    count_col: str,
    besthit_cols: list[str],
    needed_features: set[str],
    additional_cols: list[str] | None = None,
    group_col: str | None = None,
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

    columns = df.collect_schema().names()
    hit_cols = [col for col in besthit_cols if col in columns]
    hit_cols += [
        col
        for col in columns
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
) -> dict[str, np.ndarray]:
    """Hit counts per needed (uppercase) Gene ID, one value per sample.

    Only Gene IDs with a hit get an array: the Evaluator reads a missing one as
    all zero, and at the gene level most Gene IDs have no hit among a genome's genes.
    """
    df = df.drop_nulls(subset=sample_col)
    sample_index = {s: i for i, s in enumerate(sample_names)}
    present_map: dict[str, np.ndarray] = {}
    for (hit,), sub in df.group_by("hit"):
        if hit not in needed_features:
            continue
        arr = present_map[hit] = np.zeros(len(sample_names))
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
        annotations: pl.DataFrame | pl.LazyFrame | None,
    ):
        self.samples = samples
        self.present_map = present_map
        self.annotations = annotations
        self.columns = (
            set(annotations.collect_schema().names())
            if annotations is not None
            else set()
        )
        self.sample_col = sample_col
        self._no_hits = np.zeros(len(samples))
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
                return self.present_map.get(value.upper(), self._no_hits)
            case Not(part=part):
                return np.logical_not(self.evaluate(part))
            case Combination(parts=parts):
                stacked = np.stack([self.evaluate(p) for p in parts], axis=1)
                return (
                    np.all(stacked, axis=1)
                    if expr.needs_all
                    else np.any(stacked, axis=1)
                )
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
            raise RuleError(
                f"{expr} reads annotation rows, but no annotations were given"
            )
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
        counts = np.stack(
            [np.asarray(self.evaluate(p), dtype=float) for p in parts], axis=1
        )
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
                if value not in self.columns:
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


def evaluate(
    compiled: CompiledRules,
    hits: pl.DataFrame,
    annotations: pl.DataFrame | pl.LazyFrame,
    level: EvaluationLevel,
    genomes_or_samples: list[str],
    *,
    group_col: str | None = None,
    label_col: str | None = None,
    additional_cols: list[str] | None = None,
    abundance_col: str | None = None,
    gene_abundances: pl.DataFrame | None = None,
) -> dict[str, pl.DataFrame]:
    """Evaluate every rule of `compiled` at an Evaluation Level.

    Args:
        hits: one row per (gene, Gene ID it was annotated with), as made by
            `prepare_present_map_df`: `query_id`, `genome`, `hit` (uppercase),
            plus `sample` and `abundance_col` when a mapping file is used.
        annotations: the annotation rows, for rules that read rows. Has a
            `genome` column, a `sample` column at the sample level, and
            `query_id` at the gene level.
        level: "genome" or "sample": one result per genome or sample, the rule's
            Presence (or Coverage), with abundance when `abundance_col` is given.
            "gene": each gene is evaluated, and the result is one row per rule
            with a `gene_count` struct: per genome (one field each), the number
            of its genes where the rule holds.
        genomes_or_samples: the genomes (at the gene and genome levels) or
            samples to report, in output order.
        group_col, label_col: the sheet's Group and label columns; default to
            the ones `compiled` found.
        abundance_col: column of `hits` holding the abundance to sum per rule
            (genome and sample levels only).
        gene_abundances: gene level only: `query_id` or `genome`, then one column
            per sample. Adds an `abundance` struct column, one field per sample:
            the summed abundance of the genes where the rule holds, over all
            genomes.

    Returns:
        One frame per Group, keyed by the Group's value. Rules left out by their
        `include_when` are absent; a Group with no rules left is absent.
    """
    if level not in EVALUATION_LEVELS:
        raise ValueError(f"level must be one of {EVALUATION_LEVELS}, got {level!r}")
    group_col = group_col or compiled.group_col
    label_col = label_col or compiled.label_col
    if group_col is None or group_col not in compiled.df.columns:
        raise RuleError(
            f"Rule sheet has no Group column {group_col!r}: expected one of"
            f" {[DEFAULT_GROUPBY_COLUMN, *BACKUP_GROUPBY_COLUMNS]}"
        )
    if level == "gene":
        if abundance_col:
            raise ValueError("abundance_col is not supported at the gene level")
        check_gene_level(compiled)
        return _evaluate_genes(
            compiled,
            hits,
            annotations,
            genomes_or_samples,
            group_col,
            label_col,
            additional_cols,
            gene_abundances,
        )
    if gene_abundances is not None:
        raise ValueError("gene_abundances is only supported at the gene level")
    present_map = build_present_map(
        hits, level, compiled.needed_features, genomes_or_samples
    )
    return _evaluate_units(
        compiled,
        genomes_or_samples,
        present_map,
        annotations,
        group_col=group_col,
        label_col=label_col,
        sample_col=level,
        additional_cols=additional_cols,
        anno_df=hits if abundance_col else None,
        value_col=abundance_col,
    )


def check_gene_level(compiled: CompiledRules) -> None:
    """Coverage is only reported at the genome or sample level, so a rule asking
    for it can't be evaluated at the gene level."""
    for name, options in compiled.rule_options.items():
        if coverage_mode(options):
            raise RuleError(
                f"Rule '{name}': rule_options coverage is only reported at the"
                " genome or sample level, not the gene level (summarize)"
            )


def _rules_by_group(
    compiled: CompiledRules, group_col: str, label_col: str, columns: set[str]
):
    """(group, rule names) in sheet order, leaving out rules whose
    `include_when` fails for these annotation `columns`."""
    named = compiled.df.drop_nulls(subset=[label_col])
    missing = named.filter(pl.col(group_col).is_null())[label_col]
    if len(missing):
        raise RuleError(
            "Group column, and backup group columns contain no entries! Make sure each"
            f" row for the group column contains values. Rule names with null group"
            f" columns: {', '.join(missing)}"
        )
    for (group,), frame in named.group_by(group_col, maintain_order=True):
        yield (
            group,
            [
                rn
                for rn in frame[label_col]
                if should_include(compiled.rule_options.get(rn, {}), columns)
            ],
        )


def _with_rule_columns(
    df: pl.DataFrame,
    compiled: CompiledRules,
    group_col: str,
    label_col: str,
    additional_cols: list[str] | None,
) -> pl.DataFrame:
    """Add each rule's Group (and `additional_cols`) from the sheet."""
    return df.join(
        compiled.df.select(label_col, group_col, *(additional_cols or [])),
        on=label_col,
    )


def _evaluate_genes(
    compiled: CompiledRules,
    hits: pl.DataFrame,
    annotations: pl.DataFrame | pl.LazyFrame,
    genomes: list[str],
    group_col: str,
    label_col: str,
    additional_cols: list[str] | None,
    abundances: pl.DataFrame | None,
) -> dict[str, pl.DataFrame]:
    """Gene Counts: evaluate each gene, then count per genome the genes where
    the rule holds; with `abundances`, also sum those genes' abundance per
    sample. One genome at a time, so the per-gene arrays stay small (all genes
    at once would need genes x Gene IDs values)."""
    annotations = annotations.lazy()
    columns = set(annotations.collect_schema().names())
    groups = list(_rules_by_group(compiled, group_col, label_col, columns))
    genes = annotations.select("query_id", "genome")
    samples = abundances.columns[1:] if abundances is not None else []
    if abundances is not None:
        # Keyed by query_id or genome; a gene the mapping file lacks is 0.
        genes = genes.join(
            abundances.lazy(),
            on=abundances.columns[0],
            how="left",
            maintain_order="left",
        ).with_columns(pl.col(samples).fill_null(0))
    genes = genes.collect()
    genes_by_genome = genes.partition_by("genome", as_dict=True)
    # Hits may be shared by several sheets; keep only this sheet's Gene IDs.
    hits = hits.filter(pl.col("hit").is_in(compiled.needed_features))
    hits_by_genome = hits.partition_by("genome", as_dict=True)
    no_hits = hits.clear()

    counts: dict[str, list[int]] = {rn: [] for _, names in groups for rn in names}
    # Per rule, abundance per sample summed over genomes as they are evaluated:
    # a vector per (genome, rule) would need genomes x rules x samples values.
    sums = {rn: np.zeros(len(samples)) for rn in counts}
    for genome in genomes:
        genome_genes = genes_by_genome.get((genome,), genes.clear())
        gene_names = genome_genes["query_id"].to_list()
        if samples:  # genes x samples
            gene_abundance = genome_genes.select(
                pl.col(samples).cast(pl.Float64)
            ).to_numpy()
        present_map = build_present_map(
            hits_by_genome.get((genome,), no_hits),
            "query_id",
            compiled.needed_features,
            gene_names,
        )
        # shortcut: rules that read rows filter the annotations once per genome; group them once if that gets slow.
        ev = Evaluator(
            gene_names,
            present_map,
            "query_id",
            annotations.filter(pl.col("genome") == genome),
        )
        for rn, genome_counts in counts.items():
            with _in_rule(rn):
                present = ev.evaluate(compiled.rules[rn]).astype(bool)
            genome_counts.append(count := int(present.sum()))
            if samples and count:  # most rules hold for no gene of a genome
                sums[rn] += gene_abundance[present].sum(axis=0)

    dfs = {}
    for group, names in groups:
        if not names:
            continue
        # Structs, so genome and sample names can't collide with other columns.
        df = pl.DataFrame({label_col: names}).with_columns(
            pl.DataFrame(
                [counts[rn] for rn in names],
                schema={g: pl.Int64 for g in genomes},
                orient="row",
            ).to_struct("gene_count")
        )
        if samples:
            df = df.with_columns(
                pl.DataFrame(
                    [sums[rn] for rn in names], schema=samples, orient="row"
                ).to_struct("abundance")
            )
        dfs[group] = _with_rule_columns(
            df, compiled, group_col, label_col, additional_cols
        )
    return dfs


def _evaluate_units(
    compiled: CompiledRules,
    samples: list[str],
    present_map: dict[str, np.ndarray],
    annotations: pl.DataFrame | pl.LazyFrame | None,
    group_col: str | None = "group",
    label_col: str | None = "name",
    sample_col: str | None = None,
    additional_cols: list[str] | None = None,
    anno_df: pl.DataFrame | None = None,
    value_col: str | None = None,
) -> dict[str, pl.DataFrame]:
    """Presence or Coverage per genome or sample, with abundance from
    `anno_df[value_col]` when given."""
    ev = Evaluator(
        samples=samples,
        present_map=present_map,
        sample_col=sample_col,
        annotations=annotations,
    )
    dfs = {}
    for group, names in _rules_by_group(compiled, group_col, label_col, ev.columns):
        frames = []
        for rn in names:
            options = compiled.rule_options.get(rn, {})
            expr = compiled.rules[rn]
            mode = coverage_mode(options)
            with _in_rule(rn):
                if mode:
                    out = ev.coverage(expr, mode)
                else:
                    out = pl.DataFrame(
                        {
                            "present": ev.evaluate(expr).astype(bool),
                            sample_col: ev.samples,
                        }
                    )
            out = out.with_columns(pl.lit(rn).alias(label_col))
            # anno_df/present_map hold uppercased ids (CompiledRules.needed_features)
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
                out = out.join(
                    # Sum the abundance of this rule's hits per sample,
                    anno_df.filter(pl.col("hit").is_in(features))
                    .group_by(sample_col)
                    .agg(pl.col(value_col).sum()),
                    on=sample_col,
                    how="left",
                    validate="1:1",
                ).with_columns(
                    # scaled by the rule's presence or coverage.
                    pl.col(value_col).fill_null(strategy="zero") * pl.col(check_col)
                )
            frames.append(out)

        if frames:
            dfs[group] = _with_rule_columns(
                pl.concat(frames), compiled, group_col, label_col, additional_cols
            )
    return dfs
