"""Public API for writing rule functions.

Built-in rule functions and user functions (loaded with ``--functions FILE``)
are written the same way: a Python function decorated with ``@rule_function``,
with a type hint on every parameter and on the return value. The hints decide
how each rule argument is read and what the function may be used for. See
docs/rule-functions-spec.md for the rule language.

Example::

    from dram_rules.rules_api import Group, PerSample, rule_function

    @rule_function
    def majority(group: Group) -> PerSample:
        return group.present.mean(axis=1) > 0.5
"""

from __future__ import annotations

import difflib
import enum
import importlib.util
import inspect
import operator
import typing
from contextvars import ContextVar
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Iterable, Literal

import numpy as np
import polars as pl

__all__ = [
    "ArgKind",
    "Column",
    "FunctionRegistry",
    "Group",
    "Op",
    "PerSample",
    "RowFilter",
    "Rows",
    "RuleError",
    "RuleFunction",
    "RuleFunctionError",
    "ValueKind",
    "default_functions",
    "rule_function",
]


class RuleError(Exception):
    """A problem with a rule, found while loading or evaluating rule sheets."""


class RuleFunctionError(RuleError):
    """A problem with the definition of a rule function."""


# ---------------------------------------------------------------------------
# Type hints for parameters and return values
# ---------------------------------------------------------------------------


class PerSample:
    """One value per sample: a numpy array of shape ``(n_samples,)``, holding
    hit counts or True/False. As an argument, a sub-rule evaluated per sample."""


class RowFilter:
    """True/False per annotation row: a boolean ``pl.Expr``. Functions returning
    it must be named ``filter_...`` and can only be used before ``->``."""


class Column(str):
    """The name of an annotation column (passed as a plain ``str``). Rules fail
    with a clear error if the annotations don't have it."""


@dataclass(frozen=True)
class Op:
    """A comparison operator from a rule (``>=``, ``>``, ``<=``, ``<``, ``==``,
    ``!=``). Call it like a function: ``op(values, threshold)``. Works on numpy
    arrays and Polars expressions."""

    symbol: str

    def __post_init__(self):
        if self.symbol not in _OPERATORS:
            raise RuleError(
                f"Unknown comparison operator {self.symbol!r};"
                f" expected one of {' '.join(_OPERATORS)}"
            )

    def __call__(self, left, right):
        return _OPERATORS[self.symbol](left, right)

    def __str__(self) -> str:
        return self.symbol


_OPERATORS: dict[str, Callable[[Any, Any], Any]] = {
    ">=": operator.ge,
    ">": operator.gt,
    "<=": operator.le,
    "<": operator.lt,
    "==": operator.eq,
    "!=": operator.ne,
}


@dataclass(frozen=True, eq=False)
class Group:
    """The parts of a gene group such as ``[a, b, c]`` or ``[a | b | c]``.

    Attributes:
        counts: ``(n_samples, n_parts)`` hit counts. A compound part like
            ``[x & y]`` is 1 or 0.
        labels: text of each part, for messages.
    """

    counts: np.ndarray
    labels: tuple[str, ...]

    @property
    def present(self) -> np.ndarray:
        """``(n_samples, n_parts)`` True where the part is present."""
        return self.counts > 0


@dataclass(eq=False)
class Rows:
    """The annotation rows a function summarizes: all rows, or the rows left
    after the row filters before ``->``. Filled in automatically: a parameter
    hinted ``Rows`` is not a rule argument.

    Attributes:
        df: the annotation rows.
        sample_col: name of the sample column in ``df``.
        samples: sample names, in output order.
    """

    df: pl.DataFrame
    sample_col: str
    samples: list[str]

    def per_sample(self, agg: pl.Expr, fill: Any = 0) -> np.ndarray:
        """Aggregate ``agg`` over each sample's rows. Returns an array in
        ``samples`` order. Samples with no rows (or a null result) get ``fill``,
        e.g. ``rows.per_sample(pl.len(), fill=0)`` counts rows per sample."""
        values = self.df.group_by(self.sample_col).agg(agg.alias("_value"))
        samples = pl.Series(
            self.sample_col, self.samples, dtype=self.df.schema[self.sample_col]
        )
        return (
            samples.to_frame()
            .join(values, on=self.sample_col, how="left", maintain_order="left")
            .get_column("_value")
            .fill_null(fill)
            .to_numpy()
        )

    def filter(self, predicate: pl.Expr) -> Rows:
        """The rows where `predicate` (a row filter) is true."""
        return Rows(self.df.filter(predicate), self.sample_col, self.samples)


# ---------------------------------------------------------------------------
# Function definitions
# ---------------------------------------------------------------------------


class ValueKind(enum.Enum):
    """What a rule expression evaluates to."""

    PER_SAMPLE = "per-sample value"
    ROW_FILTER = "row filter"


class ArgKind(enum.Enum):
    """How a rule argument is read, from its parameter's type hint."""

    INT = "integer"
    FLOAT = "number"
    BOOL = "true/false"
    STR = "text"
    CHOICE = "choice"
    COLUMN = "column name"
    OP = "comparison operator"
    PER_SAMPLE = "per-sample value"
    GROUP = "gene group"
    ROW_FILTER = "row filter"

    @property
    def is_expression(self) -> bool:
        """True for arguments that are sub-rules rather than literals."""
        return self in _EXPRESSION_KINDS


_EXPRESSION_KINDS = {ArgKind.PER_SAMPLE, ArgKind.GROUP, ArgKind.ROW_FILTER}

_HINT_KINDS = {
    int: ArgKind.INT,
    float: ArgKind.FLOAT,
    bool: ArgKind.BOOL,
    str: ArgKind.STR,
    Column: ArgKind.COLUMN,
    Op: ArgKind.OP,
    PerSample: ArgKind.PER_SAMPLE,
    Group: ArgKind.GROUP,
    RowFilter: ArgKind.ROW_FILTER,
}

# The return hint decides where a function can be used in a rule.
_RETURN_KINDS = {PerSample: ValueKind.PER_SAMPLE, RowFilter: ValueKind.ROW_FILTER}

# Names that are part of the rule language, not functions.
RESERVED_NAMES = {"not"}

REQUIRED = inspect.Parameter.empty


@dataclass(frozen=True)
class Param:
    """One rule argument of a function."""

    name: str
    kind: ArgKind
    choices: tuple[str, ...] = ()
    default: Any = REQUIRED
    variadic: bool = False
    keyword_only: bool = False

    @property
    def required(self) -> bool:
        return self.default is REQUIRED and not self.variadic

    def describe(self) -> str:
        kind = "one of " + "|".join(self.choices) if self.choices else self.kind.value
        prefix = "*" if self.variadic else ""
        default = "" if self.default is REQUIRED else f" = {self.default}"
        return f"{prefix}{self.name}: {kind}{default}"


@dataclass(frozen=True, eq=False)
class RuleFunction:
    """A function usable in rules, built from a ``@rule_function`` definition."""

    name: str
    fn: Callable[..., Any]
    params: tuple[Param, ...]  # rule arguments, in order (the `Rows` parameter excluded)
    rows_param: str | None  # name of the `Rows` parameter, if any
    returns: ValueKind  # from the return hint
    source: str  # "file:line", for messages
    # How to pass arguments to `fn`: these go positionally, a `*args` parameter
    # (if any) is expanded, and everything else goes by keyword.
    _positional: tuple[str, ...]
    _varargs: str | None

    def __call__(self, values: dict[str, Any]) -> Any:
        """Call the Python function with argument values by parameter name
        (a tuple for `*args`)."""
        args = [values[name] for name in self._positional]
        if self._varargs:
            args.extend(values[self._varargs])
        skip = {*self._positional, self._varargs}
        return self.fn(*args, **{k: v for k, v in values.items() if k not in skip})

    def signature(self) -> str:
        args = ", ".join(p.describe() for p in self.params)
        return f"{self.name}({args}) -> {self.returns.value}"


def rule_function(fn: Callable | None = None, *, name: str | None = None):
    """Make a Python function usable in rules.

    Use as ``@rule_function`` or ``@rule_function(name="rule_name")``. This is
    the one place every rule function is registered: into the registry of the
    ``--functions`` file being loaded, or else into `default_functions` (the
    built-ins, and functions in modules imported normally).

    The function itself is returned unchanged (no wrapper, so it keeps its
    name, docstring and signature) and can still be called directly in tests.
    """

    def decorate(func: Callable) -> Callable:
        function = _build_rule_function(func, name or func.__name__)
        _registering_into.get().add(function)
        return func

    return decorate(fn) if fn is not None else decorate


def _source(func: Callable) -> str:
    try:
        file = inspect.getsourcefile(func) or "?"
        return f"{file}:{func.__code__.co_firstlineno}"
    except (TypeError, AttributeError):
        return repr(func)


def _build_rule_function(func: Callable, name: str) -> RuleFunction:
    where = f"{_source(func)}: rule function '{name}'"
    if name in RESERVED_NAMES:
        raise RuleFunctionError(f"{where}: '{name}' is part of the rule language, not a function")
    try:
        hints = typing.get_type_hints(func)
    except Exception as e:
        raise RuleFunctionError(f"{where}: can't read its type hints: {e}") from e

    params: list[Param] = []
    rows_param = varargs = None
    positional: list[str] = []
    signature = inspect.signature(func).parameters.values()
    has_varargs = any(p.kind is inspect.Parameter.VAR_POSITIONAL for p in signature)
    for p in signature:
        if p.kind is inspect.Parameter.VAR_KEYWORD:
            raise RuleFunctionError(f"{where}: **{p.name} isn't supported")
        if p.name not in hints:
            raise RuleFunctionError(f"{where}: parameter '{p.name}' needs a type hint")
        # Parameters before a `*args` must be passed positionally; the rest by keyword.
        if p.kind is inspect.Parameter.VAR_POSITIONAL:
            varargs = p.name
        elif p.kind is inspect.Parameter.POSITIONAL_ONLY or (
            has_varargs and p.kind is inspect.Parameter.POSITIONAL_OR_KEYWORD
        ):
            positional.append(p.name)
        if hints[p.name] is Rows:
            if rows_param or p.kind is inspect.Parameter.VAR_POSITIONAL:
                raise RuleFunctionError(f"{where}: only one plain 'Rows' parameter is allowed")
            rows_param = p.name
            continue
        params.append(_param(p, hints[p.name], where))

    returns = _RETURN_KINDS.get(hints.get("return"))
    if returns is None:
        raise RuleFunctionError(f"{where}: return type hint must be PerSample or RowFilter")
    # The `filter_` prefix lets readers of a rule sheet see which calls are row
    # filters, so it must agree with the return hint.
    is_filter = returns is ValueKind.ROW_FILTER
    if name.startswith("filter_") != is_filter:
        raise RuleFunctionError(
            f"{where}: functions returning RowFilter must be named 'filter_...',"
            " and only those may be"
            + (f" (rename it to 'filter_{name}')" if is_filter else "")
        )
    return RuleFunction(
        name, func, tuple(params), rows_param, returns, _source(func),
        tuple(positional), varargs,
    )


def _param(p: inspect.Parameter, hint: Any, where: str) -> Param:
    default = p.default
    variadic = p.kind is inspect.Parameter.VAR_POSITIONAL
    keyword_only = p.kind is inspect.Parameter.KEYWORD_ONLY
    if typing.get_origin(hint) is Literal:
        choices = typing.get_args(hint)
        if not all(isinstance(c, str) for c in choices):
            raise RuleFunctionError(f"{where}: Literal choices for '{p.name}' must be strings")
        return Param(p.name, ArgKind.CHOICE, tuple(choices), default, variadic, keyword_only)
    if hint in _HINT_KINDS:
        kind = _HINT_KINDS[hint]
    else:
        supported = ", ".join(sorted(h.__name__ for h in _HINT_KINDS)) + ", Literal[...], Rows"
        raise RuleFunctionError(
            f"{where}: unsupported type hint for '{p.name}': {hint!r}. Use one of: {supported}"
        )
    if kind.is_expression and default is not REQUIRED:
        raise RuleFunctionError(f"{where}: '{p.name}' is a sub-rule and can't have a default")
    return Param(p.name, kind, (), default, variadic, keyword_only)


# ---------------------------------------------------------------------------
# Registry
# ---------------------------------------------------------------------------


class FunctionRegistry:
    """The functions available to rules, by name."""

    def __init__(self, functions: Iterable[RuleFunction] = ()):
        self._functions: dict[str, RuleFunction] = {}
        for function in functions:
            self.add(function)

    def add(self, function: RuleFunction) -> None:
        existing = self._functions.get(function.name)
        # The same definition again (e.g. a module reloaded) just replaces it.
        if existing is not None and existing.source != function.source:
            raise RuleFunctionError(
                f"{function.source}: rule function '{function.name}' is already"
                f" defined at {existing.source}"
            )
        self._functions[function.name] = function

    def get(self, name: str) -> RuleFunction:
        try:
            return self._functions[name]
        except KeyError:
            close = difflib.get_close_matches(name, self._functions, n=3)
            hint = f" Did you mean {' or '.join(close)}?" if close else ""
            raise RuleError(
                f"Unknown function '{name}'.{hint} User functions are loaded with"
                " --functions FILE."
            ) from None

    def __contains__(self, name: str) -> bool:
        return name in self._functions

    def __iter__(self):
        return iter(self._functions.values())

    def __len__(self) -> int:
        return len(self._functions)

    def with_files(self, paths: Iterable[str | Path]) -> FunctionRegistry:
        """A copy of this registry with the functions from each Python file added."""
        registry = FunctionRegistry(self)
        for path in paths:
            for function in load_functions_file(path):
                registry.add(function)
        return registry


# Built-in functions, and any defined in modules imported the normal way.
default_functions = FunctionRegistry()

# Where `@rule_function` registers. `load_functions_file` points it at a fresh
# registry while a user's file runs, so user functions never leak between runs.
_registering_into: ContextVar[FunctionRegistry] = ContextVar(
    "_registering_into", default=default_functions
)


def load_functions_file(path: str | Path) -> FunctionRegistry:
    """Run a Python file and return the ``@rule_function`` functions it defines."""
    path = Path(path)
    if not path.is_file():
        raise RuleFunctionError(f"Functions file not found: {path}")
    module_name = f"dram_rules_user_functions_{abs(hash(path.resolve()))}"
    spec = importlib.util.spec_from_file_location(module_name, path)
    module = importlib.util.module_from_spec(spec)
    registry = FunctionRegistry()
    token = _registering_into.set(registry)
    try:
        spec.loader.exec_module(module)
    except RuleError:
        raise
    except Exception as e:
        raise RuleFunctionError(f"Error loading functions file {path}: {e}") from e
    finally:
        _registering_into.reset(token)
    if not registry:
        raise RuleFunctionError(f"No @rule_function functions found in {path}")
    return registry
