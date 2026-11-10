# Rule functions: language and API spec

Status: agreed design (see [Decisions](#decisions)). Not implemented yet.

## Goals

- One definition per function: its signature drives validation, gene collection,
  evaluation and documentation.
- Built-in and user-defined functions are written the same way.
- Rules are checked when the rule sheets load, with errors that name the rule.
- Evaluation stays vectorized (numpy / Polars); no per-sample Python loops.

Since no rule sheets exist outside this repo, old function names and signatures
are replaced, not deprecated. The bundled sheets are migrated by script, and the
full outputs for every rule system are compared before and after.

## What stays the same

Gene IDs (`K00001`, `` `odd name` ``, `hmm@fefe_group_a_hydrogenase`), aliases
(`@Name`), `&`, `|`, steps `[a, b, c]`, grouping with `[...]` or `(...)`, and
`rule_options` / `include_when`.

## Kinds of values

Every expression evaluates to one of these. Function signatures are written in
terms of them.

| Kind | Meaning | Python value |
|---|---|---|
| `PerSample` | one value per sample: hit count or True/False | `np.ndarray`, shape `(n_samples,)` |
| `Group` | the parts of a gene group, per sample | `Group` object (below) |
| `RowFilter` | True/False per annotation row (gene) | `pl.Expr` (boolean) |

- A gene ID is `PerSample` (its hit count). Wherever True/False is needed, a
  count above 0 means present.
- `&`, `|` and steps combine `PerSample` values. `&` and `|` also combine
  `RowFilter`s.
- A `Group` argument accepts `[a, b, c]`, `[a | b | c]`, an alias for either, or
  a single gene. Each part becomes one column: a gene part keeps its counts, and
  a compound part (`[x & y]`) is 1 or 0.

## Call syntax

```
name(positional, ..., keyword=value, ...)
```

- Function calls need no prefix: `name(` can never be a gene ID.
- Keyword arguments are new. Positional arguments come first.
- Literal arguments: integers, decimals, `"quoted strings"`, bare words (for
  choices like `COUNTS` and for column names), and comparison operators
  `>=  >  <=  <  ==  !=` (replacing `ge`, `gt`, ...).
- A parameter's type decides how its argument is read, so `score`, `COUNTS` and
  `K00001` are never ambiguous.

## Pipes

```
row_filter -> row_filter -> summary
```

Everything left of the last `->` must be a `RowFilter`. The filters narrow the
annotation rows, and the last call summarizes the remaining rows per sample.
Any function with a `Rows` parameter can be that summary. Without a pipe, it
sees all rows.

Summaries see *empty* samples (no rows left) as well: counts and sums are 0
there, and the comparison still applies, so `count_rows(<, 3)` is True for a
sample with no rows.

### Row filters vs per-sample values

Both use the same call syntax. They are told apart by where they may appear:

```
filter_contains(kegg_description, "hydrogenase") & filter_compare(kegg_bitScore, >=, 100) -> count_rows(>=, 2)
└──────────────────────────── row filters ──────────────────────────────────────────────┘    └ per-sample ┘
```

- **Left of `->`:** only `RowFilter`s, combined with `&`, `|` and `not`.
- **After the last `->`, and everywhere outside a pipe:** only `PerSample` values.

To make the difference visible when reading a rule, **`filter_` is reserved for
row filters.** A function returns `RowFilter` exactly when its name starts with
`filter_`, for built-in and user functions alike. The one exception is the
built-in `not`, which returns the same kind as its argument.

Aliases can't be checked this way, so by convention an alias that holds a
filter is named `...Filter` (e.g. `@NitrateReductaseFilter`).

This is enforced by two checks:

1. **When a function is defined** (`@rule_function`, for built-ins and for each
   `--functions` file): the name must start with `filter_` if and only if the
   return hint is `RowFilter`. Otherwise registration fails with an error
   pointing at the function, e.g. `my_functions.py:12: 'high_score' returns
   RowFilter, so it must be named 'filter_...'`. This makes names trustworthy.
2. **When rule sheets load** (after aliases are expanded): each piece of a rule
   is given its kind from the registry's return types, not from its name. Gene
   IDs are `PerSample`; `not(x)` takes the kind of `x`; `&` and `|` need parts of
   one kind; then the position rules above are checked. Checking kinds rather
   than names also covers aliases, `not(...)` and `&`/`|` combinations, which a
   name check would miss. Example error: `Rule 'hemeCount': count_rows(...)
   is per-sample and can't be used before '->'`.

The grammar itself doesn't distinguish the two: any `name(...)` is a call, so it
doesn't need to know user function names in advance.

Because the distinction is by kind, a gene ID could later also act as a row
filter ("rows annotated with K00001"), e.g. `K00001 -> sum_column(kegg_bitScore, >=, 100)`.
That is not planned now.

## Built-in functions

| Function | Signature | Returns |
|---|---|---|
| `not` | `not(x: PerSample \| RowFilter)` | same kind as `x` |
| `percent` | `percent(n: int, group: Group)` | `PerSample`: at least n% of parts present |
| `at_least` | `at_least(n: int, group: Group, mode: PRESENCE \| COUNTS = PRESENCE)` | `PerSample`: at least n parts present, or n total hits |
| `tax` | `tax(label: str, regex: bool = false)` + rows | `PerSample`: any row's taxonomy contains `label` |
| `filter_contains` | `filter_contains(col: Column, text: str, regex: bool = false)` | `RowFilter` |
| `filter_compare` | `filter_compare(col: Column, op: Op, value: float)` | `RowFilter` |
| `count_rows` | `count_rows(op: Op, n: float)` + rows | `PerSample`: number of rows `op` n |
| `sum_column` | `sum_column(col: Column, op: Op, value: float)` + rows | `PerSample`: column total `op` value |

Changes from today:
- `path_steps` and `path_subunits` are no longer functions. They describe how a
  rule is reported, not what it means, so they move to `rule_options`. See
  [Coverage output](#coverage-output).
- `at_least`: `mode` becomes an optional keyword, so `at_least(4, PRESENCE, [...])`
  becomes `at_least(4, [...])`, and `at_least(2, COUNTS, [...])` becomes
  `at_least(2, [...], mode=COUNTS)`. This touches 88 rules (migrated by script).
- `column_count_values(col, val_op, val, count_op, count)` is removed. It is a
  filter plus a count:
  ```
  # before
  @NotUnironDescription -> column_count_values(heme_regulatory_motif_count,ge,4,ge,3)
  # after
  @NotNitrateReductaseFilter -> filter_compare(heme_regulatory_motif_count, >=, 4) -> count_rows(>=, 3)
  ```
  The filter aliases `@UnironDescription` / `@NotUnironDescription` are renamed
  `@NitrateReductaseFilter` / `@NotNitrateReductaseFilter`, following the alias
  naming convention above.
- `column_sum_values` is renamed `sum_column`.
- `filter_contains` and `tax` match plain text (anywhere in the cell) by
  default. Today they treat their text as a regex, so `.` in EC numbers or `(`
  in descriptions can match something unintended or fail. Pass `regex=true` to
  match a regex, e.g. `tax("g__Methano(sarcina|coccoides)", regex=true)`.
- `tax` honours pipe filters. Today `filter_contains(...) -> tax(...)` crashes.

Missing columns: a `Column` argument that the annotations don't have raises a
`RuleError` naming the rule and column. To skip the rule instead, give it
`include_when` (see [include_when](#include_when) below).

`tax` is the exception. It works like a gene ID match on a level of the
taxonomy, so a missing `taxonomy` column is treated like missing gene IDs: no
sample matches (all False). Rules that are only meaningful with taxonomy should
say so with `include_when`.

## Writing a function

```python
# my_functions.py
import numpy as np
import polars as pl
from dram_rules.rules_api import Column, Group, Op, PerSample, Rows, rule_function


@rule_function
def majority(group: Group) -> PerSample:
    """More than half of the group's parts are present."""
    return group.present.mean(axis=1) > 0.5


@rule_function
def mean_column(rows: Rows, col: Column, op: Op, value: float) -> PerSample:
    """Per-sample mean of `col` over the (filtered) rows, compared with `value`."""
    means = rows.per_sample(pl.col(col).cast(pl.Float64).mean(), fill=0.0)
    return op(means, value)
```

Run with `--functions my_functions.py` (repeatable), then use them in rules:
`majority(@SomeComplex)`, `filter_contains(kegg_description, "hydrogenase") -> mean_column(kegg_bitScore, >=, 100)`.

### Rules for a function

- Decorate it with `@rule_function`. The rule name is the Python name, unless
  you pass `@rule_function(name="...")`.
- Every parameter needs a type hint, and the hint decides how the rule argument
  is read:

  | Hint | Rule argument |
  |---|---|
  | `int`, `float`, `bool`, `str` | literal |
  | `Literal["A", "B"]` | one of these bare words |
  | `Column` | a column name (`str`); must exist in the annotations |
  | `Op` | a comparison operator; called as `op(a, b)` |
  | `PerSample` | a sub-rule, evaluated per sample |
  | `Group` | a gene group (see above) |
  | `RowFilter` | a filter sub-rule (`pl.Expr`) |
  | `Rows` | not an argument: the annotation rows, filled in automatically |

- Parameters with defaults are optional and can be given by keyword. `*args`
  is allowed (e.g. `*steps: PerSample`).
- The return hint is `PerSample` or `RowFilter`. The result is checked (e.g. a
  `PerSample` must have one value per sample).
- A function returning `RowFilter` must be named `filter_...`, and only those
  may be (see [Row filters vs per-sample values](#row-filters-vs-per-sample-values)).
- Gene IDs inside `Group` and `PerSample` arguments are collected
  automatically, so the annotations are searched for them.
- A name that clashes with a built-in or another loaded function is an error
  when the functions file loads.

### `Group`

| Attribute | Meaning |
|---|---|
| `counts` | `(n_samples, n_parts)` hit counts (compound parts are 1 or 0) |
| `present` | `counts > 0` |
| `labels` | text of each part, for messages |

### `Rows`

| Member | Meaning |
|---|---|
| `df` | `pl.DataFrame` of annotation rows left after pipe filters |
| `sample_col` | name of the sample column in `df` |
| `samples` | sample names, in output order |
| `per_sample(agg, fill)` | aggregate expression per sample, as an array in `samples` order; samples with no rows get `fill` |

## When things are checked

- **When a function is defined** (built-ins at import, user functions when
  their `--functions` file loads): missing or unsupported type hints, a
  `filter_` name that doesn't match the return type, a name clash.
- **When rule sheets load:** unknown function (with "did you mean" suggestions),
  wrong number or names of arguments, literals of the wrong type, a row filter
  or per-sample value in the wrong place (see
  [Row filters vs per-sample values](#row-filters-vs-per-sample-values)), `&` or
  `|` mixing the two, a `coverage` option that doesn't fit its rule, rules in
  one heatmap group with different `coverage` settings.
- **When evaluating:** missing `Column`s, and a function's return value.

Errors name the rule, for example: `Rule 'hemeCount': count_rows() got 'ge'
for op; expected one of >=, >, <=, <, ==, !=`.

## Coverage output

A rule is normally reported as present/absent. `rule_options` `coverage` reports
how much of it is present instead (today's `path_steps` / `path_subunits`):

| `coverage` | Parts counted | Output columns |
|---|---|---|
| `"steps"` | the rule's top-level steps (`a, b, c`) | `steps`, `steps_present`, `coverage_percentage` |
| `"subunits"` | the rule's top-level `&` parts (`a & b & c`) | `subunits`, `subunits_present`, `coverage_percentage` |

```
# before
M00001    path_steps([K00844 | K12407 | ...], [K01810 | ...], ...)
# after
M00001    [K00844 | K12407 | ...], [K01810 | ...], ...    {"coverage": "steps"}
```

The rule itself stays an ordinary expression, so it still evaluates to
present/absent: as an alias, in traits, or in a heatmap group reporting
presence. This removes the separate kind of value these functions returned,
and the `allow_visualize_functions` switch.

All rules in one heatmap group must share the same `coverage` setting. This is
checked when the sheets load, instead of failing during the heatmap build.

## include_when

`rule_options` `{"include_when": {"all_columns_present": [...], "all_columns_absent": [...]}}`
leaves the rule out of all outputs when its conditions aren't met. It is not
the same as "False". A separate option for "output as False" may come later.

## Decisions

Settled on 2026-11-10:

1. `tax` without a `taxonomy` column gives False for every sample, like missing
   gene IDs. Rules that need taxonomy use `include_when`.
2. `filter_contains` and `tax` match plain text by default, with `regex=true`
   to match a regex.
3. `at_least` takes `mode` as an optional keyword, with `PRESENCE` as the default.
4. `coverage` goes in `rule_options`. It is JSON so that it can grow to hold
   settings like this.
5. User functions are loaded with `--functions FILE` (repeatable).
6. Row filters and per-sample functions share one call syntax. They are told
   apart by position: row filters only left of `->`, per-sample values
   everywhere else. This is checked by kind when rule sheets load.
7. The `filter_` prefix is reserved for functions returning `RowFilter`, built-in
   or user-defined, and is enforced when a function is defined. `not` is the
   only exception. Aliases holding a filter are named `...Filter` by convention
   (not enforced).
