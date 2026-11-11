# Writing Custom Rule Functions

If the [built-in functions](rule_functions.md) don't express what you need, you
can write your own in Python and use them in rule sheets like any built-in. The
built-in functions are written exactly the same way: see
[`functions.py`](https://github.com/BortonWrightonLabs/dram-rules/blob/main/src/dram_rules/rule_parser/functions.py)
for more examples.

You need to be comfortable with [numpy](https://numpy.org) and
[Polars](https://pola.rs). Functions run in the same environment as dram-viz,
so they can use any package dram-viz already depends on.

## A First Example

Save this as `my_functions.py`:

```python
import polars as pl

from dram_rules.rules_api import Column, Group, Op, PerSample, RowFilter, Rows, rule_function


@rule_function
def majority(group: Group) -> PerSample:
    """More than half of the group's parts are present."""
    return group.present.mean(axis=1) > 0.5


@rule_function
def filter_bitscore(col: Column, min_score: float = 60) -> RowFilter:
    """Rows whose score column is at least `min_score`."""
    return pl.col(col).cast(pl.Float64) >= min_score


@rule_function
def mean_column(rows: Rows, col: Column, op: Op, value: float) -> PerSample:
    """The mean of `col` over each sample's rows compares to `value`."""
    means = rows.per_sample(pl.col(col).cast(pl.Float64).mean(), fill=0.0)
    return op(means, value)
```

Then use the functions in a rule sheet:

```text
name	rule	group
mostly-complete	majority([K00370, K00371, K00374])	Custom
strong-hydrogenase-hits	filter_contains(kegg_description, "hydrogenase") -> filter_bitscore(kegg_bitScore, min_score=100) -> count_rows(>=, 2)	Custom
high-mean-score	filter_contains(kegg_description, "cytochrome") -> mean_column(kegg_bitScore, >=, 150)	Custom
```

and load the file with `--functions`:

```bash
dram_viz --annotations raw-annotations.tsv --rules_tsv my_rules.tsv \
    --functions my_functions.py --output_dir out --save_dataframes
```

`--functions` can be given more than once to load several files.

### In the DRAM Pipeline

DRAM does not have a dedicated parameter for custom functions yet. The
`TRAITS` and `PRODUCT_HEATMAP` processes pass extra arguments to dram-viz from
the process's `ext.args`, so a custom config can add `--functions`:

```groovy
// functions.config
process {
    withName: 'TRAITS' {
        ext.args = '--functions /absolute/path/to/my_functions.py'
    }
    withName: 'PRODUCT_HEATMAP' {
        ext.args = '--functions /absolute/path/to/my_functions.py'
    }
}
```

```bash
nextflow run BortonWrightonLabs/DRAM -r dev ... -c functions.config
```

Use an absolute path that is readable inside the container (with Apptainer,
paths under your home directory usually are). Rules that call your functions
must be in the sheets you pass with `--trait_rules_tsv` or `--viz_rules_tsv`.

## How a Function Is Defined

Decorate a function with `@rule_function`. Its Python name becomes the name
used in rules, unless you give another with `@rule_function(name="...")`. The
function is returned unchanged, so you can still call it directly, for example
in tests.

Every parameter, and the return value, needs a **type hint**. The hints are
how dram-viz knows what each rule argument is, and they are checked when the
file loads.

### What a Function Returns

The return hint decides where the function can be used in a rule:

| Return hint | Python value to return | Used |
| --- | --- | --- |
| `PerSample` | A numpy array with one number or True/False per sample, in the order of `rows.samples` (or of the group's rows). | Anywhere a gene identifier can be used. |
| `RowFilter` | A Polars expression giving True/False per annotation row, such as `pl.col(col) > 1`. | Only before `->`. |

Functions returning `RowFilter` must be named `filter_...`, and only those may
be. This lets readers of a rule sheet see which calls select rows. A
`PerSample` function with a `Rows` parameter (below) can also be used after
`->`, like `count_rows`.

`not(...)`, `&`, and `|` are part of the rule language and work on both kinds,
so your functions don't need to handle negation or combination.

### Argument Types

Each parameter's hint decides how the matching rule argument is read:

| Hint | Rule argument | Python value |
| --- | --- | --- |
| `int`, `float` | A number, e.g. `50` or `0.5` | `int`, `float` |
| `bool` | `true` or `false` (also `yes`/`no`) | `bool` |
| `str` | Text: `"quoted text"` or a single bare word | `str` |
| `Literal["A", "B"]` | One of these words, in any case, e.g. `mode=b` | the matching string, e.g. `"B"` |
| `Column` | An annotation column name | `str`; the run stops if the annotations don't have it |
| `Op` | A comparison: `>=`, `>`, `<=`, `<`, `==`, `!=` | an `Op`; call it as `op(values, threshold)` |
| `PerSample` | Any rule expression, e.g. `K00370 & K00371` or an alias | numpy array, one value per sample |
| `Group` | A gene group, e.g. `[K00370, K00371]` | `Group` (below) |
| `RowFilter` | A row filter expression | Polars expression |
| `Rows` | Nothing: filled in automatically | `Rows` (below) |

- Parameters with a default are optional, and can be given by keyword in rules:
  `filter_bitscore(kegg_bitScore, min_score=100)`.
- A `*args` parameter takes any number of extra arguments, e.g.
  `def any_of(*parts: PerSample) -> PerSample`.
- `PerSample`, `Group`, and `RowFilter` parameters can't have defaults, and
  `**kwargs` isn't supported.
- At most one parameter may be hinted `Rows`. A function with a `Rows`
  parameter reads annotation rows; one without can't be used after `->`.
- Gene identifiers inside `PerSample` and `Group` arguments are collected
  automatically, so the annotations are searched for them. Literal arguments,
  such as column names, are not.

### `Group`

The parts of a gene group, as in [Gene Groups](rule_functions.md#gene-groups):

| Member | Meaning |
| --- | --- |
| `counts` | `(n_samples, n_parts)` array of hit counts. A part that is an expression, such as `[A & B]`, is 1 or 0. |
| `present` | `counts > 0`, as True/False. |
| `labels` | The rule text of each part, for messages. |

### `Rows`

The annotation rows a function summarizes: all rows, or the rows left after the
filters before `->`.

| Member | Meaning |
| --- | --- |
| `df` | Polars DataFrame of the annotation rows. |
| `sample_col` | Name of the sample column in `df`. |
| `samples` | Sample names, in output order. |
| `per_sample(agg, fill=0)` | Evaluates the Polars aggregation `agg` for each sample. Returns a numpy array in `samples` order; samples with no rows, or a null result, get `fill`. |
| `filter(predicate)` | The rows where the Polars expression `predicate` is true, as a new `Rows`. |

For example, `rows.per_sample(pl.len())` counts rows per sample, and
`rows.per_sample(pl.col(col).max(), fill=0.0)` gives each sample's largest
value.

### `Op` and `Column`

`Op` is a comparison from the rule. Call it with the values and the threshold:
`op(means, value)` works on numpy arrays and on Polars expressions. `str(op)`
gives the symbol, such as `">="`.

`Column` arguments are plain strings. Before your function runs, dram-viz checks
that the annotations have the column, and stops with an error suggesting
`include_when` if they don't.

## What Is Checked

When the functions file loads:

- every parameter and the return value have supported type hints;
- the `filter_` prefix matches the `RowFilter` return hint;
- no name clashes with a built-in function or another loaded function, and no
  function is named `not`;
- the file defines at least one `@rule_function`.

Errors point to the file and line of the function:

```text
my_functions.py:12: rule function 'high': functions returning RowFilter must be named 'filter_...', and only those may be (rename it to 'filter_high')
```

When the rule sheets load, every call to your function is checked like a
built-in: the function exists, the arguments match its parameters, and it is
used in a valid place. While rules are evaluated, every result is checked
against the return hint, and a mismatch says what was returned and how to fix
the function:

```text
Rule 'r': high() has the return hint PerSample (my_functions.py:3) but returned a Polars expression, which is a row filter. If it is a filter, give it the return hint RowFilter and name it 'filter_high'; otherwise return one value per sample (an array of length 2).
Rule 'r': anyof() has the return hint PerSample (my_functions.py:2) but returned a single value (np.True_); expected one value per sample (an array of length 2).
```

An exception raised inside your function is reported with the rule, the call,
and the function's file and line, and the original exception is kept as its
cause, so the full traceback shows where it happened.

### Don't Change Arguments in Place

`PerSample` and `Group` values are shared between rules, so they are given to
your function as read-only arrays. Changing one in place, such as
`x[x > 1] = 1`, raises an error. Compute a new array instead, for example
`np.minimum(x, 1)`.

## Testing Functions

A decorated function is still a plain Python function, so it can be tested
directly. `Group` and `Rows` can be built by hand:

```python
import numpy as np
import polars as pl

from dram_rules.rules_api import Group, Op, Rows
from my_functions import majority, mean_column


def test_majority():
    # Two samples, three parts each.
    group = Group(counts=np.array([[1, 1, 0], [1, 0, 0]]), labels=("A", "B", "C"))
    assert majority(group).tolist() == [True, False]


def test_mean_column():
    df = pl.DataFrame({"genome": ["g1", "g1", "g2"], "score": [100.0, 200.0, 50.0]})
    rows = Rows(df, sample_col="genome", samples=["g1", "g2", "g3"])
    assert mean_column(rows, "score", Op(">="), 100).tolist() == [True, False, False]
```

`g3` has no rows, so its mean is the `fill` value, 0.0.

To check the whole path from rule text to result, run `dram_viz` on a small
annotations file with samples that should and shouldn't pass, as described in
[Testing Rule Syntax and Results](rules_parser.md#testing-rule-syntax-and-results).

## Security

A functions file is ordinary Python and runs with your permissions when it
loads. Only load functions files you wrote or trust.
