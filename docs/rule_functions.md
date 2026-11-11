# Rule Functions

These functions are available in every rule sheet. The rule language itself
(`&`, `|`, `not(...)`, steps, aliases, and pipes) is described in
[Rule Expressions and Rule Sheets](rules_parser.md). To add functions of your
own, see [Writing Custom Rule Functions](custom_rule_functions.md).

## Overview

| Function | Kind | Purpose |
| --- | --- | --- |
| `percent(n, group)` | per-sample | At least `n`% of a group's parts are present. |
| `at_least(n, group, mode=PRESENCE)` | per-sample | At least `n` parts present, or at least `n` matching genes. |
| `tax(label, regex=false)` | per-sample, reads rows | A sample's taxonomy contains `label`. |
| `count_rows(op, n)` | per-sample, reads rows | The number of annotation rows compares to `n`. |
| `sum_column(col, op, value)` | per-sample, reads rows | A column's total compares to `value`. |
| `filter_contains(col, text, regex=false)` | row filter | Rows whose `col` contains `text`. |
| `filter_compare(col, op, value)` | row filter | Rows whose numeric `col` compares to `value`. |

- **Per-sample** functions give one value per sample, like a gene identifier,
  and can be used anywhere in a rule.
- Functions that **read rows** summarize each sample's annotation rows. Used
  after `->`, they see only the rows the filters before it kept; used on their
  own, they see all rows.
- **Row filters** (`filter_...`) select annotation rows. They can only be used
  before `->`, combined with `&`, `|`, and `not(...)`.

In the signatures below, `op` is a comparison: `>=`, `>`, `<=`, `<`, `==`, or
`!=`. `col` is the name of an annotation column. Arguments with a default, such
as `regex=false`, are optional and given by keyword.

## Gene Groups

`percent` and `at_least` take a **gene group**: a bracketed list whose parts
are separated by commas or by `|`. Each part is counted separately:

```text
[K00399, K00400, K00401]       three parts
[K00399 | K00400 | K00401]     three parts
[K00399, [K00401 & K00402]]    two parts; the second is present only when both genes are
@CytochromeCOxidase            the parts of the alias's list
K00399                         one part
```

A part that is a single gene counts that gene's matching genes in the sample (its
hit count). A part that is an expression, such as `[K00401 & K00402]`, counts
as 1 when the expression is true and 0 otherwise.

## Per-Sample Functions

### `percent(n, group)`

True when at least `n` percent of the group's parts are present. `n` is an
integer percentage.

```text
percent(50, @CytochromeCOxidase)
percent(60, [A & B, C | D, E])
```

In the second example the group has three parts: `A & B`, `C | D` and `E`.

### `at_least(n, group, mode=PRESENCE)`

True when the group's total is at least `n`. The `mode` decides what is
counted:

| Mode | What is counted | When to use it |
| --- | --- | --- |
| `PRESENCE` (default) | The number of parts that are present. Each part contributes at most one, even when its genes occur multiple times. | When the rule requires a minimum number of distinct genes, subunits, or pathway steps. |
| `COUNTS` | The total number of matching genes across the parts. Repeated hits contribute to the total. | When copy number, or the combined number of matching genes, is biologically meaningful. |

For example, suppose a genome contains two copies of `A`, no copies of `B`, and
one copy of `C`. Two parts are present, but there are three matching genes in
total. With a threshold of three, the modes differ:

```text
at_least(3, [A, B, C])
at_least(3, [A, B, C], mode=COUNTS)
```

The first expression is false because only `A` and `C` are present. The second
is true because their counts sum to three. A part that is an expression, such
as `[A | B]`, contributes at most one in both modes. `mode` is not case
sensitive.

### `tax(label, regex=false)`

True when any of the sample's annotation rows has a `taxonomy` value
containing `label`. `label` is matched as plain text anywhere in the value; with
`regex=true` it is a regular expression.

```text
tax(g__Methanosarcina)
tax("g__Methano(sarcina|coccoides)", regex=true)
filter_contains(kegg_description, "methyl-coenzyme M reductase") -> tax(p__Halobacteriota)
```

Without a `taxonomy` column, `tax` is false for every sample, the same as genes
that were not found. Give a rule `rule_options`
`{"include_when":{"all_columns_present":["taxonomy"]}}` when it should be
omitted instead.

### `count_rows(op, n)`

True when the number of the sample's annotation rows compares to `n` with `op`.
It is normally used after row filters:

```text
filter_contains(kegg_description, "nitrate reductase") -> count_rows(>=, 3)
```

A sample with no rows left has a count of 0, so `count_rows(<, 1)` is true for
it.

### `sum_column(col, op, value)`

True when the total of the numeric column `col` over the sample's annotation
rows compares to `value` with `op`. Empty cells are ignored, and a sample with
no rows has a total of 0.

```text
sum_column(heme_regulatory_motif_count, >=, 10)
filter_contains(kegg_description, "cytochrome") -> sum_column(heme_regulatory_motif_count, >=, 10)
```

## Row Filters

### `filter_contains(col, text, regex=false)`

Selects rows whose `col` contains `text`. `text` is matched as plain text
anywhere in the value; with `regex=true` it is a regular expression. Rows with
an empty value are not selected.

```text
filter_contains(kegg_description, "nitrate reductase") -> count_rows(>=, 1)
filter_contains(kegg_EC, "1.7.2.1") -> count_rows(>=, 1)
filter_contains(kegg_description, "nitr(ate|ite) reductase", regex=true) -> count_rows(>=, 1)
```

Plain-text matching means characters such as `.`, `(` and `[` match
themselves, which is usually what EC numbers and descriptions need.

### `filter_compare(col, op, value)`

Selects rows whose numeric `col` compares to `value` with `op`. The column must
hold numbers.

```text
filter_compare(heme_regulatory_motif_count, >=, 4) -> count_rows(>=, 3)
```

## Combining Filters

Row filters combine with `&`, `|`, and `not(...)` before the pipe, and
consecutive filters in a chain are applied one after another. These two rules
are equivalent:

```text
not(filter_contains(kegg_description, "nitrate reductase")) & filter_compare(heme_regulatory_motif_count, >=, 4) -> count_rows(>=, 3)
not(filter_contains(kegg_description, "nitrate reductase")) -> filter_compare(heme_regulatory_motif_count, >=, 4) -> count_rows(>=, 3)
```

## Missing Columns

`filter_contains`, `filter_compare`, and `sum_column` stop the run with an
error when the annotations don't have `col`. If the rule should be skipped for
such inputs instead, give it `rule_options`
`{"include_when":{"all_columns_present":["<col>"]}}`. `tax` and `count_rows`
never stop the run for a missing column.
