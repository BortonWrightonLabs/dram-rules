# Rule Expressions and Rule Sheets

DRAM uses a small expression language to turn gene annotations and annotation
metadata into traits, summaries, and visualizations. This page describes the
rule-sheet format, the expression language, aliases, and rule options. The
functions you can call in rules are listed in
[Rule Functions](rule_functions.md), and
[Writing Custom Rule Functions](custom_rule_functions.md) explains how to add
your own in Python.

## Rule Sheets in DRAM

Rules are stored in tab-separated value (TSV) files. DRAM uses them when
evaluating traits with `--traits` and when building heatmaps with `--visualize`.
Custom sheets can be supplied with `--trait_rules_tsv` and `--viz_rules_tsv`.
Traits are evaluated through dram-viz and produce both a traits table and a
heatmap, so trait sheets can use the same options as visualization sheets.

DRAM automatically loads its common rule sheet for both traits and
visualizations. The common sheet contains reusable alias definitions, so a
custom rule sheet can reference a shared definition such as `@DNRA` without
copying its full expression. The maintained
[rule sheets and common rules sheet](https://github.com/BortonWrightonLabs/dram-rules/tree/main/src/dram_rules/data)
(`rules_common.tsv`) are useful examples when creating a custom sheet.

### Rule-Sheet Columns

| Column | Required | Purpose |
| --- | --- | --- |
| `name` | Yes | Output label for an evaluated rule. Leave it empty on a helper-only alias row. |
| `rule` | Yes | Rule expression to parse and evaluate. |
| `alias` | No | Reusable name assigned to the expression in the same row. |
| `rule_options` | No | Strict JSON with per-rule settings: when the rule is included, and how it is reported. See [Rule Options](#rule-options). |
| `long_name` | Visualization only | Longer label shown in visualization metadata. |
| `group`, `topic_ecosystem`, or `module` | Yes, one of them | Groups related rules into a heatmap. The first of these columns present is used; empty cells in it are filled from the others. Common alias sheets don't need one. |
| Other metadata | No | Columns such as `module`, or `notes` may organize downstream output. |

A row with a `name` is emitted as a result. A row with an `alias` defines a
reusable expression. A row may contain both, or it may leave `name` empty to
define a helper alias that is not emitted as a result.

Blank lines are allowed in rule sheets. Whitespace inside an expression is
ignored, so long rules can be formatted freely. A rule may also span physical
lines when the complete TSV field is enclosed in double quotes. For
portability, prefer one rule per physical line unless a long expression needs
formatting.

### Passing a Custom Sheet to DRAM

Use `--trait_rules_tsv` with `--traits`, or `--viz_rules_tsv` with
`--visualize`. Each custom sheet replaces the default sheet for that step;
the common alias sheet is still loaded separately. You can supply an existing
DRAM annotations TSV with `--annotations` to evaluate new rules without
repeating gene calling or database searches.

For example, save this small sheet as `my_rules.tsv`, using actual tabs
between columns:

```text
name	rule	long_name	group
nitrate-reducer	[K00370 & K00371] | [K02567 & K02568]	Nitrate reduction	Nitrogen
nitrite-reducer	K00362 | K03385	Nitrite reduction	Nitrogen
```

Evaluate it as a custom traits sheet:

```bash
nextflow run BortonWrightonLabs/DRAM -r dev \
    --annotations /path/to/raw-annotations.tsv \
    --traits --trait_rules_tsv /path/to/my_rules.tsv \
    --outdir /path/to/rule-test-traits \
    -profile apptainer
```

Or evaluate it as a custom visualization sheet:

```bash
nextflow run BortonWrightonLabs/DRAM -r dev \
    --annotations /path/to/raw-annotations.tsv \
    --visualize --viz_rules_tsv /path/to/my_rules.tsv \
    --outdir /path/to/rule-test-visualization \
    -profile apptainer
```

Replace the paths and select your installed runtime. For SLURM, add the
`slurm` profile and your site config, for example
`-profile apptainer,slurm -c cluster.config`. For a local checkout, replace
`BortonWrightonLabs/DRAM -r dev` with `/path/to/DRAM`. See the
[usage guide](usage.md) for custom configs and parameters files. Avoid
`full_mode` for these focused tests because it enables other pipeline steps.
When using `--viz_rules_tsv`, leave `--viz_rules_system` and `--sum_ecos`
unset so ecosystem selections do not also select a built-in visualization
system.

If your sheet needs a custom common alias sheet, add
`--trait_common_rules_tsv /path/to/my_common_rules.tsv` for traits, or
`--viz_common_rules_tsv /path/to/my_common_rules.tsv` for visualization.
These options replace the corresponding default common sheet. Include every
alias definition your custom rules depend on, including any built-in aliases
you still use. The [common alias sheet example](#common-alias-sheet-example)
below shows the format.

### Running dram-viz Directly

Outside the pipeline, the `dram_viz` command evaluates a sheet against an
annotations file:

```bash
dram_viz --annotations raw-annotations.tsv --rules_tsv my_rules.tsv \
    --output_dir rule-test --save_dataframes
```

`--common_rules_tsv` replaces the common alias sheet, `--functions` loads
[custom rule functions](custom_rule_functions.md), and `--save_needed_features`
writes the gene identifiers the rules look for. Run `dram_viz --help` for all
options.

### Testing Rule Syntax and Results

Every sheet is fully checked when it is loaded, before any rule is evaluated:
syntax, aliases, function names, function arguments, and where each function
may be used. Most mistakes therefore stop the run right away with a message
naming the rule; see [Troubleshooting](#troubleshooting).

The commands above then evaluate the sheet as part of the `TRAITS` or
`PRODUCT_HEATMAP` task. Start with a small, valid DRAM annotations TSV
containing the columns and representative samples your rules need.

When a task fails, inspect the error shown by Nextflow and the failed task's
`.command.err` and `.command.out` files in its `work/` directory. The launch
directory's `.nextflow.log` includes the task failure and work-directory path.
Fix the reported rule, then repeat the command with `-resume` to reuse eligible
completed tasks while rerunning tasks affected by the edited sheet.

A successful run confirms that the rules can be parsed and run against that
input. It does not establish that they express the intended biology. Check the
generated heatmaps and tables: visualization heatmaps and tables are written
under `VISUALIZE/`; traits write their heatmap (`heatmap_traits.html`) under
`VISUALIZE/` and their tables under `SUMMARIZE/`. Include samples expected to
satisfy each rule and samples expected to fail it. Missing feature hits can
yield a false result even when syntax is valid, and rules excluded by
`rule_options` are omitted from results.

### Common Alias Sheet Example

```text
name	alias	rule
	NitrateReduction	[K00370 & K00371] | [K02567 & K02568]
	NitriteReduction	K00362 | K03385
	Methanogen	K00399
```

### Traits Sheet Example

```text
name	rule	rule_options	topic_ecosystem
nitrate-reducer	@NitrateReduction	{}	traits
taxonomy-supported-methanogen	@Methanogen	{"include_when":{"all_columns_present":["taxonomy"]}}	traits
methanogen-without-taxonomy	@Methanogen	{"include_when":{"all_columns_absent":["taxonomy"]}}	traits
```

### Visualization Sheet Example

```text
name	rule	rule_options	long_name	group
nitrogen-reduction	@NitrateReduction, @NitriteReduction	{"coverage":"steps"}	Nitrate and nitrite reduction	Nitrogen
ferredoxin-complex	K00196 & K00197 & K00198	{"coverage":"subunits"}	Ferredoxin oxidoreductase complex	Carbon
```

### Multiline Rule Example

```text
name	rule	rule_options	group
carbon-pathway	"[K00192 | K00195],
[K00193 & K00197 & K00194]"	{"coverage":"steps"}	Carbon
```

## Identifiers and Aliases

### Gene Identifiers

A bare identifier is treated as a gene or feature identifier:

```text
K00370
nirK
EC:1.7.2.1
```

Matching is case-insensitive, and an `EC:` prefix is optional on either side:
`EC:1.7.2.1` in a rule matches `1.7.2.1` in the annotations, and the reverse.
Identifiers containing spaces or other characters outside the normal identifier
syntax must be enclosed in backticks:

```text
`multi word gene identifier`
```

Inside a function call, a bare word may instead be read as a column name, a
number, or a choice such as `COUNTS`, depending on the function's argument. For
example, in `filter_contains(kegg_description, "nitrate")`, `kegg_description`
is a column name, not a gene. Each function's arguments are listed in
[Rule Functions](rule_functions.md).

### Alias References

Reference an alias by prefixing its name with `@`:

```text
@NitrateReduction & @NitriteReduction
```

Alias names containing spaces or other special characters can be enclosed in
backticks:

```text
@`Nitrate reduction`
```

Aliases are expanded recursively before evaluation. An alias can hold any
expression, including steps, function calls, row filters, and pipes. The parser
reports an error when:

- an `@Alias` reference has no definition;
- aliases form a cycle; or
- a bare gene identifier has the same name as a defined alias.

The last check prevents a misspelled alias reference from silently being
interpreted as a gene. If an alias named `NitrateReduction` exists, use
`@NitrateReduction`; a bare `NitrateReduction` is rejected. Gene identifiers
and alias names must therefore be unique within the combined custom and common
rule sheets.

Every alias is checked on its own when the sheets load, even if no rule uses
it. Error messages quote aliases by name, as written in the sheet, and an error
inside an alias says which alias (and which aliases around it) it is in; see
[Reading Error Messages](#reading-error-messages).

### Qualified Identifiers

The grammar accepts `qualifier@identifier` for source- or namespace-qualified
feature names used by maintained rule sheets. Currently this is a convenience
feature that makes it clearer which database an identifier comes from; it does
not change rule evaluation:

```text
Fe@MtrA
Fe@MtrB_TIGR03509
```

The position of `@` is significant: `Fe@MtrA` is a qualified gene identifier,
while `@MtrA` is an alias reference.

## Expression Grammar

### Boolean Operators

| Operator | Meaning |
| --- | --- |
| `&` | Logical AND |
| `\|` | Logical OR |
| `not(...)` | Logical NOT |

Chains containing only one operator are valid:

```text
A & B & C
A | B | C
```

DRAM does not assign implicit precedence when `&` and `|` are mixed. Use square
brackets to state the intended grouping:

```text
[A & B] | C
A | [B & C]
```

Expressions such as `A | B & C` and `A & B | C` are invalid. Square brackets or
parentheses group rule expressions: `[A | B] & C` and `(A | B) & C` are the
same.

`not(...)` negates its argument: `not([K02591 | K02586])` is true for samples
with neither gene. It is part of the language, like `&` and `|`, so it also
works on [row filters](#row-filters-and-pipes).

### Step Expressions

Commas separate steps in a pathway or ordered set of evaluations:

```text
[A & B], [C | D], E
```

Boolean operators describe alternatives or required features within a step;
commas describe separate steps. Evaluated on its own, a step expression is true
when every step is present. Step expressions are commonly passed to `percent`
and `at_least`, and reported as coverage with
`{"coverage":"steps"}` (see [Rule Options](#rule-options)).

### Function Calls

A function is called with its name followed by arguments in parentheses. A
name followed by `(` is always a function call, so functions need no prefix:

```text
percent(50, [K00399, K00400, K00401])
at_least(2, [Cyc2_repCluster1, Cyc2_repCluster2], mode=COUNTS)
```

- Arguments are separated by commas. Optional arguments are given by keyword,
  as `name=value`, after all positional arguments.
- Text is written in double quotes: `"nitrate reductase"`. Single quotes are
  not supported. A single word without spaces may also be written bare, as in
  `tax(g__Methanosarcina)`.
- Comparisons are written `>=`, `>`, `<=`, `<`, `==`, or `!=`.
- True/false arguments accept `true` or `false` (also `yes`/`no`), in any case.
- A function argument that takes a gene group or sub-rule accepts any
  expression, including aliases.

The functions themselves, their arguments, and examples are listed in
[Rule Functions](rule_functions.md).

### Row Filters and Pipes

Most expressions give one value per sample: a gene is present in a genome or
it is not. Some functions instead work on the **annotation rows** (one row per
gene call) of each sample: they select rows by a column value, count rows, or
sum a column. These are combined with the pipe operator, `->`:

```text
filter_contains(kegg_description, "nitrate reductase") -> count_rows(>=, 3)
```

Read this as: take the annotation rows whose `kegg_description` contains
"nitrate reductase", then, for each sample, check whether at least 3 rows are
left.

- Everything before the last `->` must be a **row filter**: a function whose
  name starts with `filter_`, or `&`, `|`, and `not(...)` of row filters, or an
  alias holding one. Filters in a chain are applied one after another.
- The call after the last `->` must **summarize rows per sample**, such as
  `count_rows`, `sum_column`, or `tax`. These functions can also be used
  without a pipe, in which case they see all of a sample's rows.
- A pipe binds more loosely than `&` and `|`. To combine a pipe with other
  terms, put it in brackets:

```text
not(filter_contains(kegg_description, "nitrate reductase"))
    -> filter_compare(heme_regulatory_motif_count, >=, 4)
    -> count_rows(>=, 3)

[filter_contains(kegg_description, "hydrogenase") -> count_rows(>=, 2)] & K00001

filter_contains(kegg_description, "cytochrome") & filter_compare(kegg_bitScore, >=, 100)
    -> count_rows(>=, 1)
```

A row filter can't be used on its own as a rule, and a per-sample expression
can't be used before `->`; both are reported when the sheets load. By
convention, aliases that hold a row filter are named `...Filter`, such as
`@NitrateReductaseFilter`, so their role is clear wherever they are used.

## Rule Options

The optional `rule_options` column contains strict JSON. JSON object keys and
string values must use double quotes. Leave the cell empty or use `{}` when no
options are needed. Several options can be combined in one object:

```json
{"coverage":"steps","include_when":{"all_columns_present":["taxonomy"]}}
```

### `include_when`: Rule Availability

Rule availability can depend on the annotation-file schema:

```json
{"include_when":{"all_columns_present":["taxonomy"]}}
```

The rule is included only when every listed column is present.

```json
{"include_when":{"all_columns_absent":["taxonomy"]}}
```

The rule is included only when every listed column is absent. Both conditions
may be combined when they refer to different columns:

```json
{"include_when":{"all_columns_present":["taxonomy"],"all_columns_absent":["dbcan_id"]}}
```

Presence refers to a column existing in the annotation TSV, regardless of
whether individual values are empty. A rule that does not satisfy its
`include_when` condition is omitted from the output; it is not evaluated as a
false result. This applies to traits and visualizations alike.

Functions that read an annotation column, such as `filter_contains` or
`sum_column`, stop the run with an error if the column is missing. Give such
rules `include_when` with `all_columns_present` when they should be skipped for
inputs without that column; the error message suggests the exact option.

### `coverage`: Reporting Steps or Subunits

By default a rule is reported as present or absent. With `coverage`, the
heatmap and tables instead report how much of the rule is present:

| Value | Counts | Output columns |
| --- | --- | --- |
| `"steps"` | The rule's top-level steps, separated by commas (`A, B, C`). | `steps`, `steps_present`, `coverage_percentage` |
| `"subunits"` | The rule's top-level `&` parts (`A & B & C`). | `subunits`, `subunits_present`, `coverage_percentage` |

```text
name	rule	rule_options	group
M00001	[K00844 | K12407 | K00845], [K01810 | K06859], K01803	{"coverage":"steps"}	Module
M00150	K00244 & K00245 & K00246 & K00247	{"coverage":"subunits"}	Complex II
```

A step or subunit counts as present when its expression is true, however many
genes match it. A rule of another shape counts as a single part. The rule
itself is an ordinary expression, so the same rule or alias can also be used
elsewhere as present/absent.

All rules in one heatmap group must use the same `coverage` setting, because a
heatmap shows one kind of value. Mixing settings in a group, or an unknown
`coverage` value, is reported when the sheets load.

## Changes From Earlier Rule Syntax

Rule sheets written for earlier releases need these updates. The bundled sheets
have already been converted.

| Earlier syntax | Current syntax |
| --- | --- |
| `path_steps(A, B, C)` | `A, B, C` with `rule_options` `{"coverage":"steps"}` |
| `path_subunits(A & B & C)` | `A & B & C` with `rule_options` `{"coverage":"subunits"}` |
| `at_least(4, PRESENCE, [A, B, C])` | `at_least(4, [A, B, C])` |
| `at_least(2, COUNTS, [A, B, C])` | `at_least(2, [A, B, C], mode=COUNTS)` |
| `column_count_values(col, ge, 4, ge, 3)` | `filter_compare(col, >=, 4) -> count_rows(>=, 3)` |
| `column_sum_values(col, ge, 10)` | `sum_column(col, >=, 10)` |
| `gt`, `ge`, `lt`, `le`, `eq`, `ne` | `>`, `>=`, `<`, `<=`, `==`, `!=` |
| `filter_contains` and `tax` text was a regular expression | Plain text; add `regex=true` for a regular expression |

Two behaviour changes may affect results:

- `at_least(..., mode=COUNTS)` now counts every matching gene. Earlier releases
  counted each listed part at most once, so `COUNTS` behaved like `PRESENCE`.
- `rule_options` `include_when` now also applies to traits. Earlier releases
  ignored it for traits.

## Troubleshooting

### Reading Error Messages

Errors found when the sheets load name the rule or alias and quote the rule as
it is written in the sheet. An error inside an alias is prefixed with the
alias, and with the aliases around it, outermost first, so the quoted text can
be found in that alias's own definition:

```text
Rule 'hemeCount': 'K00001 & @NitrateReductaseFilter' mixes row filters and per-sample values with '&'
Alias '@NotNitrateReductaseFilter': in @NitrateReductaseFilter: filter_contains() argument 'regex' expects true or false, got maybe (expected filter_contains(col: column name, text: text, regex: true/false = False) -> row filter)
```

An error in a function's arguments ends with the function's expected
arguments, as in the second example.

### Ambiguous Boolean Operators

If an expression contains both `&` and `|`, group each mixed portion with
square brackets. Also check that every opening bracket has a closing bracket.

```text
Possible ambiguous use of '&' or '|' without surrounding brackets `[ ]`. ...
These should be written as `[A | B] & C` or `A | [B & C]`: A | B & C.
```

### Alias Errors

For an undefined alias, check spelling and capitalization and confirm that the
alias exists in either the custom or common sheet. For a bare identifier that
matches an alias, add the `@` prefix if the alias was intended or rename the
alias if the identifier is genuinely a gene ID. Break alias cycles by replacing
at least one recursive reference with its underlying expression.

```text
Rule 'r': Undefined alias: @Nope
Alias '@X': Cycle detected: Y -> X -> Y
Rule 'r': Gene ID(s) Methanogen match alias names but are written without '@'. Gene IDs and alias names must differ; did you mean @Methanogen?
```

### Unknown Functions and Wrong Arguments

Unknown function names suggest close matches. Functions, comparison words and
argument orders from an earlier release are reported with what to write
instead (see also
[Changes From Earlier Rule Syntax](#changes-from-earlier-rule-syntax)):

```text
Rule 'r': Unknown function 'at_lest'. Did you mean at_least? User functions are loaded with --functions FILE.
Rule 'r': path_steps() is no longer part of the rule language; write the steps as the rule itself (A, B, C) and give it rule_options {"coverage":"steps"}
Rule 'r': at_least(): PRESENCE is a value for 'mode', which is optional and comes later; write it as mode=PRESENCE after the other arguments (...)
Rule 'r': filter_compare() argument 'op' expects a comparison operator (>=, >, <=, <, ==, !=), got ge; use >= (...)
```

### Row Filters in the Wrong Place

```text
Rule 'r': filter_contains(kegg_description, "x") is a row filter; row filters can only be used before '->'
Rule 'r': K00001 is per-sample and can't be used before '->'; only row filters (filter_...) can
```

A row filter needs a summary after it, such as `-> count_rows(>=, 1)`. Only
row filters can come before `->`. See [Row Filters and Pipes](#row-filters-and-pipes).

### Malformed `rule_options`

Use strict, single-line JSON with double-quoted keys and values. Check matching
braces and brackets. For example, this value is missing its final brace and is
invalid:

```text
{"include_when":{"all_columns_present":["taxonomy"]}
```

The error names the rule whose `rule_options` could not be read.

### Missing Annotation Columns

Functions that read an annotation column require it. A missing column stops
the run with an error naming the rule and column. Use
`rule_options.include_when.all_columns_present` when a rule should exist only
for inputs containing those columns.

`tax` is the exception: without a `taxonomy` column it is false for every
sample, the same as genes that were not found. Give taxonomy-based rules
`include_when` when they are only meaningful with taxonomy.

### Inspecting Required Features

When running dram-viz directly, add `--save_needed_features` to write
`needed_features.txt`. This lists the gene identifiers collected from the
expanded rules and can help diagnose spelling, capitalization, qualification,
and alias-expansion problems. Column names, numbers, and other literal function
arguments are not included.
