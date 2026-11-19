# DRAM Visualization Library

This directory contains the visualization code for DRAM. The visualization code is written in Python and uses the [Bokeh](https://bokeh.org) and [Panel](https://panel.holoviz.org) libraries. The visualization code is used to generate the figures
and dashboards for the [DRAM v2 gene annotation tool](https://github.com/WrightonLabCSU/DRAM).

## Nextflow Integration

The DRAM Visualization Library can be integrated into the DRAM Nextflow pipeline to generate figures in the form of html files. The html files can be viewed in a web browser and can be used to explore the results of the DRAM gene annotation software. This can be incorporated into a larger DRAM Nextflow pipeline or run as a standalone DRAM Nextflow
pipeline by by running the `nextflow run` (To install DRAM with nextflow see: https://github.com/WrightonLabCSU/DRAM/tree/dev) command:


```bash
nextflow run DRAM --product --annotations <path/to/annotations.tsv> --outdir <path/to/output/directory/>
```

more options can be found by running:

```bash
nextflow run DRAM --product --help
```

DRAM's processes call this library's older `dram_viz` options. They switch to the `dram-rules` command described below when DRAM bumps its `rules` submodule.
## Standalone Usage

The DRAM Visualization Library can also be used as a standalone Python package to generate figures and dashboards. First you will need to install the DRAM Visualization Library as a standalone package.

### Installation

The DRAM Visualization Library can also be used as a standalone Python package to generate figures and dashboards. To install the DRAM Visualization Library stable release, in whatever environment you are using, run:

```
conda install conda-forge::dram-viz
```

or with pip

```bash
pip install git+https://github.com/WrightonLabCSU/dram-viz.git
```

This will install the DRAM Visualization Library and all of its dependencies from the main branch of the GitHub repository.

#### Development Installation

To work on the DRAM Visualization Library, clone the repository and cd into it. With [pixi](https://pixi.sh) (recommended), the environments are defined in `pyproject.toml` and the package is installed in editable mode, so your changes take effect without reinstalling:

```bash
pixi run test                 # run the tests (in the dev environment)
pixi run lint                 # run ruff
pixi run viz --help           # run dram-rules --viz
pixi shell -e dev             # or open a shell in the dev environment
```

Without pixi, create your own environment and install the package in editable mode with the development dependencies:

```bash
pip install -e '.[dev]'
```

Note that an environment always runs the copy of the code it was installed from. To run a different checkout (for example a `git worktree` of an older commit) with the same environment, put that checkout first on the path: `PYTHONPATH=<checkout>/src dram-rules ...`.

We use [pre-commit](https://pre-commit.com) to manage our pre-commit hooks (e.g., linting, formatting that run commit to keep code formatting in check). To install the pre-commit hooks, run:

```bash
pre-commit install
```

now pre-commit will run automatically on git commit, but you can run it manually by running:

```bash
pre-commit run --all-files
```

PRs into `main` are squash merged, so the PR title is the commit title, and the `pr-title` check requires a [Conventional Commit](https://www.conventionalcommits.org/) one. The changelog is built from these commits by [git-cliff](https://git-cliff.org) (2.14.1+, `cliff.toml`):

- `rules:` is for changes to the built-in rule sheets (`src/dram_rules/data/`), such as `rules: update antibiotic rules for gut`. They get their own Rules section.
- `style:` is for formatting only and stays out of the changelog. A formatting commit with another type can be listed in `.git-blame-ignore-revs`, which git-cliff and `git blame` both skip.
- `fix(internal):` is a fix users can't see (code, not rule results or the CLI). It goes under an Internal subheading. If a change could alter outputs, make it a plain `fix:`, and split mixed changes into separate PRs.
- A fix for a bug that hasn't been released yet ends its message with a `Fixes: <sha>` trailer naming the commit that introduced it. The changelog leaves it out when that commit is in the same release, since users never saw the bug.

### Usage

One command, `dram-rules`, runs any mix of three refinements on a DRAM annotations TSV. `dram_viz` and `python -m dram_rules` run the same command. The annotations are loaded once and shared by every refinement. Every rule sheet is checked before any refinement starts, so a mistake in a sheet stops the run before the slow work.

| Refinement | Flags | What it gives, per genome |
|---|---|---|
| visualize | `--visualize`, `--viz`, `-v` | Heatmaps and tables of each rule's Presence or Coverage |
| traits | `--traits`, `-t` | The same, for the built-in traits sheet |
| summarize | `--summarize`, `--sum`, `-s` | Summary Collections: each rule's Gene Count |

Choose at least one refinement. For example, to visualize the built-in gut sheet and score traits:

```bash
dram-rules --annotations <path/to/annotations.tsv> --output_dir <path/to/output/directory/> --viz --traits
```

#### Choosing Rule Sheets

- `--ecos` picks Ecosystem Sheets, comma separated: one per file in `src/dram_rules/data/ecosystems/`, named after it (such as `gut` or `bgc`; `dram-rules --help` lists them all), or `all` for every one of them. A sheet named twice, as in `all,gut`, runs once. Visualize and summarize share it: visualize draws each sheet, and summarize [derives its summary](#derived-summarize-sheets) from the same sheet. When it is left out, visualize uses `default` (the gut sheet) and summarize uses all of them. `legacy` (the older all-in-one sheet) and `default` are for visualize only: summarize skips them with a warning.
- `--topics` picks built-in topic Summarize Sheets for summarize: such as `energy` or `env` (`dram-rules --help` lists them all), or `all` for every one of them. `default`, the default, is also every one.
- `--viz_custom`, `--traits_custom` and `--sum_custom` each add your own sheets, comma separated, to their refinement. They add to the built-in sheets rather than replace them. Pass `--ecos none` or `--topics none` to run only your own sheets. With `--sum`, each `--viz_custom` sheet also gets a derived summary.
- `--descriptions` adds your own [Description Tables](#derived-summarize-sheets), comma separated.
- `--custom_label_column`/`-lc`, `--custom_alias_column`/`-acc` and `--custom_group_column`/`-gc` name the label, alias and Group columns of every custom sheet, when those columns aren't named as the [Rules File](#rules-file) section describes. Built-in sheets always use the usual names.
- `--common_rules_tsv`/`-c` replaces the Common Alias Sheet (aliases every sheet can use), and `--functions` loads [your own rule functions](docs/custom_rule_functions.md).

For example, to visualize the gut sheet and your own, and summarize two topics, the gut ecosystem and your own Summarize Sheets:

```bash
dram-rules --annotations <path/to/annotations.tsv> --output_dir <path/to/output/directory/> \
    --viz --ecos gut --viz_custom <my_rules.tsv> \
    --sum --topics cell,env --sum_custom <sheet1.tsv>,<sheet2.tsv> [--trna_path <trnas.tsv>]
```

`--ecos gut` applies to both refinements here: visualize draws the gut sheet and summarize derives its summary from it. `--trna_path` takes the tRNAscan output from DRAM annotate. Each tRNA becomes a gene of its genome for every refinement, so rules can name tRNAs (as `` `Ala (AGC)` ``) in any sheet. With a mapping file keyed by gene, tRNAs belong to no sample, because the mapping file doesn't name them; keyed by genome, they take their genome's abundance.

`--mapping`/`-m` adds abundances per sample to visualize and traits. Its first column holds gene names (matching the annotations' `query_id`) or genome names, and each numeric column after it is a sample. Summarize ignores the mapping file for now.

Instead of a mapping file, the annotations can carry the abundances as Sample Columns: one column per sample, headed `sample:<name>` (DRAM writes them when given a mapping file; `join_mapping` in `dram_rules.processing.process_annotations` attaches one in Python). They are read as genome abundances when every gene of each genome has the same values, so a genome's tRNAs still count in its samples, and as gene abundances otherwise. Sample Columns are never searched for Gene IDs. When `--mapping` is given too, it is used and the Sample Columns are ignored, with a warning.

Genomes are named by the annotations' `input_fasta` column, or else `fasta` or `sample`. `--fasta_column`/`-f` names another column.

#### Outputs

```
<output_dir>/visualize/<sheet>/   heatmap HTML, wide (df_genome.tsv) and vertical (df_vertical_genome.tsv) tables
<output_dir>/traits/              the same, for the built-in traits sheet
<output_dir>/traits/default/      the built-in traits sheet, when --traits_custom is given
<output_dir>/traits/<sheet>/      each --traits_custom sheet
<output_dir>/summarize/           topic_summary, ecosystem_summary, custom_summary and custom_derived_summary, each as .xlsx and .tsv
```

`<sheet>` is the built-in name (such as `gut`, or `default` when `--ecos` is left out) or the custom file's name without `.tsv`. With a mapping file, visualize also writes the sample tables (`df_sample.tsv`, `df_vertical_sample.tsv`). `--save_needed_features` writes each visualize sheet's Gene IDs to `needed_features.txt` in its folder. A Summary Collection is only written when it has sheets.

Every refinement is evaluated before anything is written. If evaluation fails, the run writes no outputs at all.

#### Dashboard

To serve a live dashboard, run one visualize or traits sheet with `--dashboard`:

```bash
dram-rules --annotations <path/to/annotations.tsv> --output_dir <path/to/output/directory/> --viz --dashboard
```

Then open http://localhost:5006 (`--port` picks another port). The run still writes its outputs. `--no_save` writes nothing, and then `--output_dir` can be left out; it can't be used with `--summarize`. `--dashboard` needs exactly one visualize or traits sheet in the run, so any other number is an error.

Use:

```bash
dram-rules --help
```

to see all options.

### SSH Tunneling

If you are using the DRAM Visualization Library as a standalone Python package, you can run the dashboard on a remote server and use SSH tunneling to view the dashboard on your local machine. This will allow you to avoid downloading large data files to your local machine. To do this, first launch the dashboard on the remote server by ssh'ing into the server, navigating to the DRAM visualization directory, and running the above dashboard command. Then, on your local machine, run the following command:

```bash
ssh -NfL localhost:5006:localhost:5006 <username>@<remote-server>
```
and navigate to http://localhost:5006 to view the dashboard.

When you are finished viewing the dashboard, you should kill the process on the remote server by hitting `Ctrl+C` and then locally closing the SSH tunnel by running:

```bash
kill $(lsof -ti:5006)
```

### Rules File

The DRAM Visualization Library uses a rules file to generate the figures. This rules file is a TSV rules file similar to the DRAM traits rules. Information on the general rules parsing can be found [on the DRAM Rules Parsing Page](https://dramit.readthedocs.io/en/latest/rules_parser.html). The rules files for the DRAM Visualization Library are in `src/dram_rules/data/` (e.g. `rules.tsv`). You can use this rule file as a template for your own rules. The rules file contains the following columns:


| name | long_name | alias | rule | rule_options | group |
|---|---|---|---|---|---|
| Name to appear for x column of heatmap | Optional name to add to heatmap hover | Optional rule alias to allow aliasing one rule line to another to break it up | The actual rule | Optional JSON settings for the rule (see below) | The Group: one heatmap per Group, shown next to each other |

The label column is the first of `name`, `label` and `gene_id` in the sheet, the alias column the first of `alias`, `aliases` and `parent`, and the Group column the first of `group`, `topic_ecosystem` and `module`. A row with no `rule` is its own rule: its label is the Gene ID to look for.

Example of a rules file (excerpt from `src/dram_rules/data/rules.tsv`):

```tsv
name	long_name	alias	rule	rule_options	group
M00422	"Acetyl-CoA pathway, CO2 => acetyl-CoA"	m422	K00192 & K00195,K00193 & K00197 & K00194	{"coverage":"steps"}	Module
M00150	"Fumarate reductase, prokaryotes"	m150	K00244 & K00245 & K00246 & K00247	{"coverage":"subunits"}	Complex II
tetrathionate => thiosulfate			K08357		Sulfur
```

By default, rules evaluate as True/False (presence or absence). With `rule_options` `{"coverage":"steps"}`, a rule is instead shown as the percentage of its steps that are present: use it for pathways whose steps are separated by commas. `{"coverage":"subunits"}` does the same for subunits separated by ANDs (`&`). All rules in one group must use the same coverage setting. `{"include_when":{"all_columns_present":[...]}}` (or `all_columns_absent`) leaves a rule out entirely when the annotations don't have (or do have) the given columns.

The full rule language is described in [docs/rules_parser.md](docs/rules_parser.md), the functions available in rules (`percent`, `at_least`, `tax`, `filter_contains`, `count_rows`, ...) in [docs/rule_functions.md](docs/rule_functions.md), and how to write your own and load them with `--functions` in [docs/custom_rule_functions.md](docs/custom_rule_functions.md).

### Summarize Sheets

A Summarize Sheet is a rule sheet whose results are Gene Counts: for each genome, how many of its genes (annotation rows) the rule finds. A rule whose parts are on different genes counts 0. The built-in topic Summarize Sheets are in `src/dram_rules/data/summarize_sheets/`. They usually have a `gene_id` column and no `rule` column, so each Gene ID is its own rule, with its Group in `topic_ecosystem`:

| gene_id | gene_description | module | header | subheader | topic_ecosystem |
|---|---|---|---|---|---|
| The label, and the rule when there is no `rule` column | Optional | Optional | Optional | Optional | The Group |

Any rule sheet works as a Summarize Sheet, with its columns found as for any [rule sheet](#rules-file). Coverage settings in `rule_options` are an error, since coverage is only reported per genome. A rule name repeated within one sheet is an error in any sheet.

The sheets of one hand-written Summary Collection (topic or custom) are written together as `<collection>.xlsx`, one worksheet per sheet sorted by Group, header, subheader, module and label, and `<collection>.tsv`. A worksheet is named after its sheet: the `--topics` name (`assim`, `cell`, `energy`, `env`) or the custom file's name without its extension. So with `--sum`, the names of the `--sum_custom` files, and of the `--viz_custom` files, must make valid Excel worksheet names: at most 31 characters, none of `[]:*?/\`, and unique regardless of case. The run checks them before doing anything else. Each worksheet's columns are the label, `gene_description`, the Group (as `topic_ecosystem`, whatever the sheet calls it), `module`, `header` and `subheader`, those the sheet has, then one Gene Count column per genome. Every Summary Collection, derived ones included, orders its columns this way: the label, `gene_description`, `rule`, `topic_ecosystem`, `module`, `group`, `header`, `subheader`, then the Gene Counts. The `.tsv` has the same columns plus the worksheet name, as `sheet`, and keeps the sheets' row order. Other sheet columns (`notes`, `rule`, `alias`, ...) are dropped. If the sheets of one Collection use different label columns, the Collection labels its rows with `name`.

#### Derived Summarize Sheets

Ecosystem summaries have no hand-written sheet. Each run derives them from the visualize sheets (`--ecos` into `ecosystem_summary`, `--viz_custom` into `custom_derived_summary`), one worksheet per sheet, so the summary always lists the Gene IDs the heatmaps use. Each Gene ID any rule looks up (aliases expanded, negated ones included) is one row, in this order: by module, then Group, each in sheet order, then by where the Gene ID first appears.

| gene_id | gene_description | rule | module | group | Gene Counts |
|---|---|---|---|---|---|
| As written in the sheet | See below | Every rule using it, joined with ` \| ` | Their modules | Their Groups | One column per genome: its genes annotated with the Gene ID |

The `.tsv` also has the worksheet name, as `sheet`.

Descriptions come from the Description Tables. In order, each winning over the ones before: `dbcan.fam-activities.tsv` (an unedited dbCAN release, by CAZy family), the built-in topic Summarize Sheets (`gene_id` and `gene_description`; a `|`-joined `gene_id` describes each of its Gene IDs, and the first sheet wins where they disagree), and `gene_descriptions.tsv` (hand curated). The two tables are in `src/dram_rules/data/definitions/`. A Gene ID is looked up whole, regardless of case. A CAZy Family;EC Gene ID such as `GH18;EC:3.2.1.14` is described from the sheet instead, as `<module>; <group>; <name>` of each rule using it, but its family must be in a table. `--descriptions` adds your own tables, with `gene_id` and `description` columns; they win over the built-in ones.

The built-in sheets are strict: a Gene ID or family without a description stops the run, listing every one, before anything is written. Your `--viz_custom` sheets are not: missing descriptions are left empty, with a warning.

When you add a Gene ID to a built-in sheet, or add a sheet to `src/dram_rules/data/ecosystems/`, describe each new Gene ID in `gene_descriptions.tsv`. `test_every_built_in_gene_id_has_a_description` fails until you do, and lists what is missing. Taxa named in `tax(...)` are not Gene IDs and need no description.

Rules can be defined across multiple lines if they are enclosed in double quotes (`"`), they do not have to be enclosed in quotes if contained on a single line. White space can be used for formatting.

You can convert a KEGG module definition to a DRAM rule by first substituting the commas (`,`) for pipes (`|`), and then the pluses (`+`) for ampersands (`&`), then substitue the spaces (` `) for commas (`,`). Substitute parenthesis for brackers (`(`) for (`[`) and (`)`) for (`]`). For example, the KEGG module definition for M00422 is:

```
(K00844,K12407,K00845,K25026,K00886,K08074,K00918) (K01810,K06859,K13810,K15916) (K00850,K16370,K21071,K24182,K00918) (K01623,K01624,K11645,K16305,K16306) K01803 ((K00134,K00150) K00927,K11389) (K01834,K15633,K15634,K15635) (K01689,K27394) (K00873,K12406)
```

becomes:

```
[K00844|K12407|K00845|K25026|K00886|K08074|K00918],[K01810|K06859|K13810|K15916] [K00850|K16370|K21071|K24182|K00918],[K01623|K01624|K11645|K16305|K16306],K01803 [[K00134|K00150],K00927|K11389],[K01834|K15633|K15634|K15635],[K01689|K27394] [K00873|K12406]
```

#### No Implicit Boolean Precedence

Sometimes DRAM's rules parsing can be stricter about binary operator grouping than KEGG module definitions. This is mostly to prevent confusion on the order of ANDs and ORs with custom rules. For example, `A | B | C & D` is not a valid rule because though most parsing languages (include KEGG module definitions) would parse that as `A | B | [C & D]`, it can and has caused confusion. In DRAM's rules parsing, you would need to add brackets to make the grouping explicit: `A | B | [C & D]`. So, when converting KEGG module definitions to DRAM rules, you may need to add brackets to make the grouping explicit sometimes. See the [DRAM Rules Parsing Page](https://dramit.readthedocs.io/en/latest/rules_parser.html#core-design-principle-no-implicit-boolean-precedence) for more information on the rules parsing syntax and how to write rules.

DRAM Parsing will also warn you and error if your rules file is not grouped properly.
