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
pixi run viz --help           # run the dram_viz command
pixi shell -e dev             # or open a shell in the dev environment
```

Without pixi, create your own environment and install the package in editable mode with the development dependencies:

```bash
pip install -e '.[dev]'
```

Note that an environment always runs the copy of the code it was installed from. To run a different checkout (for example a `git worktree` of an older commit) with the same environment, put that checkout first on the path: `PYTHONPATH=<checkout>/src dram_viz ...`.

We use [pre-commit](https://pre-commit.com) to manage our pre-commit hooks (e.g., linting, formatting that run commit to keep code formatting in check). To install the pre-commit hooks, run:

```bash
pre-commit install
```

now pre-commit will run automatically on git commit, but you can run it manually by running:

```bash
pre-commit run --all-files
```

### Usage

To generate a figure like the one shown above [TODO: insert figure], run the following command:

```bash
dram_viz --annotations <path/to/annotations.tsv> --output_dir <path/to/output/directory/>
```

`python -m dram_rules ...` runs the same command.

To launch a dashboard, run the following command:

```bash
dram_viz --annotations <path/to/annotations.tsv> --output_dir <path/to/output/directory/> --dashboard
```
This should open your default web browser and display the dashboard. If the dashboard does not open automatically, you can navigate to http://localhost:5006 to view the dashboard.

There are a number of other options available such as alternative rules to generate other dashbaords (`--rules_system ag` for example), your own rule functions (`--functions my_functions.py`, see `docs/rule-functions-spec.md`), ability to add abundance data through a mapping file and others. Use:

```bash
dram_viz --help
```

to see all options

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
| Name to appear for x column of heatmap | Optional name to add to heatmap hover | Optional rule alias to allow aliasing one rule line to another to break it up | The actual rule | Optional JSON settings for the rule (see below) | Optional grouping variable to add multiple heatmaps next to each other |

Example of a rules file (excerpt from `src/dram_rules/data/rules.tsv`):

```tsv
name	long_name	alias	rule	rule_options	group
M00422	"Acetyl-CoA pathway, CO2 => acetyl-CoA"	m422	K00192 & K00195,K00193 & K00197 & K00194	{"coverage":"steps"}	Module
M00150	"Fumarate reductase, prokaryotes"	m150	K00244 & K00245 & K00246 & K00247	{"coverage":"subunits"}	Complex II
tetrathionate => thiosulfate			K08357		Sulfur
```

By default, rules evaluate as True/False (presence or absence). With `rule_options` `{"coverage":"steps"}`, a rule is instead shown as the percentage of its steps that are present: use it for pathways whose steps are separated by commas. `{"coverage":"subunits"}` does the same for subunits separated by ANDs (`&`). All rules in one group must use the same coverage setting. `{"include_when":{"all_columns_present":[...]}}` (or `all_columns_absent`) leaves a rule out entirely when the annotations don't have (or do have) the given columns.

The functions available in rules (`percent`, `at_least`, `tax`, `filter_contains`, `count_rows`, ...), and how to write your own and load them with `--functions`, are described in [docs/rule-functions-spec.md](docs/rule-functions-spec.md).

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
