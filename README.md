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

To install the DRAM Visualization Library as a development package, clone the repository, cd into the repository, create your environment, and run:

```bash
pip install -e '.[dev]'
```

This will install the DRAM Visualization Library and all of its dependencies from the main branch of the GitHub repository, as well as the development dependencies.

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
python -m dram_viz --annotations <path/to/annotations.tsv> --outdir <path/to/output/directory/>
```

To launch a dashboard, run the following command:

```bash
python -m dram_viz --annotations <path/to/annotations.tsv> --outdir <path/to/output/directory/> --dashboard
```
This should open your default web browser and display the dashboard. If the dashboard does not open automatically, you can navigate to http://localhost:5006 to view the dashboard.

There are a number of other options available such as alternative rules to generate other dashbaords (`--rule_sysem ag` for example), ability to add abundance data through a mapping file and others. Use:

```bash
python -m dram_viz --help
```

to see all options

### Gene count summaries from visualization rules

Generate one TSV per rule system, with one row per required gene ID and one integer count column per genome:

```bash
python -m dram_rules.make_summarize \
  --annotations raw-annotations.tsv \
  --rules-system ag \
  --all-metadata \
  --output-dir summaries
```

This writes `summaries/rules_ag_summary.tsv`. Use `--rules-tsv` instead of `--rules-system` for a custom rule sheet. Aliases are expanded, and each distinct `query_id` is counted once per gene ID per genome across annotation databases. Features absent from a genome have count zero.

By default, metadata includes `topic_ecosystem` when present and the originating rule `name`. `--all-metadata` also includes columns such as `module`, `group`, `notes`, and `rule_options`. Multiple associated values are deduplicated, sorted, and joined with commas within each column. The rule expressions and alias definitions are not output as metadata.

Descriptions come from the three TSVs in `src/dram_rules/data/definitions`: `summarize_genes.tsv`, `fam-substrate-mapping.tsv`, and `dbcan.fam-activities.tsv`. The substrate mapping supplies `Family;EC:EC_Number` IDs with `Name` descriptions; rows missing a family or EC number are skipped. The activities file supplies bare family IDs from its first column and descriptions from its second. All definition IDs are matched case-insensitively. Any duplicate ID, including identical records or duplicates across files, stops the command and lists the duplicate IDs and their sources.

Use `--definitions-dir` to select another directory containing these three files. A custom `--lookup-tsv` with `gene_id` and `gene_description` columns adds or overrides descriptions after the built-in lookup passes duplicate validation. Duplicate IDs in the custom table also raise an error. Empty descriptions, punctuation-only descriptions containing colons or semicolons, and descriptions equal to their gene ID are treated as missing. EC labels are normalized to `EC:1.2.3.4`.

Missing descriptions remain blank unless `--strict` is enabled; strict mode requires a description for every feature, including CAZy features, even when the feature has no annotation hits. No summary is written if validation fails. `--save-converted-sheet` additionally writes the feature sheet without counts, and `--save-lookup` writes `master_definitions.tsv` with the final lookup. Genome names are taken from `input_fasta`, `genome`, `fasta`, or `sample`, or from the column selected with `--fasta-column`.

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

The DRAM Visualization Library uses a rules file to generate the figures. This rules file is a TSV rules file similar to the DRAM traits rules. Information on the general rules parsing can be found [on the DRAM Rules Parsing Page](https://dramit.readthedocs.io/en/latest/rules_parser.html). The rules file for the DRAM Visualization Library is located in at `dram_viz/data/rules.tsv`. You can use this rule file as a template for your own rules. The rules file contains the following columns:


| name | long_name | alia | rule | group |
|---|---|---|---|---|
| Name to appear for x column of heatmap | Optional name to add to heatmap hover | Optional rule alias to allow alliasing one rule line to another to break it up | The actual rule | Optional grouping variable to add multiple heatmaps next to each other |

Example of a rules file (excerpt from `dram_viz/data/rules.tsv`):

```tsv
name	long_name	alias	rule	group
M00422	"Acetyl-CoA pathway, CO2 => acetyl-CoA"	m422	"path_steps(K00192 & K00195,K00193 & K00197 & K00194)"	Module
M00150	"Fumarate reductase, prokaryotes"	m150	path_subunits(K00244 & K00245 & K00246 & K00247)	Complex II
tetrathionate => thiosulfate			K08357	Sulfur
```

By default, rules evuate as True/False (Presence or Absence), but the rules can also be set to evaluate as a percentage of steps present in the rule by adding either `path_steps()` or `path_subunits()` around a rule. `path_steps()` is utilized for when there are mutliple steps in a pathway seperated by a comma and you want to capture the number of steps and what percentage of steps are present. `path_subunits()` is utilized for when there are multiple subunits in a pathway seperated by ANDs (`&`) and you want to capture how many subunits are present and what percentage. You can not mix and match `path_steps()`, `path_subunits()`, and bare True/False rules in the same group, but all other rule parsing syntax is the same as the DRAM rules parsing syntax (see the DRAM Rules Parsing Page for more information).

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
