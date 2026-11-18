# Changelog

All notable changes to this project will be documented in this file.

## 0.4.0 - 2026-11-18

### Bug Fixes

- Make duplicate checking have a flag ([194e4de](https://github.com/BortonWrightonLabs/dram-rules/commit/194e4dec3691edf0d63108d1bf0b9f20d3534847))


- Pass all default columns for rules parsing outside base meta columns ([e391836](https://github.com/BortonWrightonLabs/dram-rules/commit/e39183636b88c449bc797e5222cbb96823efef46))


- Match EC nums with/without EC: prefix ([51a953f](https://github.com/BortonWrightonLabs/dram-rules/commit/51a953fb16a19629c09f80d707f0ea931a279c45))


- Raise error if no group or backup group defined on row ([894d9ac](https://github.com/BortonWrightonLabs/dram-rules/commit/894d9accb8a75d5985eb930618ab302e5cc3e83a))

  If no group entry or backup group entry defined, raise error instead of silently eating the row


- Check what group col is before filling in backup group columns ([63a5397](https://github.com/BortonWrightonLabs/dram-rules/commit/63a5397863a5bacbdf809be425b13fb9414adae0))


- Backup group check fix part 2, electric boogaloo ([1b098eb](https://github.com/BortonWrightonLabs/dram-rules/commit/1b098ebf7f0d6a2065af6961738e9d7cfbfe79cb))


- Add missing group values to ag rules ([977b9d8](https://github.com/BortonWrightonLabs/dram-rules/commit/977b9d8db814cd284f7b3208713a21cf85660445))


- Evaluator bugs in not(), at_least COUNTS, missing columns and caching ([474c660](https://github.com/BortonWrightonLabs/dram-rules/commit/474c660443027603b90c7fbbc27c07bed8e11a9e))

  - not() after a filter in a pipe raised TypeError (not_() got an unexpected
    'masks' argument). not() now always negates its own argument; inside a
    pipe a negated filter is ANDed with the filters before it.
  - at_least(n, COUNTS, ...) never counted: eval_cycle cast hit counts to
    bool, so COUNTS behaved like PRESENCE. Counts are kept; path_steps still
    counts steps, not hits. The mode is now case-insensitive. This makes
    Cyc2 (rules_common.tsv) work, e.g. A. ferrooxidans now scores Iron
    Oxidation = true.
  - filter_contains on a missing column was silently ignored, so the pipe
    ran unfiltered. It now raises a RuleError like filter_compare, and the
    message shows the rule_options to skip the rule instead.
  - evaluate_rules (the traits path) now honours rule_options include_when,
    parsed once into CompiledRules.rule_options; it previously ignored it.
    probable-iron-reduction in traits_rules.tsv now requires
    kegg_description, so traits runs without KEGG skip it instead of
    failing. Existing taxonomy include_when options on traits now apply too.
  - tax/column_count_values/column_sum_values crashed with a ShapeError when
    a sample had no annotation rows; that sample is now False.
  - get_val_cols built an AssertionError without raising it.
  - functools.cache on Evaluator methods was shared by the whole class and
    kept every Evaluator (and its annotations) alive; caches are now
    per instance.


- Build taxonomy tree by lineage path ([c373cab](https://github.com/BortonWrightonLabs/dram-rules/commit/c373cab2bf3a08c0da904ebac3b18dd427a8739b))

  Partial lineages (e.g. 'Unclassified Bacteria') no longer recurse forever,
  and empty ranks shared across lineages (f__, g__) no longer merge their
  children into each other's branches.


- Mixed-case gene ids get their abundance with --mapping ([4b5d2d4](https://github.com/BortonWrightonLabs/dram-rules/commit/4b5d2d4c02949cb608a976ebd2bcda269c6b2152))

  The abundance join compared uppercased hits against the rules' original
  spelling, so ids like McrA summed to 0. Also fix --save_dataframes with
  --mapping, which assumed every table had a genome column. Adds an
  end-to-end CLI test.


- Plot genomes without a full taxonomy lineage ([3a52594](https://github.com/BortonWrightonLabs/dram-rules/commit/3a525949e8601679a7001bd4990a05849417977f))

  The taxonomy filter only matched 7-rank lineages, so genomes with a
  partial or missing taxonomy never appeared. Match any leaf of the tree
  and always keep genomes that have no taxonomy.


- Keep decimal abundances in mapping files ([0116573](https://github.com/BortonWrightonLabs/dram-rules/commit/01165730f868e30a6d6d2710ae052148dd508d55))

  Column types were guessed from the first 100 rows, and ignore_errors turned
  later decimals in an integer-looking column into null, then 0. Infer from
  all rows and let malformed files fail loudly.


- Rules that read annotation rows work with --mapping ([6496e39](https://github.com/BortonWrightonLabs/dram-rules/commit/6496e390fab0fe23b1005c4714dbf21645aac320))

  Per-sample evaluation passed the raw annotations, which have no sample
  column, so count_rows, sum_column and tax() crashed. Give it the annotation
  rows of each sample: every gene with a non-zero abundance in that sample.


- Accept a label column other than name ([975a24d](https://github.com/BortonWrightonLabs/dram-rules/commit/975a24d0f8d765f57f0c56bf0ac913614d287086))

  Rule loading honoured --label_column (and the 'label' backup), but
  evaluation and the heatmaps read 'name', so every group came out empty and
  the run crashed. Rename the chosen column to 'name' once, on load.


- Plot genomes whose lineage is a prefix of another's ([6356d05](https://github.com/BortonWrightonLabs/dram-rules/commit/6356d054eef5b46f056d66cbe61323c3c860cf0c))

  The taxonomy filter only matched leaf nodes, so a genome at an inner node
  (d__B;p__C next to d__B;p__C;c__X) never appeared. Node ids are whole
  lineages, so match each genome against every selected node.


- Round abundances in the wide tables' row key ([2f39b37](https://github.com/BortonWrightonLabs/dram-rules/commit/2f39b37e2efafd5f86314510119859fcf09eec4d))

  mean_sample_abundance and sample_abundance are part of the wide tables' row
  key. Means that are equal on paper but were summed in a different order differ
  in their last bits, which split a rule's rows and changed the row count from
  run to run. They are now rounded to 10 significant figures (ABUNDANCE_SIG_FIGS)
  in the wide tables only; the vertical tables keep the full values.

  On the merged drep bins data with its mapping file, three runs now give
  identical wide tables (eng_sys genome: 3,131 rows every run).


- Deduplicate the built-in Summarize Sheets so plain --sum runs ([8fdc3c6](https://github.com/BortonWrightonLabs/dram-rules/commit/8fdc3c657172952c0580887afa84bf64ccf8bec4))

  The built-in Summarize Sheets repeated labels and failed the duplicate
  check, and distill_bgc/distill_eng_sys wrote the Gene ID Cyc2, which the
  default Common Alias Sheet defines as an alias. Deduplicate the sheets and
  write @Cyc2, so `dram-rules --sum` runs with the built-in sheets and the
  default -c.

  Tests drop the fake built-in sheets fixture and use the real sheets: the
  CLI tests use the default -c, the topic/ecosystem tests check the real
  sheets' Groups, and tests of custom sheets alone skip the built-ins.



### Build

- Pixi tasks, macOS platforms, dependency bounds; fix README commands ([f46ead1](https://github.com/BortonWrightonLabs/dram-rules/commit/f46ead126fc0cf181e49d02fcbf19c530be2e212))

  - pixi tasks: `pixi run test`, `pixi run lint` (dev environment) and
    `pixi run viz` for the dram_viz command
  - add osx-arm64 and osx-64 to the pixi platforms (pixi.lock updated)
  - lower bounds for pyarrow and lark, which had none
  - drop the Python 3.15 classifier, which requires-python excludes
  - README: `python -m dram_viz` doesn't exist; use `dram_viz` (or
    `python -m dram_rules`), and `--output_dir`, not `--outdir`. Add pixi
    development instructions and a note on running another checkout.
  - README rules section: document rule_options coverage/include_when in
    place of path_steps/path_subunits, link the rule functions spec, and
    fix the rules file path


- Pixi task to install pre-commit hooks; move agent ignores to global ([9b8965f](https://github.com/BortonWrightonLabs/dram-rules/commit/9b8965f46235ebc24cfce68f3de68b7c8516e57a))



### Documentation

- Add rule functions language and API spec ([adee24a](https://github.com/BortonWrightonLabs/dram-rules/commit/adee24abeb6ae8f159f48a632f8f2ca5f9070af8))


- User docs for rules and rule functions; migration hints ([2231b90](https://github.com/BortonWrightonLabs/dram-rules/commit/2231b905b6180a63d73d7a2b87403c1e08cfb385))

  Add user documentation, written to also replace DRAM2's rules_parser.md:
  - docs/rules_parser.md: rule sheets, the rule language (function calls,
    not(), row filters and pipes), rule_options include_when and
    coverage, changes from earlier syntax, and troubleshooting
  - docs/rule_functions.md: reference for the built-in functions
  - docs/custom_rule_functions.md: writing, loading and testing your own
    rule functions

  Remove the design spec, now superseded by these docs.

  Rules using removed syntax now say what to write instead: path_steps,
  path_subunits, column_count_values and column_sum_values name their
  replacements, and a choice word passed before an optional choice
  parameter (e.g. at_least(2, PRESENCE, [...])) suggests mode=PRESENCE.


- Describe the dram-rules command, its steps, sheets, outputs and Gene Counts ([bd05845](https://github.com/BortonWrightonLabs/dram-rules/commit/bd058451c7bc5d1019f9fd409eebadb6fd2d51d3))

  - README: the three steps and their flags, --ecos/--topics/custom sheets,
    -lc/-acc/-gc, --mapping, genome column lookup, the output layout, the
    dashboard rules, and Summarize Sheets and Summary Collections.
  - rules_parser.md, custom_rule_functions.md: dram-rules replaces the removed
    dram_viz flags; pipeline sections note DRAM switches to dram-rules when it
    bumps the submodule.
  - rule_functions.md: what functions see when summarize evaluates each gene.


- Call visualize, traits and summarize refinements, not steps ([e5337cf](https://github.com/BortonWrightonLabs/dram-rules/commit/e5337cf3e2f537ed533d2fe2ac93e3a19e2c5eb1))

  A Step is one comma-separated part of a Rule; DRAM2 now calls these
  three Refinements of its Refine stage. Wording only.



### Features

- Add accession hmm parsing ([381d6cd](https://github.com/BortonWrightonLabs/dram-rules/commit/381d6cd88c8db71bc24334f7e3a54806884d2b3e))


- Add gene_ids to heatmap hover ([7deefbd](https://github.com/BortonWrightonLabs/dram-rules/commit/7deefbdca28887af9f96977bec31a6d58ac6e6cb))

  For each rule, catch which gene_ids are present in that rule,
  in a comma seperated list, and add that to the hover for the heatmap.

  Have to adjust some outputs for the dfs for the pivots.


- Allow use of brack or parenthesis for grouping ([22f5980](https://github.com/BortonWrightonLabs/dram-rules/commit/22f598044933b6ecf5b8cd8313854dd7a2828b08))


- Show rule name in heatmap ([e5d42af](https://github.com/BortonWrightonLabs/dram-rules/commit/e5d42af44aecf9eef4bee9b01f12662bafa9acf4))


- Rule function registry, typed rule functions and user functions ([a8e84b9](https://github.com/BortonWrightonLabs/dram-rules/commit/a8e84b946f1662ef79d569ecb14d0052b45063d2))

  Rule functions are now defined once, as plain Python functions decorated
  with @rule_function (dram_rules.rules_api). Type hints on the parameters
  decide how each rule argument is read, and the return hint (PerSample or
  RowFilter) decides where the function can be used. Built-ins live in
  rule_parser/functions.py and use the same API, and users can add their
  own with --functions FILE (repeatable). See docs/rule-functions-spec.md.

  Compiling rules (RuleCompiler) matches each call to its function,
  converts its arguments, and checks placement when the sheets load: row
  filters only before `->`, per-sample values everywhere else. Gene IDs
  are collected from sub-rule arguments automatically. Errors name the
  rule, quote aliases as written, and give the chain of aliases an error
  is in. While evaluating, every function result is checked against its
  return hint, with a message saying what was returned and how to fix it.
  User functions get read-only views of shared arrays.

  Language changes (bundled sheets migrated, outputs unchanged):
  - keyword arguments and comparison operators (>=, >, <=, <, ==, !=)
  - at_least(n, group, mode=PRESENCE|COUNTS), mode optional
  - column_count_values -> filter_compare(...) -> count_rows(...);
    column_sum_values -> sum_column
  - filter_contains and tax match plain text unless regex=true; tax now
    honours pipe filters
  - not() is part of the language, like & and |, and works on row filters
    and per-sample values
  - path_steps/path_subunits replaced by rule_options {"coverage": ...};
    mixed coverage settings in one heatmap group are an error at load time
  - probable-iron-reduction also requires heme_regulatory_motif_count

  Also removes the unused step-order graph (trees_by_rules) and the
  networkx dependency.


- Taxonomy Label select relabels the y-axis by rank ([7449617](https://github.com/BortonWrightonLabs/dram-rules/commit/7449617dfd5dcad3ea229d1c775f4c660ad0a15c))

  With 'Show Taxonomy on Y Axis' ticked, the first genome chart is labelled
  '<rank> | <genome>' for the chosen rank (genome alone when the lineage has
  no such rank). Previously the select was ignored and genomes sharing a
  lineage collapsed onto one full-lineage row.


- One evaluate for every Evaluation Level; read every kind of rule sheet ([991aef8](https://github.com/BortonWrightonLabs/dram-rules/commit/991aef8636e1986373ea671ab1c070cd3fdf4f6c))

  - rules.evaluate replaces evaluate_cycles, evaluate_rules and
    evaluate_rules_on_anno. The genome and sample levels give Presence or
    Coverage (plus abundance), as before. The gene level gives the Gene Count
    per genome, evaluated one genome at a time so per-gene arrays stay small.
  - build_present_map only allocates Gene IDs with hits. Gene IDs without hits
    still read as zeros.
  - Rule sheets find their label column (name, label, gene_id). The rule
    column is optional, and an empty rule is its label, so Summarize Sheets
    load as they are.
  - A repeated label is always an error, and make_product names the sheet in
    the message.


- Dram-rules command with a summarize step for custom sheets ([62fb390](https://github.com/BortonWrightonLabs/dram-rules/commit/62fb390335aaa42453f04c68e8d1cf343101acc6))

  - New `dram-rules` command (cli.py; also `python -m dram_rules`) with the
    shared options and the summarize step: --summarize/--sum/-s, --sum_custom,
    --trna_path. dram_viz is unchanged.
  - make_summarize.py evaluates each custom Summarize Sheet at the gene level
    and writes the custom Summary Collection to summarize/custom_summary.xlsx
    (one worksheet per Group, sorted as distill sorts) and .tsv. Columns: the
    label, gene_description/module/header/subheader when present, then a Gene
    Count per genome. The label column is `name` when sheets disagree.
  - tRNA rows each count as their own gene, so tRNAs sharing a contig are no
    longer overcounted as in distill.
  - Coverage options at the gene level are an error naming the rule, checked
    before annotations load.
  - Adds xlsxwriter (fastexcel for tests).


- --topics and --ecos choose built-in Summarize Sheets ([d2982fb](https://github.com/BortonWrightonLabs/dram-rules/commit/d2982fbb850f67505c63fb860a2b1f4ef7607a5d))

  Each is written as its own Summary Collection, topic_summary and
  ecosystem_summary, next to custom_summary. The defaults match distill's
  (all four topics, all five ecosystems). legacy and default in --ecos are
  visualize-only and an error for summarize. "" selects no built-in sheets.


- None, like "", turns off the built-in --topics/--ecos sheets ([0e188a1](https://github.com/BortonWrightonLabs/dram-rules/commit/0e188a14eec4b719548d3eeffbc4bbacf4c54fe8))


- Dram-rules --viz runs several visualize sheets in one run ([653df50](https://github.com/BortonWrightonLabs/dram-rules/commit/653df50064a7901ddf83afe416fc035996439136))

  --ecos and --viz_custom choose the sheets; each gets
  <output_dir>/visualize/<sheet>/ with its heatmaps and tables, which are
  now always written. Hits, the mapping file and the taxonomy load once
  for every sheet. make_product is now the visualize step module, and
  dram_viz/dram.viz point at dram-rules. A sheet with no Group column
  fails at compile time, naming the sheet.


- Dram-rules --traits, --dashboard/--port and --no_save ([0712161](https://github.com/BortonWrightonLabs/dram-rules/commit/07121619b41efd5f52a37013175a4de91bb373b4))

  --traits evaluates the built-in traits sheet as visualize does, into
  <output_dir>/traits/; --traits_custom adds sheets beside it in
  traits/default/ and traits/<stem>/. --dashboard serves the run's one
  visualize or traits sheet (more is an error before any evaluation) and
  still writes outputs; --no_save writes nothing and makes --output_dir
  optional. A Dashboard with no output_dir no longer saves on creation.


- --port without --dashboard is an error ([e3b5c3e](https://github.com/BortonWrightonLabs/dram-rules/commit/e3b5c3e4db5067b64d87498313adb21f4276ed89))


- -lc/-acc/-gc override custom sheets' columns; every step evaluates on one shared Hits load before anything is written ([2eeb38c](https://github.com/BortonWrightonLabs/dram-rules/commit/2eeb38c1ab9c66b159d03e3c4e57aeab29b9f922))

  - --custom_label_column/-lc, --custom_alias_column/-acc, --custom_group_column/-gc
    set the columns of every custom sheet (viz, traits, sum); built-in sheets keep
    the lookup. A custom sheet without an override column is an error naming it.
  - Hits are built once in the CLI for every Gene ID of every step; tRNA rows get
    their own Hits, added for summarize only.
  - visualize now only evaluates; save_visualized/serve write and serve after
    every step has evaluated, so an evaluation failure leaves no outputs. The
    dashboard is served last, after summarize is written.
  - --mapping with --sum logs that summarize ignores it.
  - make_summarize.add_trna is now read_trna (the tRNA rows alone).


- Derive ecosystem summaries from the visualize rule sheets ([#74](https://github.com/BortonWrightonLabs/dram-rules/pull/74)) ([1a74997](https://github.com/BortonWrightonLabs/dram-rules/commit/1a749973eb281c0bea67cabe64badfdc8814b6f6))

  Ecosystem summaries are now built each run from the visualize ecosystem
  sheets instead of the hand-kept distill_*.tsv sheets, which had drifted
  from them: one row per Gene ID the sheet's rules look up (aliases
  expanded, negated ones included), with the rules, modules and Groups
  using it and a description from the Description Tables.

  - Ecosystem sheets move to data/ecosystems/<eco>.tsv; the --ecos names
    are read from that folder and shared by visualize and summarize. The
    distill_*.tsv sheets are deleted.
  - Description Tables, each winning over the ones before: the dbCAN
    fam-activities table (kept byte for byte, so the whitespace hooks skip
    it), the built-in topic Summarize Sheets, the hand-curated
    gene_descriptions.tsv, then --descriptions custom tables. A Family;EC
    Gene ID is described from the sheet as <module>; <group>; <name>.
  - Built-in sheets are strict: a Gene ID or CAZy family without a
    description stops the run before anything is written, and a test
    guards it. --viz_custom sheets are not strict and, with --sum, write
    custom_derived_summary.
  - gene_descriptions.tsv describes every Gene ID the built-in sheets use,
    including CAZy's Note text for the families the 2022 dbCAN table
    lacks (CE19, CE20, GH174, GH176, GH177, GH179-GH184, GH187, PL43,
    PL44).
  - aquatic: the hydrogenase rules use the FeFe_group*/NiFe_group* Gene
    IDs instead of the misnamed hmm@ ones.
  - Topic and --sum_custom Summarize Sheets are unchanged.

  BREAKING CHANGE: ecosystem_summary's label column is now gene_id, rule
  and group are new, and pathway, header and subheader are gone. The
  built-in ecosystem sheets moved from data/rules_<eco>.tsv to
  data/ecosystems/<eco>.tsv.


- Add all to --ecos and --topics ([b1e4f2c](https://github.com/BortonWrightonLabs/dram-rules/commit/b1e4f2c8467a46797fc3f273d9c8fcd860739fe9))

  --ecos all runs every Ecosystem Sheet (the keys of ECOSYSTEMS, never legacy
  or default) and --topics all every topic Summarize Sheet (TOPICS). all is
  expanded in split_choices, so a sheet added to data/ecosystems/ is picked up
  with no CLI edit, and names given twice run once (all,gut runs gut once).
  The help texts list the names from the same constants, and a test keeps
  built-in sheet names off all, none, default and legacy.

  Summarize now skips the visualize-only --ecos names (legacy and default)
  with one warning on the dram.rules logger instead of stopping the run, so
  --ecos gut,legacy --viz --sum draws both and summarizes gut.

  A summarize run whose --viz_custom sheets share a file name now stops with a
  usage error before any sheet compiles. custom_derived_summary names its
  worksheets by file stem, so a/x.tsv and b/x.tsv used to keep only the last
  one, silently. The check is the one visualize's sheet folders already used,
  now shared as check_unique.



### Miscellaneous Tasks

- Gitignore local agent-skill docs and issue tracker ([b8218cb](https://github.com/BortonWrightonLabs/dram-rules/commit/b8218cb902e5bf420881b921ca861996437788f9))


- Measure coverage of the dram_rules package ([3e76945](https://github.com/BortonWrightonLabs/dram-rules/commit/3e76945db7533d614390490802e23aa419a5f3fb))

  --cov=dram_viz named the distribution, not the importable package, so
  coverage recorded no files.


- Ruff 0.17.5 everywhere, coverage of changed lines, networkx for golden checks ([6c0d958](https://github.com/BortonWrightonLabs/dram-rules/commit/6c0d958d82372d092b47029dc5585d884bc35870))

  - Bump ruff in pre-commit to v0.17.5 and pin the dev env to the same
    version, so local lint matches CI. Fix the 13 errors the newer ruff
    found (implicit Optional, dict() calls); no behaviour change.
  - CI and `pixi run -e dev cover` fail when a line changed since main never
    runs in a test (diff-cover, --fail-under=100). Tests for the uncovered
    lines: a NULL tRNA file, rules left out by include_when in summarize,
    and evaluate's error paths.
  - networkx in the dev env, so DRAM2's distill runs for golden checks.


- Pre-commit rejects workstation home-folder paths in files and commit messages ([8a86d32](https://github.com/BortonWrightonLabs/dram-rules/commit/8a86d32b24c020a37ccb2e79869d0de5439d2d2d))

  Two local pygrep hooks match paths under a user's home folder (Linux home or
  macOS Users): one over staged text files, one at the commit-msg stage.
  default_install_hook_types now includes commit-msg, so pre-commit install sets
  up both; existing clones need to run it again.


- Require Conventional Commit PR titles and add git-cliff config ([#73](https://github.com/BortonWrightonLabs/dram-rules/pull/73)) ([4bf906d](https://github.com/BortonWrightonLabs/dram-rules/commit/4bf906d7ea31a48bcd4b6cdf3cbc827ec59b7fdb))

  Squash merges into main take their title from the PR title (or the single
  commit's title), so the pr-title workflow checks both. cliff.toml groups the
  changelog by conventional type and renders each commit body from the raw
  message, so paragraphs starting `Word:` (e.g. "Tests:") no longer cut the body
  short.



### Performance

- Speed up present-map building and rule loading; fix pfam/cazy parsing ([09acbf4](https://github.com/BortonWrightonLabs/dram-rules/commit/09acbf41059355e1f846fb3e983e284538debd24))

  prepare_present_map_df processed every annotation column in sequence on
  the full frame (cast, explode, uppercase), then unpivoted ~47M rows before
  the join discarded almost all of them. Each column is now its own lazy
  branch, concatenated and collected once. ~6.5s -> ~0.6s on the 63-genome
  test file, cutting the bgc run from ~10.7s to ~4.9s wall. Output is
  identical for all bundled rule systems.

  load_rules now collects once and parses each unique rule string once
  (previously the whole plan, Lark parse included, ran four times), and the
  parser is cached per process. Rule strings stay as a String column rather
  than pl.Object.

  Column parsers are now a readable table (ID_PARSERS, renamed from
  ID_EXPR_DICT) of named parsers, with SUFFIX_PARSERS and DEFAULT_PARSER as
  fallbacks so unknown database columns are still searched. Uppercasing and
  "EC:" handling happen once in _match_key.

  Fixes:
  - pfam_hits/pfam_id and cazy_hits never matched rules: str.extract_all
    returns the whole match, ignoring capture groups.
  - from_rules mutated the shared BACKUP_GROUPBY_COLUMNS default, so later
    calls in the same process lost topic_ecosystem as a backup group column.


- Read mapping files fast and keep annotations lazy ([549a10c](https://github.com/BortonWrightonLabs/dram-rules/commit/549a10c707358ef7afb3c2c043e7ae295f820a0d))

  - Read mapping sample columns as Float64 with types inferred from the
    first 10k rows instead of the whole file (58 s of a 77 s run).
  - Keep only non-zero sample abundances; join by-genome mappings on
    genome instead of copying abundances onto every gene.
  - Scan annotations lazily; Rows.df collects them on first use, so runs
    whose rules never read rows never load them.
  - Move mapping loading into processing.load_mapping; add scan_tsv.

  Big annotations + mapping: 1:20 -> 0:17, peak RSS 14.5 -> 4.3 GB.
  By-genome mapping: OOM-killed at 40.9 GB -> 0:14, 4.0 GB.


- Sheets whose rules read rows share one in-memory copy of the annotations ([825b315](https://github.com/BortonWrightonLabs/dram-rules/commit/825b3159a2e431df8fd6c1de2fe812a0a252a906))

  Before, each sheet whose rules read annotation rows (tax(), count_rows,
  filters, also through aliases) read the whole annotations file again. Now the
  CLI collects the annotations once, after building the Hits, and only when some
  sheet's rule reads rows (CompiledRules.reads_rows). Visualize, traits, the
  sample level and summarize all use that copy.

  The per-sample join with a mapping file is still rebuilt for each sheet, from
  memory: it has a row per gene per sample, so holding it for the whole run could
  cost more memory than rebuilding it.

  On raw-anno-63-plus-with-taxa.tsv:
  - --viz --ecos ag,bgc,eng_sys,gut --traits: 32.4 s, 6.3 GB -> 30.4 s, 3.5 GB peak
  - --viz --ecos gut --traits --sum (7 sheets): 104-106 s, 4.0-5.1 GB -> 100.5-101 s, 3.5 GB
  - --viz --ecos ag (no rule reads rows): unchanged, 5.8 s, 2.4 GB
  All output tables are identical.

  Also adds a test for get_val_cols' multiple-value-column error; the inserted
  lines put it in the diff, so the coverage gate now checks it.



### Refactor

- Remove dead dashboard and pipeline code ([e4cbac0](https://github.com/BortonWrightonLabs/dram-rules/commit/e4cbac01a52a9a71d7c70a07921ad3ff46dcdf4d))

  Drop the hidden Sort By widget and its options, the unused x_cols path of
  heatmap(), unused parameters of join_present_map_df_to_mapping_df, and
  build_tax_tree_selected_recurse (the dashboard already selects every node
  by default).


- Review cleanups ([7c818b7](https://github.com/BortonWrightonLabs/dram-rules/commit/7c818b7d2afb3b5a34368f85b29a5803f96cd015))

  - One sample_abundances() helper builds the long non-zero abundance table
    for both the hit join and the per-sample annotation rows.
  - Rename get_sorted_dfs -> sort_rows and the save loop's sample_col -> y_col.
  - The --label_column CLI test also checks the heatmap is written.


- One Rule Sheet loader resolves label, alias and Group columns ([77a29c7](https://github.com/BortonWrightonLabs/dram-rules/commit/77a29c761b36892ea8bd0d7ae918fcb94ed769fd))

  CompiledRules.from_rules looks up the label, alias and Group columns with
  their backups for every kind of sheet and checks the column overrides.
  Summarize Sheets now read aliases from an aliases or parent column, and a
  sheet without a Group column fails when it loads. long_name copies the
  label after whitespace is stripped.


- Hits loads the annotations once and evaluates every level ([417791f](https://github.com/BortonWrightonLabs/dram-rules/commit/417791f9010a7a383068eaeba32585c1853fe2ff))

  Hits.load owns the genome column lookup, the Hits build for every sheet,
  tRNA rows, the in-memory copy for rules that read rows and the mapping
  joins; hits.evaluate(sheet, level) is all a refinement calls. tRNA rows
  now feed every refinement, never the Meta or taxonomy rows.


- Summarize() returns each Summary Collection without click ([ccc35dd](https://github.com/BortonWrightonLabs/dram-rules/commit/ccc35dd372b593288b3b95c946e9aeca26287d56))

  make_summarize.summarize(collections, hits) builds every collection's
  table; write_collection only writes. Summarize tests call it directly,
  and the evaluate-failure test uses a sheet that really fails.



### Styling

- Modern type hints and imports in the rule modules ([11fed7a](https://github.com/BortonWrightonLabs/dram-rules/commit/11fed7a3181741895c3d2403f095af8eea65acc7))

  Apply ruff's safe fixes (dict/list/tuple/set instead of typing.Dict etc.,
  collections.abc imports) to the rule parser modules and parenthesize an
  implicit string concatenation in a test. No behaviour change.


- Ruff 0.17.5 formats the Python example in custom_rule_functions.md ([44fc7f3](https://github.com/BortonWrightonLabs/dram-rules/commit/44fc7f304f6ef5545c76b74d776845aab4305e27))



### Testing

- Pin Sample Abundance summed once per Hit ([b7d714c](https://github.com/BortonWrightonLabs/dram-rules/commit/b7d714ca3f657f4998a76443b7bfe1a1be8b6ed7))

  A gene with two of a rule's Gene IDs counts twice, and a per-genome
  mapping gives every gene the genome's abundance (ADR 0003).



### Package

- Migrate rule parser to inside repo ([1318e83](https://github.com/BortonWrightonLabs/dram-rules/commit/1318e836dadb28737696a2d7ed2747b09d568389))

  rule parser was a subtree from another repo, now it is just a
  subdirectory, with the intention of renaming this repo
  something like dram-rules or something and having the
  option to produce vizusualization from the rules
  but also produce other rule parsing things


- Restructure package for conda and submodule usage ([799b9c8](https://github.com/BortonWrightonLabs/dram-rules/commit/799b9c8e78262ff9cd9e084123249d4d5eb4ce16))


- Update and lint common rules ([043337a](https://github.com/BortonWrightonLabs/dram-rules/commit/043337a7868ed67ee10e2f557768af9c0870f96f))


- Oct 21 rules update ([108c09b](https://github.com/BortonWrightonLabs/dram-rules/commit/108c09be4213f34657be9d48165c0b1506eaa5c8))


- Update rules ([51c2e2b](https://github.com/BortonWrightonLabs/dram-rules/commit/51c2e2bfb54146aef17935bc9aa4c2fd7eebf1f3))



### Rules

- Update common rules with new dram_db hmms ([2e37190](https://github.com/BortonWrightonLabs/dram-rules/commit/2e37190e4488625d53909a79833887cd1f316399))


- Update trait rules ([1c8eb14](https://github.com/BortonWrightonLabs/dram-rules/commit/1c8eb1434afeae5529c7e939812d823bcd8cb85e))


- Update common rules with new iron rules for traits ([ebae146](https://github.com/BortonWrightonLabs/dram-rules/commit/ebae146a3cd89e9c4565cb0ac548537b606f921c))


- Update common rules with new iron rules for traits ([cf5a84a](https://github.com/BortonWrightonLabs/dram-rules/commit/cf5a84a8309fd4af10ca1623792a20464976d5a5))


- Add SoxXA, SoxYZ and SoxB aliases to the common rules ([d4dfcc2](https://github.com/BortonWrightonLabs/dram-rules/commit/d4dfcc2ea8a31bbe4beae75bcd8ba07741dbf608))

  The aquatic sheet's Sulfur-Oxidation-Sox rule uses them; without them
  --ecos aquatic stopped on an undefined alias.
