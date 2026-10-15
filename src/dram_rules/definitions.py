from __future__ import annotations

DEFAULT_LABEL_COLUMN = "name"
BACKUP_LABEL_COLUMNS = ["label"]
DEFAULT_ALIAS_COLUMN = "alias"
BACKUP_ALIAS_COLUMNS = ["aliases", "parent"]
DEFAULT_FASTA_COLUMN = "input_fasta"
BACKUP_FASTA_COLUMNS = ["fasta", "sample"]
DEFAULT_GROUPBY_COLUMN = "group"
BACKUP_GROUPBY_COLUMNS = ["topic_ecosystem"]

# regex portion not used right now, but could be useful in the future
# captures in named pandas columns that remove the prefix string (e.g. "d__")
TAXONOMY_RANKS_REGEX = {
    "domain": r"(?:;?d__)(?P<domain>.*?)",
    "phylum": r"(?:;?p__)(?P<phylum>.*?)",
    "class": r"(?:;?c__)(?P<class>.*?)",
    "order": r"(?:;?o__)(?P<order>.*?)",
    "family": r"(?:;?f__)(?P<family>.*?)",
    "genus": r"(?:;?g__)(?P<genus>.*?)",
    "species": r"(?:;?s__)(?P<species>.*)",
}
NO_TAXONOMY_RANKS = len(TAXONOMY_RANKS_REGEX)
