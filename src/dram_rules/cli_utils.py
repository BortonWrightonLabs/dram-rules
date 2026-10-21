import polars as pl


def get_column_name(
    name: str,
    backup_names: list[str],
    name_of_data: str,
    df: pl.LazyFrame | pl.DataFrame = None,
    cols: list[str] = None,
    raise_error: bool = True,
) -> str:
    assert df is not None or cols is not None, (
        "Must provide either a dataframe or a list of columns to search for the column name"
    )
    if cols is None:
        if isinstance(df, pl.LazyFrame):
            cols = df.collect_schema().names()
        else:
            cols = df.columns
    if name in cols:
        return name
    for backup_name in backup_names:
        if backup_name in cols:
            return backup_name
    if raise_error:
        raise ValueError(
            f"None of the following column names were found in the {name_of_data}: {[name] + backup_names}"
        )
    else:
        return name
