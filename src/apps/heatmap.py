from __future__ import annotations

from math import pi
from pathlib import Path
from typing import Optional

import panel as pn
import param
import polars as pl
from bokeh.embed import file_html
from bokeh.layouts import row
from bokeh.models import Div
from bokeh.palettes import BuGn, Cividis256
from bokeh.plotting import figure
from bokeh.resources import INLINE
from bokeh.transform import factor_cmap, linear_cmap, log_cmap
from panel_material_ui import Tree

from src.definitions import NO_TAXONOMY_RANKS, TAXONOMY_RANKS_REGEX
from src.rule_parser.rules import get_val_col

pn.extension("tabulator", "katex", template="bootstrap")

PALETTE_CATEGORICAL = BuGn
PALETTE_CONTINUOUS = Cividis256

HEATMAP_CELL_HEIGHT = 15
HEATMAP_CELL_WIDTH = 15


def heatmap(
    df: pl.DataFrame,
    y_col,
    tooltip_cols,
    title="",
    rect_kw=None,
    c_min=0,
    c_max=1,
    c_col: str = None,
    x_col: str = None,
    x_cols: list[str] = None,
    log_colors: bool = False,
    **fig_kwargs,
):
    """
    Make a heatmap from a dataframe

    Parameters
    ----------
    df : pl.DataFrame
        The dataframe to make the heatmap from
    x_col : str
        The column to use for the x-axis
    y_col : str
        The column to use for the y-axis
    c_col : str
        The column to use for the color
    tooltip_cols : list
        A list of tuples of columns to use for the tooltips
    title : str
        The title of the plot
    rect_kw : dict
        Keyword arguments for the rect
    c_min : float
        The minimum value for the color
    c_max : float
        The maximum value for the color

    Returns
    -------
    Plot
        The heatmap plot
    """

    if x_cols:
        df = df.unpivot(
            index=y_col, on=x_cols, variable_name="x_col", value_name="value"
        ).unique(maintain_order=True)
        x_col = "x_col"
        c_col = "value"
        tooltip_cols = [y_col, "value"]

    rect_kw = rect_kw or {}

    tooltips = []
    for col in tooltip_cols:
        if isinstance(col, tuple):
            tooltips.append(col)
        else:
            tooltips.append((col.replace("_", " ").title(), f"@{col}"))

    p = figure(
        frame_width=HEATMAP_CELL_WIDTH * df[x_col].n_unique(),
        frame_height=HEATMAP_CELL_WIDTH * df[y_col].n_unique(),
        x_range=df[x_col].unique(maintain_order=True).to_list(),
        y_range=df[y_col].unique(maintain_order=True).to_list(),
        tools="hover",
        toolbar_location=None,
        tooltips=tooltips,
        # title="\n".join(char for char in title),
        title=title,
        title_location="right",
        **fig_kwargs,
    )

    color_dtype = df.schema[c_col]
    if color_dtype.is_numeric():
        palette = tuple(reversed(PALETTE_CONTINUOUS))
        if log_colors:
            fill_color = log_cmap(c_col, palette=palette, low=c_min, high=c_max)
        else:
            fill_color = linear_cmap(c_col, palette=palette, low=c_min, high=c_max)
    else:
        df = df.with_columns(pl.col(c_col).cast(pl.String))
        factors = (
            ["false", "true"]
            if color_dtype == pl.Boolean
            else df[c_col].drop_nulls().unique().sort().to_list()
        )
        max_factors = max(PALETTE_CATEGORICAL.keys())
        palette = (
            PALETTE_CATEGORICAL[max(len(factors), 3)]
            if len(factors) <= max_factors
            else PALETTE_CONTINUOUS
        )
        fill_color = factor_cmap(
            c_col, palette=tuple(reversed(palette)), factors=factors
        )
    p.rect(
        x=x_col,
        y=y_col,
        width=0.9,
        height=0.9,
        source=df,
        fill_alpha=0.9,
        color=fill_color,
        **rect_kw,
    )

    p.title.text_font_size = "8pt"

    p.grid.grid_line_color = None
    p.axis.axis_line_color = None
    p.axis.major_tick_line_color = None
    p.axis.major_label_text_font_size = "12px"
    p.xaxis.major_label_orientation = pi / 2

    return p


class Dashboard(pn.viewable.Viewer):
    """
    A class representing a dashboard for visualizing data.

    Parameters:
    -----------
    module_df : pl.DataFrame
        The module dataframe.
    etc_df : pl.DataFrame
        The etc dataframe.
    function_df : pl.DataFrame
        The function dataframe.
    tax_tree_data : optional
        The taxonomy tree data.
    """

    min_coverage = param.Number(default=0, bounds=(0, 1), label="Minimum Coverage")
    view = param.ClassSelector(class_=pn.template.FastListTemplate)
    plot_view = param.ClassSelector(class_=pn.Row)

    def __init__(
        self,
        dfs: dict[str, dict[str, pl.DataFrame]],
        taxanomy_tree_data: Optional[list[dict]] = None,
        selected_tax_tree: Optional[list[str]] = None,
        output_dir: str | Path = None,
    ):
        super().__init__()
        self.dfs = dfs
        self.taxonomy_tree_data = taxanomy_tree_data
        self._output_dir = Path(output_dir) if output_dir else Path.cwd()
        self.plot_view = pn.Row()
        self.download_button = pn.widgets.Button(
            label="Download Heatmap", color="primary"
        )
        self.download_button.on_click(self.download_heatmap)

        self.redraw_button = pn.widgets.Button(label="Redraw", color="primary")
        self.reset_button = pn.widgets.Button(label="Reset Filters", color="warning")
        self.reset_button.on_click(self.reset_filters)

        self.redraw_button.on_click(self.update_plot)

        self.tax_axis_filter = pn.widgets.Checkbox(
            label="Show Taxonomy on Y Axis", value=False
        )
        self.tax_axis_rank = pn.widgets.Select(
            label="Taxonomy Label",
            options=list(TAXONOMY_RANKS_REGEX),
            visible=False,
            value="genus",
        )
        self.show_tax_box = pn.Column(self.tax_axis_filter, self.tax_axis_rank)
        pn.bind(self.reveal_tax_axis_rank_selector, self.tax_axis_filter, watch=True)

        sort_options = {k: [k] for k in self.dfs}
        self.taxonomy_filter = None
        self._taxonomy_paths = {}
        if self.taxonomy_tree_data is not None:

            def index_items(items, parent_path=()):
                for index, item in enumerate(items):
                    path = (*parent_path, index)
                    self._taxonomy_paths[item["id"]] = path
                    index_items(item.get("items", []), path)

            index_items(self.taxonomy_tree_data)
            self.taxonomy_filter = Tree(
                items=self.taxonomy_tree_data,
                checkboxes=True,
                multi_select=True,
                propagate_to_child=True,
                propagate_to_parent=True,
                expanded=[],
            )
            self.taxonomy_filter.active = (
                list(self._taxonomy_paths.values())
                if selected_tax_tree is None
                else [
                    self._taxonomy_paths[node_id]
                    for node_id in selected_tax_tree
                    if node_id in self._taxonomy_paths
                ]
            )
            sort_options["genome"].extend(list(TAXONOMY_RANKS_REGEX.keys()))
            # sort_options = ["genome", *list(TAXONOMY_RANKS_REGEX.keys())]

        if "Meta" in self.dfs["genome"]:
            sort_options["genome"].extend(
                [col for col in self.dfs["genome"]["Meta"].columns if col != "genome"]
            )

        column_options = {k: ["presence/coverage"] for k in self.dfs}
        abundance_cols = ["sample_abundance", "mean_sample_abundance"]
        for k, dfs_dict in self.dfs.items():
            for col in abundance_cols:
                df = dfs_dict[list(dfs_dict.keys())[-1]]
                if col in df.columns:
                    column_options[k].append(col)
                    break

        self.column_options = pn.widgets.NestedSelect(
            label="Y and C columns",
            options=column_options,
            levels=["Y Column", "Color Column"],
        )
        # self.c_col = pn.widgets.MultiChoice(name="Color Column", options=sort_options["genome"])
        self.sort_by = pn.widgets.MultiChoice(
            label="Sort By", options=sort_options[self.y_col]
        )

        self._init_view()
        self.download_heatmap()

    def __panel__(self, *args, **kwargs):
        return self.view

    def _init_view(self):
        additional_sidebar = []
        self.update_plot()

        if self.taxonomy_tree_data is not None:
            additional_sidebar.append(self.show_tax_box)
            additional_sidebar.append("## Taxonomy Filter (genome heatmaps)")
            additional_sidebar.append(self.taxonomy_filter)

        self.view = pn.template.FastListTemplate(
            title="DRAM Product Visualization",
            # main=[self.plot_view],
            main=[
                pn.Tabs(
                    ("Heatmap", self.plot_view),
                    *[
                        (
                            f"{group} df",
                            pn.widgets.Tabulator(df.to_pandas(), page_size=50),
                        )
                        for group, df in self.dfs[self.y_col].items()
                    ],
                )
            ],
            sidebar=[
                pn.Row(self.redraw_button, self.reset_button),
                self.download_button,
                pn.Row(self.column_options),
                # self.sort_by,
                self.param.min_coverage,
                *additional_sidebar,
            ],
        )

    @property
    def y_col(self):
        return self.column_options.value["Y Column"]

    @y_col.setter
    def y_col(self, value):
        self.column_options.value["Y Column"] = value

    @property
    def color_col(self):
        return self.column_options.value["Color Column"]

    @color_col.setter
    def color_col(self, value):
        self.column_options.value["Color Column"] = value

    def update_plot(self, event=None):
        """
        Make the product plot
        """
        None if not self.tax_axis_filter.value else self.tax_axis_rank.value
        charts = []
        y_mode = self.y_col
        c_mode = self.color_col
        for group, df in self.dfs[y_mode].items():
            df = self.filter_by_taxonomy(df)
            if df.is_empty():
                continue
            tooltip_cols = df.columns
            kw = {"y_col": y_mode}
            if (
                not charts
                and y_mode == "genome"
                and self.tax_axis_filter.value
                and "taxonomy" in df.columns
            ):
                kw["y_col"] = "taxonomy"
            if charts:
                kw["y_axis_location"] = None
            if "presence/coverage" in c_mode:
                c_col = get_val_col(df)
                if c_col == "coverage_percentage" and self.min_coverage > 0:
                    df = df.with_columns(
                        pl.when(pl.col(c_col) < self.min_coverage)
                        .then(0)
                        .otherwise(pl.col(c_col))
                        .alias(c_col)
                    )
            elif "abundance" in c_mode:
                c_col = c_mode
                kw["c_max"] = df[c_col].max()
                if kw["c_max"] == 0:
                    kw["c_max"] = 1

            df = self.get_sorted_dfs(df, by=self.sort_by.value)

            hm = heatmap(
                df,
                x_col="name",
                c_col=c_col,
                tooltip_cols=tooltip_cols,
                title=group,
                **kw,
            )
            # insert metadata at the beginning so it shows as the first hm, regardless of the order of the dfs dict
            # if group == "Metadata":
            #     charts.insert(0, hm)
            # else:
            charts.append(hm)

        # module_df, etc_df, function_df = self.get_sorted_dfs(module_df, etc_df, function_df, by=self.sort_by.value)

        self.plot_view[:] = charts or [
            pn.pane.Markdown("No genomes match the selected taxonomy.")
        ]

    def reset_filters(self, event=None):
        """
        Resets the filters applied to the heatmap.

        Parameters:
        - event (optional): The event that triggered the reset. Defaults to None.
        """

        self.min_coverage = self.param.min_coverage.default

        if self.taxonomy_filter is not None:
            self.taxonomy_filter.active = list(self._taxonomy_paths.values())
            self.tax_axis_filter.value = False
            self.tax_axis_rank.visible = False
            self.tax_axis_rank.value = "genus"

        self.sort_by.value = []

    def filter_by_taxonomy(self, df: pl.DataFrame) -> pl.DataFrame:
        """
        Filter the dataframes by selected taxonomy ranks
        """
        if (
            self.taxonomy_filter is None
            or self.y_col != "genome"
            or "taxonomy" not in df.columns
        ):
            return df
        selected = self.taxonomy_filter.value
        leaves = [
            node["id"].replace("; ", ";")
            for node in selected
            if len(node["id"].split(";")) == NO_TAXONOMY_RANKS
        ]
        df = df.filter(pl.col("taxonomy").is_in(leaves))

        return df

    def reveal_tax_axis_rank_selector(
        self, event=None, tax_axis_filter_value: bool = None
    ):
        """
        Set the taxonomy filter
        """
        if isinstance(tax_axis_filter_value, bool) and tax_axis_filter_value:
            self.tax_axis_filter.value = tax_axis_filter_value
        if self.tax_axis_filter.value:
            self.tax_axis_rank.visible = True
            return
        self.tax_axis_rank.visible = False

    def get_sorted_dfs(self, df, by=None):
        """
        Sort the dataframes by taxonomy
        """
        by = by or self.y_col
        return df.sort(by=by, maintain_order=True, nulls_last=True)

    def download_heatmap(self, event=None, output_dir=None):
        """
        Save the heatmap to a file
        """
        output_dir = output_dir or self._output_dir
        plots = [
            Div(text=pane.object) if isinstance(pane, pn.pane.Markdown) else pane.object
            for pane in self.plot_view
        ]
        output_path = (
            output_dir
            / f"heatmap_{self.y_col.replace('/', '-')}_{self.color_col.replace('/', '_')}.html"
        )
        output_path.write_text(
            file_html(
                row(*plots),
                INLINE,
                title="DRAM Product Visualization",
                _always_new=True,
            ),
            encoding="utf-8",
        )
