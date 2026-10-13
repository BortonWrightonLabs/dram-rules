from __future__ import annotations

from math import pi
from pathlib import Path
from typing import Optional

import pandas as pd
import panel as pn
import param
from bokeh.palettes import BuGn, Cividis256
from bokeh.plotting import figure
from bokeh.resources import INLINE
from bokeh.transform import factor_cmap, linear_cmap, log_cmap
from pandas.api.types import is_numeric_dtype
from panel_jstree import Tree

from dram_viz.definitions import NO_TAXONOMY_RANKS, TAXONOMY_RANKS_REGEX
from dram_viz.rule_parser.src.rules import get_val_col

pn.extension("tabulator", "katex", template="bootstrap")

PALETTE_CATEGORICAL = BuGn
PALETTE_CONTINUOUS = Cividis256

HEATMAP_CELL_HEIGHT = 15
HEATMAP_CELL_WIDTH = 15


def heatmap(
    df,
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
    df : pd.DataFrame
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
        df = pd.melt(
            df, id_vars=y_col, value_vars=x_cols, var_name="x_col"
        ).drop_duplicates()
        x_col = "x_col"
        c_col = "value"
        tooltip_cols = [y_col, "value"]

    rect_kw = rect_kw or {}
    # df = df.sort_values(by=[y_col], ascending=False)

    tooltips = []
    for col in tooltip_cols:
        if isinstance(col, tuple):
            tooltips.append(col)
        else:
            tooltips.append((col.replace("_", " ").title(), f"@{col}"))

    p = figure(
        frame_width=HEATMAP_CELL_WIDTH * len(df[x_col].unique()),
        frame_height=HEATMAP_CELL_WIDTH * len(df[y_col].unique()),
        x_range=list(df[x_col].unique()),
        y_range=list(df[y_col].unique()),
        tools="hover",
        toolbar_location=None,
        tooltips=tooltips,
        # title="\n".join(char for char in title),
        title=title,
        title_location="right",
        **fig_kwargs,
    )

    # if x_col:
    # if df[c_col].dtype == float:
    if is_numeric_dtype(df[c_col]) and df[c_col].dtype != bool:
        palette = tuple(reversed(PALETTE_CONTINUOUS))
        if log_colors:
            fill_color = log_cmap(c_col, palette=palette, low=c_min, high=c_max)
        else:
            fill_color = linear_cmap(c_col, palette=palette, low=c_min, high=c_max)
    else:
        # if the column is boolean, we want to treat it as categorical
        # but we need to convert it to string first so that the factor_cmap works correctly
        # And we need to make sure both false and true are included as factors
        # even if they aren't both present in the data
        if df[c_col].dtype == bool:
            df[c_col] = df[c_col].astype(str)
            factors = [str(False), str(True)]
        else:
            df[c_col] = df[c_col].astype(str)
            factors = sorted(df[c_col].unique())
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
    module_df : pd.DataFrame
        The module dataframe.
    etc_df : pd.DataFrame
        The etc dataframe.
    function_df : pd.DataFrame
        The function dataframe.
    tax_tree_data : optional
        The taxonomy tree data.
    """

    min_coverage = param.Number(default=0, bounds=(0, 1), label="Minimum Coverage")
    view = param.ClassSelector(class_=pn.template.FastListTemplate)
    plot_view = param.ClassSelector(class_=pn.Row)

    def __init__(
        self,
        dfs: dict[str, dict[str, pd.DataFrame]],
        taxanomy_tree_data: Optional[pd.DataFrame] = None,
        selected_tax_tree: Optional[list[str]] = None,
        output_dir: str | Path = None,
    ):
        super().__init__()
        self.dfs = dfs
        self.taxonomy_tree_data = taxanomy_tree_data
        self._output_dir = Path(output_dir) if output_dir else Path.cwd()
        self.plot_view = pn.Row()
        self.download_button = pn.widgets.Button(
            name="Download Heatmap", button_type="primary"
        )
        self.download_button.on_click(self.download_heatmap)

        self.redraw_button = pn.widgets.Button(name="Redraw", button_type="primary")
        self.reset_button = pn.widgets.Button(
            name="Reset Filters", button_type="warning"
        )
        self.reset_button.on_click(self.reset_filters)

        self.redraw_button.on_click(self.update_plot)

        self.tax_axis_filter = pn.widgets.Checkbox(
            name="Show Taxonomy on Y Axis", value=False
        )
        self.tax_axis_rank = pn.widgets.Select(
            name="Taxonomy Label",
            options=list(TAXONOMY_RANKS_REGEX),
            visible=False,
            value="genus",
        )
        self.show_tax_box = pn.Column(self.tax_axis_filter, self.tax_axis_rank)
        pn.bind(self.reveal_tax_axis_rank_selector, self.tax_axis_filter, watch=True)

        sort_options = {k: [k] for k in self.dfs}
        self.taxonomy_filter = None
        if self.taxonomy_tree_data is not None:
            self.taxonomy_filter = Tree(
                data=self.taxonomy_tree_data, show_icons=False, cascade=True
            )

            # hack to make sure the taxonomy filter tree.value is set since it isn't set on the first load
            # TODO: remove maybe when this is put into panel
            self.taxonomy_filter.value = selected_tax_tree or []
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
            name="Y and C columns",
            options=column_options,
            levels=["Y Column", "Color Column"],
        )
        # self.c_col = pn.widgets.MultiChoice(name="Color Column", options=sort_options["genome"])
        self.sort_by = pn.widgets.MultiChoice(
            name="Sort By", options=sort_options[self.y_col]
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
            additional_sidebar.append("## Taxonomy Filter")
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
        for i, (group, df) in enumerate(self.dfs[y_mode].items()):
            df = df.to_pandas()
            tooltip_cols = df.columns.tolist()
            kw = {"y_col": y_mode}
            if i == 0 and self.tax_axis_filter.value:
                kw["y_col"] = "taxonomy"
            if i != 0:
                kw["y_axis_location"] = None
            if "presence/coverage" in c_mode:
                c_col = get_val_col(df)
                if c_col == "coverage_percentage" and self.min_coverage > 0:
                    df.loc[df[c_col] < self.min_coverage, c_col] = 0
            elif "abundance" in c_mode:
                c_col = c_mode
                kw["c_max"] = df[c_col].max()
                if kw["c_max"] == 0:
                    kw["c_max"] = 1

            df = self.filter_by_taxonomy(df)
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

        self.plot_view[:] = charts

    def reset_filters(self, event=None):
        """
        Resets the filters applied to the heatmap.

        Parameters:
        - event (optional): The event that triggered the reset. Defaults to None.
        """

        self.min_coverage = self.param.min_coverage.default

        if self.taxonomy_filter is not None:
            self.taxonomy_filter.value = [
                node["id"] for node in self.taxonomy_filter.flat_tree
            ]
            self.tax_axis_filter.value = False
            self.tax_axis_rank.visible = False
            self.tax_axis_rank.value = "genus"

        self.sort_by.value = []

    def filter_by_taxonomy(self, df):
        """
        Filter the dataframes by selected taxonomy ranks
        """
        if self.taxonomy_filter is None:
            return df
        selected = self.taxonomy_filter.value
        # leaves = [node for node in selected if len(node.split(";")) == NO_TAXONOMY_RANKS]
        # maybe we don't need this replace, but leaving in for now to be sure we match the data
        leaves = [
            node.replace("; ", ";")
            for node in selected
            if len(node.split(";")) == NO_TAXONOMY_RANKS
        ]
        df = df.loc[df["taxonomy"].isin(leaves)]

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
        return df.sort_values(by=by, kind="stable")

    def download_heatmap(self, event=None, output_dir=None):
        """
        Save the heatmap to a file
        """
        output_dir = output_dir or self._output_dir
        self.plot_view.save(
            output_dir
            / f"heatmap_{self.y_col.replace('/', '-')}_{self.color_col.replace('/', '_')}.html",
            resources=INLINE,
        )
