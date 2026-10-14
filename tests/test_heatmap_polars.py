import unittest
from unittest.mock import patch

import polars as pl
from bokeh.models import CategoricalColorMapper, LinearColorMapper, LogColorMapper

from dram_viz.apps.heatmap import heatmap


class PolarsHeatmapTests(unittest.TestCase):
    def plot(self, values, **kwargs):
        frame = pl.DataFrame(
            {"name": ["z", "a", "z"], "genome": ["g2", "g1", "g2"], "value": values}
        )
        original = frame.clone()
        with patch.object(
            pl.DataFrame, "to_pandas", side_effect=AssertionError("pandas conversion")
        ):
            plot = heatmap(
                frame,
                y_col="genome",
                x_col="name",
                c_col="value",
                tooltip_cols=frame.columns,
                **kwargs,
            )
        self.assertTrue(frame.equals(original))
        self.assertEqual(plot.x_range.factors, ["z", "a"])
        self.assertEqual(plot.y_range.factors, ["g2", "g1"])
        renderer = plot.renderers[0]
        return renderer.data_source.data, renderer.glyph.fill_color.transform

    def test_numeric_and_log_colors(self):
        data, mapper = self.plot([0.1, 0.2, 0.3])
        self.assertIsInstance(mapper, LinearColorMapper)
        self.assertEqual(list(data["value"]), [0.1, 0.2, 0.3])
        _, mapper = self.plot([1, 2, 3], log_colors=True, c_min=1, c_max=3)
        self.assertIsInstance(mapper, LogColorMapper)

    def test_boolean_labels_and_missing_values(self):
        data, mapper = self.plot([True, False, None])
        self.assertIsInstance(mapper, CategoricalColorMapper)
        self.assertEqual(list(data["value"]), ["true", "false", None])
        self.assertEqual(mapper.factors, ["false", "true"])
        _, mapper = self.plot([True, True, True])
        self.assertEqual(mapper.factors, ["false", "true"])

    def test_categorical_labels_and_missing_values(self):
        data, mapper = self.plot(["beta", None, "alpha"])
        self.assertEqual(list(data["value"]), ["beta", None, "alpha"])
        self.assertEqual(mapper.factors, ["alpha", "beta"])

    def test_nulls_are_distinct_from_literal_none_strings(self):
        data, mapper = self.plot(["None", None, "alpha"])
        self.assertEqual(list(data["value"]), ["None", None, "alpha"])
        self.assertEqual(mapper.factors, ["None", "alpha"])

    def test_all_null_colors(self):
        data, mapper = self.plot([None, None, None])
        self.assertEqual(list(data["value"]), [None, None, None])
        self.assertEqual(mapper.factors, [])

    def test_unpivot_preserves_order_and_deduplicates(self):
        frame = pl.DataFrame(
            {"genome": ["g2", "g1", "g2"], "B": [2, 1, 2], "A": [4, 3, 4]}
        )
        with patch.object(
            pl.DataFrame, "to_pandas", side_effect=AssertionError("pandas conversion")
        ):
            plot = heatmap(frame, y_col="genome", x_cols=["B", "A"], tooltip_cols=[])
        self.assertEqual(plot.x_range.factors, ["B", "A"])
        self.assertEqual(plot.y_range.factors, ["g2", "g1"])
        data = plot.renderers[0].data_source.data
        self.assertEqual(list(data["value"]), [2, 1, 4, 3])
        self.assertEqual(list(data["genome"]), ["g2", "g1", "g2", "g1"])


if __name__ == "__main__":
    unittest.main()
