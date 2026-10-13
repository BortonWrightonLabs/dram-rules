import tempfile
import unittest
from pathlib import Path

import panel as pn
import polars as pl
from bokeh.document import Document
from panel.io.state import set_curdoc
from panel_material_ui import Tree

from dram_viz.apps.heatmap import Dashboard
from dram_viz.processing.process_annotations import (
    build_tax_edge_df_pl,
    build_tax_tree_selected_recurse,
    build_tree_pl,
)


class TaxonomyFilterTests(unittest.TestCase):
    def setUp(self):
        self.output = tempfile.TemporaryDirectory()
        self.addCleanup(self.output.cleanup)
        self.output_dir = Path(self.output.name)
        self.taxonomies = [
            "d__Bacteria;p__A;c__A;o__A;f__A;g__A;s__one",
            "d__Bacteria;p__A;c__A;o__A;f__A;g__A;s__two",
            "d__Bacteria;p__B;c__B;o__B;f__B;g__B;s__three",
        ]
        self.frame = pl.DataFrame(
            {
                "genome": ["g1", "g2", "g3", "g4"],
                "taxonomy": [
                    self.taxonomies[0],
                    self.taxonomies[0],
                    self.taxonomies[1],
                    self.taxonomies[2],
                ],
                "name": ["Trait"] * 4,
                "present": [True, False, True, True],
                "mean_sample_abundance": [1.0, 2.0, 3.0, 4.0],
            }
        )
        edges, _ = build_tax_edge_df_pl(self.frame.select("genome", "taxonomy"))
        self.tree = build_tree_pl(
            edges, id_cb=lambda source, child, parent_id: f"{parent_id};{child}"
        )
        self.dfs = {
            "genome": {
                "All": self.frame,
                "Second": self.frame.filter(pl.col("genome") == "g4"),
            },
            "sample": {
                "All": pl.DataFrame(
                    {
                        "sample": ["sample1"],
                        "name": ["Trait"],
                        "present": [True],
                        "sample_abundance": [10.0],
                    }
                )
            },
        }

    def dashboard(self, **kwargs):
        return Dashboard(
            dfs=self.dfs,
            taxanomy_tree_data=self.tree,
            output_dir=self.output_dir,
            **kwargs,
        )

    def plotted_genomes(self, dashboard):
        return [
            list(pane.object.renderers[0].data_source.data["genome"])
            for pane in dashboard.plot_view
        ]

    def test_node_schema_and_lineage_ids(self):
        def check(nodes):
            for node in nodes:
                self.assertEqual(set(node), {"id", "label", "items"})
                self.assertEqual(node["label"], node["id"].split(";")[-1])
                check(node["items"])

        check(self.tree)
        ids = build_tax_tree_selected_recurse(self.tree)
        self.assertEqual(len(ids), len(set(ids)))
        self.assertTrue(set(self.taxonomies).issubset(ids))

    def test_initial_selection(self):
        dashboard = self.dashboard()
        self.assertIsInstance(dashboard.taxonomy_filter, Tree)
        self.assertEqual(
            {item["id"] for item in dashboard.taxonomy_filter.value},
            set(build_tax_tree_selected_recurse(self.tree)),
        )
        self.assertEqual(dashboard.taxonomy_filter.expanded, [])
        self.assertTrue(dashboard.taxonomy_filter.checkboxes)
        self.assertTrue(dashboard.taxonomy_filter.propagate_to_child)
        self.assertTrue(dashboard.taxonomy_filter.propagate_to_parent)
        self.assertEqual(
            self.plotted_genomes(dashboard), [["g1", "g2", "g3", "g4"], ["g4"]]
        )

    def test_explicit_empty_and_partial_initial_selection(self):
        empty = self.dashboard(selected_tax_tree=[])
        self.assertEqual(empty.taxonomy_filter.value, [])
        self.assertIsInstance(empty.plot_view[0], pn.pane.Markdown)
        self.assertEqual(
            empty.plot_view[0].object, "No genomes match the selected taxonomy."
        )
        partial = self.dashboard(selected_tax_tree=[self.taxonomies[0]])
        self.assertEqual(self.plotted_genomes(partial), [["g1", "g2"]])

    def test_unchecked_leaf_is_not_restored_by_selected_parent(self):
        dashboard = self.dashboard()
        dashboard.taxonomy_filter.active = [
            path
            for node_id, path in dashboard._taxonomy_paths.items()
            if node_id != self.taxonomies[0]
        ]
        self.assertEqual(self.plotted_genomes(dashboard)[0], ["g1", "g2", "g3", "g4"])
        dashboard.update_plot()
        self.assertEqual(self.plotted_genomes(dashboard), [["g3", "g4"], ["g4"]])
        self.assertTrue(self.dfs["genome"]["All"].equals(self.frame))

    def test_branch_exclusion_skips_empty_groups(self):
        dashboard = self.dashboard()
        branch = "d__Bacteria;p__B"
        dashboard.taxonomy_filter.active = [
            path
            for node_id, path in dashboard._taxonomy_paths.items()
            if node_id != branch and not node_id.startswith(branch + ";")
        ]
        dashboard.update_plot()
        self.assertEqual(self.plotted_genomes(dashboard), [["g1", "g2", "g3"]])
        dashboard.taxonomy_filter.active = [
            dashboard._taxonomy_paths[self.taxonomies[2]]
        ]
        dashboard.update_plot()
        self.assertEqual(self.plotted_genomes(dashboard), [["g4"], ["g4"]])

    def test_reset_requires_redraw_and_restores_defaults(self):
        dashboard = self.dashboard()
        dashboard.taxonomy_filter.active = []
        dashboard.update_plot()
        dashboard.min_coverage = 0.5
        dashboard.tax_axis_filter.value = True
        dashboard.tax_axis_rank.value = "species"
        dashboard.sort_by.value = ["genome"]
        dashboard.reset_filters()
        self.assertEqual(dashboard.min_coverage, 0)
        self.assertFalse(dashboard.tax_axis_filter.value)
        self.assertFalse(dashboard.tax_axis_rank.visible)
        self.assertEqual(dashboard.tax_axis_rank.value, "genus")
        self.assertEqual(dashboard.sort_by.value, [])
        self.assertEqual(
            len(dashboard.taxonomy_filter.active), len(dashboard._taxonomy_paths)
        )
        self.assertIsInstance(dashboard.plot_view[0], pn.pane.Markdown)
        dashboard.update_plot()
        self.assertEqual(self.plotted_genomes(dashboard)[0], ["g1", "g2", "g3", "g4"])

    def test_first_nonempty_group_keeps_y_axis(self):
        self.dfs["genome"]["All"] = self.frame.filter(pl.col("genome") == "g1")
        dashboard = self.dashboard(selected_tax_tree=[self.taxonomies[2]])
        self.assertEqual(self.plotted_genomes(dashboard), [["g4"]])
        self.assertTrue(dashboard.plot_view[0].object.yaxis)

    def test_lineage_whitespace_normalization(self):
        edges, _ = build_tax_edge_df_pl(self.frame.select("genome", "taxonomy"))
        self.tree = build_tree_pl(
            edges, id_cb=lambda source, child, parent_id: f"{parent_id}; {child}"
        )
        dashboard = self.dashboard()
        self.assertEqual(self.plotted_genomes(dashboard)[0], ["g1", "g2", "g3", "g4"])

    def test_sample_mode_ignores_taxonomy_selection(self):
        dashboard = self.dashboard(selected_tax_tree=[])
        dashboard.column_options.value = {
            "Y Column": "sample",
            "Color Column": "sample_abundance",
        }
        dashboard.tax_axis_filter.value = True
        dashboard.update_plot()
        source = dashboard.plot_view[0].object.renderers[0].data_source.data
        self.assertEqual(list(source["sample"]), ["sample1"])
        self.assertEqual(list(source["sample_abundance"]), [10.0])
        with_taxonomy = self.frame.to_pandas()
        self.assertIs(dashboard.filter_by_taxonomy(with_taxonomy), with_taxonomy)

    def test_without_taxonomy(self):
        frame = self.frame.drop("taxonomy")
        dashboard = Dashboard(
            dfs={"genome": {"All": frame}}, output_dir=self.output_dir
        )
        self.assertIsNone(dashboard.taxonomy_filter)
        self.assertEqual(self.plotted_genomes(dashboard)[0], ["g1", "g2", "g3", "g4"])
        with_tree = self.dashboard(selected_tax_tree=[])
        self.assertIs(with_tree.filter_by_taxonomy(frame), frame)

    def test_html_export_and_restoring_selection(self):
        dashboard = self.dashboard(selected_tax_tree=[self.taxonomies[0]])
        exported = self.output_dir / "heatmap_genome_presence_coverage.html"
        html = exported.read_text()
        self.assertIn('"g1"', html)
        self.assertNotIn('"g3"', html)
        dashboard.taxonomy_filter.active = []
        dashboard.update_plot()
        dashboard.download_heatmap()
        self.assertIn("No genomes match the selected taxonomy.", exported.read_text())
        dashboard.taxonomy_filter.active = [
            dashboard._taxonomy_paths[self.taxonomies[2]]
        ]
        dashboard.redraw_button.clicks += 1
        dashboard.download_heatmap()
        self.assertEqual(self.plotted_genomes(dashboard), [["g4"], ["g4"]])
        self.assertIn('"g4"', exported.read_text())

    def test_html_export_preserves_live_document(self):
        document = Document()
        with set_curdoc(document):
            dashboard = self.dashboard()
            dashboard.plot_view.server_doc(document, location=False)
            dashboard.download_heatmap()
            dashboard.taxonomy_filter.active = []
            dashboard.update_plot()
            dashboard.download_heatmap()
            dashboard.reset_filters()
            dashboard.update_plot()
            dashboard.download_heatmap()
            self.assertEqual(
                self.plotted_genomes(dashboard)[0], ["g1", "g2", "g3", "g4"]
            )
            self.assertTrue(document.roots)
            self.assertTrue(
                all(model.document is document for model in document.models)
            )


if __name__ == "__main__":
    unittest.main()
