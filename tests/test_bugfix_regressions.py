"""Regression tests for bugs found during the codebase review.

Each test pins the behaviour of a specific fix so the crash cannot silently
return. See the review plan for the corresponding bug descriptions.
"""

import unittest

from lipidmaps.data.data_manager import DataManager
from lipidmaps.data.ingestion.csv_reader import CSVIngestion
from lipidmaps.data.models.sample import (
    LipidDataset,
    QuantifiedLipid,
    SampleMetadata,
)


class TestStructurePropertyNoneSafety(unittest.TestCase):
    """`QuantifiedLipid.structure` must not crash when standardized_name is None."""

    def test_unrecognized_non_fa_lipid_returns_none_without_error(self):
        lipid = QuantifiedLipid(input_name="TotallyUnknownThing", standardized_name=None)
        # Must not raise AttributeError on None.lower()
        self.assertIsNone(lipid.structure)

    def test_recognized_lipid_still_parses(self):
        lipid = QuantifiedLipid(input_name="PC 34:1", standardized_name="PC 34:1")
        self.assertIsNotNone(lipid.structure)
        self.assertEqual(lipid.structure.headgroup, "PC")


class TestZscoreMissingValues(unittest.TestCase):
    """`QuantifiedLipid.zscore()` must tolerate None/missing values."""

    def test_zscore_with_none_value(self):
        lipid = QuantifiedLipid(
            input_name="PC 34:1",
            values={"s1": 1.0, "s2": None, "s3": 3.0},
        )
        scores = lipid.zscore()  # must not raise
        self.assertEqual(scores["s2"], 0.0)
        # Present values [1.0, 3.0]: mean 2.0, sample std (ddof=1) sqrt(2)
        self.assertAlmostEqual(scores["s1"], -1.0 / (2 ** 0.5))
        self.assertAlmostEqual(scores["s3"], 1.0 / (2 ** 0.5))

    def test_zscore_all_missing(self):
        lipid = QuantifiedLipid(input_name="PC 34:1", values={"s1": None, "s2": None})
        self.assertEqual(lipid.zscore(), {"s1": 0.0, "s2": 0.0})


class TestGroupStatisticsMissingValues(unittest.TestCase):
    """`DataManager.get_group_statistics` must skip None values, not crash."""

    def test_group_statistics_skips_none(self):
        samples = [
            SampleMetadata(sample_name="c1", group="control"),
            SampleMetadata(sample_name="c2", group="control"),
            SampleMetadata(sample_name="t1", group="treated"),
        ]
        lipids = [
            QuantifiedLipid(
                input_name="PC 34:1",
                values={"c1": 2.0, "c2": None, "t1": 5.0},
            ),
        ]
        dataset = LipidDataset(samples=samples, lipids=lipids)
        manager = DataManager(dataset=dataset)

        stats = manager.get_group_statistics()  # must not raise
        self.assertAlmostEqual(stats["control"]["mean_values"]["PC 34:1"], 2.0)
        # only one non-None value in control -> std defined as 0.0
        self.assertEqual(stats["control"]["std_values"]["PC 34:1"], 0.0)
        self.assertAlmostEqual(stats["treated"]["mean_values"]["PC 34:1"], 5.0)


class TestHasLabelsHeaderOnly(unittest.TestCase):
    """`has_labels` must not raise IndexError on a header-only CSV."""

    def test_header_only_csv_with_labels(self, tmp_path=None):
        import tempfile
        from pathlib import Path

        with tempfile.TemporaryDirectory() as d:
            csv_path = Path(d) / "header_only.csv"
            csv_path.write_text("lipid,s1,s2\n")  # header, no data rows

            ingestion = CSVIngestion(has_labels=True)
            raw = ingestion.read_standard_csv(csv_path)  # must not raise
            self.assertEqual(len(raw.rows), 0)


if __name__ == "__main__":
    unittest.main()
