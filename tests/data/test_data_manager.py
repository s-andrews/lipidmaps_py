import unittest
from pathlib import Path

from lipidmaps.data.data_manager import DataManager


class TestDataManager(unittest.TestCase):
    def setUp(self):
        self.data_manager = DataManager()

    def test_process_csv_populates_dataset(self):
        csv_path = Path(__file__).parent / "inputs" / "small_demo.csv"
        dataset = self.data_manager.process_csv(csv_path)
        # dataset should be populated
        self.assertIsNotNone(dataset)
        # samples should match header columns minus the NAME column
        with csv_path.open() as fh:
            header = fh.readline().strip().split(",")
        expected_samples = len(header) - 1
        self.assertEqual(len(dataset.samples), expected_samples)
        # there should be at least one lipid row
        self.assertGreater(len(dataset.lipids), 0)


if __name__ == "__main__":
    unittest.main()
