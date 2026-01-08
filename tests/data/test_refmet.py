import sys
import os
from types import SimpleNamespace
import logging
import unittest
from lipidmaps.data.models.refmet import RefMet, RefMetResult
from lipidmaps.data.models.sample import QuantifiedLipid
import csv


sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
current_dir = os.path.dirname(os.path.abspath(__file__))
test_file = os.path.join(current_dir, "inputs", "large_demo.csv")
# We should also display the result when we get the annotation of this file

logging.basicConfig(
    level=logging.DEBUG, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)


class TestRefMetValidation(unittest.TestCase):
    def validate_file_structure(
        self, csv_path, required_columns=None, quant_columns=None
    ):
        """Validate file structure: required columns, field count, numeric quant columns."""
        errors = []
        required_columns = required_columns or ["sample_name", "lm_id"]
        quant_columns = quant_columns or ["quant1", "quant2"]
        with open(csv_path, newline="", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            header = reader.fieldnames
            # Check required columns
            for col in required_columns + quant_columns:
                if col not in header:
                    errors.append(f"Missing required column: {col}")
            for i, row in enumerate(reader, 1):
                # Check field count
                if len(row) != len(header):
                    errors.append(f"Row {i} field count mismatch: {row}")
                # Check quant columns are numeric
                for qc in quant_columns:
                    val = row.get(qc)
                    if val is not None:
                        try:
                            float(val)
                        except (ValueError, TypeError):
                            errors.append(
                                f"Row {i} non-numeric quant column {qc}: {val}"
                            )
        return errors

    def test_file_structure_positive(self):
        """Test that a well-formed file passes structure validation."""
        csv_path = os.path.join(current_dir, "inputs", "file_structure_positive.csv")
        errors = self.validate_file_structure(csv_path)
        if errors:
            print("\nFile structure errors (positive):", errors)
        self.assertEqual(
            len(errors), 0, f"Positive file structure test failed: {errors}"
        )

    def test_file_structure_negative(self):
        """Test that a malformed file fails structure validation and reports errors."""
        csv_path = os.path.join(current_dir, "inputs", "file_structure_negative.csv")
        errors = self.validate_file_structure(csv_path)
        print("\nFile structure errors (negative):", errors)
        self.assertGreater(
            len(errors),
            0,
            "Negative file structure test did not report errors as expected.",
        )

    """Test RefMet validation and annotation methods."""

    def test_validate_metabolite_names(self):
        """Test RefMet.validate_metabolite_names() returns correct RefMetResult objects."""
        metabolite_names = ["PC(16:0/18:1)", "Cholesterol", "TAG(16:0/18:1/18:2)"]

        results = RefMet.validate_metabolite_names(metabolite_names)
        print(f"\n{results}")
        # Verify return type and length
        self.assertIsInstance(results, list)
        self.assertEqual(len(results), 3)

        # Verify all items are RefMetResult objects
        for result in results:
            self.assertIsInstance(result, RefMetResult)
            self.assertIsNotNone(result.input_name)
            self.assertIsNotNone(result.standardized_name)

        # Check specific result for PC
        pc_result = results[0]
        self.assertEqual(pc_result.input_name, "PC(16:0/18:1)")
        self.assertEqual(pc_result.standardized_name, "PC 16:0/18:1")
        self.assertEqual(pc_result.sub_class, "PC")
        self.assertIsNotNone(pc_result.formula)

        logger.info(f"PC result: {pc_result.model_dump()}")

    def test_attach_results_to_samples(self):
        """Test RefMet.attach_results_to_samples() properly attaches results."""
        # Create sample objects
        samples = [
            SimpleNamespace(sample_name="PC(16:0/18:1)"),
            SimpleNamespace(sample_name="Cholesterol"),
        ]

        # Get validation results
        names = [s.sample_name for s in samples]
        results = RefMet.validate_metabolite_names(names)

        # Attach results
        RefMet.attach_results_to_samples(samples, results)

        # Verify attachment
        for sample in samples:
            self.assertTrue(hasattr(sample, "refmet_result"))
            self.assertIsInstance(sample.refmet_result, dict)
            self.assertIn("input_name", sample.refmet_result)
            self.assertIn("standardized_name", sample.refmet_result)
            self.assertEqual(sample.refmet_result["input_name"], sample.sample_name)

    def test_attach_results_length_mismatch(self):
        """Test that attach_results_to_samples raises error on length mismatch."""
        samples = [SimpleNamespace(sample_name="PC(16:0/18:1)")]
        results = [
            RefMetResult(input_name="PC(16:0/18:1)", standardized_name="PC 16:0/18:1"),
            RefMetResult(input_name="Cholesterol", standardized_name="Cholesterol"),
        ]

        with self.assertRaises(ValueError) as context:
            RefMet.attach_results_to_samples(samples, results)

        self.assertIn("Length mismatch", str(context.exception))

    def test_annotate_samples(self):
        """Test RefMet.annotate_samples() convenience method."""
        # Create test samples
        test_samples = [
            QuantifiedLipid(input_name="PC(16:0/18:1)", values={"sample1": 1234.5})
        ]
        temp_objs = [SimpleNamespace(sample_name=q.input_name) for q in test_samples]

        logger.info("Calling RefMet.annotate_samples()...")
        results = RefMet.annotate_samples(temp_objs)

        # Verify return value is a list of RefMetResult objects
        self.assertIsInstance(results, list)
        self.assertEqual(len(results), 1)
        self.assertIsInstance(results[0], RefMetResult)

        # Check the first result
        result = results[0]
        self.assertEqual(result.input_name, "PC(16:0/18:1)")
        self.assertEqual(result.standardized_name, "PC 16:0/18:1")
        self.assertEqual(result.sub_class, "PC")

        # Verify samples have refmet_result attached
        for obj in temp_objs:
            self.assertTrue(hasattr(obj, "refmet_result"))
            self.assertIsInstance(obj.refmet_result, dict)
            self.assertEqual(obj.refmet_result["input_name"], "PC(16:0/18:1)")

    def test_get_lm_ids(self):
        """Test RefMet.get_lm_ids() extracts LM IDs from samples."""
        # Create samples with refmet results
        samples = [
            SimpleNamespace(
                sample_name="Sample1",
                refmet_result={
                    "lm_id": "LMGP01050001",
                    "standardized_name": "PC 16:0/18:1",
                },
            ),
            SimpleNamespace(
                sample_name="Sample2",
                refmet_result={"lm_id": None, "standardized_name": "Unknown"},
            ),
            SimpleNamespace(
                sample_name="Sample3",
                refmet_result={
                    "lm_id": "LMST01010001",
                    "standardized_name": "Cholesterol",
                },
            ),
        ]

        lm_ids = RefMet.get_lm_ids(samples)

        # Verify we got the LM IDs
        self.assertIsInstance(lm_ids, list)
        self.assertIn("LMGP01050001", lm_ids)
        self.assertIn("LMST01010001", lm_ids)
        # Should be unique
        self.assertEqual(len(set(lm_ids)), len(lm_ids))

    def test_get_unmatched_results(self):
        """Test RefMet.get_unmatched_results() returns standardized names without LM IDs."""
        samples = [
            SimpleNamespace(
                sample_name="Sample1",
                refmet_result={
                    "lm_id": "LMGP01050001",
                    "standardized_name": "PC 16:0/18:1",
                },
            ),
            SimpleNamespace(
                sample_name="Sample2",
                refmet_result={"lm_id": None, "standardized_name": "Unknown Lipid"},
            ),
            SimpleNamespace(
                sample_name="Sample3",
                refmet_result={"lm_id": None, "standardized_name": "Mystery Compound"},
            ),
        ]

        unmatched = RefMet.get_unmatched_results(samples)

        # Verify unmatched names
        self.assertIsInstance(unmatched, list)
        self.assertIn("Unknown Lipid", unmatched)
        self.assertIn("Mystery Compound", unmatched)
        self.assertNotIn("PC 16:0/18:1", unmatched)
        # Should be unique
        self.assertEqual(len(set(unmatched)), len(unmatched))

    def test_refmet_positive_csv_matches_lm_id(self):
        """Test that sample_name in refmet_positive.csv returns the expected lm_id from RefMet."""
        positive_csv = os.path.join(current_dir, "inputs", "refmet_positive.csv")
        mismatches = []
        with open(positive_csv, newline="", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            rows = [
                row for row in reader if row.get("sample_name") and row.get("lm_id")
            ]
        sample_names = [row["sample_name"] for row in rows]
        expected_lm_ids = [row["lm_id"] for row in rows]
        results = RefMet.validate_metabolite_names(sample_names)
        for name, expected_lm_id, result in zip(
            sample_names, expected_lm_ids, results
        ):
            found_lm_id = getattr(result, "lm_id", None)
            if found_lm_id != expected_lm_id:
                mismatches.append((name, expected_lm_id, found_lm_id))
        if mismatches:
            print("\nRefMet LMID mismatches:")
            for name, expected, found in mismatches:
                print(
                    f"Sample: {name}, Expected LMID: {expected}, Found LMID: {found}"
                )
        self.assertEqual(
            len(mismatches),
            0,
            "Some sample_names did not return expected LMIDs. See output above.",
        )


def annotate_samples_and_export_csv():
    """Extract sample names from large_demo.csv, annotate with RefMet, and output results to CSV."""
    input_csv = os.path.join(current_dir, "inputs", "large_demo.csv")
    output_csv = os.path.join(current_dir, "large_demo_refmet_results.csv")
    # Read sample names
    with open(input_csv, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        sample_names = [row["sample"] for row in reader if row.get("sample")]
    # Annotate using RefMet
    results = RefMet.validate_metabolite_names(sample_names)

    # Write results to CSV
    with open(output_csv, "w", newline="", encoding="utf-8") as f:
        fieldnames = [
            "sample",
            "lm_id",
            "standardized_name",
            "refmet_id",
            "formula",
            "exact_mass",
        ]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for name, result in zip(sample_names, results):
            writer.writerow(
                {
                    "sample": name,
                    "lm_id": getattr(result, "lm_id", ""),
                    "standardized_name": getattr(result, "standardized_name", ""),
                    "refmet_id": getattr(result, "refmet_id", ""),
                    "formula": getattr(result, "formula", ""),
                    "exact_mass": getattr(result, "exact_mass", ""),
                }
            )
    print(f"Results written to {output_csv}")


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "export":
        annotate_samples_and_export_csv()
    else:
        unittest.main()
