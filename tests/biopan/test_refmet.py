import sys
import os
from types import SimpleNamespace
import logging
import unittest
from lipidmaps.biopan.models.refmet import RefMet, RefMetResult
from lipidmaps.biopan.models.sample import QuantifiedLipid


sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
current_dir = os.path.dirname(os.path.abspath(__file__))
biopan_test_file = os.path.join(current_dir, "inputs", "biopan_small_demo.csv")


logging.basicConfig(
    level=logging.DEBUG, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)


class TestRefMetValidation(unittest.TestCase):
    """Test RefMet validation and annotation methods."""

    def test_validate_metabolite_names(self):
        """Test RefMet.validate_metabolite_names() returns correct RefMetResult objects."""
        metabolite_names = ["PC(16:0/18:1)", "Cholesterol", "TAG(16:0/18:1/18:2)"]
        
        results = RefMet.validate_metabolite_names(metabolite_names)
        
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
            SimpleNamespace(sample_name="Cholesterol")
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
            RefMetResult(input_name="Cholesterol", standardized_name="Cholesterol")
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
                refmet_result={"lm_id": "LMGP01050001", "standardized_name": "PC 16:0/18:1"}
            ),
            SimpleNamespace(
                sample_name="Sample2",
                refmet_result={"lm_id": None, "standardized_name": "Unknown"}
            ),
            SimpleNamespace(
                sample_name="Sample3",
                refmet_result={"lm_id": "LMST01010001", "standardized_name": "Cholesterol"}
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
                refmet_result={"lm_id": "LMGP01050001", "standardized_name": "PC 16:0/18:1"}
            ),
            SimpleNamespace(
                sample_name="Sample2",
                refmet_result={"lm_id": None, "standardized_name": "Unknown Lipid"}
            ),
            SimpleNamespace(
                sample_name="Sample3",
                refmet_result={"lm_id": None, "standardized_name": "Mystery Compound"}
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
