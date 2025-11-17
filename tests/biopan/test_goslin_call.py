import unittest
from lipidmaps.biopan.models.goslin import Goslin, GoslinResult


class Sample:
    def __init__(self, name):
        self.sample_name = name


class TestGoslinAnnotation(unittest.TestCase):
    def test_annotate_samples(self):
        """Test Goslin annotation of lipid samples."""
        samps = [Sample("Cholesterol"), Sample("PC(16:0/18:1)")]
        results = Goslin.annotate_samples(samps)
        
        # Verify return value is a list of GoslinResult objects
        self.assertIsInstance(results, list)
        self.assertEqual(len(results), 2)
        
        for result in results:
            self.assertIsInstance(result, GoslinResult)
            self.assertIsNotNone(result.input_name)
            self.assertIsNotNone(result.normalized_name)
        
        # Verify samples have goslin_result attribute added
        for s in samps:
            self.assertTrue(hasattr(s, "goslin_result"))
            self.assertIsInstance(s.goslin_result, dict)
            self.assertIn("input_name", s.goslin_result)
            self.assertIn("normalized_name", s.goslin_result)
            
    def test_cholesterol_annotation(self):
        """Test specific annotation for Cholesterol."""
        samps = [Sample("Cholesterol")]
        Goslin.annotate_samples(samps)
        
        # Verify sample was mutated with goslin_result
        self.assertTrue(hasattr(samps[0], "goslin_result"))
        self.assertEqual(samps[0].sample_name, "Cholesterol")
        self.assertEqual(samps[0].goslin_result["input_name"], "Cholesterol")
        self.assertIsNotNone(samps[0].goslin_result.get("normalized_name"))
        
    def test_phosphatidylcholine_annotation(self):
        """Test specific annotation for PC lipid."""
        samps = [Sample("PC(16:0/18:1)")]
        Goslin.annotate_samples(samps)
        
        # Verify sample was mutated with goslin_result
        self.assertTrue(hasattr(samps[0], "goslin_result"))
        self.assertEqual(samps[0].sample_name, "PC(16:0/18:1)")
        self.assertEqual(samps[0].goslin_result["input_name"], "PC(16:0/18:1)")
        self.assertIsNotNone(samps[0].goslin_result.get("normalized_name"))


if __name__ == "__main__":
    unittest.main()
