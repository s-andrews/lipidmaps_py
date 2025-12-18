import unittest
from lipidmaps.data.data_manager import DataManager
from lipidmaps.data.models.sample import QuantifiedLipid, LipidDataset, SampleMetadata
from lipidmaps.data.models.reaction import Reaction

class TestAnnotateLipidsWithReactions(unittest.TestCase):
    def setUp(self):
        # Create mock lipids
        self.lipids = [
            QuantifiedLipid(input_name="PC(16:0/18:1)", values={"S1": 1.0}, lm_id="LMGP01010001"),
            QuantifiedLipid(input_name="LPC(16:0)", values={"S1": 2.0}, lm_id="LMGP02010001"),
        ]
        self.samples = [SampleMetadata(sample_id="S1", group="Control")]
        self.dataset = LipidDataset(samples=self.samples, lipids=self.lipids)
        self.manager = DataManager()
        self.manager.dataset = self.dataset

        # Create mock reactions
        self.reactions = [
            Reaction(
                reaction_id="R1",
                reaction_name="PC to LPC",
                reactants=[{"lm_id": "LMGP01010001", "input_name": "PC(16:0/18:1)"}],
                products=[{"lm_id": "LMGP02010001", "input_name": "LPC(16:0)"}],
                type="class-level",
                pathway_id=None,
                enzyme_id=None,
            ),
            Reaction(
                reaction_id="R2",
                reaction_name="LPC to PC",
                reactants=[{"lm_id": "LMGP02010001", "input_name": "LPC(16:0)"}],
                products=[{"lm_id": "LMGP01010001", "input_name": "PC(16:0/18:1)"}],
                type="class-level",
                pathway_id=None,
                enzyme_id=None,
            ),
        ]

    def test_annotate_lipids_with_reactions(self):
        self.manager.annotate_lipids_with_reactions(self.reactions)
        # Check that each lipid has reactions
        for lipid in self.lipids:
            self.assertIsInstance(lipid.reactions, list)
            self.assertGreaterEqual(len(lipid.reactions), 1)
            # Check that the reaction_id is present in the reactions
            reaction_ids = [r["reaction_id"] for r in lipid.reactions]
            self.assertTrue(any(rid in ["R1", "R2"] for rid in reaction_ids))

if __name__ == "__main__":
    unittest.main()
