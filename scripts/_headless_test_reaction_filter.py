# Headless test for reaction filtering using headgroup rules
import sys
sys.path.insert(0, 'src')

from lipidmaps.data.utils.lipid_reaction_rules import reactions_possible_in_dataset
from lipidmaps.data.models.sample import LipidDataset, QuantifiedLipid, SampleMetadata

# Create dataset with PC lipids
samples = [SampleMetadata(sample_name='S1', group='A'), SampleMetadata(sample_name='S2', group='A')]
lipids = [
    QuantifiedLipid(input_name='PC(16:0/18:1)', values={'S1': 10, 'S2': 12}),
    QuantifiedLipid(input_name='PE(16:0/18:1)', values={'S1': 5, 'S2': 6}),
]
ds = LipidDataset(samples=samples, lipids=lipids)

# Fake reaction objects (dict-like) with compound_name fields
reactions = [
    {
        'reaction_id': 'R1',
        'reaction_name': 'PC_to_PA',
        'reactants': [{'compound_name': 'PC(16:0/18:1)'}],
        'products': [{'compound_name': 'PA(16:0/18:1)'}],
    },
    {
        'reaction_id': 'R2',
        'reaction_name': 'TAG_maker',
        'reactants': [{'compound_name': 'DAG(18:1/18:1)'}],
        'products': [{'compound_name': 'TAG(18:1/18:1/18:1)'}],
    },
    {
        'reaction_id': 'R3',
        'reaction_name': 'PE_to_LPE',
        'reactants': [{'compound_name': 'PE(16:0/18:1)'}],
        'products': [{'compound_name': 'LPE(16:0)'}],
    },
]

possible = reactions_possible_in_dataset(ds, reactions)
print('Possible reactions IDs:', [r.get('reaction_id') if isinstance(r, dict) else getattr(r, 'reaction_id', None) for r in possible])
