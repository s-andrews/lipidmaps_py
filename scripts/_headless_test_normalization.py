# Headless normalization test
import sys
from pprint import pprint
sys.path.insert(0, 'src')

from lipidmaps.data.models.sample import LipidDataset, QuantifiedLipid, SampleMetadata
from lipidmaps.data.quantitation import QuantitationAnalyzer

# Build synthetic dataset
samples = [SampleMetadata(sample_name=f'Sample{i}', group='A') for i in (1,2,3)]
lipids = [
    QuantifiedLipid(input_name='Lipid1', values={'Sample1': 10.0, 'Sample2': 20.0, 'Sample3': 30.0}),
    QuantifiedLipid(input_name='Lipid2', values={'Sample1': 5.0, 'Sample2': 15.0, 'Sample3': 25.0}),
    QuantifiedLipid(input_name='Lipid3', values={'Sample1': 2.0, 'Sample2': 4.0, 'Sample3': 8.0}),
]

ds = LipidDataset(samples=samples, lipids=lipids)

# Instantiate analyzer (use keyword arg to satisfy pydantic v2)
try:
    analyzer = QuantitationAnalyzer(dataset=ds)
except TypeError:
    analyzer = QuantitationAnalyzer(dataset=ds)  # best-effort

# Run total-lipid normalization (scale to 1e6)
try:
    norm = analyzer.normalize_total_lipid(scale_factor=1e6)
except Exception as e:
    print('normalize_total_lipid failed:', e)
    raise

method_key = 'total_lipid:scale=1e6'
# store per-lipid NormalizedResult via set_normalized
for ln, vals in norm.items():
    lipid = next((l for l in ds.lipids if l.input_name == ln), None)
    if lipid is not None:
        try:
            lipid.set_normalized(method_key, vals=vals)
        except Exception:
            # best-effort fallback: attach raw dict
            setattr(lipid, 'normalized', getattr(lipid, 'normalized', {}))
            lipid.normalized[method_key] = vals

print('Normalized mapping (first 3 items):')
pprint({k: v for k, v in list(norm.items())[:3]})

# Build DataFrame via dataset helper
try:
    df = ds.normalized_dataframe(method_key)
    print('\nNormalized DataFrame:')
    print(df.to_string())
except Exception as e:
    print('normalized_dataframe failed:', e)

print('\nPer-lipid stored normalized keys:')
for l in ds.lipids:
    keys = getattr(l, 'normalized', {}).keys()
    print(l.input_name, list(keys))
