````markdown
# LIPID MAPS Python API

Welcome — Examples for common workflows using the `lipidmaps` package.

Quick import

```python
from lipidmaps import process_csv
ds = process_csv("data/myfile.csv", validate_data=True, use_refmet=True)
```

Inspect dataset

```python
print(ds.list_sample_names()[:5])
print(ds.list_lipid_names()[:5])
```

Fetch reactions (by LM ID)

```python
reactions = ds.fetch_reactions_by_lm_id(reaction_type="species-level", only_lipid_components=False)
print(f"Fetched {len(reactions)} reactions")
lipids_in_rx = ds.get_lipids_for_reaction(reactions[0])
```

Normalization (minimal)

```python
from lipidmaps.data.quantitation import QuantitationAnalyzer
qa = QuantitationAnalyzer(dataset=ds)
method_key = "total_lipid:scale=1e6"
norm_map = qa.normalize_total_lipid(scale_factor=1e6)
for l in ds.lipids:
    if l.input_name in norm_map:
        l.set_normalized(method_key, norm_map[l.input_name])
df_norm = ds.normalized_dataframe(method_key)
print(df_norm.head())
```

More examples and detailed API usage are available in the repository `README.md` and the docs sidebar.
````
