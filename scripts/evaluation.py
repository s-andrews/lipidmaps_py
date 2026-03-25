from lipidmaps import process_csv
dataset = process_csv("../tests/data/inputs/demo/BioPan_large_demo.csv", validate_data=True, use_refmet=True, use_headgroups=True, taxonomy_group="mammalia")
# print(dataset.list_headgroups())
# print(dataset.list_generic_lm_ids())