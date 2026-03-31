from lipidmaps import process_csv
dataset = process_csv("/lipidmaps/lipidmaps_py/tests/data/inputs/demo/BioPan_large_demo.csv", 
                      validate_data=True, 
                      use_refmet=True, 
                      use_headgroups=True, 
                      taxonomy_group="mammalia")
# print(dataset.list_headgroups())
# print(dataset.list_generic_lm_ids())

# for lipid in dataset.lipids:
#     if "FA" in lipid.input_name:
#         print(lipid.structure)