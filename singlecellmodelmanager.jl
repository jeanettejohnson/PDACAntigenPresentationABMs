xml_path = icCellsPath("CAF", "annulus", 1, "inner_radius")

 dv = DiscreteVariation(xml_path, [0.0; 200.0; 0.0; 200.0; 0.0; 200.0])


dv_caf_inner = DiscreteVariation(xml_path, [0.0; 200.0; 0.0; 200.0; 0.0; 200.0])

xml_path = icCellsPath("CAF", "annulus", 1, "outer_radius")


dv_caf_outer = DiscreteVariation(xml_path, [200.0; 400.0; 200.0; 400.0; 200.0; 400.0])

cv = CoVariation(dv_caf_inner, dv_caf_outer)


dv_caf_outer = _ # length 6
dv_caf_inner = _
dv_imm_outer = _
dv_imm_inner = _
dv_tum_outer = _
dv_tum_inner = _

cv = CoVariation(dv_caf_outer, dv_caf_inner, dv_imm_outer, dv_imm_inner, dv_tum_outer, dv_tum_inner)

# for each patient, get counts of these
dv_caf_count = [500, 200, 300]
dv_imm_count = _
dv_tum_count = _

cv_counts = CoVariation(dv_caf_count, dv_imm_count, dv_tum_count)

run(inputs, cv, cv_counts; n_replicates=3)

ic_cell_folder_names = ["tissue_1", "tissue_2", "tissue_3"] # folders in data/inputs/ics/cells/
for f in ic_cell_folder_names
    inputs = InputFolders(config_folder, custom_code_folder;
        ic_cells_folder=f, 
        rulesets_collection_folder = "rulesets_collection_1")
    out = run(inputs, cv; n_replicates=3)
end