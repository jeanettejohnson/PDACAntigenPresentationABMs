ENV["PHYSICELL_CPP"] = "g++-15"

using PhysiCellModelManager
using CSV, DataFrames

initializeModelManager(
    joinpath(@__DIR__, "PhysiCell"),
    joinpath(@__DIR__, "data")
)


df = CSV.read(expanduser("~/assignmentsummary_HTAN_singlecell.csv"), DataFrame)

inputs = InputFolders(
    "antigen_presentation_htan_singlecell",   # config
    "antigen_presentation_htan_singlecell";   # custom_code
    rulesets_collection = "antigen_presentation_htan_singlecell",
    ic_cell = "antigen_presentation_htan_singlecell"
)

for row in eachrow(df)
    sample    = row.sample_id
    caf_count  = round(Int, row.CAF)
    cd4_count  = round(Int, row.CD4_T)
    cd8_count  = round(Int, row.CD8_T) + round(Int, row.CD8_T_cytotoxic) + round(Int, row.Proliferating_T)
    treg_count = round(Int, row.Treg)
    apcaf_count = round(Int, row.apCAF)

    epithelial_count = round(Int, row.Pattern2_Pattern7) + round(Int, row.Pattern2)
    epithelial_class1_count = round(Int, row.Pattern2_Pattern7_class_1) + round(Int, row.Pattern2_class_1)
    epithelial_class1_class2_count = round(Int, row.Pattern2_Pattern7_class_1_class_2) + round(Int, row.Pattern2_class_1_class_2)
    epithelial_class2_count = round(Int, row.Pattern2_Pattern7_class_2) + round(Int, row.Pattern2_class_2)

    mesenchymal_count = round(Int, row.Pattern7)
    mesenchymal_class1_count = round(Int, row.Pattern7_class_1)
    mesenchymal_class1_class2_count = round(Int, row.Pattern7_class_1_class_2)
    mesenchymal_class2_count = round(Int, row.Pattern7_class_2)

    pdac_unspecified_count = round(Int, row.PDAC_unclassified)

    dv_caf  = DiscreteVariation(icCellsPath("CAF",                       "annulus", 1, "number"), caf_count)
    dv_cd4  = DiscreteVariation(icCellsPath("CD4_Tcell",                 "annulus", 1, "number"), cd4_count)
    dv_cd8  = DiscreteVariation(icCellsPath("CD8_Tcell",                 "annulus", 1, "number"), cd8_count)
    dv_treg = DiscreteVariation(icCellsPath("Treg",                      "annulus", 1, "number"), treg_count)
    dv_apcaf = DiscreteVariation(icCellsPath("apCAF",                    "annulus", 1, "number"), apcaf_count)
    dv_epithelial = DiscreteVariation(icCellsPath("epithelial_tumor",                       "annulus", 1, "number"), epithelial_count)
    dv_epithelial_class1 = DiscreteVariation(icCellsPath("epithelial_tumor_class1",         "annulus", 1, "number"), epithelial_class1_count)
    dv_epithelial_class1_class2 = DiscreteVariation(icCellsPath("epithelial_tumor_class1_class2", "annulus", 1, "number"), epithelial_class1_class2_count)
    dv_epithelial_class2 = DiscreteVariation(icCellsPath("epithelial_tumor_class2",         "annulus", 1, "number"), epithelial_class2_count)
    dv_mesenchymal = DiscreteVariation(icCellsPath("mesenchymal_tumor",                       "annulus", 1, "number"), mesenchymal_count)
    dv_mesenchymal_class1 = DiscreteVariation(icCellsPath("mesenchymal_tumor_class1",         "annulus", 1, "number"), mesenchymal_class1_count)
    dv_mesenchymal_class1_class2 = DiscreteVariation(icCellsPath("mesenchymal_tumor_class1_class2", "annulus", 1, "number"), mesenchymal_class1_class2_count)
    dv_mesenchymal_class2 = DiscreteVariation(icCellsPath("mesenchymal_tumor_class2",         "annulus", 1, "number"), mesenchymal_class2_count)
    dv_pdac_unspecified = DiscreteVariation(icCellsPath("PDAC_unclassified",             "annulus", 1, "number"), pdac_unspecified_count)

    cv = CoVariation(dv_caf, dv_cd4, dv_cd8, dv_treg, dv_apcaf, dv_epithelial, dv_epithelial_class1, dv_epithelial_class1_class2, dv_epithelial_class2, dv_mesenchymal, dv_mesenchymal_class1, dv_mesenchymal_class1_class2, dv_mesenchymal_class2, dv_pdac_unspecified)

    println("Running $sample  CAF=$caf_count  CD4=$cd4_count  CD8=$cd8_count  Treg=$treg_count Epithelial=$epithelial_count  Mesenchymal=$mesenchymal_count  PDAC_unspecified=$pdac_unspecified_count")
    monad = createTrial(inputs, cv; n_replicates=1, use_previous=false)
    PhysiCellModelManager.run(monad)
end
