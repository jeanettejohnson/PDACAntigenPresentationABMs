# HPC-compatible copy of ../run_htan_singlecell_tme_geometries.jl -- submits
# to SLURM instead of running locally. See hpc_setup.jl for the job options.

ENV["PHYSICELL_CPP"] = "g++"

using PhysiCellModelManager
using CSV, DataFrames

initializeModelManager(
    joinpath(@__DIR__, "..", "PhysiCell"),
    joinpath(@__DIR__, "..", "data")
)

include(joinpath(@__DIR__, "hpc_setup.jl"))
include(joinpath(@__DIR__, "caf_mhc2.jl"))

df = CSV.read(joinpath(@__DIR__, "..", "assignmentsummary_HTAN_singlecell.csv"), DataFrame)

# Samples to run, as sample IDs or ID prefixes ("HT056" is all of that
# patient's samples). Empty runs every sample; set it for a test run. A later
# full run reuses the runs already made (use_previous=true).
const SUBSET = String[]

# CAF contact -> MHC-II induction rate, per minute (see caf_mhc2.jl). 0 is the
# baseline and leaves the runs as they are; 2.3e-4 enables it at about one
# conversion per 3 days of contact. Enabled runs are separate monads, so they
# never replace the baseline ones.
const CAF_MHC2_RATE = 0.0
const CAF_MHC2_VARIATIONS = cafMhc2Variations(CAF_MHC2_RATE)
CAF_MHC2_RATE > 0 && println("CAF-contact MHC-II induction on: $(CAF_MHC2_RATE) /min")
checkBaseRulesets("antigen_presentation_htan_singlecell")

if !isempty(SUBSET)
    df = df[[any(p -> startswith(s, p), SUBSET) for s in df.sample_id], :]
    isempty(df) && error("SUBSET $SUBSET matched none of the samples.")
    println("SUBSET $SUBSET -> $(nrow(df)) samples: ", join(df.sample_id, ", "))
end

inputs = InputFolders(
    "antigen_presentation_htan_singlecell",   # config
    "antigen_presentation_htan_singlecell";   # custom_code
    rulesets_collection = "antigen_presentation_htan_singlecell",
    ic_cell = "antigen_presentation_htan_singlecell"
)

# Six spatial geometry configurations (c1–c6): every way of splitting the three
# groups -- tumor (all tumor types and PDAC_unclassified), immune (CD4, CD8,
# CD8_exhausted, Treg) and stroma (CAF, apCAF) -- between the inner disk
# (0–200 μm) and the outer ring (200–400 μm). The two placements with every group
# in one region are left out; the well-mixed runs cover that.
#
# Ordered by tumor, then immune, then stroma, inner before outer. A label names
# the inner region, then the outer, each listing its groups in that order.
#
#     label                     tumor  immune  stroma   2026-09-25 runs
#     c1_tumor-immune_stroma    inner  inner   outer    (new)
#     c2_tumor-stroma_immune    inner  outer   inner    c1
#     c3_tumor_immune-stroma    inner  outer   outer    c2
#     c4_immune-stroma_tumor    outer  inner   inner    c3
#     c5_immune_tumor-stroma    outer  inner   outer    c4
#     c6_stroma_tumor-immune    outer  outer   inner    c5
#
# The 2026-09-25 runs numbered five of these c1–c5 in another order (last
# column), so a c-number does not mean the same layout in the two sets.
# spatial_config_index is the position in this list (1–6).
const GEOMETRY_CONFIGS = [
    (label="c1_tumor-immune_stroma", tum_inner=  0.0, tum_outer=200.0, imm_inner=  0.0, imm_outer=200.0, caf_inner=200.0, caf_outer=400.0),
    (label="c2_tumor-stroma_immune", tum_inner=  0.0, tum_outer=200.0, imm_inner=200.0, imm_outer=400.0, caf_inner=  0.0, caf_outer=200.0),
    (label="c3_tumor_immune-stroma", tum_inner=  0.0, tum_outer=200.0, imm_inner=200.0, imm_outer=400.0, caf_inner=200.0, caf_outer=400.0),
    (label="c4_immune-stroma_tumor", tum_inner=200.0, tum_outer=400.0, imm_inner=  0.0, imm_outer=200.0, caf_inner=  0.0, caf_outer=200.0),
    (label="c5_immune_tumor-stroma", tum_inner=200.0, tum_outer=400.0, imm_inner=  0.0, imm_outer=200.0, caf_inner=200.0, caf_outer=400.0),
    (label="c6_stroma_tumor-immune", tum_inner=200.0, tum_outer=400.0, imm_inner=200.0, imm_outer=400.0, caf_inner=  0.0, caf_outer=200.0),
]

# Build every (sample, geometry) monad up front (no jobs submitted yet), then
# run them all together in one Trial so the worker pool can submit up to
# setNumberOfParallelSims concurrently instead of waiting for each job to
# finish before starting the next.
monads = []
for row in eachrow(df)
    for config in GEOMETRY_CONFIGS
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

        # ── Geometry label — makes this (sample, config) combination unique in the DB ──
        config_index = findfirst(c -> c.label == config.label, GEOMETRY_CONFIGS)
        dv_config_tag = DiscreteVariation(configPath("user_parameter", "spatial_config_index"), config_index)

        dv_ecm  = DiscreteVariation(configPath("ecm", "initial_condition"), 1)

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

        # Scalar radius values from the current geometry config
        imm_inner = config.imm_inner
        imm_outer = config.imm_outer
        caf_inner = config.caf_inner
        caf_outer = config.caf_outer
        tum_inner = config.tum_inner
        tum_outer = config.tum_outer

        # ── Immune group: T cell subsets ──────────────────────────────────────────
        dv_cd4_inner  = DiscreteVariation(icCellsPath("CD4_Tcell",      "annulus", 1, "inner_radius"), imm_inner)
        dv_cd4_outer  = DiscreteVariation(icCellsPath("CD4_Tcell",      "annulus", 1, "outer_radius"), imm_outer)
        dv_cd8_inner  = DiscreteVariation(icCellsPath("CD8_Tcell",      "annulus", 1, "inner_radius"), imm_inner)
        dv_cd8_outer  = DiscreteVariation(icCellsPath("CD8_Tcell",      "annulus", 1, "outer_radius"), imm_outer)
        dv_treg_inner = DiscreteVariation(icCellsPath("Treg",           "annulus", 1, "inner_radius"), imm_inner)
        dv_treg_outer = DiscreteVariation(icCellsPath("Treg",           "annulus", 1, "outer_radius"), imm_outer)
        dv_exh_inner  = DiscreteVariation(icCellsPath("CD8_exhausted",  "annulus", 1, "inner_radius"), imm_inner)
        dv_exh_outer  = DiscreteVariation(icCellsPath("CD8_exhausted",  "annulus", 1, "outer_radius"), imm_outer)

        cv_imm_positions = CoVariation(
            dv_cd4_inner,  dv_cd4_outer,
            dv_cd8_inner,  dv_cd8_outer,
            dv_treg_inner, dv_treg_outer,
            dv_exh_inner,  dv_exh_outer,
        )

        # ── CAF group: CAF + apCAF ────────────────────────────────────────────────
        dv_caf_inner   = DiscreteVariation(icCellsPath("CAF",   "annulus", 1, "inner_radius"), caf_inner)
        dv_caf_outer   = DiscreteVariation(icCellsPath("CAF",   "annulus", 1, "outer_radius"), caf_outer)
        dv_apcaf_inner = DiscreteVariation(icCellsPath("apCAF", "annulus", 1, "inner_radius"), caf_inner)
        dv_apcaf_outer = DiscreteVariation(icCellsPath("apCAF", "annulus", 1, "outer_radius"), caf_outer)

        cv_caf_positions = CoVariation(
            dv_caf_inner,   dv_caf_outer,
            dv_apcaf_inner, dv_apcaf_outer,
        )

        # ── Tumor group: all tumor subtypes ───────────────────────────────────────
        dv_epi_inner          = DiscreteVariation(icCellsPath("epithelial_tumor",                "annulus", 1, "inner_radius"), tum_inner)
        dv_epi_outer          = DiscreteVariation(icCellsPath("epithelial_tumor",                "annulus", 1, "outer_radius"), tum_outer)
        dv_epi_c1_inner       = DiscreteVariation(icCellsPath("epithelial_tumor_class1",         "annulus", 1, "inner_radius"), tum_inner)
        dv_epi_c1_outer       = DiscreteVariation(icCellsPath("epithelial_tumor_class1",         "annulus", 1, "outer_radius"), tum_outer)
        dv_epi_c1c2_inner     = DiscreteVariation(icCellsPath("epithelial_tumor_class1_class2",  "annulus", 1, "inner_radius"), tum_inner)
        dv_epi_c1c2_outer     = DiscreteVariation(icCellsPath("epithelial_tumor_class1_class2",  "annulus", 1, "outer_radius"), tum_outer)
        dv_epi_c2_inner       = DiscreteVariation(icCellsPath("epithelial_tumor_class2",         "annulus", 1, "inner_radius"), tum_inner)
        dv_epi_c2_outer       = DiscreteVariation(icCellsPath("epithelial_tumor_class2",         "annulus", 1, "outer_radius"), tum_outer)
        dv_mes_inner          = DiscreteVariation(icCellsPath("mesenchymal_tumor",               "annulus", 1, "inner_radius"), tum_inner)
        dv_mes_outer          = DiscreteVariation(icCellsPath("mesenchymal_tumor",               "annulus", 1, "outer_radius"), tum_outer)
        dv_mes_c1_inner       = DiscreteVariation(icCellsPath("mesenchymal_tumor_class1",        "annulus", 1, "inner_radius"), tum_inner)
        dv_mes_c1_outer       = DiscreteVariation(icCellsPath("mesenchymal_tumor_class1",        "annulus", 1, "outer_radius"), tum_outer)
        dv_mes_c1c2_inner     = DiscreteVariation(icCellsPath("mesenchymal_tumor_class1_class2", "annulus", 1, "inner_radius"), tum_inner)
        dv_mes_c1c2_outer     = DiscreteVariation(icCellsPath("mesenchymal_tumor_class1_class2", "annulus", 1, "outer_radius"), tum_outer)
        dv_mes_c2_inner       = DiscreteVariation(icCellsPath("mesenchymal_tumor_class2",        "annulus", 1, "inner_radius"), tum_inner)
        dv_mes_c2_outer       = DiscreteVariation(icCellsPath("mesenchymal_tumor_class2",        "annulus", 1, "outer_radius"), tum_outer)
        dv_pdac_inner         = DiscreteVariation(icCellsPath("PDAC_unclassified",               "annulus", 1, "inner_radius"), tum_inner)
        dv_pdac_outer         = DiscreteVariation(icCellsPath("PDAC_unclassified",               "annulus", 1, "outer_radius"), tum_outer)

        cv_tum_positions = CoVariation(
            dv_epi_inner,      dv_epi_outer,
            dv_epi_c1_inner,   dv_epi_c1_outer,
            dv_epi_c1c2_inner, dv_epi_c1c2_outer,
            dv_epi_c2_inner,   dv_epi_c2_outer,
            dv_mes_inner,      dv_mes_outer,
            dv_mes_c1_inner,   dv_mes_c1_outer,
            dv_mes_c1c2_inner, dv_mes_c1c2_outer,
            dv_mes_c2_inner,   dv_mes_c2_outer,
            dv_pdac_inner,     dv_pdac_outer,
        )

        # ── Counts + config tag (all scalar, one value per simulation) ──────────────
        cv_counts = CoVariation(
            dv_config_tag,
            dv_ecm,
            dv_caf, dv_cd4, dv_cd8, dv_treg, dv_apcaf,
            dv_epithelial, dv_epithelial_class1, dv_epithelial_class1_class2, dv_epithelial_class2,
            dv_mesenchymal, dv_mesenchymal_class1, dv_mesenchymal_class1_class2, dv_mesenchymal_class2,
            dv_pdac_unspecified,
        )

        println("Queuing $(config.label) $sample  CAF=$caf_count  CD4=$cd4_count  CD8=$cd8_count  Treg=$treg_count  Epithelial=$epithelial_count  Mesenchymal=$mesenchymal_count  PDAC_unspecified=$pdac_unspecified_count")
        flush(stdout)
        push!(monads, createTrial(inputs, cv_counts, cv_imm_positions, cv_caf_positions, cv_tum_positions, CAF_MHC2_VARIATIONS...; n_replicates=1, use_previous=true))
    end
end

trial = createTrial(monads)
PhysiCellModelManager.run(trial)
