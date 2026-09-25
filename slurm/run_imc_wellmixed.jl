# PCMM runner for IMC well-mixed: the IMC spatial model with the cells placed
# at random. Same config, rules, custom code and per-ROI cell volumes as
# run_imc_spatial.jl; only the initial condition differs:
#
#   - cells: <ROI>_wellmixed_r<k>, one random layout per replicate, built by
#     prep_imc_spatial/make_imc_wellmixed_ics.py -- the ROI's biological cells
#     in a disk with the ROI's free area, no overlaps, measured volumes kept
#   - ECM: <ROI>_wellmixed, uniform at the mean of the ROI's ECM image
#   - domain: +/-800 um, as HTAN well-mixed
#
# Each layout is its own IC folder, so PCMM treats it as its own monad rather
# than a replicate of the others: PCMM counts runs as replicates only when they
# share input folders. To add replicates, raise n_layouts in the prep script,
# run it, commit the new layouts, raise N_LAYOUTS here and run this again --
# PCMM reuses the runs already made.
#
# The previous version, which ran these samples on the HTAN config and ruleset,
# is slurm/archive/run_imc_wellmixed_htan_config.jl.

ENV["PHYSICELL_CPP"] = "g++"

using PhysiCellModelManager
using CSV, DataFrames

initializeModelManager(
    joinpath(@__DIR__, "..", "PhysiCell"),
    joinpath(@__DIR__, "..", "data")
)

include(joinpath(@__DIR__, "hpc_setup.jl"))
include(joinpath(@__DIR__, "caf_mhc2.jl"))

# Random layouts per ROI to run; must not exceed the n_layouts the prep script
# built.
const N_LAYOUTS = 1

# PhysiCell-seed replicates of each layout (the configs use
# random_seed=system_clock). 1 means each layout runs once.
const N_REPLICATES_PER_LAYOUT = 1

# ROIs to run, as ROI keys or prefixes ("JHH368" is its 4 ROIs). Empty runs
# all 48; set it for a test run. A later full run reuses the runs already made
# (use_previous=true).
const SUBSET = String[]

# CAF contact -> MHC-II induction rate, per minute (see caf_mhc2.jl). 0 is the
# baseline and leaves the runs as they are; 2.3e-4 enables it at about one
# conversion per 3 days of contact. Enabled runs are separate monads, so they
# never replace the baseline ones.
const CAF_MHC2_RATE = 0.0
const CAF_MHC2_VARIATIONS = cafMhc2Variations(CAF_MHC2_RATE)
CAF_MHC2_RATE > 0 && println("CAF-contact MHC-II induction on: $(CAF_MHC2_RATE) /min")
checkBaseRulesets("antigen_presentation")

const PROJ = "antigen_presentation"
const DOMAIN_HALF_WIDTH = 800.0

const SPEC_PATH = joinpath(@__DIR__, "..", "prep_imc_spatial", "imc_spatial_roi_specs.csv")
const IC_ROOT = joinpath(@__DIR__, "..", "data", "inputs", "ics")

df = CSV.read(SPEC_PATH, DataFrame)

if !isempty(SUBSET)
    df = df[[any(p -> startswith(roi, p), SUBSET) for roi in df.roi], :]
    isempty(df) && error("SUBSET $SUBSET matched none of the ROIs in $(basename(SPEC_PATH)).")
    println("SUBSET $SUBSET -> $(nrow(df)) of 48 ROIs: ", join(df.roi, ", "))
end

# Fail before queuing anything if a layout or ECM field was never built.
missing_inputs = String[]
for roi in df.roi, k in 1:N_LAYOUTS
    cells = joinpath(IC_ROOT, "cells", "$(roi)_wellmixed_r$(k)", "cells.csv")
    isfile(cells) || push!(missing_inputs, cells)
end
for roi in df.roi
    ecm = joinpath(IC_ROOT, "substrates", "$(roi)_wellmixed", "substrates.csv")
    isfile(ecm) || push!(missing_inputs, ecm)
end
if !isempty(missing_inputs)
    error("""
          Missing IMC well-mixed inputs ($(length(missing_inputs)), first: $(first(missing_inputs))).
          Build them with prep_imc_spatial/make_imc_wellmixed_ics.py (n_layouts >= N_LAYOUTS = $N_LAYOUTS).
          """)
end

# The varying cell types come from the spec table, as in run_imc_spatial.jl.
volume_columns = filter(n -> occursin(r"^vol_.+_total$", n), names(df))
cell_types = [match(r"^vol_(.+)_total$", c).captures[1] for c in volume_columns]

samplings = Sampling[]
for row in eachrow(df), k in 1:N_LAYOUTS
    roi = String(row.roi)

    inputs = InputFolders(
        PROJ,                                     # config
        PROJ;                                     # custom_code
        rulesets_collection = PROJ,
        ic_cell             = "$(roi)_wellmixed_r$(k)",
        ic_substrate        = "$(roi)_wellmixed",
    )

    dvs = DiscreteVariation[]
    append!(dvs, domainVariations(x_min=-DOMAIN_HALF_WIDTH, x_max=DOMAIN_HALF_WIDTH,
                                  y_min=-DOMAIN_HALF_WIDTH, y_max=DOMAIN_HALF_WIDTH))
    for ct in cell_types
        push!(dvs, DiscreteVariation(configPath(ct, "total"),   row[Symbol("vol_$(ct)_total")]))
        push!(dvs, DiscreteVariation(configPath(ct, "nuclear"), row[Symbol("vol_$(ct)_nuclear")]))
    end

    # One parameter set per ROI: pass the vector so the values move together
    # (see run_imc_spatial.jl).
    cv = CoVariation(dvs)

    trial_piece = createTrial(inputs, cv, CAF_MHC2_VARIATIONS...; n_replicates=N_REPLICATES_PER_LAYOUT, use_previous=true)
    push!(samplings, Sampling(trial_piece; n_replicates=N_REPLICATES_PER_LAYOUT, use_previous=true))

    println("Queuing $roi layout $k  domain +/-$(DOMAIN_HALF_WIDTH)")
    flush(stdout)
end

trial = Trial(samplings)
PhysiCellModelManager.run(trial)
