# PCMM-based IMC spatial runner.
#
# Replaces the earlier hand-rolled version that bypassed PhysiCellModelManager
# and drove `sbatch` directly. Every per-ROI difference in the 48
# PhysiCell_settings_<ROI>.xml files is expressible through PCMM:
#
#   - cell positions / substrate field -> ic_cell and ic_substrate input
#     folders, one per ROI, which PCMM passes to PhysiCell as -i and -s
#   - 4 domain bounds + 26 cell volumes -> config DiscreteVariations, bundled
#     per ROI into a CoVariation so they resolve to one config_variation_id
#   - save/folder -> irrelevant; PCMM owns output paths and passes -o
#
# Because each ROI carries its own ic_cell/ic_substrate folders, each ROI forms
# its own Sampling (Sampling asserts all its monads share input folders) and the
# Trial groups all 48. Adding a within-ROI parameter sweep later would populate
# each Sampling with several monads without changing this shape.
#
# Run prep_imc_spatial/setup_imc_spatial_pcmm.py first -- it builds the PCMM input folders and
# writes imc_spatial_roi_specs.csv, both derived from the 48 ROI configs.

ENV["PHYSICELL_CPP"] = "g++"

using PhysiCellModelManager
using CSV, DataFrames

initializeModelManager(
    joinpath(@__DIR__, "..", "PhysiCell"),
    joinpath(@__DIR__, "..", "data")
)

include(joinpath(@__DIR__, "hpc_setup.jl"))
include(joinpath(@__DIR__, "caf_mhc2.jl"))

# No per-script setJobOptions override: the base config is written with
# <omp_num_threads>1</omp_num_threads> by prep_imc_spatial/setup_imc_spatial_pcmm.py, so the
# shared hpc_setup.jl request (1 CPU, 2G) is already correct. Measured headroom
# is large -- a JHH368ROI4 run (3937 cells -> 7978 agents over 10080 min) peaked
# at 52 MB RSS, against the 2G request.

# Replicates per ROI. Each is an independent seed: the configs use
# random_seed=system_clock, so replicates differ only by randomness.
const N_REPLICATES = 1

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

const SPEC_PATH = joinpath(@__DIR__, "..", "prep_imc_spatial", "imc_spatial_roi_specs.csv")

# The spec table is tracked in git, so this should not fire on a clean clone.
# Fail early and specifically rather than letting PCMM report a missing folder.
if !isfile(SPEC_PATH)
    error("""
          Missing $(basename(SPEC_PATH)).
          This file is tracked and ships with the repo, so it should already be
          present -- restore it with:
              git checkout -- prep_imc_spatial/imc_spatial_roi_specs.csv
          Regenerate it only if the initial conditions were rebuilt:
              python3 prep_imc_spatial/setup_imc_spatial_pcmm.py
          """)
end

df = CSV.read(SPEC_PATH, DataFrame)

if !isempty(SUBSET)
    df = df[[any(p -> startswith(roi, p), SUBSET) for roi in df.roi], :]
    isempty(df) && error("SUBSET $SUBSET matched none of the ROIs in $(basename(SPEC_PATH)).")
    println("SUBSET $SUBSET -> $(nrow(df)) of 48 ROIs: ", join(df.roi, ", "))
end

# Derive the varying cell types from the spec table rather than restating them,
# so this stays in step with prep_imc_spatial/setup_imc_spatial_pcmm.py.
volume_columns = filter(n -> occursin(r"^vol_.+_total$", n), names(df))
cell_types = [match(r"^vol_(.+)_total$", c).captures[1] for c in volume_columns]

println("Building $(nrow(df)) ROI samplings ($(length(cell_types)) varying cell volumes each)")
flush(stdout)

samplings = Sampling[]
for row in eachrow(df)
    # CSV.jl yields InlineStrings (String15) for short string columns; the
    # InputFolders constructor requires plain String.
    roi = String(row.roi)

    # All locations are shared except the two ICs, which are per-ROI.
    inputs = InputFolders(
        PROJ,                             # config
        PROJ;                             # custom_code
        rulesets_collection = PROJ,
        ic_cell             = roi,
        ic_substrate        = roi,
    )

    dvs = DiscreteVariation[]
    append!(dvs, domainVariations(x_min=row.x_min, x_max=row.x_max,
                                  y_min=row.y_min, y_max=row.y_max))
    for ct in cell_types
        push!(dvs, DiscreteVariation(configPath(ct, "total"),   row[Symbol("vol_$(ct)_total")]))
        push!(dvs, DiscreteVariation(configPath(ct, "nuclear"), row[Symbol("vol_$(ct)_nuclear")]))
    end

    # CoVariation binds the 30 values so they move together. Passed as
    # independent variations they would expand into a full-factorial grid
    # instead of a single parameter set for this ROI.
    #
    # Pass the vector, not `dvs...`: the vararg constructor requires a single
    # element type, and these mix DiscreteVariation{Float64} (domain bounds)
    # with DiscreteVariation{Int64} (volumes).
    cv = CoVariation(dvs)

    # createTrial returns a Simulation when n_replicates == 1 and a Monad
    # otherwise; Sampling accepts either and normalizes to a Monad internally.
    trial_piece = createTrial(inputs, cv, CAF_MHC2_VARIATIONS...; n_replicates=N_REPLICATES, use_previous=true)
    sampling = Sampling(trial_piece; n_replicates=N_REPLICATES, use_previous=true)

    println("Queuing $roi  domain x[$(row.x_min),$(row.x_max)] y[$(row.y_min),$(row.y_max)]")
    flush(stdout)
    push!(samplings, sampling)
end

trial = Trial(samplings)
PhysiCellModelManager.run(trial)
