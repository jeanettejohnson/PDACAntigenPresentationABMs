# HPC-compatible copy of ../test_run_wellmixed_imc.jl -- submits to SLURM
# instead of running locally. Test version: small sample set, no reuse of
# previous results. See hpc_setup.jl for the job options.

ENV["PHYSICELL_CPP"] = "g++"

using PhysiCellModelManager
using CSV, DataFrames

initializeModelManager(
    joinpath(@__DIR__, "..", "PhysiCell"),
    joinpath(@__DIR__, "..", "data")
)

include(joinpath(@__DIR__, "hpc_setup.jl"))

df = CSV.read(joinpath(@__DIR__, "..", "assignmentsummary_JHH_IMC_test.csv"), DataFrame)

inputs = InputFolders(
    "antigen_presentation_htan_singlecell",   # config
    "antigen_presentation_htan_singlecell";   # custom_code
    rulesets_collection = "antigen_presentation_htan_singlecell",
    ic_cell = "antigen_presentation_htan_singlecell"
)

# All columns except these bookkeeping ones are cell-type count columns; each
# name must match a <cell_patches name=...> entry in the well-mixed cells.xml.
bookkeeping = ["sample_id", "total"]
cell_types = filter(n -> !(n in bookkeeping), names(df))

# Build every sample's monad up front (no jobs submitted yet), then run them
# all together in one Trial so the worker pool can submit up to
# setNumberOfParallelSims concurrently instead of waiting for each sample's
# job to finish before starting the next.
monads = []
for row in eachrow(df)
    sample = row.sample_id

    # ECM initial condition on, matching the original well-mixed driver.
    dvs = DiscreteVariation[DiscreteVariation(configPath("ecm", "initial_condition"), 1)]

    for ct in cell_types
        n = round(Int, row[ct])
        push!(dvs, DiscreteVariation(icCellsPath(ct, "annulus", 1, "number"), n))
    end

    cv = CoVariation(dvs...)

    println("Queuing $sample  ($(round(Int, row.total)) cells across $(length(cell_types)) types)")
    flush(stdout)
    push!(monads, createTrial(inputs, cv; n_replicates=1, use_previous=false))
end

trial = createTrial(monads)
PhysiCellModelManager.run(trial)
