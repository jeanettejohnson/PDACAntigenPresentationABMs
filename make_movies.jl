ENV["PHYSICELL_CPP"] = "g++-15"

using PhysiCellModelManager

initializeModelManager(
    joinpath(@__DIR__, "PhysiCell"),
    joinpath(@__DIR__, "data")
)

completed_ids = [10, 11, 12, 13]

for sim_id in completed_ids
    println("Making movie for simulation $sim_id...")
    try
        makeMovie(sim_id)
        println("  Done: simulation $sim_id")
    catch e
        println("  ERROR for simulation $sim_id: $e")
    end
end
