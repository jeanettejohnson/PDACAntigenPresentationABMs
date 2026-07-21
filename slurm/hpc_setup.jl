# Shared PCMM/HPC setup for the PhysiCellModelManager-based scripts in this
# directory. `include` this after `initializeModelManager(...)`. Submits
# simulations to SLURM instead of running them locally, capped at 50
# concurrent jobs.

include(joinpath(@__DIR__, "slurm_common.jl"))

useHPC(true)
setJobOptions(Dict(
    "cpus-per-task" => JOB_CPUS_PER_TASK,
    "mem" => JOB_MEM,
    "time" => JOB_TIME,
    "account" => slurmAccount(),
))
setNumberOfParallelSims(50)
