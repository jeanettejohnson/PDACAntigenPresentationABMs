# HPC-compatible copy of ../run_imc_spatial_configs.jl -- submits one SLURM
# job per ROI instead of running the compiled executable directly, capped at
# 50 concurrent via a bounded worker pool.
#
# Unlike the other scripts in this directory, ../run_imc_spatial_configs.jl
# bypasses PhysiCellModelManager entirely (see its header comment for why), so
# there's no Trial/Monad/useHPC mechanism to hook into here. This hand-rolls
# the same queue + bounded-worker-pool pattern ModelManager's own runner uses
# internally (see ModelManager/src/runner.jl), just built around raw `sbatch`
# calls instead.

ENV["PHYSICELL_CPP"] = "g++"

include(joinpath(@__DIR__, "slurm_common.jl"))

const PROJECT_ROOT = joinpath(@__DIR__, "..")
const PHYSICELL_DIR = joinpath(PROJECT_ROOT, "PhysiCell")
const CONFIG_DIR = joinpath(PHYSICELL_DIR, "config")
const EXECUTABLE = joinpath(PHYSICELL_DIR, "project")
const PROJ_NAME = "antigen_presentation"
const MAX_CONCURRENT = 50

# Only real ROI configs, e.g. PhysiCell_settings_JHH317ROI1.xml.
# Excludes the base PhysiCell_settings.xml and leftover sample-project
# templates (PDAC.xml, coculture.xml, therapy.xml, biorobots.xml, ...) that
# also live in this config dir.
const ROI_PATTERN = r"^PhysiCell_settings_(JHH\w+ROI\d+)\.xml$"

function discover_rois()
    rois = Tuple{String, String}[]  # (roi_name, xml_path)
    for f in sort(readdir(CONFIG_DIR))
        m = match(ROI_PATTERN, f)
        m === nothing && continue
        push!(rois, (m.captures[1], joinpath(CONFIG_DIR, f)))
    end
    return rois
end

function build_project()
    println("Building $PROJ_NAME (CC=$(ENV["PHYSICELL_CPP"]))...")
    flush(stdout)
    cd(PHYSICELL_DIR) do
        run(`make clean`)
        run(`make load PROJ=$PROJ_NAME`)
        run(`make -j$(Sys.CPU_THREADS)`)
    end
end

"""
    prepCmdForWrap(cmd::Cmd)

Strip surrounding backticks from the string representation of `cmd`, matching
ModelManager's own `prepareHPCCommand` so a multi-token command becomes the
single string sbatch's --wrap expects.
"""
prepCmdForWrap(cmd::Cmd) = strip(string(cmd), '`')

function submitROI(roi_name::String, xml_path::String; skip_existing::Bool=true)
    output_dir = joinpath(PHYSICELL_DIR, "outputs", roi_name)

    if skip_existing && isfile(joinpath(output_dir, "final.xml"))
        println("SKIP $roi_name (final.xml already present)")
        flush(stdout)
        return :skipped
    end

    mkpath(output_dir)
    out_path = joinpath(output_dir, "run.out")
    err_path = joinpath(output_dir, "run.err")
    xml_relpath = relpath(xml_path, PHYSICELL_DIR)

    println("Submitting $roi_name  (config: $xml_relpath)")
    flush(stdout)
    wrap_cmd = prepCmdForWrap(`$EXECUTABLE $xml_relpath`)
    flags = [
        "--wrap=$wrap_cmd",
        "--wait",
        "--account=$(slurmAccount())",
        "--time=$JOB_TIME",
        "--cpus-per-task=$JOB_CPUS_PER_TASK",
        "--mem=$JOB_MEM",
        "--job-name=ROI_$roi_name",
        "--output=$out_path",
        "--error=$err_path",
        "--chdir=$PHYSICELL_DIR",
    ]
    try
        run(`sbatch $flags`)
        return :ok
    catch err
        println("FAILED $roi_name: $err")
        flush(stdout)
        return :failed
    end
end

"""
    submitAllROIs(rois; skip_existing, max_concurrent) -> Dict

Submit one `sbatch --wait` job per ROI through a bounded worker pool (at most
`max_concurrent` jobs in flight at once), collecting each result as it
completes. Mirrors ModelManager's queue + worker-pool pattern in runner.jl.
"""
function submitAllROIs(rois; skip_existing::Bool=true, max_concurrent::Int=MAX_CONCURRENT)
    tasks = [
        @task (roi_name, submitROI(roi_name, xml_path; skip_existing=skip_existing))
        for (roi_name, xml_path) in rois
    ]

    queue = Channel{Task}(length(tasks))
    results_ch = Channel{Tuple{String,Symbol}}(length(tasks))
    @async begin
        for t in tasks
            put!(queue, t)
        end
        close(queue)
    end
    for _ in 1:max_concurrent
        @async for t in queue
            schedule(t)
            put!(results_ch, fetch(t))
        end
    end

    results = Dict(:ok => String[], :skipped => String[], :failed => String[])
    for _ in 1:length(tasks)
        roi_name, status = take!(results_ch)
        push!(results[status], roi_name)
    end
    return results
end

function main(; skip_existing::Bool=true, rebuild::Bool=true)
    rebuild && build_project()

    rois = discover_rois()
    println("Found $(length(rois)) ROI configs.\n")
    flush(stdout)

    results = submitAllROIs(rois; skip_existing=skip_existing)

    println("\nDone. ok=$(length(results[:ok]))  skipped=$(length(results[:skipped]))  failed=$(length(results[:failed]))")
    isempty(results[:failed]) || println("Failed ROIs: ", join(results[:failed], ", "))
    flush(stdout)
    return results
end

main()
