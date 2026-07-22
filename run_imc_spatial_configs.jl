# run_imc_spatial_configs.jl
# -----------------------------------------------------------------------------
# Runs every per-ROI IMC simulation from its real spatial layout, instead of
# the well-mixed annulus scatter used by run_wellmixed_imc.jl.
#
# This targets PhysiCell/user_projects/antigen_presentation, NOT the pcmm
# custom code (antigen_presentation_htan_singlecell) used by the
# run_wellmixed_*.jl / run_htan_singlecell_tme_geometries.jl scripts. These
# are two separate custom-code lineages: user_projects/antigen_presentation is
# a standard PhysiCell user project (its own main.cpp/custom_modules), and its
# config/ already holds one PhysiCell_settings_<ROI>.xml per ROI (via
# generate_roi_configs.py) with the ROI's real cell positions and ECM baked
# in directly. Nothing here is expressible as a DiscreteVariation, so this
# script bypasses PhysiCellModelManager entirely: it builds via the standard
# PhysiCell `make load PROJ=...` workflow, then invokes the compiled
# executable once per ROI.
#
# NOTE: this project's main.cpp takes the config path as a plain positional
# argument (argv[1]), not a -c flag -- unlike the pcmm custom code's main.cpp.
# -----------------------------------------------------------------------------

# g++ compiler is machine dependent -- plain "g++" on macOS is an Apple
# clang alias without real OpenMP support; use Homebrew's real GNU g++.
ENV["PHYSICELL_CPP"] = "g++-16"

const PROJECT_ROOT = @__DIR__
const PHYSICELL_DIR = joinpath(PROJECT_ROOT, "PhysiCell")
const PROJ_NAME = "antigen_presentation"
# NOTE: discover from user_projects/antigen_presentation/config, NOT
# PhysiCell/config -- the latter is a local "deploy" copy that is NOT
# git-tracked (confirmed via `git ls-files`) and will be empty/stale on any
# fresh clone (e.g. an HPC cluster). user_projects/antigen_presentation/config
# is the tracked, canonical source; `make load` copies it into PhysiCell/config
# at build time, so we only need to read the tracked copy here.
const CONFIG_DIR = joinpath(PHYSICELL_DIR, "user_projects", PROJ_NAME, "config")
const EXECUTABLE = joinpath(PHYSICELL_DIR, "project")

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
    cd(PHYSICELL_DIR) do
        run(`make clean`)
        run(`make load PROJ=$PROJ_NAME`)
        run(`make -j$(Sys.CPU_THREADS)`)
    end
end

function run_roi(roi_name::String, xml_path::String; skip_existing::Bool=true)
    output_dir = joinpath(PHYSICELL_DIR, "outputs", roi_name)

    if skip_existing && isfile(joinpath(output_dir, "final.xml"))
        println("SKIP $roi_name (final.xml already present)")
        return :skipped
    end

    mkpath(output_dir)
    log_path = joinpath(output_dir, "run.log")
    xml_relpath = relpath(xml_path, PHYSICELL_DIR)

    println("Running $roi_name  (config: $xml_relpath)")
    try
        cd(PHYSICELL_DIR) do
            open(log_path, "w") do log
                run(pipeline(`$EXECUTABLE $xml_relpath`, stdout=log, stderr=log))
            end
        end
        return :ok
    catch err
        println("FAILED $roi_name: $err")
        return :failed
    end
end

function main(; skip_existing::Bool=true, rebuild::Bool=true)
    rebuild && build_project()

    rois = discover_rois()
    println("Found $(length(rois)) ROI configs.\n")

    results = Dict(:ok => String[], :skipped => String[], :failed => String[])
    for (roi_name, xml_path) in rois
        status = run_roi(roi_name, xml_path; skip_existing=skip_existing)
        push!(results[status], roi_name)
    end

    println("\nDone. ok=$(length(results[:ok]))  skipped=$(length(results[:skipped]))  failed=$(length(results[:failed]))")
    isempty(results[:failed]) || println("Failed ROIs: ", join(results[:failed], ", "))
    return results
end

main()
