ENV["PHYSICELL_CPP"] = "g++"

using PhysiCellModelManager

initializeModelManager(
    joinpath(@__DIR__, "PhysiCell"),
    joinpath(@__DIR__, "data")
)

# ── Batch movie generation for PCMM simulation outputs ────────────────────────
#
# PCMM writes each simulation's SVG snapshots to
# data/outputs/simulations/<id>/output/. This converts those to JPEGs and then
# to an mp4, matching what PCMM's makeMovie does internally.
#
# Movies are written to movies/<id>_<ROI>.mp4 -- collected in one place rather
# than buried in each simulation's output folder, and labelled with the ROI so
# the numeric simulation id is not the only handle.
#
# Requirements: ImageMagick (magick) and ffmpeg, both provided by the
# physicell-sim-260606 conda environment.
#
# Usage: julia make_movies.jl
#        julia make_movies.jl 2 5      # specific simulation ids only
#

SIMULATIONS_DIR = joinpath(@__DIR__, "data", "outputs", "simulations")
MOVIES_DIR      = joinpath(@__DIR__, "movies")
FRAMERATE       = 24

"""Map simulation id -> ic_cell folder name (the ROI), for labelling movies."""
function roi_labels()
    labels = Dict{String,String}()
    try
        df = PhysiCellModelManager.constructSelectQuery(
            "simulations", "", "simulation_id, ic_cell_id") |> PhysiCellModelManager.queryToDataFrame
        ics = PhysiCellModelManager.constructSelectQuery(
            "ic_cells", "", "ic_cell_id, folder_name") |> PhysiCellModelManager.queryToDataFrame
        lookup = Dict(r.ic_cell_id => r.folder_name for r in eachrow(ics))
        for r in eachrow(df)
            haskey(lookup, r.ic_cell_id) && (labels[string(r.simulation_id)] = lookup[r.ic_cell_id])
        end
    catch e
        @warn "Could not read ROI labels from the database; naming movies by id only." exception=e
    end
    return labels
end

const ROI_LABELS = roi_labels()

function make_movie_for(sim_dir::String)
    sim_id = basename(sim_dir)
    # PCMM nests the PhysiCell output one level down.
    output_dir = joinpath(sim_dir, "output")
    if !isdir(output_dir)
        println("  SKIP $sim_id: no output/ folder")
        return false
    end

    svgs = filter(f -> startswith(f, "snapshot") && endswith(f, ".svg"),
                  readdir(output_dir))
    if isempty(svgs)
        println("  SKIP $sim_id: no snapshot SVGs found")
        return false
    end

    label = get(ROI_LABELS, sim_id, "")
    name  = isempty(label) ? sim_id : "$(sim_id)_$(label)"
    mkpath(MOVIES_DIR)
    mp4_path = joinpath(MOVIES_DIR, "$name.mp4")
    if isfile(mp4_path)
        println("  SKIP $sim_id: $(basename(mp4_path)) already exists")
        return false
    end

    # Read dimensions from first SVG to set resize target
    first_svg = joinpath(output_dir, sort(svgs)[1])
    h_str = readchomp(`magick identify -format "%h" $first_svg`)
    w_str = readchomp(`magick identify -format "%w" $first_svg`)
    h = 2 * div(parse(Int, h_str), 2)
    w = 2 * div(parse(Int, w_str), 2)
    resize = "$(w)x$(h)!"

    # Convert all snapshot SVGs to JPEG
    println("  Converting $(length(svgs)) SVGs to JPEG ($(w)×$(h))...")
    run(`magick mogrify -format jpg -resize $resize $(joinpath(output_dir, "s*.svg"))`)

    # Compile JPEGs into mp4
    println("  Compiling mp4...")
    run(`ffmpeg -r $FRAMERATE -f image2
        -i $(joinpath(output_dir, "snapshot%08d.jpg"))
        -vcodec libx264 -pix_fmt yuv420p -strict -2 -tune animation -crf 15 -acodec none
        $mp4_path`)

    # Clean up JPEGs
    for f in readdir(output_dir)
        if endswith(f, ".jpg")
            rm(joinpath(output_dir, f))
        end
    end

    println("  ✓ $sim_id → movies/$(basename(mp4_path))")
    return true
end

# Determine which simulations to process
sim_dirs = if !isempty(ARGS)
    [joinpath(SIMULATIONS_DIR, a) for a in ARGS]
elseif isdir(SIMULATIONS_DIR)
    # numeric sort so 2 comes before 10
    ids = filter(d -> isdir(joinpath(SIMULATIONS_DIR, d)) && all(isdigit, d),
                 readdir(SIMULATIONS_DIR))
    [joinpath(SIMULATIONS_DIR, d) for d in sort(ids; by=x -> parse(Int, x))]
else
    String[]
end

if isempty(sim_dirs)
    println("No simulation outputs found in $SIMULATIONS_DIR")
    exit(0)
end

println("Processing $(length(sim_dirs)) simulation(s) from $SIMULATIONS_DIR\n")
n_done = 0
n_skip = 0
errors  = String[]

for dir in sim_dirs
    println("── $(basename(dir)) ──────────────────────────────────────────")
    try
        result = make_movie_for(dir)
        result ? (n_done += 1) : (n_skip += 1)
    catch e
        println("  ERROR: $e")
        push!(errors, basename(dir))
    end
    println()
end

println("Done. $n_done movie(s) generated, $n_skip skipped.")
if !isempty(errors)
    println("Errors: $(join(errors, ", "))")
    exit(1)
end
