ENV["PHYSICELL_CPP"] = "g++-15"

using PhysiCellModelManager

initializeModelManager(
    joinpath(@__DIR__, "PhysiCell"),
    joinpath(@__DIR__, "data")
)

# ── Batch movie generation for all per-ROI output folders ─────────────────────
#
# Each simulation writes its SVG snapshots to PhysiCell/outputs/<ROI>/.
# This script converts those SVGs to JPEGs and then to an mp4, matching what
# PCMM's makeMovie does internally (make jpeg + make movie).
#
# Requirements: ImageMagick (magick) and ffmpeg must be on PATH.
#
# Usage: julia make_movies.jl
#        julia make_movies.jl JHH317ROI1 JHH317ROI2   # specific ROIs only
#

OUTPUTS_DIR = joinpath(@__DIR__, "PhysiCell", "outputs")
FRAMERATE   = 24

function make_movie_for(output_dir::String)
    svgs = filter(f -> startswith(f, "snapshot") && endswith(f, ".svg"),
                  readdir(output_dir))
    if isempty(svgs)
        println("  SKIP $(basename(output_dir)): no snapshot SVGs found")
        return false
    end

    mp4_path = joinpath(output_dir, "out.mp4")
    if isfile(mp4_path)
        println("  SKIP $(basename(output_dir)): out.mp4 already exists")
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

    println("  ✓ $(basename(output_dir)) → out.mp4")
    return true
end

# Determine which ROIs to process
roi_dirs = if !isempty(ARGS)
    [joinpath(OUTPUTS_DIR, a) for a in ARGS]
else
    filter(isdir, [joinpath(OUTPUTS_DIR, d) for d in readdir(OUTPUTS_DIR)])
end

if isempty(roi_dirs)
    println("No output folders found in $OUTPUTS_DIR")
    exit(0)
end

println("Processing $(length(roi_dirs)) output folder(s) in $OUTPUTS_DIR\n")
n_done = 0
n_skip = 0
errors  = String[]

for dir in roi_dirs
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
