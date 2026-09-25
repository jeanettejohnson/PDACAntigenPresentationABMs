# CAF contact -> MHC-II induction, shared by the four PCMM runners.
#
# Both rulesets carry eight rows for it (README, "CAF contact induces MHC-II"):
# a tumour cell touching a CAF or apCAF gains class II -- class1 ->
# class1_class2 and no class -> class2 -- in both lineages, at a fixed rate and
# with no way back. The rows' rate is 0 in the ruleset, so the baseline runs are
# unchanged. A runner enables them by setting its CAF_MHC2_RATE above 0, which
# adds the rate as a rules variation: enabled runs are then separate monads and
# never replace the baseline ones.
#
# `include` this after `initializeModelManager(...)`.

#: The four behaviours the rows drive. The CAF and apCAF rows of a behaviour
#: share one max_response, so one value sets both.
const CAF_MHC2_BEHAVIOURS = (
    ("epithelial_tumor_class1",  "transform to epithelial_tumor_class1_class2"),
    ("mesenchymal_tumor_class1", "transform to mesenchymal_tumor_class1_class2"),
    ("epithelial_tumor",         "transform to epithelial_tumor_class2"),
    ("mesenchymal_tumor",        "transform to mesenchymal_tumor_class2"),
)

"""
    cafMhc2Variations(rate)

The rules variation that sets the CAF-contact MHC-II rate to `rate` (1/min), as
a vector to splat into `createTrial`. Empty for the baseline, `rate == 0`.
"""
function cafMhc2Variations(rate::Real)
    rate >= 0 || error("CAF_MHC2_RATE must be >= 0, got $rate")
    rate == 0 && return AbstractVariation[]
    dvs = [DiscreteVariation(rulePath(cell_type, behaviour, "increasing_signals", "max_response"), Float64(rate))
           for (cell_type, behaviour) in CAF_MHC2_BEHAVIOURS]
    return AbstractVariation[CoVariation(dvs)]
end

"""
    checkBaseRulesets(collection)

Stop if `base_rulesets.xml` of the rulesets collection is older than its
`base_rulesets.csv`. PCMM builds the XML from the CSV only when the XML is
missing, so an XML left by an earlier run would silently hide later rule edits.
"""
function checkBaseRulesets(collection::AbstractString)
    folder = joinpath(@__DIR__, "..", "data", "inputs", "rulesets_collections", collection)
    csv = joinpath(folder, "base_rulesets.csv")
    xml = joinpath(folder, "base_rulesets.xml")
    if isfile(xml) && mtime(xml) < mtime(csv)
        error("""
              $(xml) is older than $(basename(csv)), so PCMM would run the old rules.
              It builds the XML from the CSV only when the XML is missing. Run from a
              fresh clone, or delete the XML (it is generated and gitignored) if no
              runs in this clone need their old rules.
              """)
    end
    return nothing
end
