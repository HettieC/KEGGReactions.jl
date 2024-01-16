module KEGGReactions

using HTTP, DocStringExtensions, Scratch, Serialization

const CACHE_DIRS = ["reaction", "reaction_metabolites"]
CACHE_LOCATION = ""

function __init__()
    global CACHE_LOCATION = @get_scratch!("kegg_data")

    for dir in CACHE_DIRS
        !isdir(joinpath(CACHE_LOCATION, dir)) && mkdir(joinpath(CACHE_LOCATION, dir))
    end

end

include("Types.jl")
include("Cache.jl")
include("Utils.jl")



export get_kegg_cmpd
export get_kegg_rxn
export clear_cache!

end # module KEGGReactions
